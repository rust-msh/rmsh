use glam::DVec3;
use rcad_kernel::BRep;
use rcad_kernel::PCurve;
use rcad_kernel::geom::*;

use crate::inttools::pcurve_derive::fallback_pcurve_by_projection;

/// Populates `brep.geom` with analytic geometry for a box BRep.
///
/// After this call, every edge has a `Curve3::Line` and every face has a `Surface3::Plane`.
/// Precondition: brep was created by `BRep::from_primitive(Box{..})`.
pub fn populate_box_geom(brep: &mut BRep) {
    let geom = &mut brep.geom;
    geom.curves.clear();
    geom.edge_curve.clear();
    geom.edge_curve_range.clear();
    geom.edge_degenerated.clear();
    geom.surfaces.clear();
    geom.face_surface.clear();

    // Edges → Line3
    for edge in &brep.edges {
        let p0 = brep.vertices[edge.start].point;
        let p1 = brep.vertices[edge.end].point;
        let delta = p1 - p0;
        let len = delta.length();
        let dir = if len > 1e-12 { delta / len } else { DVec3::X };
        let curve_idx = geom.curves.len();
        geom.curves.push(Curve3::Line(Line3 {
            origin: p0,
            direction: dir,
        }));
        geom.edge_curve.push(Some(curve_idx));
        // t_range: project endpoints onto the line
        let t0 = 0.0_f64;
        let t1 = (p1 - p0).dot(dir); // = len
        geom.edge_curve_range.push(Some([t0, t1]));
        geom.edge_degenerated.push(len <= 1e-12);
    }

    // Faces → Plane
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                // Use the first wire vertex rather than face.triangles (triangles
                // are rendering metadata and must not be used in modeling code).
                let origin = face
                    .outer_wire
                    .edges
                    .first()
                    .and_then(|we| brep.edges.get(we.idx))
                    .map(|e| brep.vertices[e.start].point)
                    .unwrap_or(DVec3::ZERO);
                let surf_idx = geom.surfaces.len();
                geom.surfaces.push(Surface3::Plane(Plane {
                    origin,
                    normal: face.normal,
                }));
                geom.face_surface.push(Some(surf_idx));
            }
        }
    }
}

/// Populate `edge_pcurves` for edges adjacent to curved faces that currently
/// lack a PCurve entry.
///
/// After a boolean operation, intersection edges on curved surfaces (cylinder,
/// sphere, cone, torus) often have no PCurve.  This function uses
/// [`fallback_pcurve_by_projection`] to derive a 2D parameter-space curve on
/// each adjacent curved surface and stores it in `brep.geom.edge_pcurves`.
///
/// Call this after [`boolean_op`] when downstream code needs PCurves
/// (e.g. parametric queries, STEP export of trimmed surfaces).
pub fn populate_boolean_result_pcurves(brep: &mut BRep) {
    // Collect all (edge_idx, face_idx) pairs where the face has a curved surface.
    // Use solids[0].shells[0] like the rest of the algorithms library.
    let face_indices: Vec<(usize, usize)> = {
        let Some(shell) = brep.solids.first().and_then(|s| s.shells.first()) else {
            return;
        };
        shell
            .faces
            .iter()
            .enumerate()
            .flat_map(|(fi, face)| {
                face.outer_wire
                    .edges
                    .iter()
                    .map(move |we| (we.idx, fi))
                    .chain(face.inner_wires.iter().flat_map(move |iw| {
                        iw.edges.iter().map(move |we| (we.idx, fi))
                    }))
            })
            .collect()
    };

    // Ensure edge_pcurves is large enough.
    let n_edges = brep.edges.len();
    if brep.geom.edge_pcurves.len() < n_edges {
        brep.geom.edge_pcurves.resize(n_edges, Vec::new());
    }

    for (edge_idx, face_idx) in face_indices {
        // Look up surface for this face.
        let surf_idx = match brep
            .geom
            .face_surface
            .get(face_idx)
            .and_then(|&si| si)
        {
            Some(si) => si,
            None => continue,
        };
        let surface = match brep.geom.surfaces.get(surf_idx) {
            Some(s) => s.clone(),
            None => continue,
        };
        // Only fill for curved surfaces.
        if matches!(surface, Surface3::Plane(_)) {
            continue;
        }

        // Check if a PCurve for this surface already exists on this edge.
        let already_has = brep.geom.edge_pcurves[edge_idx]
            .iter()
            .any(|pc| pc.surface_idx == surf_idx);
        if already_has {
            continue;
        }

        // Need a 3D curve or at least the two endpoint vertices to project.
        let (curve_opt, t_range_opt): (Option<Curve3>, Option<[f64; 2]>) = match brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|&ci| ci)
            .and_then(|ci| brep.geom.curves.get(ci))
        {
            Some(c) => {
                let r = brep
                    .geom
                    .edge_curve_range
                    .get(edge_idx)
                    .and_then(|r| *r)
                    .unwrap_or_else(|| c.default_domain());
                (Some(c.clone()), Some(r))
            }
            None => (None, None),
        };

        // Derive the 2D PCurve.
        let pcurve2d = if let (Some(curve), Some(t_range)) = (curve_opt, t_range_opt) {
            fallback_pcurve_by_projection(&curve, &t_range, &surface)
        } else {
            // No analytic curve: sample from vertex endpoints as a straight line.
            let Some(edge) = brep.edges.get(edge_idx) else {
                continue;
            };
            let Some(p0) = brep.vertices.get(edge.start).map(|v| v.point) else {
                continue;
            };
            let Some(p1) = brep.vertices.get(edge.end).map(|v| v.point) else {
                continue;
            };
            if (p1 - p0).length_squared() < 1e-20 {
                continue; // degenerate
            }
            // Project a polyline of 17 equally-spaced points between the endpoints.
            let polyline: Vec<_> = (0..17)
                .map(|i| p0 + (p1 - p0) * (i as f64 / 16.0))
                .collect();
            match crate::inttools::pcurve_derive::polyline_pcurve_by_projection(
                &polyline, &surface,
            ) {
                Some(c2d) => c2d,
                None => continue,
            }
        };

        // Store in geom.
        let curve2d_idx = brep.geom.curve2ds.len();
        brep.geom.curve2ds.push(pcurve2d);
        brep.geom.edge_pcurves[edge_idx].push(PCurve {
            surface_idx: surf_idx,
            curve2d_idx,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;

    #[test]
    fn box_geom_populated() {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 2.0,
            height: 3.0,
            depth: 4.0,
        });
        populate_box_geom(&mut brep);

        assert_eq!(brep.geom.edge_curve.len(), 12);
        assert!(brep.geom.edge_curve.iter().all(|c| c.is_some()));
        assert_eq!(brep.geom.face_surface.len(), 6);
        assert!(brep.geom.face_surface.iter().all(|s| s.is_some()));

        // All curves should be lines
        for c in &brep.geom.curves {
            assert!(matches!(c, Curve3::Line(_)));
        }
        // All surfaces should be planes
        for s in &brep.geom.surfaces {
            assert!(matches!(s, Surface3::Plane(_)));
        }
    }
}
