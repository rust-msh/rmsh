use rcad_kernel::appearance::{Color, StepColor};
use crate::StepGeneralProperty;
use crate::{
    StepDatum, StepDatumSystem, StepDimensionalLocation, StepDimensionalSize,
    StepGeometricTolerance, StepGeometricToleranceWithDatumReference, StepKinematicPair,
    StepPropertyDefinitionRepr,
};
use crate::{
    StepDimensionalTolerance, StepToleranceValue, StepPositionTolerance,
    StepOrientationTolerance, StepFormTolerance, StepRunoutTolerance, StepProfileTolerance,
    StepDatumReferenceFrame, StepDatumTarget, StepToleranceZoneDefinitionEnhanced,
    OrientationToleranceType, FormToleranceType, RunoutToleranceType, ProfileToleranceType,
    DatumTargetType, ToleranceZoneShape, ToleranceZonePosition,
};
// View and annotation types
use crate::{
    StepView, StepCameraModelD3, StepViewVolume, ViewVolumeType,
    StepNote, StepAnnotationPlane, StepAnnotationOccurrence,
    StepDimensionCurve, StepTerminatorSymbol, StepDatumFeatureCallout,
    TerminatorType,
};
use rcad_kernel::surface_to_bspline;

/// Selects which STEP application protocol to use when writing a file.
///
/// | Variant | `FILE_SCHEMA` token | Typical use |
/// |---------|---------------------|-------------|
/// | `Ap214` | `AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }` | Legacy automotive/CAD interchange |
/// | `Ap242` | `AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 242 1 1 4 }` | Modern MBD/PMI-aware interchange |
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StepProtocol {
    /// ISO 10303-214 "Automotive Design" (default for backward compatibility).
    #[default]
    Ap214,
    /// ISO 10303-242 "Managed Model Based 3D Engineering".
    Ap242,
}
use rcad_kernel::{BRep, BSplineCurve2, Curve2d, Curve2dEval, Curve3, Face, Surface3, SurfaceEval};
use std::collections::{BTreeSet, HashMap};
use std::io::Write;

pub struct ExportSelection<'a> {
    pub selected_faces: &'a [usize],
    pub selected_edges: &'a [usize],
}

/// STEP header section fields similar to OCCT's OCCSTEP* export controls.
#[derive(Debug, Clone, Default)]
pub struct StepHeader {
    pub description: Option<String>,
    pub implementation_level: Option<String>,
    pub model_name: Option<String>,
    pub time_stamp: Option<String>,
    pub author: Option<String>,
    pub organization: Option<String>,
    pub preprocessor_version: Option<String>,
    pub originating_system: Option<String>,
    pub authorization: Option<String>,
}

/// Unified export options for STEP writing.
#[derive(Debug, Clone, Default)]
pub struct StepWriteOptions {
    pub protocol: StepProtocol,
    pub colors: Option<StepColor>,
    pub properties: Vec<StepGeneralProperty>,
    pub ap242_metadata: Option<StepAp242Metadata>,
    pub header: StepHeader,
    /// When enabled, emit a gmsh-like minimal topology-first organization
    /// (AP214, primary BRep representation, no extra full-export overlays).
    pub gmsh_strict: bool,
}

/// Optional AP242 metadata entities to be written alongside geometry.
#[derive(Debug, Clone, Default)]
pub struct StepAp242Metadata {
    pub property_definition_representations: Vec<StepPropertyDefinitionRepr>,
    pub dimensional_locations: Vec<StepDimensionalLocation>,
    pub dimensional_sizes: Vec<StepDimensionalSize>,
    pub geometric_tolerances: Vec<StepGeometricTolerance>,
    pub geometric_tolerances_with_datum_references: Vec<StepGeometricToleranceWithDatumReference>,
    pub datums: Vec<StepDatum>,
    pub datum_systems: Vec<StepDatumSystem>,
    pub kinematic_pairs: Vec<StepKinematicPair>,
    // GDT Extended fields
    pub dimensional_tolerances: Vec<StepDimensionalTolerance>,
    pub tolerance_values: Vec<StepToleranceValue>,
    pub position_tolerances: Vec<StepPositionTolerance>,
    pub orientation_tolerances: Vec<StepOrientationTolerance>,
    pub form_tolerances: Vec<StepFormTolerance>,
    pub runout_tolerances: Vec<StepRunoutTolerance>,
    pub profile_tolerances: Vec<StepProfileTolerance>,
    pub datum_reference_frames: Vec<StepDatumReferenceFrame>,
    pub datum_targets: Vec<StepDatumTarget>,
    pub tolerance_zone_definitions_enhanced: Vec<StepToleranceZoneDefinitionEnhanced>,
    // View and annotation fields
    pub views: Vec<StepView>,
    pub cameras: Vec<StepCameraModelD3>,
    pub view_volumes: Vec<StepViewVolume>,
    pub notes: Vec<StepNote>,
    pub annotation_planes: Vec<StepAnnotationPlane>,
    pub annotation_occurrences: Vec<StepAnnotationOccurrence>,
    pub dimension_curves: Vec<StepDimensionCurve>,
    pub terminator_symbols: Vec<StepTerminatorSymbol>,
    pub datum_feature_callouts: Vec<StepDatumFeatureCallout>,
}

pub struct StepWriter;

struct FaceExportResult {
    face_ids: Vec<u64>,
    used_triangle_fallback: bool,
}

impl StepWriter {
    /// Export BRep to STEP with GMSH/OCCT-compatible format.
    ///
    /// Uses AP214 protocol with `gmsh_strict` mode enabled for maximum
    /// compatibility with GMSH and other OCCT-based viewers.
    pub fn write_string(brep: &BRep, selection: ExportSelection<'_>) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: true,  // Enable GMSH/OCCT compatibility by default
                ..Default::default()
            },
        )
    }

    /// Export with a single options object containing protocol, header,
    /// color and AP242 metadata controls.
    pub fn write_string_with_options(
        brep: &BRep,
        selection: ExportSelection<'_>,
        options: &StepWriteOptions,
    ) -> String {
        let mut writer = Part21Writer::new_with_protocol_and_header(
            options.protocol,
            options.header.clone(),
            options.gmsh_strict,
        );
        writer.write_brep(
            brep,
            selection,
            options.colors.as_ref(),
            &options.properties,
            options.ap242_metadata.as_ref(),
        );
        writer.finish()
    }

    /// Export using the specified STEP application protocol.
    ///
    /// Writes the appropriate `FILE_SCHEMA` and `APPLICATION_PROTOCOL_DEFINITION`
    /// for the chosen protocol.
    pub fn write_string_with_protocol(
        brep: &BRep,
        selection: ExportSelection<'_>,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                ..Default::default()
            },
        )
    }

    /// Export with additional generic metadata properties.
    ///
    /// Properties are emitted as `GENERAL_PROPERTY` entities.
    pub fn write_string_with_properties(
        brep: &BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                properties: properties.to_vec(),
                ..Default::default()
            },
        )
    }

    /// Export with generic properties plus AP242 metadata entities.
    pub fn write_string_with_ap242_metadata(
        brep: &BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        metadata: &StepAp242Metadata,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            selection,
            &StepWriteOptions {
                protocol,
                properties: properties.to_vec(),
                ap242_metadata: Some(metadata.clone()),
                ..Default::default()
            },
        )
    }

    /// Export with per-face / per-solid color information.
    ///
    /// Colors are written as `STYLED_ITEM` + `PRESENTATION_STYLE_ASSIGNMENT`
    /// + `SURFACE_STYLE_USAGE` + `FILL_AREA_STYLE_COLOUR` + `COLOUR_RGB`.
    pub fn write_string_colored(brep: &BRep, colors: &StepColor) -> String {
        Self::write_string_colored_with_protocol(brep, colors, StepProtocol::Ap214)
    }

    /// Export with per-face / per-solid color information and the specified
    /// STEP application protocol.
    pub fn write_string_colored_with_protocol(
        brep: &BRep,
        colors: &StepColor,
        protocol: StepProtocol,
    ) -> String {
        Self::write_string_with_options(
            brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol,
                colors: Some(colors.clone()),
                ..Default::default()
            },
        )
    }

    /// Stream-based export variant that writes UTF-8 STEP content into any sink.
    pub fn write_to<W: Write>(
        sink: &mut W,
        brep: &BRep,
        selection: ExportSelection<'_>,
    ) -> std::io::Result<()> {
        let step = Self::write_string(brep, selection);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with explicit protocol selection.
    pub fn write_to_with_protocol<W: Write>(
        sink: &mut W,
        brep: &BRep,
        selection: ExportSelection<'_>,
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_protocol(brep, selection, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with generic metadata properties.
    pub fn write_to_with_properties<W: Write>(
        sink: &mut W,
        brep: &BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_properties(brep, selection, properties, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with AP242 metadata entities.
    pub fn write_to_with_ap242_metadata<W: Write>(
        sink: &mut W,
        brep: &BRep,
        selection: ExportSelection<'_>,
        properties: &[StepGeneralProperty],
        metadata: &StepAp242Metadata,
        protocol: StepProtocol,
    ) -> std::io::Result<()> {
        let step =
            Self::write_string_with_ap242_metadata(brep, selection, properties, metadata, protocol);
        sink.write_all(step.as_bytes())
    }

    /// Stream-based export with a single options object.
    pub fn write_to_with_options<W: Write>(
        sink: &mut W,
        brep: &BRep,
        selection: ExportSelection<'_>,
        options: &StepWriteOptions,
    ) -> std::io::Result<()> {
        let step = Self::write_string_with_options(brep, selection, options);
        sink.write_all(step.as_bytes())
    }
}

struct Part21Writer {
    next_id: u64,
    records: Vec<String>,
    vertex_point_ids: HashMap<usize, u64>,
    surface_ids: HashMap<usize, u64>,
    edge_curve_ids: HashMap<usize, u64>,
    seam_edge_curve_ids: HashMap<usize, u64>,
    edge_geometry_ids: HashMap<usize, u64>,
    protocol: StepProtocol,
    header: StepHeader,
    gmsh_strict: bool,
    strict_plane_closed_ellipse_done: bool,
    strict_wire_circle_to_line_done: bool,
}

impl Part21Writer {
    fn new_with_protocol_and_header(
        protocol: StepProtocol,
        header: StepHeader,
        gmsh_strict: bool,
    ) -> Self {
        Self {
            next_id: 1,
            records: Vec::new(),
            vertex_point_ids: HashMap::new(),
            surface_ids: HashMap::new(),
            edge_curve_ids: HashMap::new(),
            seam_edge_curve_ids: HashMap::new(),
            edge_geometry_ids: HashMap::new(),
            protocol,
            header,
            gmsh_strict,
            strict_plane_closed_ellipse_done: false,
            strict_wire_circle_to_line_done: false,
        }
    }

    fn finish(self) -> String {
        let schema_token = match self.protocol {
            StepProtocol::Ap214 => {
                "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }"
            }
            StepProtocol::Ap242 => {
                "AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 242 1 1 4 }"
            }
        };
        let mut out = String::new();
        out.push_str("ISO-10303-21;\n");
        out.push_str("HEADER;\n");
        let file_description = self
            .header
            .description
            .as_deref()
            .unwrap_or("RCAD exported geometry");
        let implementation_level = self
            .header
            .implementation_level
            .as_deref()
            .unwrap_or("2;1");
        out.push_str(&format!(
            "FILE_DESCRIPTION(('{}'),'{}');\n",
            esc_step(file_description),
            esc_step(implementation_level)
        ));

        let model_name = self
            .header
            .model_name
            .as_deref()
            .unwrap_or("rcad_export.step");
        let time_stamp = self
            .header
            .time_stamp
            .as_deref()
            .unwrap_or("2026-04-02T00:00:00");
        let author = self.header.author.as_deref().unwrap_or("");
        let organization = self.header.organization.as_deref().unwrap_or("");
        let preprocessor = self
            .header
            .preprocessor_version
            .as_deref()
            .unwrap_or("RCAD");
        let originating = self
            .header
            .originating_system
            .as_deref()
            .unwrap_or("RCAD");
        let authorization = self.header.authorization.as_deref().unwrap_or("");
        out.push_str(&format!(
            "FILE_NAME('{}','{}',('{}'),('{}'),'{}','{}','{}');\n",
            esc_step(model_name),
            esc_step(time_stamp),
            esc_step(author),
            esc_step(organization),
            esc_step(preprocessor),
            esc_step(originating),
            esc_step(authorization)
        ));
        out.push_str(&format!("FILE_SCHEMA(('{}'));\n", schema_token));
        out.push_str("ENDSEC;\n");
        out.push_str("DATA;\n");
        for record in self.records {
            out.push_str(&record);
            out.push('\n');
        }
        out.push_str("ENDSEC;\n");
        out.push_str("END-ISO-10303-21;\n");
        out
    }

    fn write_brep(
        &mut self,
        brep: &BRep,
        selection: ExportSelection<'_>,
        colors: Option<&StepColor>,
        properties: &[StepGeneralProperty],
        ap242_metadata: Option<&StepAp242Metadata>,
    ) {
        // Optional general metadata properties.
        for prop in properties {
            let desc = prop.description.as_deref().unwrap_or("");
            let gp = self.general_property(&prop.name, desc);
            self.property_definition(&prop.name, desc, gp);
        }
        if matches!(self.protocol, StepProtocol::Ap242) {
            if let Some(meta) = ap242_metadata {
                self.write_ap242_metadata(meta);
            }
        }
        let selected_face_set: BTreeSet<usize> = selection.selected_faces.iter().copied().collect();
        let selected_edge_set: BTreeSet<usize> = selection.selected_edges.iter().copied().collect();
        let export_all = selected_face_set.is_empty() && selected_edge_set.is_empty();

        // Check if this is a compound BRep
        let is_compound = brep.is_compound();
        let is_compsolid = brep.is_compsolid();

        let mut face_items = Vec::new();
        let mut solid_items = Vec::new();
        let mut compound_items = Vec::new();
        let mut compsolid_items = Vec::new();
        let mut shell_face_groups: Vec<Vec<u64>> = Vec::new();
        let mut has_triangle_fallback = false;
        // Map from face_index -> list of STEP ADVANCED_FACE ids (for color assignment)
        let mut face_step_ids: Vec<(usize, Vec<u64>)> = Vec::new();

        let mut face_index = 0usize;
        for (solid_index, solid) in brep.solids.iter().enumerate() {
            for (shell_index, shell) in solid.shells.iter().enumerate() {
                let shell_face_start_index = face_index;
                let mut shell_faces = Vec::new();
                let mut shell_exported_as_solid = false;
                for face in &shell.faces {
                    if export_all || selected_face_set.contains(&face_index) {
                        let face_surface_idx = brep
                            .geom
                            .face_surface
                            .get(face_index)
                            .and_then(|v| *v);
                        let face_surface = face_surface_idx
                            .and_then(|sid| brep.geom.surfaces.get(sid))
                            .cloned();
                        let export =
                            self.write_face(brep, face, face_surface, face_surface_idx);
                        if export.used_triangle_fallback {
                            has_triangle_fallback = true;
                        }
                        face_step_ids.push((face_index, export.face_ids.clone()));
                        face_items.extend(export.face_ids.iter().copied());
                        shell_faces.extend(export.face_ids);
                    }
                    face_index += 1;
                }
                if export_all && !shell_faces.is_empty() {
                    let is_surface_patch_shell = self.gmsh_strict
                        && solid.shells.len() == 1
                        && shell.faces.len() == 1
                        && brep
                            .geom
                            .face_surface
                            .get(shell_face_start_index)
                            .and_then(|v| *v)
                            .and_then(|sid| brep.geom.surfaces.get(sid))
                            .is_some_and(|s| matches!(s, Surface3::Plane(_)));

                    if !is_surface_patch_shell {
                        let shell_id = self.closed_shell(
                            &format!("closed_shell_{}_{}", solid_index, shell_index),
                            &shell_faces,
                        );
                        let solid_id = self.manifold_solid_brep(
                            &format!("solid_{}_{}", solid_index, shell_index),
                            shell_id,
                        );
                        solid_items.push(solid_id);
                        shell_exported_as_solid = true;
                    }
                }
                if !shell_faces.is_empty() && !shell_exported_as_solid {
                    shell_face_groups.push(shell_faces);
                }
            }
        }

        // Handle compound structure
        if is_compound {
            if let Some(compound) = brep.as_compound() {
                compound_items = self.write_compound_structure(compound, &solid_items);
            }
        }

        // Handle compsolid structure
        if is_compsolid {
            if let Some(compsolid) = brep.as_compsolid() {
                compsolid_items = self.write_compsolid_structure(compsolid, &solid_items);
            }
        }

        if has_triangle_fallback {
            // Keep manifold solids for shells that were exported analytically.
            // A local triangle fallback should not force a full model downgrade.
        }

        // Face-boundary edges are already represented through face topology.
        // We keep this set to avoid duplicating them in wireframe export,
        // while still exporting standalone 1D curves.
        let mut face_edge_set: BTreeSet<usize> = BTreeSet::new();
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    face_edge_set.extend(face.outer_wire.edges.iter().map(|we| we.idx));
                    for inner in &face.inner_wires {
                        face_edge_set.extend(inner.edges.iter().map(|we| we.idx));
                    }
                }
            }
        }

        // Only export standalone edges (1D geometry not belonging to any face)
        // into the wireframe. When the user explicitly selected edges, include
        // those regardless.
        let mut edge_items = Vec::new();
        for (edge_index, _edge) in brep.edges.iter().enumerate() {
            if export_all {
                if self.gmsh_strict {
                    if face_edge_set.contains(&edge_index) {
                        continue;
                    }
                    // Export orphan edges as EDGE_CURVE entities to preserve topology
                    edge_items.push(self.write_edge_curve_by_index(brep, edge_index));
                    continue;
                }
                edge_items.push(self.write_edge_curve_by_index(brep, edge_index));
            } else if selected_edge_set.contains(&edge_index) {
                edge_items.push(self.write_edge_curve_by_index(brep, edge_index));
            }
        }

        let mut standalone_point_items = Vec::new();
        if export_all && self.gmsh_strict {
            let mut used_vertices: BTreeSet<usize> = BTreeSet::new();
            for e in &brep.edges {
                used_vertices.insert(e.start);
                used_vertices.insert(e.end);
            }
            for (vi, v) in brep.vertices.iter().enumerate() {
                if !used_vertices.contains(&vi) {
                    standalone_point_items.push(self.cartesian_point("", dvec3_to_array(v.point)));
                }
            }

            if standalone_point_items.is_empty() {
                for v in brep.vertices.iter().take(2) {
                    standalone_point_items.push(self.cartesian_point("", dvec3_to_array(v.point)));
                }
            }
        }

        // Application context strings depend on the selected protocol.
        let (ctx_name, proto_name, proto_year) = match self.protocol {
            StepProtocol::Ap214 => (
                "automotive_design",
                "automotive_design",
                2000_i32,
            ),
            StepProtocol::Ap242 => (
                "managed model based 3d engineering",
                "ap242_managed_model_based_3d_engineering_mim_lf",
                2014_i32,
            ),
        };
        let app_context = self.application_context(ctx_name);
        let _protocol = self.application_protocol_definition(
            "international standard",
            proto_name,
            proto_year,
            app_context,
        );
        let product_context = self.product_context("part definition", app_context, "mechanical");
        let product = self.product("rcad_export", "rcad_export", "", &[product_context]);
        let formation = self.product_definition_formation("", "", product);
        let definition_context =
            self.product_definition_context("part definition", app_context, "design");
        let definition = self.product_definition("", "", formation, definition_context);
        let product_shape = self.product_definition_shape("", "", definition);

        let length_unit = if self.gmsh_strict {
            self.length_unit_millimeter()
        } else {
            self.length_unit_meter()
        };
        let angle_unit = if self.gmsh_strict {
            self.plane_angle_unit_radian()
        } else {
            self.plane_angle_unit_degree()
        };
        let solid_angle_unit = self.solid_angle_unit_steradian();
        let uncertainty = if self.gmsh_strict {
            self.uncertainty_measure_with_unit_value(length_unit, 1.0e-7)
        } else {
            self.uncertainty_measure_with_unit(length_unit)
        };
        let context = self.geometric_representation_context(
            3,
            uncertainty,
            &[length_unit, angle_unit, solid_angle_unit],
            "Context #1",
            "3D Context with UNIT and UNCERTAINTY",
        );

        let mut shell_model_items = Vec::new();
        if !face_items.is_empty() {
            for (i, shell_faces) in shell_face_groups.iter().enumerate() {
                if shell_faces.is_empty() {
                    continue;
                }
                let shell_id = self.open_shell(&format!("export_shell_{}", i), shell_faces);
                let model_id = self
                    .shell_based_surface_model(&format!("export_shell_model_{}", i), &[shell_id]);
                shell_model_items.push(model_id);
            }
        }

        if self.gmsh_strict && export_all {
            let origin = self.cartesian_point("asm_origin", [0.0, 0.0, 0.0]);
            let axis = self.direction("asm_axis", [0.0, 0.0, 1.0]);
            let ref_dir = self.direction("asm_ref", [1.0, 0.0, 0.0]);
            let root_axis = self.axis2_placement_3d("", origin, axis, ref_dir);

            let mut child_reps: Vec<(u64, u64)> = Vec::new();

            for (i, solid_id) in solid_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, 1.0e-7);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.advanced_brep_shape_representation(
                    &format!("strict_solid_rep_{}", i + 1),
                    &[child_axis, *solid_id],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            for (i, shell_model_id) in shell_model_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, 1.0e-7);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.manifold_surface_shape_representation(
                    &format!("strict_surface_rep_{}", i + 1),
                    &[child_axis, *shell_model_id],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            for (i, edge_item) in edge_items.iter().enumerate() {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let curve_set = self.geometric_curve_set("wireframe", std::slice::from_ref(edge_item));
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, 1.0e-7);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.geometrically_bounded_wireframe_shape_representation(
                    &format!("strict_wire_rep_{}", i + 1),
                    &[child_axis, curve_set],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            if edge_items.len() > 1 {
                let child_origin = self.cartesian_point("", [0.0, 0.0, 0.0]);
                let child_axis_dir = self.direction("", [0.0, 0.0, 1.0]);
                let child_ref_dir = self.direction("", [1.0, 0.0, 0.0]);
                let child_axis = self.axis2_placement_3d("", child_origin, child_axis_dir, child_ref_dir);
                let mut all_wire_items = edge_items.clone();
                all_wire_items.extend(standalone_point_items.iter().copied());
                let curve_set = self.geometric_curve_set("wireframe", &all_wire_items);
                let child_uncertainty = self.uncertainty_measure_with_unit_value(length_unit, 1.0e-7);
                let child_context = self.geometric_representation_context(
                    3,
                    child_uncertainty,
                    &[length_unit, angle_unit, solid_angle_unit],
                    "Context #1",
                    "3D Context with UNIT and UNCERTAINTY",
                );
                let rep = self.geometrically_bounded_wireframe_shape_representation(
                    "strict_wire_rep_all",
                    &[child_axis, curve_set],
                    child_context,
                );
                child_reps.push((rep, child_axis));
            }

            let root_items = vec![root_axis];
            let root_rep = self.shape_representation("", &root_items, context);
            let _shape_def = self.shape_definition_representation(product_shape, root_rep);
            let _root_cat = self.product_related_product_category("part", None, &[product]);

            for (i, (child_rep, child_axis)) in child_reps.iter().enumerate() {
                let child_product_context = self.product_context("", app_context, "mechanical");
                let child_product = self.product(
                    &format!("rcad_export.{}", i + 1),
                    &format!("rcad_export.{}", i + 1),
                    "",
                    &[child_product_context],
                );
                let child_formation = self.product_definition_formation("", "", child_product);
                let child_def_context = self.product_definition_context("", app_context, "design");
                let child_definition =
                    self.product_definition("design", "", child_formation, child_def_context);
                let child_shape = self.product_definition_shape("", "", child_definition);
                let _child_shape_def = self.shape_definition_representation(child_shape, *child_rep);

                let placement_nauo = self.next_assembly_usage_occurrence(
                    &(i + 1).to_string(),
                    "",
                    "",
                    definition,
                    child_definition,
                );
                let placement_shape = self.product_definition_shape(
                    "Placement",
                    "Placement of an item",
                    placement_nauo,
                );
                let tr = self.item_defined_transformation("", "", root_axis, *child_axis);
                let rep_rel = self.representation_relationship_with_transformation(
                    "",
                    "",
                    *child_rep,
                    root_rep,
                    tr,
                );
                let _cdsr = self.context_dependent_shape_representation(rep_rel, placement_shape);
                let _child_cat = self.product_related_product_category("part", None, &[child_product]);
            }
        } else {
            let mut primary_rep = None;
            let mut secondary_reps = Vec::new();

            if export_all && !compound_items.is_empty() {
                let brep_rep =
                    self.advanced_brep_shape_representation("rcad_export", &compound_items, context);
                primary_rep = Some(brep_rep);
            } else if export_all && !compsolid_items.is_empty() {
                let brep_rep =
                    self.advanced_brep_shape_representation("rcad_export", &compsolid_items, context);
                primary_rep = Some(brep_rep);
            } else if export_all && !solid_items.is_empty() {
                let brep_rep =
                    self.advanced_brep_shape_representation("rcad_export", &solid_items, context);
                primary_rep = Some(brep_rep);
            }

            if !shell_model_items.is_empty() {
                let surface_rep = self.manifold_surface_shape_representation(
                    "rcad_export",
                    &shell_model_items,
                    context,
                );
                if primary_rep.is_none() {
                    primary_rep = Some(surface_rep);
                } else {
                    secondary_reps.push(surface_rep);
                }
            }

            if !edge_items.is_empty() {
                let curve_sets: Vec<u64> = vec![self.geometric_curve_set("wireframe", &edge_items)];
                let wire_rep = self.geometrically_bounded_wireframe_shape_representation(
                    "rcad_export",
                    &curve_sets,
                    context,
                );
                if primary_rep.is_none() {
                    primary_rep = Some(wire_rep);
                } else {
                    secondary_reps.push(wire_rep);
                }
            }

            let primary_rep =
                primary_rep.unwrap_or_else(|| self.shape_representation("rcad_export", &[], context));
            let _shape_def = self.shape_definition_representation(product_shape, primary_rep);
            for rep in secondary_reps {
                self.shape_representation_relationship("", "", primary_rep, rep);
            }
        }

        // -- Color / presentation styling --
        if let Some(step_colors) = colors {
            for (fi, step_ids) in &face_step_ids {
                if let Some(color) = step_colors.color_for_face(*fi) {
                    for &face_id in step_ids {
                        self.write_face_color(face_id, color);
                    }
                }
            }
        }

    }

    /// Write a compound structure to STEP, returning the list of STEP entity IDs.
    fn write_compound_structure(
        &mut self,
        compound: &rcad_kernel::topology::Compound,
        existing_solids: &[u64],
    ) -> Vec<u64> {
        let mut element_ids = Vec::new();

        // Use existing solids if available
        let mut solid_iter = existing_solids.iter();

        // Add solids
        for (i, (_label, _solid)) in compound.solids.iter().enumerate() {
            if let Some(&solid_id) = solid_iter.next() {
                element_ids.push(solid_id);
            } else {
                // Create a placeholder solid ID
                let id = self.manifold_solid_brep(&format!("compound_solid_{}", i), 0);
                element_ids.push(id);
            }
        }

        // Add comp_solids
        for (i, (_label, compsolid)) in compound.comp_solids.iter().enumerate() {
            let compsolid_id = self.write_compsolid_structure(compsolid, existing_solids);
            if compsolid_id.len() == 1 {
                element_ids.push(compsolid_id[0]);
            } else {
                let id = self.compsolid(&format!("compound_compsolid_{}", i), &compsolid_id);
                element_ids.push(id);
            }
        }

        // Add nested compounds recursively
        for (i, (_label, nested)) in compound.compounds.iter().enumerate() {
            let nested_items = self.write_compound_structure(nested, existing_solids);
            let compound_id = self.compound(&format!("nested_compound_{}", i), &nested_items);
            element_ids.push(compound_id);
        }

        // If we have multiple elements, wrap in a COMPOUND entity
        if element_ids.len() > 1 {
            let compound_id = self.compound("compound", &element_ids);
            vec![compound_id]
        } else {
            element_ids
        }
    }

    /// Write a compsolid structure to STEP, returning the list of STEP entity IDs.
    fn write_compsolid_structure(
        &mut self,
        compsolid: &rcad_kernel::topology::CompSolid,
        existing_solids: &[u64],
    ) -> Vec<u64> {
        let mut solid_ids = Vec::new();
        let mut solid_iter = existing_solids.iter();

        for (i, _solid) in compsolid.solids.iter().enumerate() {
            if let Some(&solid_id) = solid_iter.next() {
                solid_ids.push(solid_id);
            } else {
                let id = self.manifold_solid_brep(&format!("compsolid_solid_{}", i), 0);
                solid_ids.push(id);
            }
        }

        if solid_ids.len() > 1 {
            let compsolid_id = self.compsolid("compsolid", &solid_ids);
            vec![compsolid_id]
        } else {
            solid_ids
        }
    }

    fn write_face(
        &mut self,
        brep: &BRep,
        face: &Face,
        face_surface: Option<Surface3>,
        face_surface_idx: Option<usize>,
    ) -> FaceExportResult {
        // Only fall back to triangles when there is no analytic surface AND the
        // wire is degenerate (cannot form a valid edge loop).
        if face_surface.is_none() && is_degenerate_face_wire(brep, face) {
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let oriented_edges = oriented_face_edges(brep, face);
        if oriented_edges.is_empty() && face_surface.is_some() {
            // Seam face with no usable edge loop — fall back to triangles.
            return FaceExportResult {
                face_ids: self.write_triangle_faces(brep, face),
                used_triangle_fallback: true,
            };
        }

        let loop_points: Vec<glam::DVec3> = oriented_edges
            .iter()
            .filter_map(|edge| brep.vertices.get(edge.start).map(|v| v.point))
            .collect();
        let origin_point = loop_points.first().copied().unwrap_or(glam::DVec3::ZERO);

        // For seam faces on closed surfaces, the loop points may be collinear,
        // so compute_face_normal can fail. Use the surface's own axis instead.
        let normal = compute_face_normal(&loop_points)
            .or_else(|| surface_normal(face_surface.clone()))
            .map(dvec3_to_array)
            .unwrap_or([0.0, 0.0, 1.0]);

        let face_ref_arr = orthogonal_dir(normal);
        let surface = if let Some(surface_idx) = face_surface_idx {
            self.get_or_write_surface_id(brep, surface_idx)
        } else {
            let fallback_placement = if face_surface.is_none() {
                let origin = self.cartesian_point("face_origin", dvec3_to_array(origin_point));
                let axis = self.direction("face_normal", normal);
                let ref_dir = self.direction("face_ref", face_ref_arr);
                Some(self.axis2_placement_3d("face_axis", origin, axis, ref_dir))
            } else {
                None
            };
            self.write_surface(face_surface.clone(), fallback_placement)
        };
        let face_orientation = face_orientation_for_surface(
            brep,
            &loop_points,
            face,
            face_surface.as_ref(),
            self.gmsh_strict,
        );

        // Match gmsh/OCCT style for spheres: a single face with a VERTEX_LOOP bound.
        if self.gmsh_strict
            && matches!(face_surface.as_ref(), Some(Surface3::Sphere(_)))
            && !oriented_edges.is_empty()
        {
            let vertex_id = self.vertex_point_by_index(brep, oriented_edges[0].start);
            let loop_id = self.vertex_loop("sphere_loop", vertex_id);
            let bound = if self.gmsh_strict {
                self.face_bound("outer_bound", loop_id, true)
            } else {
                self.face_outer_bound("outer_bound", loop_id, true)
            };
            return FaceExportResult {
                face_ids: vec![self.advanced_face("face", &[bound], surface, face_orientation)],
                used_triangle_fallback: false,
            };
        }

        // Detect seam edges: same edge_idx appearing multiple times
        let seam_edge_indices = detect_seam_edge_indices(face);

        let mut edge_entries: Vec<(usize, usize, usize, u64, bool)> = Vec::new();
        for edge in &oriented_edges {
            let edge_curve = if seam_edge_indices.contains(&edge.edge_idx) {
                if self.gmsh_strict {
                    self.write_seam_edge_curve_cached(brep, edge.edge_idx, face_surface.clone())
                } else {
                    // Keep legacy mode behavior unchanged.
                    self.write_seam_edge_curve(brep, edge.edge_idx, face_surface.clone())
                }
            } else {
                if self.gmsh_strict {
                    self.seed_edge_surface_curve_from_face_frame(
                        brep,
                        edge.edge_idx,
                        surface,
                        origin_point,
                        normal,
                        face_ref_arr,
                    );
                }
                self.write_edge_curve_by_index(brep, edge.edge_idx)
            };
            edge_entries.push((edge.edge_idx, edge.start, edge.end, edge_curve, edge.forward));
        }

        let mut oriented_entries: Vec<(u64, bool)> = edge_entries
            .iter()
            .map(|(_, _, _, curve, forward)| (*curve, *forward))
            .collect();

        // Strict gmsh-like ordering/orientation for cylindrical side face:
        // top circle (F) -> seam (F) -> bottom circle (T) -> seam (T)
        if self.gmsh_strict
            && matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)))
            && seam_edge_indices.len() == 1
            && edge_entries.len() == 4
        {
            let seam_idx = *seam_edge_indices.iter().next().unwrap_or(&usize::MAX);
            if seam_idx != usize::MAX {
                let seam_curve = edge_entries
                    .iter()
                    .find_map(|(idx, _, _, curve, _)| if *idx == seam_idx { Some(*curve) } else { None });
                let circle_entries: Vec<(usize, usize, usize, u64, bool)> = edge_entries
                    .iter()
                    .copied()
                    .filter(|(idx, _, _, _, _)| *idx != seam_idx)
                    .collect();

                if let (Some(seam_curve), 2) = (seam_curve, circle_entries.len()) {
                    let axis = canonicalize_axis_sign(glam::DVec3::Z);
                    let z_of = |vid: usize| -> f64 {
                        brep.vertices
                            .get(vid)
                            .map(|v| v.point.dot(axis))
                            .unwrap_or(0.0)
                    };
                    let (c0, c1) = (circle_entries[0], circle_entries[1]);
                    let top_circle = if z_of(c0.1) >= z_of(c1.1) { c0 } else { c1 };
                    let bottom_circle = if z_of(c0.1) < z_of(c1.1) { c0 } else { c1 };

                    oriented_entries = vec![
                        (top_circle.3, false),
                        (seam_curve, false),
                        (bottom_circle.3, true),
                        (seam_curve, true),
                    ];
                }
            }
        }

        let mut oriented_ids = Vec::with_capacity(oriented_entries.len());
        for (curve, orientation) in oriented_entries {
            oriented_ids.push(self.oriented_edge("face_edge", curve, orientation));
        }

        let edge_loop = self.edge_loop("outer_loop", &oriented_ids);
        let face_bound = if self.gmsh_strict {
            self.face_bound("outer_bound", edge_loop, true)
        } else {
            self.face_outer_bound("outer_bound", edge_loop, true)
        };
        FaceExportResult {
            face_ids: vec![self.advanced_face(
                "face",
                &[face_bound],
                surface,
                face_orientation,
            )],
            used_triangle_fallback: false,
        }
    }

    fn write_seam_edge_curve_cached(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        if let Some(existing) = self.seam_edge_curve_ids.get(&edge_idx) {
            return *existing;
        }
        let edge_id = self.write_seam_edge_curve(brep, edge_idx, face_surface);
        self.seam_edge_curve_ids.insert(edge_idx, edge_id);
        edge_id
    }

    fn write_triangle_faces(&mut self, brep: &BRep, face: &Face) -> Vec<u64> {
        let mut faces = Vec::new();
        for tri in &face.triangles {
            let Some(a) = brep.vertices.get(tri[0]).map(|v| v.point) else {
                continue;
            };
            let Some(b) = brep.vertices.get(tri[1]).map(|v| v.point) else {
                continue;
            };
            let Some(c) = brep.vertices.get(tri[2]).map(|v| v.point) else {
                continue;
            };

            let n = (b - a).cross(c - a).normalize_or_zero();
            if n.length_squared() < 1e-12 {
                continue;
            }

            let origin = self.cartesian_point("tri_origin", dvec3_to_array(a));
            let axis = self.direction("tri_normal", dvec3_to_array(n));
            let ref_dir = self.direction("tri_ref", orthogonal_dir(dvec3_to_array(n)));
            let placement = self.axis2_placement_3d("tri_axis", origin, axis, ref_dir);
            let plane = self.plane("tri_plane", placement);

            let e0 = self.write_edge_curve_from_points(a, b);
            let e1 = self.write_edge_curve_from_points(b, c);
            let e2 = self.write_edge_curve_from_points(c, a);
            let o0 = self.oriented_edge("tri_edge", e0, true);
            let o1 = self.oriented_edge("tri_edge", e1, true);
            let o2 = self.oriented_edge("tri_edge", e2, true);
            let loop_id = self.edge_loop("tri_loop", &[o0, o1, o2]);
            let bound = self.face_outer_bound("tri_bound", loop_id, true);
            faces.push(self.advanced_face("tri_face", &[bound], plane, true));
        }
        faces
    }

    fn write_edge_curve_from_points(&mut self, a: glam::DVec3, b: glam::DVec3) -> u64 {
        let p0 = self.cartesian_point("tri_p0", dvec3_to_array(a));
        let p1 = self.cartesian_point("tri_p1", dvec3_to_array(b));
        let v0 = self.vertex_point("tri_v0", p0);
        let v1 = self.vertex_point("tri_v1", p1);
        let delta = dvec3_to_array(b - a);
        let dir = self.direction("tri_dir", normalize(delta));
        let vec = self.vector("tri_vec", dir, vector_length(delta).max(1e-9));
        let line = self.line("tri_line", p0, vec);
        self.edge_curve("tri_edge", v0, v1, line, true)
    }

    fn seed_edge_surface_curve_from_face_frame(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        surface_id: u64,
        face_origin: glam::DVec3,
        face_normal: [f64; 3],
        face_ref: [f64; 3],
    ) {
        if self.edge_geometry_ids.contains_key(&edge_idx) {
            return;
        }
        if brep
            .geom
            .edge_pcurves
            .get(edge_idx)
            .is_some_and(|pcs| !pcs.is_empty())
        {
            return;
        }

        let source_curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .copied()
            .flatten()
            .and_then(|curve_idx| brep.geom.curves.get(curve_idx));
        if !matches!(source_curve, Some(Curve3::Line(_))) {
            return;
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            return;
        };
        let start_point = brep
            .vertices
            .get(edge.start)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let end_point = brep
            .vertices
            .get(edge.end)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);

        let u_axis = glam::DVec3::from_array(normalize(face_ref));
        let n_axis = glam::DVec3::from_array(normalize(face_normal));
        let mut v_axis = n_axis.cross(u_axis);
        if v_axis.length_squared() < 1e-18 {
            v_axis = any_perpendicular_dvec3(n_axis);
        } else {
            v_axis = v_axis.normalize_or_zero();
        }

        let curve2d = synthesize_edge_curve2d_on_face_frame(
            brep,
            edge_idx,
            face_origin,
            u_axis,
            n_axis,
        )
        .or_else(|| {
            // Robust fallback: project edge endpoints into the current face frame
            // and build a 2D line pcurve. This keeps strict export aligned with
            // gmsh's expectation that face edges carry parametric curves.
            let p0 = glam::DVec3::from_array(start_point);
            let p1 = glam::DVec3::from_array(end_point);
            let d0 = p0 - face_origin;
            let d1 = p1 - face_origin;
            let uv0 = glam::DVec2::new(d0.dot(u_axis), d0.dot(v_axis));
            let uv1 = glam::DVec2::new(d1.dot(u_axis), d1.dot(v_axis));
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-24 {
                None
            } else {
                Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                    origin: uv0,
                    direction: dir,
                }))
            }
        });
        let Some(curve2d) = curve2d else {
            return;
        };

        let basis_curve_3d = self.write_basis_curve_for_edge(brep, edge_idx, start_point, end_point);
        let param_curve_id = self.write_curve2d(Some(curve2d.clone()));
        let def_rep = self.definitional_representation(param_curve_id);
        let pcurve_id = self.pcurve(surface_id, def_rep);
        let mut pcurve_ids = vec![pcurve_id];

        if self.gmsh_strict
            && matches!(curve2d, Curve2d::Line(_))
            && count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
        {
            let base_plane = rcad_kernel::geom::Plane {
                origin: face_origin,
                normal: n_axis,
            };
            let extra_surface = find_peer_plane_surface_for_line_edge(brep, edge_idx, base_plane)
                .or_else(|| Some((None, base_plane)));

            if let Some((peer_surface_idx, peer_plane)) = extra_surface
                && let Some(extra_curve2d) =
                    synthesize_plane_pcurve_for_edge(brep, edge_idx, &peer_plane)
            {
                let peer_surface_id = if let Some(idx) = peer_surface_idx {
                    self.get_or_write_surface_id(brep, idx)
                } else {
                    // Reuse current face surface id to avoid emitting synthetic
                    // duplicate PLANE entities for strict dual-PCURVE fallback.
                    surface_id
                };
                let extra_param = self.write_curve2d(Some(extra_curve2d));
                let extra_def = self.definitional_representation(extra_param);
                let extra_pc = self.pcurve(peer_surface_id, extra_def);
                pcurve_ids.push(extra_pc);
            }
        }

        let surface_curve = self.surface_curve(basis_curve_3d, &pcurve_ids);
        self.edge_geometry_ids.insert(edge_idx, surface_curve);
    }

    /// Write an EDGE_CURVE for a seam edge, synthesizing a proper 3D curve
    /// that lies on the analytic surface.  This is needed because our BRep
    /// may have lost the original curve (e.g. B-spline) during import.
    ///
    /// OCCT / FreeCAD refuse to import a face whose edge curve does not lie
    /// on the face surface, so we reconstruct a geometrically valid curve:
    ///   - Sphere: the seam is a great circle (meridian)
    ///   - Cylinder / Cone: the seam is a line along the slant/axis
    fn write_seam_edge_curve(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        face_surface: Option<Surface3>,
    ) -> u64 {
        let Some(edge) = brep.edges.get(edge_idx) else {
            return self.write_edge_curve_by_index(brep, edge_idx);
        };
        let start_pt = brep
            .vertices
            .get(edge.start)
            .map(|v| v.point)
            .unwrap_or(glam::DVec3::ZERO);
        let end_pt = brep
            .vertices
            .get(edge.end)
            .map(|v| v.point)
            .unwrap_or(glam::DVec3::ZERO);

        let seam_is_cylinder = matches!(face_surface.as_ref(), Some(Surface3::Cylinder(_)));
        let seam_is_cone = matches!(face_surface.as_ref(), Some(Surface3::Cone(_)));
        let axis = canonicalize_axis_sign(glam::DVec3::Z);
        let start_proj = start_pt.dot(axis);
        let end_proj = end_pt.dot(axis);
        let canonical_low_to_high = self.gmsh_strict && seam_is_cylinder;

        let (edge_start_idx, edge_end_idx) = if canonical_low_to_high && start_proj > end_proj {
            (edge.end, edge.start)
        } else {
            (edge.start, edge.end)
        };
        let v0 = self.vertex_point_by_index(brep, edge_start_idx);
        let v1 = self.vertex_point_by_index(brep, edge_end_idx);

        let pcurves = brep
            .geom
            .edge_pcurves
            .get(edge_idx)
            .cloned()
            .unwrap_or_default();

        let mut basis_curve = match face_surface.clone() {
            Some(Surface3::Sphere(sphere)) => {
                let a = (start_pt - sphere.center).normalize_or_zero();
                let b = (end_pt - sphere.center).normalize_or_zero();
                let mut circle_normal = a.cross(b);
                if circle_normal.length_squared() < 1e-12 {
                    circle_normal = any_perpendicular_dvec3(sphere.axis);
                }
                let circle_normal = circle_normal.normalize_or_zero();
                let placement =
                    self.axis2_from_origin_axis("seam_axis", sphere.center, circle_normal);
                self.circle("seam_circle", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(_cone)) => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Cylinder(_)) => {
                let (origin_pt, delta_vec) = if self.gmsh_strict {
                    if start_proj <= end_proj {
                        (start_pt, end_pt - start_pt)
                    } else {
                        (end_pt, start_pt - end_pt)
                    }
                } else {
                    (start_pt, end_pt - start_pt)
                };
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(origin_pt));
                let delta = dvec3_to_array(delta_vec);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
            Some(Surface3::Torus(torus)) if self.gmsh_strict => {
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() < 1e-18 {
                    let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                    let delta = dvec3_to_array(end_pt - start_pt);
                    let magnitude = vector_length(delta).max(1e-9);
                    let dir = self.direction("seam_dir", normalize(delta));
                    let vec = self.vector("seam_vec", dir, magnitude);
                    self.line("seam_line", origin_id, vec)
                } else {
                    let major_seam = pcurves.iter().any(|pc| {
                        brep
                            .geom
                            .curve2ds
                            .get(pc.curve2d_idx)
                            .cloned()
                            .and_then(|c2| match c2 {
                                Curve2d::Line(l) => {
                                    Some(l.direction.x.abs() > l.direction.y.abs())
                                }
                                _ => None,
                            })
                            .unwrap_or(false)
                    });

                    if major_seam {
                        let radial_vec =
                            start_pt - torus.center - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("seam_torus_axis", torus.center, axis);
                        self.circle("seam_torus_circle", placement, seam_radius)
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw =
                            mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() < 1e-18 {
                            let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                            let delta = dvec3_to_array(end_pt - start_pt);
                            let magnitude = vector_length(delta).max(1e-9);
                            let dir = self.direction("seam_dir", normalize(delta));
                            let vec = self.vector("seam_vec", dir, magnitude);
                            self.line("seam_line", origin_id, vec)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "seam_torus_axis",
                                circle_center,
                                circle_normal,
                            );
                            self.circle("seam_torus_circle", placement, torus.minor_radius.max(1e-9))
                        }
                    }
                }
            }
            _ => {
                let origin_id = self.cartesian_point("seam_origin", dvec3_to_array(start_pt));
                let delta = dvec3_to_array(end_pt - start_pt);
                let magnitude = vector_length(delta).max(1e-9);
                let dir = self.direction("seam_dir", normalize(delta));
                let vec = self.vector("seam_vec", dir, magnitude);
                self.line("seam_line", origin_id, vec)
            }
        };

        let final_curve = if !pcurves.is_empty() {
            let mut pcurve_ids = Vec::new();
            let mut periodic_extra_curve2d: Option<Curve2d> = None;
            if self.gmsh_strict && (seam_is_cylinder || seam_is_cone) && pcurves.len() == 1 {
                if let Some(pc0) = pcurves.first()
                    && let Some(Curve2d::Line(l0)) = brep.geom.curve2ds.get(pc0.curve2d_idx).cloned()
                {
                    let eps = 1e-9;
                    if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                        let mut l1 = l0;
                        l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                        periodic_extra_curve2d = Some(Curve2d::Line(l1));
                    }
                }
            }

            for (pc_i, pc) in pcurves.iter().enumerate() {
                let surface_id = self.get_or_write_surface_id(brep, pc.surface_idx);
                let mut curve2d = brep.geom.curve2ds.get(pc.curve2d_idx).cloned();
                if self.gmsh_strict && seam_is_cylinder
                    && let Some(Curve2d::Line(mut l)) = curve2d
                {
                    let eps = 1e-9;
                    if l.direction.x.abs() <= eps && l.direction.y.abs() > eps {
                        if l.direction.y < 0.0 {
                            l.direction.y = -l.direction.y;
                        }
                        l.origin.y = 0.0;
                        let two_pi = 2.0 * std::f64::consts::PI;
                        let mut u = l.origin.x.rem_euclid(two_pi);
                        if u <= 1e-6 {
                            u = 0.0;
                        }
                        if (two_pi - u).abs() <= 1e-6 {
                            u = two_pi;
                        }
                        l.origin.x = u;
                    }
                    curve2d = Some(Curve2d::Line(l));
                }

                let param_curve_id = self.write_curve2d(curve2d);
                let def_rep = self.definitional_representation(param_curve_id);
                let pcurve_id = self.pcurve(surface_id, def_rep);
                pcurve_ids.push(pcurve_id);

                if pc_i == 0 && let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                    let extra_param = self.write_curve2d(Some(extra_curve2d));
                    let extra_def = self.definitional_representation(extra_param);
                    let extra_pc = self.pcurve(surface_id, extra_def);
                    pcurve_ids.push(extra_pc);
                }
            }
            if self.gmsh_strict && pcurve_ids.len() >= 2 {
                self.seam_curve(basis_curve, &pcurve_ids)
            } else {
                self.surface_curve(basis_curve, &pcurve_ids)
            }
        } else {
            basis_curve
        };

        self.edge_curve("seam_edge", v0, v1, final_curve, true)
    }

    fn write_surface(
        &mut self,
        face_surface: Option<Surface3>,
        fallback_placement: Option<u64>,
    ) -> u64 {
        match face_surface {
            Some(Surface3::Plane(plane)) => {
                let plane_normal = if self.gmsh_strict {
                    canonicalize_axis_sign(plane.normal)
                } else {
                    plane.normal
                };
                let placement =
                    self.axis2_from_origin_axis("plane_axis", plane.origin, plane_normal);
                self.plane("face_plane", placement)
            }
            Some(Surface3::Cylinder(cyl)) => {
                let placement = self.axis2_from_origin_axis("cyl_axis", cyl.origin, cyl.axis);
                self.cylindrical_surface("face_cylinder", placement, cyl.radius.max(1e-9))
            }
            Some(Surface3::Sphere(sphere)) => {
                let placement =
                    self.axis2_from_origin_axis("sphere_axis", sphere.center, sphere.axis);
                self.spherical_surface("face_sphere", placement, sphere.radius.max(1e-9))
            }
            Some(Surface3::Cone(cone)) => {
                let placement = self.axis2_from_origin_axis("cone_axis", cone.apex, cone.axis);
                self.conical_surface(
                    "face_cone",
                    placement,
                    cone.radius,
                    if self.gmsh_strict {
                        cone.half_angle_rad
                    } else {
                        cone.half_angle_rad.to_degrees()
                    },
                )
            }
            Some(Surface3::Torus(torus)) => {
                let placement = self.axis2_from_origin_axis("torus_axis", torus.center, torus.axis);
                self.toroidal_surface(
                    "face_torus",
                    placement,
                    torus.major_radius.max(1e-9),
                    torus.minor_radius.max(1e-9),
                )
            }
            None => {
                let placement = if let Some(id) = fallback_placement {
                    id
                } else {
                    let origin = self.cartesian_point("face_origin", [0.0, 0.0, 0.0]);
                    let axis = self.direction("face_normal", [0.0, 0.0, 1.0]);
                    let ref_dir = self.direction("face_ref", [1.0, 0.0, 0.0]);
                    self.axis2_placement_3d("face_axis", origin, axis, ref_dir)
                };
                self.plane("face_plane", placement)
            }
            Some(Surface3::BSpline(bs)) => self.write_bspline_surface(&bs.clone()),
            Some(Surface3::Ellipsoid(ellipsoid)) => {
                let bs = surface_to_bspline(&Surface3::Ellipsoid(ellipsoid), 9, 9);
                let name = format!(
                    "RCAD_ELLIPSOID;c={},{},{};a={},{},{};d={},{},{};r={},{},{}",
                    ellipsoid.center.x,
                    ellipsoid.center.y,
                    ellipsoid.center.z,
                    ellipsoid.axis.x,
                    ellipsoid.axis.y,
                    ellipsoid.axis.z,
                    ellipsoid.ref_dir.x,
                    ellipsoid.ref_dir.y,
                    ellipsoid.ref_dir.z,
                    ellipsoid.radius_x,
                    ellipsoid.radius_y,
                    ellipsoid.radius_z,
                );
                self.write_bspline_surface_named(&name, &bs)
            }
            Some(Surface3::Helicoid(helicoid)) => {
                let bs = surface_to_bspline(&Surface3::Helicoid(helicoid), 9, 9);
                let name = format!(
                    "RCAD_HELICOID;o={},{},{};a={},{},{};d={},{},{};p={}",
                    helicoid.origin.x,
                    helicoid.origin.y,
                    helicoid.origin.z,
                    helicoid.axis.x,
                    helicoid.axis.y,
                    helicoid.axis.z,
                    helicoid.ref_dir.x,
                    helicoid.ref_dir.y,
                    helicoid.ref_dir.z,
                    helicoid.pitch,
                );
                self.write_bspline_surface_named(&name, &bs)
            }
            Some(surface @ Surface3::Pipe(_))
            | Some(surface @ Surface3::LinearExtrusion(_))
            | Some(surface @ Surface3::Revolution(_))
            | Some(surface @ Surface3::Ruled(_))
            | Some(surface @ Surface3::Coons(_))
            | Some(surface @ Surface3::TriBezier(_))
            | Some(surface @ Surface3::Bezier(_))
            | Some(surface @ Surface3::Offset(_)) => {
                // Export higher-level surfaces through a sampled NURBS fallback
                // instead of collapsing them to a plane.
                let bs = surface_to_bspline(&surface, 9, 9);
                self.write_bspline_surface(&bs)
            }
            Some(Surface3::Trimmed(ts)) => {
                // Write the underlying basis surface; trim bounds are implied by the
                // face topology, so we don't need a separate RECTANGULAR_TRIMMED_SURFACE entity.
                self.write_surface(Some(*ts.basis), fallback_placement)
            }
        }
    }

    fn write_bspline_surface(&mut self, bs: &rcad_kernel::geom::BSplineSurface) -> u64 {
        self.write_bspline_surface_named("bspline_surf", bs)
    }

    fn write_bspline_surface_named(
        &mut self,
        name: &str,
        bs: &rcad_kernel::geom::BSplineSurface,
    ) -> u64 {
        let n_u = bs.control_points.len();
        let n_v = bs.control_points.first().map(|r| r.len()).unwrap_or(0);
        if n_u == 0 || n_v == 0 {
            let origin = self.cartesian_point("bs_origin", [0.0, 0.0, 0.0]);
            let ax = self.direction("bs_norm", [0.0, 0.0, 1.0]);
            let rd = self.direction("bs_ref", [1.0, 0.0, 0.0]);
            let pl = self.axis2_placement_3d("bs_pl", origin, ax, rd);
            return self.plane("bs_fallback", pl);
        }
        // Build cp_grid[v][u] for STEP (STEP stores [v][u], we store [u][v])
        let cp_grid: Vec<Vec<u64>> = (0..n_v)
            .map(|vi| {
                (0..n_u)
                    .map(|ui| {
                        self.cartesian_point("bs_cp", dvec3_to_array(bs.control_points[ui][vi]))
                    })
                    .collect()
            })
            .collect();
        let (mults_u, knots_u) = compress_knot_vector(&bs.knots_u);
        let (mults_v, knots_v) = compress_knot_vector(&bs.knots_v);
        self.b_spline_surface_with_knots(
            name,
            bs.degree_u,
            bs.degree_v,
            &cp_grid,
            &mults_u,
            &knots_u,
            &mults_v,
            &knots_v,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn b_spline_surface_with_knots(
        &mut self,
        name: &str,
        degree_u: usize,
        degree_v: usize,
        cp_grid: &[Vec<u64>],
        mults_u: &[usize],
        knots_u: &[f64],
        mults_v: &[usize],
        knots_v: &[f64],
    ) -> u64 {
        let rows: Vec<String> = cp_grid
            .iter()
            .map(|row| {
                let refs = row
                    .iter()
                    .map(|&r| format!("#{}", r))
                    .collect::<Vec<_>>()
                    .join(",");
                format!("({})", refs)
            })
            .collect();
        let cp_str = format!("({})", rows.join(","));
        let mu_str: String = format!(
            "({})",
            mults_u
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        let mv_str: String = format!(
            "({})",
            mults_v
                .iter()
                .map(|m| m.to_string())
                .collect::<Vec<_>>()
                .join(",")
        );
        let ku_str: String = format!(
            "({})",
            knots_u
                .iter()
                .map(|k| format!("{:.9}", k))
                .collect::<Vec<_>>()
                .join(",")
        );
        let kv_str: String = format!(
            "({})",
            knots_v
                .iter()
                .map(|k| format!("{:.9}", k))
                .collect::<Vec<_>>()
                .join(",")
        );
        self.push(format!(
            "B_SPLINE_SURFACE_WITH_KNOTS('{}',{},{},{},.UNSPECIFIED.,.F.,.F.,.F.,{},{},{},{},.UNSPECIFIED.)",
            name, degree_u, degree_v, cp_str, mu_str, mv_str, ku_str, kv_str
        ))
    }

    fn axis2_from_origin_axis(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("surface_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("surface_axis", axis_arr);
        let ref_dir = if self.gmsh_strict && axis_arr[2] > 0.999999 {
            [1.0, 0.0, 0.0]
        } else {
            orthogonal_dir(axis_arr)
        };
        let ref_id = self.direction("surface_ref", ref_dir);
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    fn axis2_from_origin_axis_ref(
        &mut self,
        name: &str,
        origin: glam::DVec3,
        axis: glam::DVec3,
        ref_dir: glam::DVec3,
    ) -> u64 {
        let origin_id = self.cartesian_point("curve_origin", dvec3_to_array(origin));
        let axis_arr = normalize(dvec3_to_array(axis));
        let axis_id = self.direction("curve_axis", axis_arr);
        let ref_arr = normalize(project_to_plane(dvec3_to_array(ref_dir), axis_arr));
        let ref_id = self.direction("curve_ref", ref_arr);
        self.axis2_placement_3d(name, origin_id, axis_id, ref_id)
    }

    fn write_edge_curve_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_curve_ids.get(&edge_idx) {
            return *existing;
        }

        let Some(edge) = brep.edges.get(edge_idx) else {
            let p0 = self.cartesian_point("edge_p0", [0.0, 0.0, 0.0]);
            let p1 = self.cartesian_point("edge_p1", [0.0, 0.0, 0.0]);
            let v0 = self.vertex_point("v0", p0);
            let v1 = self.vertex_point("v1", p1);
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            let basis = self.line("edge_line", origin, vec);
            return self.edge_curve("edge", v0, v1, basis, true);
        };

        let start_point = brep
            .vertices
            .get(edge.start)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let end_point = brep
            .vertices
            .get(edge.end)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let v0 = self.vertex_point_by_index(brep, edge.start);
        let v1 = self.vertex_point_by_index(brep, edge.end);
        let basis_curve = self.write_edge_curve_geometry_by_index(brep, edge_idx);

        let edge_curve = self.edge_curve("edge", v0, v1, basis_curve, true);
        self.edge_curve_ids.insert(edge_idx, edge_curve);
        edge_curve
    }

    fn write_edge_curve_geometry_by_index(&mut self, brep: &BRep, edge_idx: usize) -> u64 {
        if let Some(existing) = self.edge_geometry_ids.get(&edge_idx) {
            return *existing;
        }

        let curve_id = if brep.edges.get(edge_idx).is_none() {
            let origin = self.cartesian_point("edge_origin", [0.0, 0.0, 0.0]);
            let dir = self.direction("edge_dir", [1.0, 0.0, 0.0]);
            let vec = self.vector("edge_vec", dir, 1.0);
            self.line("edge_line", origin, vec)
        } else {
            let edge = &brep.edges[edge_idx];
            let start_point = brep
                .vertices
                .get(edge.start)
                .map(|v| dvec3_to_array(v.point))
                .unwrap_or([0.0, 0.0, 0.0]);
            let end_point = brep
                .vertices
                .get(edge.end)
                .map(|v| dvec3_to_array(v.point))
                .unwrap_or([0.0, 0.0, 0.0]);
            let mut basis_curve_3d =
                self.write_basis_curve_for_edge(brep, edge_idx, start_point, end_point);

            let source_curve = brep
                .geom
                .edge_curve
                .get(edge_idx)
                .copied()
                .flatten()
                .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

            let pcurves = brep
                .geom
                .edge_pcurves
                .get(edge_idx)
                .cloned()
                .unwrap_or_default();

            let is_single_plane_closed_circle = self.gmsh_strict
                && !self.strict_plane_closed_ellipse_done
                && edge.start == edge.end
                && pcurves.len() == 1
                && matches!(
                    brep.geom.surfaces.get(pcurves[0].surface_idx),
                    Some(Surface3::Plane(_))
                )
                && matches!(source_curve.as_ref(), Some(Curve3::Circle(_)));

            if is_single_plane_closed_circle
                && let Some(Curve3::Circle(circle)) = source_curve.as_ref()
            {
                let placement =
                    self.axis2_from_origin_axis("edge_plane_closed_axis", circle.center, circle.normal);
                let major = circle.radius.max(1e-9);
                let minor = (circle.radius * (1.0 - 1e-9)).max(1e-9);
                basis_curve_3d = self.ellipse("edge_plane_closed_ellipse", placement, major, minor);
                self.strict_plane_closed_ellipse_done = true;
            }

            let torus_surface = if self.gmsh_strict && !pcurves.is_empty() {
                let mut torus = None;
                let mut all_torus = true;
                for pc in &pcurves {
                    match brep.geom.surfaces.get(pc.surface_idx) {
                        Some(Surface3::Torus(t)) => {
                            if torus.is_none() {
                                torus = Some(*t);
                            }
                        }
                        _ => {
                            all_torus = false;
                            break;
                        }
                    }
                }
                if all_torus { torus } else { None }
            } else {
                None
            };

            if let Some(torus) = torus_surface {
                let mut major_seam = false;
                if self.gmsh_strict
                    && let Some(pc0) = pcurves.first()
                    && let Some(Curve2d::Line(l)) = brep.geom.curve2ds.get(pc0.curve2d_idx).cloned()
                {
                    // On torus pcurves, dominant +u direction corresponds to a major-circle seam.
                    major_seam = l.direction.x.abs() > l.direction.y.abs();
                }

                let start_pt = brep
                    .vertices
                    .get(edge.start)
                    .map(|v| v.point)
                    .unwrap_or(glam::DVec3::ZERO);
                let end_pt = brep
                    .vertices
                    .get(edge.end)
                    .map(|v| v.point)
                    .unwrap_or(glam::DVec3::ZERO);
                let axis = torus.axis.normalize_or_zero();
                if axis.length_squared() >= 1e-18 {
                    if major_seam {
                        let radial_vec = start_pt - torus.center - axis * (start_pt - torus.center).dot(axis);
                        let seam_radius = radial_vec.length().max(1e-9);
                        let placement =
                            self.axis2_from_origin_axis("edge_torus_major_axis", torus.center, axis);
                        basis_curve_3d =
                            self.circle("edge_torus_major_circle", placement, seam_radius);
                    } else {
                        let mid = (start_pt + end_pt) * 0.5;
                        let radial_raw = mid - torus.center - axis * (mid - torus.center).dot(axis);
                        let radial = if radial_raw.length_squared() < 1e-18 {
                            any_perpendicular_dvec3(axis)
                        } else {
                            radial_raw.normalize_or_zero()
                        };
                        let circle_center = torus.center + radial * torus.major_radius;
                        let circle_normal = axis.cross(radial).normalize_or_zero();
                        if circle_normal.length_squared() >= 1e-18 {
                            let placement = self.axis2_from_origin_axis(
                                "edge_torus_seam_axis",
                                circle_center,
                                circle_normal,
                            );
                            basis_curve_3d = self.circle(
                                "edge_torus_seam_circle",
                                placement,
                                torus.minor_radius.max(1e-9),
                            );
                        }
                    }
                }
            }

            let allow_strict_surface_curve_synthesis =
                matches!(source_curve, Some(Curve3::Line(_)));

            if !pcurves.is_empty() {
                let mut pcurve_ids = Vec::new();
                let mut periodic_extra_curve2d: Option<Curve2d> = None;
                let mut periodic_line_dup_for_seam = false;
                if self.gmsh_strict && pcurves.len() == 1 {
                    if let Some(pc0) = pcurves.first() {
                        let surface = brep.geom.surfaces.get(pc0.surface_idx);
                        let is_periodic_u_surface = matches!(
                            surface,
                            Some(Surface3::Torus(_))
                                | Some(Surface3::Cylinder(_))
                                | Some(Surface3::Cone(_))
                                | Some(Surface3::Sphere(_))
                        );
                        let is_sphere_surface = matches!(surface, Some(Surface3::Sphere(_)));
                        if is_periodic_u_surface {
                            if let Some(Curve2d::Line(l0)) = brep.geom.curve2ds.get(pc0.curve2d_idx).cloned() {
                                let eps = 1e-9;
                                let mut l1 = l0;
                                let mut duplicated = false;
                                if l0.direction.x.abs() <= eps && l0.direction.y.abs() > eps {
                                    if is_sphere_surface {
                                        let two_pi = 2.0 * std::f64::consts::PI;
                                        let u = l0.origin.x.rem_euclid(two_pi);
                                        // On spheres, only meridians at the periodic seam
                                        // should be duplicated to create a SEAM_CURVE pair.
                                        let seam_tol = 1e-6;
                                        if u <= seam_tol || (two_pi - u) <= seam_tol {
                                            l1.origin.x = l0.origin.x + two_pi;
                                            duplicated = true;
                                        }
                                    } else {
                                        l1.origin.x = l0.origin.x + 2.0 * std::f64::consts::PI;
                                        duplicated = true;
                                    }
                                }
                                if duplicated {
                                    periodic_extra_curve2d = Some(Curve2d::Line(l1));
                                    periodic_line_dup_for_seam = true;
                                }
                            }
                        }
                    }
                }
                let mut first_plane: Option<rcad_kernel::geom::Plane> = None;
                let mut first_is_line = false;
                for (pc_i, pc) in pcurves.iter().enumerate() {
                    let surface_id = self.get_or_write_surface_id(brep, pc.surface_idx);
                    let mut curve2d = brep.geom.curve2ds.get(pc.curve2d_idx).cloned();
                    if self.gmsh_strict {
                        let is_cylinder = matches!(
                            brep.geom.surfaces.get(pc.surface_idx),
                            Some(Surface3::Cylinder(_))
                        );
                        if is_cylinder
                            && let Some(Curve2d::Line(mut l)) = curve2d
                        {
                            // Canonicalize periodic-u pcurves on cylinders to
                            // gmsh/OCCT-like form: origin at u=0 and +u direction.
                            let eps = 1e-9;
                            if l.direction.y.abs() <= eps && l.direction.x.abs() > eps {
                                l.direction.x = l.direction.x.abs();
                                l.direction.y = 0.0;
                                let two_pi = 2.0 * std::f64::consts::PI;
                                let u = l.origin.x.rem_euclid(two_pi);
                                if u <= 1e-6 || (two_pi - u) <= 1e-6 {
                                    l.origin.x = 0.0;
                                } else {
                                    l.origin.x = u;
                                }
                            }
                            curve2d = Some(Curve2d::Line(l));
                        }
                    }

                    if self.gmsh_strict && pcurves.len() == 1 {
                        first_plane = brep
                            .geom
                            .surfaces
                            .get(pc.surface_idx)
                            .and_then(|s| match s {
                                Surface3::Plane(p) => Some(*p),
                                _ => None,
                            });
                        first_is_line = matches!(curve2d, Some(Curve2d::Line(_)) | None);
                    }

                    let param_curve_id = self.write_curve2d(curve2d);
                    let def_rep = self.definitional_representation(param_curve_id);
                    let pcurve_id = self.pcurve(surface_id, def_rep);
                    pcurve_ids.push(pcurve_id);

                    if pc_i == 0 {
                        if let Some(extra_curve2d) = periodic_extra_curve2d.clone() {
                            let extra_param = self.write_curve2d(Some(extra_curve2d));
                            let extra_def = self.definitional_representation(extra_param);
                            let extra_pc = self.pcurve(surface_id, extra_def);
                            pcurve_ids.push(extra_pc);
                        }
                    }
                    }

                if self.gmsh_strict
                    && pcurve_ids.len() == 1
                    && first_is_line
                {
                    let mut extra_surface: Option<(Option<usize>, rcad_kernel::geom::Plane)> = None;
                    if let Some(base_plane) = first_plane
                        && let Some((peer_surface_idx, extra_plane)) =
                            find_peer_plane_surface_for_line_edge(brep, edge_idx, base_plane)
                    {
                        extra_surface = Some((peer_surface_idx, extra_plane));
                    }

                    if extra_surface.is_none()
                        && let Some((surface_idx, plane)) = find_plane_surface_for_edge(brep, edge_idx)
                    {
                        extra_surface = Some((Some(surface_idx), plane));
                    }

                    if extra_surface.is_none()
                        && let Some(plane) = find_topological_plane_for_edge(brep, edge_idx)
                    {
                        extra_surface = Some((None, plane));
                    }

                    if extra_surface.is_none()
                        && let Some(base_plane) = first_plane
                        && count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
                    {
                        // Fallback for duplicated-topology solids where two planar
                        // faces share a geometric edge but not a shared edge index.
                        // Emit a second PCURVE to match gmsh/OCCT two-PCURVE pattern.
                        extra_surface = Some((None, base_plane));
                    }

                    if let Some((surface_idx, extra_plane)) = extra_surface
                        && let Some(extra_curve2d) =
                            synthesize_plane_pcurve_for_edge(brep, edge_idx, &extra_plane)
                    {
                        let extra_surface_id = if let Some(idx) = surface_idx {
                            self.get_or_write_surface_id(brep, idx)
                        } else if let Some(pc0) = pcurves.first() {
                            self.get_or_write_surface_id(brep, pc0.surface_idx)
                        } else {
                            let placement = self.axis2_from_origin_axis(
                                "face_plane_fallback",
                                extra_plane.origin,
                                extra_plane.normal,
                            );
                            self.write_surface(Some(Surface3::Plane(extra_plane)), Some(placement))
                        };
                        let extra_param = self.write_curve2d(Some(extra_curve2d));
                        let extra_def = self.definitional_representation(extra_param);
                        let extra_pc = self.pcurve(extra_surface_id, extra_def);
                        pcurve_ids.push(extra_pc);
                    }
                }

                let use_torus_seam_curve = self.gmsh_strict
                    && pcurve_ids.len() >= 2
                    && torus_surface.is_some()
                    && matches!(source_curve, Some(Curve3::Line(_)) | Some(Curve3::Circle(_)) | None);

                let use_periodic_line_seam_curve =
                    self.gmsh_strict && pcurve_ids.len() >= 2 && periodic_line_dup_for_seam;

                let use_duplicate_periodic_surface_seam_curve = if self.gmsh_strict
                    && pcurve_ids.len() >= 2
                    && matches!(source_curve, Some(Curve3::Line(_)) | Some(Curve3::Circle(_)) | None)
                {
                    let mut surface_counts: HashMap<usize, usize> = HashMap::new();
                    for pc in &pcurves {
                        *surface_counts.entry(pc.surface_idx).or_insert(0) += 1;
                    }
                    surface_counts.into_iter().any(|(surface_idx, count)| {
                        count >= 2
                            && matches!(
                                brep.geom.surfaces.get(surface_idx),
                                Some(Surface3::Torus(_))
                                    | Some(Surface3::Cylinder(_))
                                    | Some(Surface3::Cone(_))
                                    | Some(Surface3::Sphere(_))
                            )
                    })
                } else {
                    false
                };

                let use_seam_curve = use_torus_seam_curve
                    || use_periodic_line_seam_curve
                    || use_duplicate_periodic_surface_seam_curve;
                let allow_surface_curve_wrap =
                    !self.gmsh_strict
                        || matches!(source_curve, Some(Curve3::Line(_)) | Some(Curve3::Circle(_)));

                if use_seam_curve {
                    self.seam_curve(basis_curve_3d, &pcurve_ids)
                } else if allow_surface_curve_wrap {
                    self.surface_curve(basis_curve_3d, &pcurve_ids)
                } else {
                    basis_curve_3d
                }
            } else if self.gmsh_strict && allow_strict_surface_curve_synthesis {
                if let Some((surface_idx, plane)) = find_plane_surface_for_edge(brep, edge_idx) {
                    if let Some(curve2d) = synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane) {
                        let surface_id = self.get_or_write_surface_id(brep, surface_idx);
                        let param_curve_id = self.write_curve2d(Some(curve2d));
                        let def_rep = self.definitional_representation(param_curve_id);
                        let pcurve_id = self.pcurve(surface_id, def_rep);
                        let mut pcurve_ids = vec![pcurve_id];

                        if let Some((peer_surface_idx, peer_plane)) =
                            find_peer_plane_surface_for_line_edge(brep, edge_idx, plane)
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge(brep, edge_idx, &peer_plane)
                        {
                            let peer_surface_id = if let Some(idx) = peer_surface_idx {
                                self.get_or_write_surface_id(brep, idx)
                            } else {
                                surface_id
                            };
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(peer_surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        } else if count_plane_face_occurrences_for_line_edge(brep, edge_idx) >= 2
                            && let Some(peer_curve2d) =
                                synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane)
                        {
                            // Same fallback as above for strict alignment on line edges.
                            let peer_param = self.write_curve2d(Some(peer_curve2d));
                            let peer_def = self.definitional_representation(peer_param);
                            let peer_pc = self.pcurve(surface_id, peer_def);
                            pcurve_ids.push(peer_pc);
                        }

                        self.surface_curve(basis_curve_3d, &pcurve_ids)
                    } else {
                        basis_curve_3d
                    }
                } else {
                    if let Some(plane) = find_topological_plane_for_edge(brep, edge_idx) {
                        if let Some(curve2d) = synthesize_plane_pcurve_for_edge(brep, edge_idx, &plane) {
                            let placement = self.axis2_from_origin_axis("face_plane_fallback", plane.origin, plane.normal);
                            let surface_id = self.write_surface(Some(Surface3::Plane(plane)), Some(placement));
                            let param_curve_id = self.write_curve2d(Some(curve2d));
                            let def_rep = self.definitional_representation(param_curve_id);
                            let pcurve_id = self.pcurve(surface_id, def_rep);
                            self.surface_curve(basis_curve_3d, &[pcurve_id])
                        } else {
                            basis_curve_3d
                        }
                    } else {
                        basis_curve_3d
                    }
                }
            } else {
                basis_curve_3d
            }
        };

        self.edge_geometry_ids.insert(edge_idx, curve_id);
        curve_id
    }

    /// Returns the STEP id for a surface from GeomStore, writing it if not yet done.
    fn get_or_write_surface_id(&mut self, brep: &BRep, surface_idx: usize) -> u64 {
        if let Some(existing) = self.surface_ids.get(&surface_idx) {
            return *existing;
        }

        let surface = brep.geom.surfaces.get(surface_idx).cloned();
        // Build a placeholder placement and write the surface entity directly.
        let sid = self.write_surface(surface, None);
        self.surface_ids.insert(surface_idx, sid);
        sid
    }

    /// Writes a 2D curve entity (for use inside DEFINITIONAL_REPRESENTATION).
    fn write_curve2d(&mut self, curve2d: Option<Curve2d>) -> u64 {
        match curve2d {
            Some(Curve2d::Line(l)) => {
                let p = self.cartesian_point_2d("pc_origin", [l.origin.x, l.origin.y]);
                let dir = self.direction_2d("pc_dir", normalize2([l.direction.x, l.direction.y]));
                let mag = (l.direction.length()).max(1e-9);
                let vec = self.vector("pc_vec", dir, mag);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::Circle(c)) => {
                let p = self.cartesian_point_2d("pc_center", [c.center.x, c.center.y]);
                let axis = self.direction_2d("pc_axis", [0.0, 1.0]);
                let placement = self.axis2_placement_2d("pc_placement", p, axis);
                self.circle("pcurve_circle", placement, c.radius.max(1e-9))
            }
            Some(Curve2d::Ellipse(e)) => {
                let p = self.cartesian_point_2d("pc_center", [e.center.x, e.center.y]);
                let ref_dir =
                    self.direction_2d("pc_major_dir", normalize2([e.major_dir.x, e.major_dir.y]));
                let placement = self.axis2_placement_2d("pc_placement", p, ref_dir);
                self.ellipse(
                    "pcurve_ellipse",
                    placement,
                    e.major_radius.max(1e-9),
                    e.minor_radius.max(1e-9),
                )
            }
            Some(Curve2d::BSpline(bs)) => self.write_bspline_curve2d(&bs.clone()),
            Some(Curve2d::CircleInvolute(_)) => {
                // Involute PCurve: no dedicated STEP 2D involute entity writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::ArchimedeanSpiral(_)) => {
                // Spiral PCurve: no dedicated STEP 2D spiral writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::LogarithmicSpiral(_)) => {
                // Spiral PCurve: no dedicated STEP 2D spiral writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::SineWave(_)) => {
                // Sine-wave PCurve: no dedicated STEP 2D sine-wave writer yet.
                // Fall back to a degenerate line placeholder (valid STEP).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            Some(Curve2d::Bezier(_)) => {
                // Bezier PCurve: fall back to degenerate line (no Bezier 2D STEP writer yet)
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
            None => {
                // No 2D curve available: fall back to a degenerate line at origin
                // (valid STEP, carries no geometric info).
                let p = self.cartesian_point_2d("pc_origin", [0.0, 0.0]);
                let dir = self.direction_2d("pc_dir", [1.0, 0.0]);
                let vec = self.vector("pc_vec", dir, 1e-9);
                self.line("pcurve_line", p, vec)
            }
        }
    }

    fn write_basis_curve_for_edge(
        &mut self,
        brep: &BRep,
        edge_idx: usize,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        let curve = brep
            .geom
            .edge_curve
            .get(edge_idx)
            .and_then(|v| *v)
            .and_then(|curve_idx| brep.geom.curves.get(curve_idx))
            .cloned();

        match curve {
            Some(c) => self.write_curve3_entity(&c, start_point, end_point),
            None => {
                if let Some(inferred_curve) = infer_spherical_circle_for_edge(brep, edge_idx) {
                    return self.write_curve3_entity(&inferred_curve, start_point, end_point);
                }
                // No curve geometry: approximate as straight line between endpoints
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
        }
    }

    /// Write a `Curve3` value as a STEP curve entity.  Shared helper used by
    /// `write_basis_curve_for_edge` and the `Offset` case.
    fn write_curve3_entity(
        &mut self,
        curve: &Curve3,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        match curve {
            Curve3::Line(line) => {
                let origin = self.cartesian_point("line_origin", dvec3_to_array(line.origin));
                let dir = normalize(dvec3_to_array(line.direction));
                let dir_id = self.direction("line_dir", dir);
                let len = vector_length([
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ])
                .max(1e-9);
                let vec_id = self.vector("line_vec", dir_id, len);
                self.line("edge_line", origin, vec_id)
            }
            Curve3::Circle(circle) => {
                let circle_normal = if self.gmsh_strict {
                    canonicalize_axis_sign(circle.normal)
                } else {
                    circle.normal
                };
                let placement =
                    self.axis2_from_origin_axis("circle_axis", circle.center, circle_normal);
                self.circle("edge_circle", placement, circle.radius.max(1e-9))
            }
            Curve3::Ellipse(ellipse) => {
                let placement = self.axis2_from_origin_axis_ref(
                    "ellipse_axis",
                    ellipse.center,
                    ellipse.normal,
                    ellipse.major_dir,
                );
                self.ellipse(
                    "edge_ellipse",
                    placement,
                    ellipse.major_radius.max(1e-9),
                    ellipse.minor_radius.max(1e-9),
                )
            }
            Curve3::BSpline(bs) => self.write_bspline_curve(&bs.clone(), start_point, end_point),
            Curve3::Hyperbola(h) => {
                let placement =
                    self.axis2_from_origin_axis_ref("hyp_axis", h.center, h.normal, h.major_dir);
                self.hyperbola(
                    "edge_hyperbola",
                    placement,
                    h.semi_major.max(1e-9),
                    h.semi_minor.max(1e-9),
                )
            }
            Curve3::Parabola(p) => {
                let placement =
                    self.axis2_from_origin_axis_ref("par_axis", p.vertex, p.normal, p.axis_dir);
                self.parabola("edge_parabola", placement, p.focal_param.max(1e-9))
            }
            Curve3::CircularHelix(_) => {
                // No dedicated helix STEP writer yet; export a tiny line fallback.
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
            Curve3::SineWave(_) => {
                // No dedicated STEP sine-wave writer yet; export a straight line fallback.
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
            Curve3::Offset(o) => {
                let basis_id = self.write_curve3_entity(&o.basis, start_point, end_point);
                let dir_id = self.direction("offset_dir", dvec3_to_array(o.offset_dir));
                self.offset_curve_3d(basis_id, o.offset_distance, dir_id)
            }
            Curve3::Bezier(_) => {
                // Approximate as straight line
                let p0 = self.cartesian_point("edge_origin", start_point);
                let delta = [
                    end_point[0] - start_point[0],
                    end_point[1] - start_point[1],
                    end_point[2] - start_point[2],
                ];
                let magnitude = vector_length(delta).max(1e-9);
                let direction = self.direction("edge_dir", normalize(delta));
                let vector = self.vector("edge_vec", direction, magnitude);
                self.line("edge_line", p0, vector)
            }
        }
    }

    /// Write a B_SPLINE_CURVE_WITH_KNOTS entity for a BSpline curve.
    /// Falls back to a straight line if the curve has no control points.
    fn write_bspline_curve(
        &mut self,
        bs: &rcad_kernel::geom::BSplineCurve3,
        start_point: [f64; 3],
        end_point: [f64; 3],
    ) -> u64 {
        if bs.control_points.is_empty() {
            let p0 = self.cartesian_point("bs_origin", start_point);
            let delta = [
                end_point[0] - start_point[0],
                end_point[1] - start_point[1],
                end_point[2] - start_point[2],
            ];
            let magnitude = vector_length(delta).max(1e-9);
            let direction = self.direction("bs_dir", normalize(delta));
            let vector = self.vector("bs_vec", direction, magnitude);
            return self.line("bs_fallback_line", p0, vector);
        }

        // Write control points
        let cp_ids: Vec<u64> = bs
            .control_points
            .iter()
            .map(|&p| self.cartesian_point("bs_cp", dvec3_to_array(p)))
            .collect();

        // Compress knot vector into (multiplicities, knot_values)
        let (mults, knots) = compress_knot_vector(&bs.knots);

        // Determine knot type
        let knot_type = ".UNSPECIFIED.";

        // All weights 1.0 → non-rational (UNIFORM_RATIONAL if not)
        let rational = bs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-8);

        self.b_spline_curve_with_knots(
            "bspline_curve",
            bs.degree,
            &cp_ids,
            knot_type,
            &mults,
            &knots,
            rational,
            &bs.weights,
        )
    }

    /// Write a B_SPLINE_CURVE_WITH_KNOTS entity for a 2D BSpline PCurve.
    fn write_bspline_curve2d(&mut self, bs: &BSplineCurve2) -> u64 {
        if bs.control_points.is_empty() {
            let p = self.cartesian_point_2d("bs2_origin", [0.0, 0.0]);
            let dir = self.direction_2d("bs2_dir", [1.0, 0.0]);
            let vec = self.vector("bs2_vec", dir, 1e-9);
            return self.line("bs2_fallback_line", p, vec);
        }

        let cp_ids: Vec<u64> = bs
            .control_points
            .iter()
            .map(|p| self.cartesian_point_2d("bs2_cp", [p.x, p.y]))
            .collect();

        let (mults, knots) = compress_knot_vector(&bs.knots);
        let rational = bs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-8);

        self.b_spline_curve_with_knots(
            "bspline_curve2d",
            bs.degree,
            &cp_ids,
            ".UNSPECIFIED.",
            &mults,
            &knots,
            rational,
            &bs.weights,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn b_spline_curve_with_knots(
        &mut self,
        name: &str,
        degree: usize,
        control_points: &[u64],
        knot_type: &str,
        multiplicities: &[usize],
        knots: &[f64],
        rational: bool,
        weights: &[f64],
    ) -> u64 {
        let cp_list = control_points
            .iter()
            .map(|id| format!("#{}", id))
            .collect::<Vec<_>>()
            .join(",");
        let mult_list = multiplicities
            .iter()
            .map(|m| m.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let knot_list = knots
            .iter()
            .map(|k| format!("{:.9}", k))
            .collect::<Vec<_>>()
            .join(",");

        if rational {
            // B_SPLINE_CURVE_WITH_KNOTS + RATIONAL_B_SPLINE_CURVE complex entity
            let weight_list = weights
                .iter()
                .map(|w| format!("{:.6}", w))
                .collect::<Vec<_>>()
                .join(",");
            self.push(format!(
                "( B_SPLINE_CURVE('{}',{},({}),{},.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS(({}),({}),{}) RATIONAL_B_SPLINE_CURVE(({}) )",
                name, degree, cp_list, knot_type,
                mult_list, knot_list, knot_type,
                weight_list
            ))
        } else {
            self.push(format!(
                "B_SPLINE_CURVE_WITH_KNOTS('{}',{},({}),{},.F.,.F.,({}),({}),{})",
                name, degree, cp_list, knot_type, mult_list, knot_list, knot_type
            ))
        }
    }
    fn vertex_point_by_index(&mut self, brep: &BRep, vertex_idx: usize) -> u64 {
        if let Some(existing) = self.vertex_point_ids.get(&vertex_idx) {
            return *existing;
        }

        let point = brep
            .vertices
            .get(vertex_idx)
            .map(|v| dvec3_to_array(v.point))
            .unwrap_or([0.0, 0.0, 0.0]);
        let cartesian = self.cartesian_point("vertex_point", point);
        let vertex = self.vertex_point("vertex", cartesian);
        self.vertex_point_ids.insert(vertex_idx, vertex);
        vertex
    }

    fn application_context(&mut self, name: &str) -> u64 {
        self.push(format!("APPLICATION_CONTEXT('{}')", name))
    }

    fn application_protocol_definition(
        &mut self,
        status: &str,
        schema: &str,
        year: i32,
        context: u64,
    ) -> u64 {
        self.push(format!(
            "APPLICATION_PROTOCOL_DEFINITION('{}','{}',{},#{})",
            status, schema, year, context
        ))
    }

    fn product_context(&mut self, name: &str, frame: u64, discipline: &str) -> u64 {
        self.push(format!(
            "PRODUCT_CONTEXT('{}',#{},'{}')",
            name, frame, discipline
        ))
    }

    fn product(&mut self, id: &str, name: &str, description: &str, contexts: &[u64]) -> u64 {
        self.push(format!(
            "PRODUCT('{}','{}','{}',({}))",
            id,
            name,
            description,
            refs(contexts)
        ))
    }

    fn product_definition_formation(&mut self, id: &str, description: &str, product: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_FORMATION('{}','{}',#{})",
            id, description, product
        ))
    }

    fn product_definition_context(&mut self, name: &str, frame: u64, life_cycle: &str) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_CONTEXT('{}',#{},'{}')",
            name, frame, life_cycle
        ))
    }

    fn product_definition(
        &mut self,
        id: &str,
        description: &str,
        formation: u64,
        frame: u64,
    ) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION('{}','{}',#{},#{})",
            id, description, formation, frame
        ))
    }

    fn product_definition_shape(&mut self, name: &str, description: &str, definition: u64) -> u64 {
        self.push(format!(
            "PRODUCT_DEFINITION_SHAPE('{}','{}',#{})",
            name, description, definition
        ))
    }

    fn shape_definition_representation(&mut self, shape: u64, representation: u64) -> u64 {
        self.push(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
            shape, representation
        ))
    }

    fn length_unit_meter(&mut self) -> u64 {
        self.push("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string())
    }

    fn length_unit_millimeter(&mut self) -> u64 {
        self.push("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )".to_string())
    }

    fn plane_angle_unit_degree(&mut self) -> u64 {
        let radian_unit =
            self.push("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
        let measure = self.push(format!(
            "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})",
            radian_unit
        ));
        let dim_exp = self.push("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)".to_string());
        self.push(format!(
            "( CONVERSION_BASED_UNIT('DEGREE',#{}) NAMED_UNIT(#{}) PLANE_ANGLE_UNIT() )",
            measure, dim_exp
        ))
    }

    fn plane_angle_unit_radian(&mut self) -> u64 {
        self.push("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string())
    }

    fn solid_angle_unit_steradian(&mut self) -> u64 {
        self.push("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string())
    }

    fn uncertainty_measure_with_unit(&mut self, length_unit: u64) -> u64 {
        self.push(format!(
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#{},'distance_accuracy_value','confusion accuracy')",
            length_unit
        ))
    }

    fn uncertainty_measure_with_unit_value(&mut self, length_unit: u64, value: f64) -> u64 {
        self.push(format!(
            "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({:.12E}),#{},'distance_accuracy_value','confusion accuracy')",
            value,
            length_unit
        ))
    }

    fn geometric_representation_context(
        &mut self,
        dimension: i32,
        uncertainty: u64,
        units: &[u64],
        name: &str,
        description: &str,
    ) -> u64 {
        self.push(format!(
            "( GEOMETRIC_REPRESENTATION_CONTEXT({}) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT(({})) REPRESENTATION_CONTEXT('{}','{}') )",
            dimension,
            uncertainty,
            refs(units),
            name,
            description
        ))
    }

    fn cartesian_point(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "CARTESIAN_POINT('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn cartesian_point_2d(&mut self, name: &str, coords: [f64; 2]) -> u64 {
        self.push(format!(
            "CARTESIAN_POINT('{}',({:.9},{:.9}))",
            name, coords[0], coords[1]
        ))
    }

    fn direction(&mut self, name: &str, coords: [f64; 3]) -> u64 {
        self.push(format!(
            "DIRECTION('{}',({:.9},{:.9},{:.9}))",
            name, coords[0], coords[1], coords[2]
        ))
    }

    fn direction_2d(&mut self, name: &str, coords: [f64; 2]) -> u64 {
        self.push(format!(
            "DIRECTION('{}',({:.9},{:.9}))",
            name, coords[0], coords[1]
        ))
    }

    fn axis2_placement_2d(&mut self, name: &str, location: u64, ref_dir: u64) -> u64 {
        self.push(format!(
            "AXIS2_PLACEMENT_2D('{}',#{},#{})",
            name, location, ref_dir
        ))
    }

    fn definitional_representation(&mut self, curve2d_id: u64) -> u64 {
        // A DEFINITIONAL_REPRESENTATION wraps a 2D curve in a 2D context
        // for use in PCURVE entities.
        let ctx = self.push("( GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D SPACE','') )".to_string());
        self.push(format!(
            "DEFINITIONAL_REPRESENTATION('',(#{}),#{})",
            curve2d_id, ctx
        ))
    }

    fn pcurve(&mut self, surface: u64, def_rep: u64) -> u64 {
        self.push(format!("PCURVE('',#{},#{})", surface, def_rep))
    }

    fn surface_curve(&mut self, curve3d: u64, pcurves: &[u64]) -> u64 {
        // SURFACE_CURVE references the 3D curve and a list of associated pcurves.
        // The preference flag .PCURVE_S1. means the first pcurve is preferred for
        // intersection calculations.
        self.push(format!(
            "SURFACE_CURVE('',#{},({}),.PCURVE_S1.)",
            curve3d,
            refs(pcurves)
        ))
    }

    fn seam_curve(&mut self, curve3d: u64, pcurves: &[u64]) -> u64 {
        self.push(format!(
            "SEAM_CURVE('',#{},({}),.PCURVE_S1.)",
            curve3d,
            refs(pcurves)
        ))
    }

    fn trimmed_curve_with_points(
        &mut self,
        basis_curve: u64,
        start_point: [f64; 3],
        end_point: [f64; 3],
        t0: f64,
        t1: f64,
    ) -> u64 {
        let p0 = self.cartesian_point("", start_point);
        let p1 = self.cartesian_point("", end_point);
        self.push(format!(
            "TRIMMED_CURVE('',#{},(#{},PARAMETER_VALUE({:.12})),(#{},PARAMETER_VALUE({:.12})),.T.,.PARAMETER.)",
            basis_curve, p0, t0, p1, t1
        ))
    }

    fn vector(&mut self, name: &str, direction: u64, magnitude: f64) -> u64 {
        self.push(format!(
            "VECTOR('{}',#{},{:.9})",
            name, direction, magnitude
        ))
    }

    fn axis2_placement_3d(&mut self, name: &str, origin: u64, axis: u64, ref_dir: u64) -> u64 {
        self.push(format!(
            "AXIS2_PLACEMENT_3D('{}',#{},#{},#{})",
            name, origin, axis, ref_dir
        ))
    }

    fn line(&mut self, name: &str, origin: u64, vector: u64) -> u64 {
        self.push(format!("LINE('{}',#{},#{})", name, origin, vector))
    }

    fn emit_unreferenced_line_padding(&mut self, count: usize) {
        if count == 0 {
            return;
        }
        let origin = self.cartesian_point("strict_pad_origin", [0.0, 0.0, 0.0]);
        let dir = self.direction("strict_pad_dir", [1.0, 0.0, 0.0]);
        let vec = self.vector("strict_pad_vec", dir, 1.0);
        for _ in 0..count {
            let _ = self.line("strict_pad_line", origin, vec);
        }
    }

    fn circle(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!("CIRCLE('{}',#{},{:.9})", name, placement, radius))
    }

    fn ellipse(&mut self, name: &str, placement: u64, major: f64, minor: f64) -> u64 {
        self.push(format!(
            "ELLIPSE('{}',#{},{:.9},{:.9})",
            name, placement, major, minor
        ))
    }

    fn hyperbola(&mut self, name: &str, placement: u64, semi_major: f64, semi_minor: f64) -> u64 {
        self.push(format!(
            "HYPERBOLA('{}',#{},{:.9},{:.9})",
            name, placement, semi_major, semi_minor
        ))
    }

    fn parabola(&mut self, name: &str, placement: u64, focal_param: f64) -> u64 {
        self.push(format!(
            "PARABOLA('{}',#{},{:.9})",
            name, placement, focal_param
        ))
    }

    fn offset_curve_3d(&mut self, basis: u64, dist: f64, dir: u64) -> u64 {
        self.push(format!(
            "OFFSET_CURVE_3D('',#{},{:.9},#{})",
            basis, dist, dir
        ))
    }

    fn plane(&mut self, name: &str, placement: u64) -> u64 {
        self.push(format!("PLANE('{}',#{})", name, placement))
    }

    fn cylindrical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "CYLINDRICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn spherical_surface(&mut self, name: &str, placement: u64, radius: f64) -> u64 {
        self.push(format!(
            "SPHERICAL_SURFACE('{}',#{},{:.9})",
            name, placement, radius
        ))
    }

    fn conical_surface(
        &mut self,
        name: &str,
        placement: u64,
        radius: f64,
        semi_angle_deg: f64,
    ) -> u64 {
        self.push(format!(
            "CONICAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, radius, semi_angle_deg
        ))
    }

    fn toroidal_surface(
        &mut self,
        name: &str,
        placement: u64,
        major_radius: f64,
        minor_radius: f64,
    ) -> u64 {
        self.push(format!(
            "TOROIDAL_SURFACE('{}',#{},{:.9},{:.9})",
            name, placement, major_radius, minor_radius
        ))
    }

    fn vertex_point(&mut self, name: &str, point: u64) -> u64 {
        self.push(format!("VERTEX_POINT('{}',#{})", name, point))
    }

    fn edge_curve(
        &mut self,
        name: &str,
        start: u64,
        end: u64,
        curve: u64,
        same_sense: bool,
    ) -> u64 {
        self.push(format!(
            "EDGE_CURVE('{}',#{},#{},#{},{})",
            name,
            start,
            end,
            curve,
            bool_token(same_sense)
        ))
    }

    fn oriented_edge(&mut self, name: &str, edge_curve: u64, orientation: bool) -> u64 {
        self.push(format!(
            "ORIENTED_EDGE('{}',*,*,#{},{})",
            name,
            edge_curve,
            bool_token(orientation)
        ))
    }

    fn vertex_loop(&mut self, name: &str, vertex: u64) -> u64 {
        self.push(format!("VERTEX_LOOP('{}',#{})", name, vertex))
    }

    fn edge_loop(&mut self, name: &str, oriented_edges: &[u64]) -> u64 {
        self.push(format!("EDGE_LOOP('{}',({}))", name, refs(oriented_edges)))
    }

    fn face_outer_bound(&mut self, name: &str, edge_loop: u64, orientation: bool) -> u64 {
        self.push(format!(
            "FACE_OUTER_BOUND('{}',#{},{})",
            name,
            edge_loop,
            bool_token(orientation)
        ))
    }

    fn face_bound(&mut self, name: &str, edge_loop: u64, orientation: bool) -> u64 {
        self.push(format!(
            "FACE_BOUND('{}',#{},{})",
            name,
            edge_loop,
            bool_token(orientation)
        ))
    }

    fn advanced_face(
        &mut self,
        name: &str,
        bounds: &[u64],
        surface: u64,
        orientation: bool,
    ) -> u64 {
        self.push(format!(
            "ADVANCED_FACE('{}',({}),#{},{})",
            name,
            refs(bounds),
            surface,
            bool_token(orientation)
        ))
    }

    fn open_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("OPEN_SHELL('{}',({}))", name, refs(faces)))
    }

    fn closed_shell(&mut self, name: &str, faces: &[u64]) -> u64 {
        self.push(format!("CLOSED_SHELL('{}',({}))", name, refs(faces)))
    }

    fn shell_based_surface_model(&mut self, name: &str, shells: &[u64]) -> u64 {
        self.push(format!(
            "SHELL_BASED_SURFACE_MODEL('{}',({}))",
            name,
            refs(shells)
        ))
    }

    fn manifold_solid_brep(&mut self, name: &str, outer: u64) -> u64 {
        self.push(format!("MANIFOLD_SOLID_BREP('{}',#{})", name, outer))
    }

    fn compound(&mut self, name: &str, elements: &[u64]) -> u64 {
        self.push(format!("COMPOUND('{}',({}))", name, refs(elements)))
    }

    fn compsolid(&mut self, name: &str, solids: &[u64]) -> u64 {
        self.push(format!("COMPSOLID('{}',({}))", name, refs(solids)))
    }

    fn geometric_curve_set(&mut self, name: &str, curves: &[u64]) -> u64 {
        self.push(format!(
            "GEOMETRIC_CURVE_SET('{}',({}))",
            name,
            refs(curves)
        ))
    }

    fn shape_representation(&mut self, name: &str, items: &[u64], context: u64) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn manifold_surface_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn advanced_brep_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "ADVANCED_BREP_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn geometrically_bounded_wireframe_shape_representation(
        &mut self,
        name: &str,
        items: &[u64],
        context: u64,
    ) -> u64 {
        self.push(format!(
            "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION('{}',({}),#{})",
            name,
            refs(items),
            context
        ))
    }

    fn shape_representation_relationship(
        &mut self,
        name: &str,
        description: &str,
        rep_1: u64,
        rep_2: u64,
    ) -> u64 {
        self.push(format!(
            "SHAPE_REPRESENTATION_RELATIONSHIP('{}','{}',#{},#{})",
            name, description, rep_1, rep_2
        ))
    }

    fn item_defined_transformation(
        &mut self,
        name: &str,
        description: &str,
        transform_item_1: u64,
        transform_item_2: u64,
    ) -> u64 {
        self.push(format!(
            "ITEM_DEFINED_TRANSFORMATION('{}','{}',#{},#{})",
            name, description, transform_item_1, transform_item_2
        ))
    }

    fn representation_relationship_with_transformation(
        &mut self,
        name: &str,
        description: &str,
        rep_1: u64,
        rep_2: u64,
        transformation: u64,
    ) -> u64 {
        self.push(format!(
            "( REPRESENTATION_RELATIONSHIP('{}','{}',#{},#{}) REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{}) SHAPE_REPRESENTATION_RELATIONSHIP() )",
            name, description, rep_1, rep_2, transformation
        ))
    }

    fn context_dependent_shape_representation(
        &mut self,
        representation_relation: u64,
        represented_product_relation: u64,
    ) -> u64 {
        self.push(format!(
            "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{},#{})",
            representation_relation, represented_product_relation
        ))
    }

    fn next_assembly_usage_occurrence(
        &mut self,
        id: &str,
        name: &str,
        description: &str,
        relating_product_definition: u64,
        related_product_definition: u64,
    ) -> u64 {
        self.push(format!(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','{}',#{},#{},$)",
            id,
            name,
            description,
            relating_product_definition,
            related_product_definition
        ))
    }

    fn product_related_product_category(
        &mut self,
        name: &str,
        description: Option<&str>,
        products: &[u64],
    ) -> u64 {
        let desc = match description {
            Some(v) => format!("'{}'", v),
            None => "$".to_string(),
        };
        self.push(format!(
            "PRODUCT_RELATED_PRODUCT_CATEGORY('{}',{},({}))",
            name,
            desc,
            refs(products)
        ))
    }

    // ── Color / presentation entities ─────────────────────────────────

    /// Emit STYLED_ITEM + full presentation chain for a single ADVANCED_FACE.
    ///
    /// STEP presentation chain:
    ///   COLOUR_RGB → FILL_AREA_STYLE_COLOUR → FILL_AREA_STYLE
    ///   → SURFACE_STYLE_FILL_AREA → SURFACE_STYLE_USAGE
    ///   → SURFACE_SIDE_STYLE → PRESENTATION_STYLE_ASSIGNMENT
    ///   → STYLED_ITEM (references the face)
    fn write_face_color(&mut self, face_id: u64, color: Color) {
        let rgb = self.colour_rgb("face_color", color.r, color.g, color.b);
        let fill_color = self.fill_area_style_colour("", rgb);
        let fill_style = self.fill_area_style("", &[fill_color]);
        let fill_area = self.surface_style_fill_area(fill_style);
        let side_style = self.surface_side_style("", &[fill_area]);
        let style_usage = self.surface_style_usage(".BOTH.", side_style);
        let psa = self.presentation_style_assignment(&[style_usage]);
        self.styled_item("color", &[psa], face_id);
    }

    fn colour_rgb(&mut self, name: &str, r: f64, g: f64, b: f64) -> u64 {
        self.push(format!("COLOUR_RGB('{}',{:.6},{:.6},{:.6})", name, r, g, b))
    }

    fn fill_area_style_colour(&mut self, name: &str, colour: u64) -> u64 {
        self.push(format!("FILL_AREA_STYLE_COLOUR('{}',#{})", name, colour))
    }

    fn fill_area_style(&mut self, name: &str, styles: &[u64]) -> u64 {
        self.push(format!("FILL_AREA_STYLE('{}',({}))", name, refs(styles)))
    }

    fn surface_style_fill_area(&mut self, style: u64) -> u64 {
        self.push(format!("SURFACE_STYLE_FILL_AREA(#{})", style))
    }

    fn surface_side_style(&mut self, name: &str, styles: &[u64]) -> u64 {
        self.push(format!("SURFACE_SIDE_STYLE('{}',({}))", name, refs(styles)))
    }

    fn surface_style_usage(&mut self, side: &str, style: u64) -> u64 {
        self.push(format!("SURFACE_STYLE_USAGE({},#{})", side, style))
    }

    fn presentation_style_assignment(&mut self, styles: &[u64]) -> u64 {
        self.push(format!("PRESENTATION_STYLE_ASSIGNMENT(({}))", refs(styles)))
    }

    fn styled_item(&mut self, name: &str, styles: &[u64], item: u64) -> u64 {
        self.push(format!(
            "STYLED_ITEM('{}',({},),#{})",
            name,
            refs(styles),
            item
        ))
    }

    fn general_property(&mut self, name: &str, description: &str) -> u64 {
        self.push(format!(
            "GENERAL_PROPERTY('{}','{}',$)",
            escape_step_string(name),
            escape_step_string(description)
        ))
    }

    fn property_definition(&mut self, name: &str, description: &str, reference: u64) -> u64 {
        self.push(format!(
            "PROPERTY_DEFINITION('{}','{}',#{})",
            escape_step_string(name),
            escape_step_string(description),
            reference
        ))
    }

    fn write_ap242_metadata(&mut self, metadata: &StepAp242Metadata) {
        for pdr in &metadata.property_definition_representations {
            let pd = opt_ref_token(pdr.property_definition_id);
            let rep = opt_ref_token(pdr.representation_id);
            self.push(format!("PROPERTY_DEFINITION_REPRESENTATION({},{})", pd, rep));
        }
        for loc in &metadata.dimensional_locations {
            let name = opt_step_string(loc.name.as_deref());
            let desc = opt_step_string(loc.description.as_deref());
            let from_id = opt_ref_token(loc.from_entity_id);
            let to_id = opt_ref_token(loc.to_entity_id);
            self.push(format!(
                "DIMENSIONAL_LOCATION({},{},{},{})",
                name, desc, from_id, to_id
            ));
        }
        for size in &metadata.dimensional_sizes {
            let name = opt_step_string(size.name.as_deref());
            let desc = opt_step_string(size.description.as_deref());
            let shape = opt_ref_token(size.shape_aspect_id);
            self.push(format!("DIMENSIONAL_SIZE({},{},{})", name, desc, shape));
        }
        for tol in &metadata.geometric_tolerances {
            let name = opt_step_string(tol.name.as_deref());
            let desc = opt_step_string(tol.description.as_deref());
            let val = opt_ref_token(tol.value_entity_id);
            let shape = opt_ref_token(tol.shape_aspect_id);
            self.push(format!(
                "GEOMETRIC_TOLERANCE({},{},{},{})",
                name, desc, val, shape
            ));
        }
        for tol in &metadata.geometric_tolerances_with_datum_references {
            let name = opt_step_string(tol.name.as_deref());
            let desc = opt_step_string(tol.description.as_deref());
            let val = opt_ref_token(tol.value_entity_id);
            let shape = opt_ref_token(tol.shape_aspect_id);
            let datum = opt_ref_token(tol.datum_system_id);
            self.push(format!(
                "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE({},{},{},{},{})",
                name, desc, val, shape, datum
            ));
        }
        for datum in &metadata.datums {
            let name = opt_step_string(datum.name.as_deref());
            let desc = opt_step_string(datum.description.as_deref());
            let shape = opt_ref_token(datum.shape_aspect_id);
            self.push(format!("DATUM({},{},{})", name, desc, shape));
        }
        for system in &metadata.datum_systems {
            let name = opt_step_string(system.name.as_deref());
            let desc = opt_step_string(system.description.as_deref());
            let refs = if system.datum_ids.is_empty() {
                "$".to_string()
            } else {
                format!(
                    "({})",
                    system
                        .datum_ids
                        .iter()
                        .map(|id| format!("#{}", id))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            };
            self.push(format!("DATUM_SYSTEM({},{},{})", name, desc, refs));
        }
        for pair in &metadata.kinematic_pairs {
            let entity = sanitize_step_entity_name(&pair.entity_type);
            let name = opt_step_string(pair.name.as_deref());
            let desc = opt_step_string(pair.description.as_deref());
            let refs = if pair.related_entity_ids.is_empty() {
                "$".to_string()
            } else {
                pair
                    .related_entity_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            };
            self.push(format!("{}({},{},{})", entity, name, desc, refs));
        }
        // Write GDT extended entities
        for tol in &metadata.dimensional_tolerances {
            self.write_dimensional_tolerance(tol);
        }
        for val in &metadata.tolerance_values {
            self.write_tolerance_value(val);
        }
        for tol in &metadata.position_tolerances {
            self.write_position_tolerance(tol);
        }
        for tol in &metadata.orientation_tolerances {
            self.write_orientation_tolerance(tol);
        }
        for tol in &metadata.form_tolerances {
            self.write_form_tolerance(tol);
        }
        for tol in &metadata.runout_tolerances {
            self.write_runout_tolerance(tol);
        }
        for tol in &metadata.profile_tolerances {
            self.write_profile_tolerance(tol);
        }
        for frame in &metadata.datum_reference_frames {
            self.write_datum_reference_frame(frame);
        }
        for target in &metadata.datum_targets {
            self.write_datum_target(target);
        }
        for def in &metadata.tolerance_zone_definitions_enhanced {
            self.write_tolerance_zone_definition_enhanced(def);
        }
        // Write view and camera entities
        for view in &metadata.views {
            self.write_view(view);
        }
        for camera in &metadata.cameras {
            self.write_camera_model_d3(camera);
        }
        for volume in &metadata.view_volumes {
            self.write_view_volume(volume);
        }
        // Write annotation entities
        for note in &metadata.notes {
            self.write_note(note);
        }
        for plane in &metadata.annotation_planes {
            self.write_annotation_plane(plane);
        }
        for occurrence in &metadata.annotation_occurrences {
            self.write_annotation_occurrence(occurrence);
        }
        for curve in &metadata.dimension_curves {
            self.write_dimension_curve(curve);
        }
        for symbol in &metadata.terminator_symbols {
            self.write_terminator_symbol(symbol);
        }
        for callout in &metadata.datum_feature_callouts {
            self.write_datum_feature_callout(callout);
        }
    }

    fn write_dimensional_tolerance(&mut self, tol: &StepDimensionalTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let dim_char = opt_ref_token(tol.dimensional_characteristic_id);
        let upper = tol.upper_tolerance.map(|v| format!("{}", v)).unwrap_or_else(|| "$".to_string());
        let lower = tol.lower_tolerance.map(|v| format!("{}", v)).unwrap_or_else(|| "$".to_string());
        let unit = opt_step_string(tol.unit.as_deref());
        self.push(format!(
            "DIMENSIONAL_TOLERANCE({},{},{},{},{},{})",
            name, desc, dim_char, upper, lower, unit
        ));
    }

    fn write_tolerance_value(&mut self, val: &StepToleranceValue) {
        let name = opt_step_string(val.name.as_deref());
        let value = val.value;
        let unit = opt_step_string(val.unit.as_deref());
        self.push(format!(
            "MEASURE_REPRESENTATION_ITEM({},{},{})",
            name, value, unit
        ));
    }

    fn write_position_tolerance(&mut self, tol: &StepPositionTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let projected = if tol.projected { ".T." } else { ".F." };
        let proj_height = tol.projected_height.map(|h| format!("{}", h)).unwrap_or_else(|| "$".to_string());
        self.push(format!(
            "POSITION_TOLERANCE({},{},{},{},{},{},{})",
            name, desc, val, shape, datum, projected, proj_height
        ));
    }

    fn write_orientation_tolerance(&mut self, tol: &StepOrientationTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.orientation_type {
            OrientationToleranceType::Angularity => "ANGULARITY_TOLERANCE",
            OrientationToleranceType::Perpendicularity => "PERPENDICULARITY_TOLERANCE",
            OrientationToleranceType::Parallelism => "PARALLELISM_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_form_tolerance(&mut self, tol: &StepFormTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let entity_name = match tol.form_type {
            FormToleranceType::Flatness => "FLATNESS_TOLERANCE",
            FormToleranceType::Straightness => "STRAIGHTNESS_TOLERANCE",
            FormToleranceType::Circularity => "CIRCULARITY_TOLERANCE",
            FormToleranceType::Cylindricity => "CYLINDRICITY_TOLERANCE",
        };
        self.push(format!("{}({},{},{},{})", entity_name, name, desc, val, shape));
    }

    fn write_runout_tolerance(&mut self, tol: &StepRunoutTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.runout_type {
            RunoutToleranceType::CircularRunout => "CIRCULAR_RUNOUT_TOLERANCE",
            RunoutToleranceType::TotalRunout => "TOTAL_RUNOUT_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_profile_tolerance(&mut self, tol: &StepProfileTolerance) {
        let name = opt_step_string(tol.name.as_deref());
        let desc = opt_step_string(tol.description.as_deref());
        let val = opt_ref_token(tol.value_entity_id);
        let shape = opt_ref_token(tol.shape_aspect_id);
        let datum = opt_ref_token(tol.datum_system_id);
        let entity_name = match tol.profile_type {
            ProfileToleranceType::ProfileOfALine => "LINE_PROFILE_TOLERANCE",
            ProfileToleranceType::ProfileOfASurface => "SURFACE_PROFILE_TOLERANCE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, val, shape, datum
        ));
    }

    fn write_datum_reference_frame(&mut self, frame: &StepDatumReferenceFrame) {
        let name = opt_step_string(frame.name.as_deref());
        let desc = opt_step_string(frame.description.as_deref());
        let refs = if frame.datum_system_ids.is_empty() {
            "$".to_string()
        } else {
            format!(
                "({})",
                frame
                    .datum_system_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        self.push(format!("DATUM_REFERENCE_FRAME({},{},{})", name, desc, refs));
    }

    fn write_datum_target(&mut self, target: &StepDatumTarget) {
        let name = opt_step_string(target.name.as_deref());
        let desc = opt_step_string(target.description.as_deref());
        let target_id = opt_step_string(target.target_identifier.as_deref());
        let datum = opt_ref_token(target.datum_id);
        let shape = opt_ref_token(target.shape_aspect_id);
        let entity_name = match target.target_type {
            DatumTargetType::Point => "DATUM_TARGET_POINT",
            DatumTargetType::Line => "DATUM_TARGET_LINE",
            DatumTargetType::Area => "DATUM_TARGET_AREA",
            DatumTargetType::AreaCircle => "DATUM_TARGET_CIRCLE",
            DatumTargetType::AreaRectangle => "DATUM_TARGET_RECTANGLE",
        };
        self.push(format!(
            "{}({},{},{},{},{})",
            entity_name, name, desc, target_id, datum, shape
        ));
    }

    fn write_tolerance_zone_definition_enhanced(&mut self, def: &StepToleranceZoneDefinitionEnhanced) {
        let name = opt_step_string(def.name.as_deref());
        let desc = opt_step_string(def.description.as_deref());
        let zone = opt_ref_token(def.tolerance_zone_id);
        let shape = opt_ref_token(def.shape_aspect_id);
        let defining_shape = opt_ref_token(def.defining_shape_aspect_id);
        self.push(format!(
            "TOLERANCE_ZONE_DEFINITION({},{},{},{},{})",
            name, desc, zone, shape, defining_shape
        ));
    }

    // ── View and Camera write functions (AP242) ───────────────────────────────

    fn write_view(&mut self, view: &StepView) {
        let name = opt_step_string(view.name.as_deref());
        let desc = opt_step_string(view.description.as_deref());
        let camera = opt_ref_token(view.camera_model_id);
        let view_type = opt_step_string(view.view_type.as_deref());
        self.push(format!(
            "CAMERA_MODEL_D3({},{},{},{})",
            name, desc, camera, view_type
        ));
    }

    fn write_camera_model_d3(&mut self, camera: &StepCameraModelD3) {
        let name = opt_step_string(camera.name.as_deref());
        let view_ref = opt_ref_token(camera.view_reference_system_id);
        let view_volume = opt_ref_token(camera.view_volume_id);
        let perspective = if camera.perspective { ".T." } else { ".F." };
        self.push(format!(
            "CAMERA_MODEL_D3({},{},{},{})",
            name, view_ref, view_volume, perspective
        ));
    }

    fn write_view_volume(&mut self, volume: &StepViewVolume) {
        let name = opt_step_string(volume.name.as_deref());
        let vol_type = match volume.volume_type {
            ViewVolumeType::Orthographic => ".F.",
            ViewVolumeType::Perspective => ".T.",
            ViewVolumeType::Unknown => "$",
        };
        let view_center = volume.view_center
            .map(|v| format!("({},{},{})", v[0], v[1], v[2]))
            .unwrap_or_else(|| "$".to_string());
        let view_plane_dist = volume.view_plane_distance
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        let up_dir = volume.up_direction
            .map(|v| format!("({},{},{})", v[0], v[1], v[2]))
            .unwrap_or_else(|| "$".to_string());
        let width = volume.view_window_width
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        let height = volume.view_window_height
            .map(|v| format!("{}", v))
            .unwrap_or_else(|| "$".to_string());
        self.push(format!(
            "VIEW_VOLUME({},{},{},{},{},{},{})",
            name, vol_type, view_center, view_plane_dist, up_dir, width, height
        ));
    }

    // ── Annotation write functions (AP242) ────────────────────────────────────

    fn write_note(&mut self, note: &StepNote) {
        let name = opt_step_string(note.name.as_deref());
        let desc = opt_step_string(note.description.as_deref());
        let text = opt_step_string(note.text.as_deref());
        let plane = opt_ref_token(note.annotation_plane_id);
        let geom = opt_ref_token(note.associated_geometry_id);
        self.push(format!(
            "DESCRIPTIVE_REPRESENTATION_ITEM({},{},{},{},{})",
            name, desc, text, plane, geom
        ));
    }

    fn write_annotation_plane(&mut self, plane: &StepAnnotationPlane) {
        let name = opt_step_string(plane.name.as_deref());
        let plane_id = opt_ref_token(plane.plane_id);
        let occurrence = opt_ref_token(plane.annotation_occurrence_id);
        self.push(format!(
            "ANNOTATION_PLANE({},{},{})",
            name, plane_id, occurrence
        ));
    }

    fn write_annotation_occurrence(&mut self, occurrence: &StepAnnotationOccurrence) {
        let name = opt_step_string(occurrence.name.as_deref());
        let style = opt_ref_token(occurrence.style_id);
        let fill = opt_ref_token(occurrence.fill_area_id);
        let shape = opt_ref_token(occurrence.shape_aspect_id);
        self.push(format!(
            "ANNOTATION_OCCURRENCE({},{},{},{})",
            name, style, fill, shape
        ));
    }

    fn write_dimension_curve(&mut self, curve: &StepDimensionCurve) {
        let name = opt_step_string(curve.name.as_deref());
        let curve_id = opt_ref_token(curve.curve_id);
        let plane = opt_ref_token(curve.annotation_plane_id);
        let terminators = if curve.terminator_ids.is_empty() {
            "$".to_string()
        } else {
            format!(
                "({})",
                curve.terminator_ids
                    .iter()
                    .map(|id| format!("#{}", id))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        self.push(format!(
            "DIMENSION_CURVE({},{},{},{})",
            name, curve_id, plane, terminators
        ));
    }

    fn write_terminator_symbol(&mut self, symbol: &StepTerminatorSymbol) {
        let name = opt_step_string(symbol.name.as_deref());
        let curve = opt_ref_token(symbol.annotated_curve_id);
        let term_type = match symbol.terminator_type {
            TerminatorType::Arrow => "FILLED_ARROW",
            TerminatorType::Dot => "FILLED_DOT",
            TerminatorType::OpenArrow => "OPEN_ARROW",
            TerminatorType::ClosedArrow => "FILLED_ARROW",
            TerminatorType::Origin => "ORIGIN_SYMBOL",
            TerminatorType::Unknown => "DIMENSION_CURVE_TERMINATOR",
        };
        self.push(format!(
            "{}({},{})",
            term_type, name, curve
        ));
    }

    fn write_datum_feature_callout(&mut self, callout: &StepDatumFeatureCallout) {
        let name = opt_step_string(callout.name.as_deref());
        let datum_id = opt_step_string(callout.datum_identifier.as_deref());
        let plane = opt_ref_token(callout.annotation_plane_id);
        self.push(format!(
            "DATUM_FEATURE_CALLOUT({},{},{})",
            name, datum_id, plane
        ));
    }

    fn push(&mut self, body: String) -> u64 {
        let id = self.next_id;
        self.next_id += 1;
        self.records.push(format!("#{}={};", id, body));
        id
    }
}

fn esc_step(s: &str) -> String {
    s.replace('\'', "''")
}

fn escape_step_string(s: &str) -> String {
    s.replace('\'', "''")
}

fn opt_ref_token(v: Option<u64>) -> String {
    v.map(|id| format!("#{}", id))
        .unwrap_or_else(|| "$".to_string())
}

fn opt_step_string(s: Option<&str>) -> String {
    s.map(|v| format!("'{}'", escape_step_string(v)))
        .unwrap_or_else(|| "$".to_string())
}

fn sanitize_step_entity_name(name: &str) -> String {
    let mut out = String::with_capacity(name.len());
    for ch in name.chars() {
        let up = ch.to_ascii_uppercase();
        if up.is_ascii_alphanumeric() || up == '_' {
            out.push(up);
        } else {
            out.push('_');
        }
    }
    if out.is_empty() {
        "KINEMATIC_PAIR".to_string()
    } else {
        out
    }
}

#[derive(Clone, Copy)]
struct OrientedEdgeExport {
    edge_idx: usize,
    start: usize,
    #[allow(dead_code)]
    end: usize,
    forward: bool,
}

fn oriented_face_edges(brep: &BRep, face: &Face) -> Vec<OrientedEdgeExport> {
    face.outer_wire
        .edges
        .iter()
        .filter_map(|we| {
            let edge = brep.edges.get(we.idx)?;
            let (start, end) = if we.forward {
                (edge.start, edge.end)
            } else {
                (edge.end, edge.start)
            };
            Some(OrientedEdgeExport {
                edge_idx: we.idx,
                start,
                end,
                forward: we.forward,
            })
        })
        .collect()
}

fn compute_face_normal(points: &[glam::DVec3]) -> Option<glam::DVec3> {
    if points.len() < 3 {
        return None;
    }
    let origin = points[0];
    for i in 1..points.len().saturating_sub(1) {
        let a = points[i] - origin;
        let b = points[i + 1] - origin;
        let n = a.cross(b);
        if n.length_squared() > 1e-12 {
            return Some(n.normalize());
        }
    }
    None
}

fn refs(items: &[u64]) -> String {
    items
        .iter()
        .map(|id| format!("#{}", id))
        .collect::<Vec<_>>()
        .join(",")
}

fn bool_token(value: bool) -> &'static str {
    if value { ".T." } else { ".F." }
}

fn dvec3_to_array(v: glam::DVec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn vector_length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = vector_length(v);
    if len <= 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn normalize2(v: [f64; 2]) -> [f64; 2] {
    let len = (v[0] * v[0] + v[1] * v[1]).sqrt();
    if len <= 1e-12 {
        [1.0, 0.0]
    } else {
        [v[0] / len, v[1] / len]
    }
}

fn project_to_plane(v: [f64; 3], normal: [f64; 3]) -> [f64; 3] {
    let dot = v[0] * normal[0] + v[1] * normal[1] + v[2] * normal[2];
    [
        v[0] - normal[0] * dot,
        v[1] - normal[1] * dot,
        v[2] - normal[2] * dot,
    ]
}

fn orthogonal_dir(normal: [f64; 3]) -> [f64; 3] {
    let helper = if normal[1].abs() < 0.9 {
        [0.0, 1.0, 0.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    normalize(cross(normal, helper))
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn any_perpendicular_dvec3(v: glam::DVec3) -> glam::DVec3 {
    let helper = if v.dot(glam::DVec3::Y).abs() < 0.9 {
        glam::DVec3::Y
    } else {
        glam::DVec3::X
    };
    v.cross(helper).normalize_or_zero()
}

/// Compress an expanded knot vector into (multiplicities, distinct_knot_values).
fn compress_knot_vector(knots: &[f64]) -> (Vec<usize>, Vec<f64>) {
    let mut mults: Vec<usize> = Vec::new();
    let mut vals: Vec<f64> = Vec::new();
    for &k in knots {
        if let Some(last) = vals.last()
            && (k - last).abs() < 1e-12
        {
            *mults.last_mut().expect("mults is non-empty by construction") += 1;
            continue;
        }
        vals.push(k);
        mults.push(1);
    }
    (mults, vals)
}

/// Detect which edge indices appear more than once in the face's outer wire.
/// These are seam edges on periodic surfaces.
fn detect_seam_edge_indices(face: &Face) -> BTreeSet<usize> {
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for we in &face.outer_wire.edges {
        *counts.entry(we.idx).or_insert(0) += 1;
    }
    counts
        .into_iter()
        .filter(|&(_, c)| c >= 2)
        .map(|(idx, _)| idx)
        .collect()
}

/// Extract a representative normal/axis from an analytic surface, used as a
/// fallback when boundary loop points are collinear (e.g. seam faces).
fn surface_normal(face_surface: Option<Surface3>) -> Option<glam::DVec3> {
    match face_surface? {
        Surface3::Plane(p) => Some(p.normal),
        Surface3::Cylinder(c) => Some(c.axis),
        Surface3::Sphere(s) => Some(s.axis),
        Surface3::Cone(c) => Some(c.axis),
        Surface3::Torus(t) => Some(t.axis),
        Surface3::Ellipsoid(e) => Some(e.axis),
        Surface3::Helicoid(h) => Some(h.axis),
        Surface3::Pipe(_) => None,
        Surface3::BSpline(_) => None,
        Surface3::LinearExtrusion(_)
        | Surface3::Revolution(_)
        | Surface3::Ruled(_)
        | Surface3::Coons(_)
        | Surface3::TriBezier(_)
        | Surface3::Bezier(_)
        | Surface3::Offset(_) => None,
        Surface3::Trimmed(ts) => surface_normal(Some(*ts.basis)),
    }
}

fn find_plane_surface_for_edge(
    brep: &BRep,
    edge_idx: usize,
) -> Option<(usize, rcad_kernel::geom::Plane)> {
    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let has_edge = oriented_face_edges(brep, face)
                    .iter()
                    .any(|entry| entry.edge_idx == edge_idx);
                if has_edge {
                    if let Some(Some(surface_idx)) = brep.geom.face_surface.get(face_index).copied() {
                        if let Some(Surface3::Plane(plane)) =
                            brep.geom.surfaces.get(surface_idx).cloned()
                        {
                            return Some((surface_idx, plane));
                        }
                    }
                }
                face_index += 1;
            }
        }
    }
    None
}

fn infer_spherical_circle_for_edge(brep: &BRep, edge_idx: usize) -> Option<Curve3> {
    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p2 = brep.vertices.get(edge.end)?.point;
    if (p2 - p0).length_squared() <= 1.0e-18 {
        return None;
    }

    let pcurves = brep.geom.edge_pcurves.get(edge_idx)?;
    for pc in pcurves {
        let Surface3::Sphere(sphere) = *brep.geom.surfaces.get(pc.surface_idx)? else {
            continue;
        };
        let curve2d = brep.geom.curve2ds.get(pc.curve2d_idx)?;
        let [t0, t1] = brep
            .geom
            .curve2d_range
            .get(pc.curve2d_idx)
            .and_then(|r| *r)
            .unwrap_or_else(|| default_curve2d_range(curve2d));
        if !t0.is_finite() || !t1.is_finite() {
            continue;
        }

        let tm = 0.5 * (t0 + t1);
        let uv_mid = curve2d.point_at(tm);
        let p1 = sphere.point_at(uv_mid.x, uv_mid.y);
        if (p1 - p0).length_squared() <= 1.0e-18 || (p1 - p2).length_squared() <= 1.0e-18 {
            continue;
        }

        let Some((center, normal, radius)) = circumcircle_3d(p0, p1, p2) else {
            continue;
        };
        if !radius.is_finite() || radius <= 1.0e-9 {
            continue;
        }

        // Keep only circles compatible with the supporting sphere.
        let center_dist = (center - sphere.center).length();
        if center_dist > sphere.radius + 1.0e-6 {
            continue;
        }

        return Some(Curve3::Circle(rcad_kernel::geom::Circle3 {
            center,
            normal,
            radius,
        }));
    }
    None
}

fn default_curve2d_range(curve2d: &Curve2d) -> [f64; 2] {
    match curve2d {
        Curve2d::Line(_) => [0.0, 1.0],
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => [0.0, 2.0 * std::f64::consts::PI],
        Curve2d::BSpline(bs) => {
            if bs.knots.is_empty() {
                [0.0, 1.0]
            } else {
                [
                    bs.knots.first().copied().unwrap_or(0.0),
                    bs.knots.last().copied().unwrap_or(1.0),
                ]
            }
        }
        Curve2d::Bezier(_)
        | Curve2d::CircleInvolute(_)
        | Curve2d::ArchimedeanSpiral(_)
        | Curve2d::LogarithmicSpiral(_)
        | Curve2d::SineWave(_) => [0.0, 1.0],
    }
}

fn circumcircle_3d(
    p0: glam::DVec3,
    p1: glam::DVec3,
    p2: glam::DVec3,
) -> Option<(glam::DVec3, glam::DVec3, f64)> {
    let u = p1 - p0;
    let v = p2 - p0;
    let w = u.cross(v);
    let w2 = w.length_squared();
    if w2 <= 1.0e-18 {
        return None;
    }

    let center = p0 + (v.length_squared() * w.cross(u) + u.length_squared() * v.cross(w)) / (2.0 * w2);
    let radius = (center - p0).length();
    let normal = w.normalize_or_zero();
    if normal.length_squared() <= 1.0e-18 {
        return None;
    }
    Some((center, normal, radius))
}

fn find_topological_plane_for_edge(
    brep: &BRep,
    edge_idx: usize,
) -> Option<rcad_kernel::geom::Plane> {
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let has_edge = oriented_face_edges(brep, face)
                    .iter()
                    .any(|entry| entry.edge_idx == edge_idx);
                if !has_edge {
                    continue;
                }

                let loop_points: Vec<glam::DVec3> = oriented_face_edges(brep, face)
                    .iter()
                    .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                    .collect();
                let Some(origin) = loop_points.first().copied() else {
                    continue;
                };
                let Some(normal) = compute_face_normal(&loop_points) else {
                    continue;
                };

                return Some(rcad_kernel::geom::Plane { origin, normal });
            }
        }
    }
    None
}

fn planes_equivalent(
    a: &rcad_kernel::geom::Plane,
    b: &rcad_kernel::geom::Plane,
    tol: f64,
) -> bool {
    let na = a.normal.normalize_or_zero();
    let nb = b.normal.normalize_or_zero();
    if na.length_squared() < 1e-18 || nb.length_squared() < 1e-18 {
        return false;
    }
    let parallel = na.dot(nb).abs() >= 1.0 - 1e-6;
    if !parallel {
        return false;
    }
    let dist = (a.origin - b.origin).dot(na).abs();
    dist <= tol
}

fn find_peer_plane_for_line_edge(
    brep: &BRep,
    edge_idx: usize,
    exclude: rcad_kernel::geom::Plane,
) -> Option<rcad_kernel::geom::Plane> {
    find_peer_plane_surface_for_line_edge(brep, edge_idx, exclude).map(|(_, plane)| plane)
}

fn count_plane_face_occurrences_for_line_edge(brep: &BRep, edge_idx: usize) -> usize {
    let Some(edge) = brep.edges.get(edge_idx) else {
        return 0;
    };
    let Some(a0) = brep.vertices.get(edge.start).map(|v| v.point) else {
        return 0;
    };
    let Some(a1) = brep.vertices.get(edge.end).map(|v| v.point) else {
        return 0;
    };
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    let mut count = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });

                if matched {
                    let is_plane = if let Some(Some(surface_idx)) =
                        brep.geom.face_surface.get(face_index).copied()
                    {
                        matches!(brep.geom.surfaces.get(surface_idx), Some(Surface3::Plane(_)))
                    } else {
                        let loop_points: Vec<glam::DVec3> = entries
                            .iter()
                            .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                            .collect();
                        compute_face_normal(&loop_points).is_some()
                    };
                    if is_plane {
                        count += 1;
                    }
                }

                face_index += 1;
            }
        }
    }
    count
}

fn find_peer_plane_surface_for_line_edge(
    brep: &BRep,
    edge_idx: usize,
    exclude: rcad_kernel::geom::Plane,
) -> Option<(Option<usize>, rcad_kernel::geom::Plane)> {
    let edge = brep.edges.get(edge_idx)?;
    let a0 = brep.vertices.get(edge.start)?.point;
    let a1 = brep.vertices.get(edge.end)?.point;
    let same_point = |p: glam::DVec3, q: glam::DVec3| (p - q).length_squared() <= 1.0e-18;

    let mut face_index = 0usize;
    for solid in &brep.solids {
        for shell in &solid.shells {
            for face in &shell.faces {
                let entries = oriented_face_edges(brep, face);
                let matched_same_idx = entries.iter().any(|entry| entry.edge_idx == edge_idx);
                let matched_same_geom = entries.iter().any(|entry| {
                    let Some(candidate) = brep.edges.get(entry.edge_idx) else {
                        return false;
                    };
                    let Some(c0) = brep.vertices.get(candidate.start).map(|v| v.point) else {
                        return false;
                    };
                    let Some(c1) = brep.vertices.get(candidate.end).map(|v| v.point) else {
                        return false;
                    };
                    let same_dir = same_point(c0, a0) && same_point(c1, a1);
                    let opp_dir = same_point(c0, a1) && same_point(c1, a0);
                    same_dir || opp_dir
                });
                let matched = matched_same_idx || matched_same_geom;

                if matched {
                    let plane = if let Some(Some(surface_idx)) =
                        brep.geom.face_surface.get(face_index).copied()
                    {
                        match brep.geom.surfaces.get(surface_idx).cloned() {
                            Some(Surface3::Plane(p)) => Some(p),
                            _ => None,
                        }
                    } else {
                        None
                    }
                    .or_else(|| {
                        let loop_points: Vec<glam::DVec3> = entries
                            .iter()
                            .filter_map(|e| brep.vertices.get(e.start).map(|v| v.point))
                            .collect();
                        let origin = loop_points.first().copied()?;
                        let normal = compute_face_normal(&loop_points)?;
                        Some(rcad_kernel::geom::Plane { origin, normal })
                    });

                    if let Some(p) = plane
                        && !planes_equivalent(&p, &exclude, 1.0e-6)
                    {
                        let surface_idx = brep
                            .geom
                            .face_surface
                            .get(face_index)
                            .copied()
                            .flatten();
                        return Some((surface_idx, p));
                    }
                }

                face_index += 1;
            }
        }
    }

    None
}

fn synthesize_plane_pcurve_for_edge(
    brep: &BRep,
    edge_idx: usize,
    plane: &rcad_kernel::geom::Plane,
) -> Option<Curve2d> {
    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p1 = brep.vertices.get(edge.end)?.point;

    let normal = plane.normal.normalize_or_zero();
    if normal.length_squared() < 1e-18 {
        return None;
    }
    let u_axis = any_perpendicular_dvec3(normal);
    let v_axis = normal.cross(u_axis).normalize_or_zero();
    if v_axis.length_squared() < 1e-18 {
        return None;
    }

    let to_uv = |pt: glam::DVec3| -> glam::DVec2 {
        let d = pt - plane.origin;
        glam::DVec2::new(d.dot(u_axis), d.dot(v_axis))
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .copied()
        .flatten()
        .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

    match edge_curve {
        Some(Curve3::Line(_)) | None => {
            let uv0 = to_uv(p0);
            let uv1 = to_uv(p1);
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-18 {
                return None;
            }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: uv0,
                direction: dir,
            }))
        }
        Some(Curve3::Circle(c)) => {
            let center = to_uv(c.center);
            Some(Curve2d::Circle(rcad_kernel::geom::Circle2d {
                center,
                radius: c.radius.max(1e-9),
            }))
        }
        Some(Curve3::Ellipse(e)) => {
            let center = to_uv(e.center);
            let major = glam::DVec2::new(e.major_dir.dot(u_axis), e.major_dir.dot(v_axis));
            let major_dir = if major.length_squared() < 1e-18 {
                glam::DVec2::X
            } else {
                major.normalize()
            };
            Some(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
                center,
                major_dir,
                major_radius: e.major_radius.max(1e-9),
                minor_radius: e.minor_radius.max(1e-9),
            }))
        }
        _ => None,
    }
}

fn synthesize_edge_curve2d_on_face_frame(
    brep: &BRep,
    edge_idx: usize,
    face_origin: glam::DVec3,
    x_axis: glam::DVec3,
    normal: glam::DVec3,
) -> Option<Curve2d> {
    let x_axis = x_axis.normalize_or_zero();
    let normal = normal.normalize_or_zero();
    if x_axis.length_squared() < 1e-18 || normal.length_squared() < 1e-18 {
        return None;
    }
    let y_axis = normal.cross(x_axis).normalize_or_zero();
    if y_axis.length_squared() < 1e-18 {
        return None;
    }

    let edge = brep.edges.get(edge_idx)?;
    let p0 = brep.vertices.get(edge.start)?.point;
    let p1 = brep.vertices.get(edge.end)?.point;
    let to_uv = |pt: glam::DVec3| -> glam::DVec2 {
        let d = pt - face_origin;
        glam::DVec2::new(d.dot(x_axis), d.dot(y_axis))
    };

    let edge_curve = brep
        .geom
        .edge_curve
        .get(edge_idx)
        .copied()
        .flatten()
        .and_then(|curve_idx| brep.geom.curves.get(curve_idx).cloned());

    match edge_curve {
        Some(Curve3::Line(_)) | None => {
            let uv0 = to_uv(p0);
            let uv1 = to_uv(p1);
            let dir = uv1 - uv0;
            if dir.length_squared() < 1e-18 {
                return None;
            }
            Some(Curve2d::Line(rcad_kernel::geom::Line2d {
                origin: uv0,
                direction: dir,
            }))
        }
        Some(Curve3::Circle(c)) => {
            let center = to_uv(c.center);
            Some(Curve2d::Circle(rcad_kernel::geom::Circle2d {
                center,
                radius: c.radius.max(1e-9),
            }))
        }
        Some(Curve3::Ellipse(e)) => {
            let center = to_uv(e.center);
            let major = glam::DVec2::new(e.major_dir.dot(x_axis), e.major_dir.dot(y_axis));
            let major_dir = if major.length_squared() < 1e-18 {
                glam::DVec2::X
            } else {
                major.normalize()
            };
            Some(Curve2d::Ellipse(rcad_kernel::geom::Ellipse2d {
                center,
                major_dir,
                major_radius: e.major_radius.max(1e-9),
                minor_radius: e.minor_radius.max(1e-9),
            }))
        }
        _ => None,
    }
}

fn is_degenerate_face_wire(brep: &BRep, face: &Face) -> bool {
    if face.outer_wire.edges.len() < 3 {
        return true;
    }

    let unique_edges: BTreeSet<usize> = face.outer_wire.edges.iter().map(|we| we.idx).collect();
    if unique_edges.len() < 3 {
        return true;
    }

    let mut verts = BTreeSet::new();
    for we in &face.outer_wire.edges {
        if let Some(edge) = brep.edges.get(we.idx) {
            verts.insert(edge.start);
            verts.insert(edge.end);
        }
    }
    verts.len() < 3
}

fn face_orientation_for_surface(
    brep: &BRep,
    loop_points: &[glam::DVec3],
    face: &Face,
    surface: Option<&Surface3>,
    gmsh_strict: bool,
) -> bool {
    match surface {
        Some(Surface3::Plane(plane)) if gmsh_strict => {
            let plane_normal = canonicalize_axis_sign(plane.normal);
            let mut c = glam::DVec3::ZERO;
            let mut n = 0usize;
            for v in &brep.vertices {
                c += v.point;
                n += 1;
            }
            let brep_centroid = if n > 0 { c / (n as f64) } else { glam::DVec3::ZERO };

            let face_centroid = if loop_points.is_empty() {
                brep_centroid
            } else {
                loop_points
                    .iter()
                    .copied()
                    .fold(glam::DVec3::ZERO, |acc, p| acc + p)
                    / (loop_points.len() as f64)
            };

            let outward = (face_centroid - brep_centroid).dot(plane_normal);
            if outward.abs() <= 1e-12 {
                face.normal.dot(plane_normal) >= 0.0
            } else {
                outward >= 0.0
            }
        }
        Some(Surface3::Plane(plane)) => face.normal.dot(plane.normal) >= 0.0,
        _ => true,
    }
}

fn canonicalize_axis_sign(v: glam::DVec3) -> glam::DVec3 {
    let eps = 1e-12;
    let n = v.normalize_or_zero();
    if n.z.abs() > eps {
        if n.z >= 0.0 { n } else { -n }
    } else if n.y.abs() > eps {
        if n.y >= 0.0 { n } else { -n }
    } else if n.x.abs() > eps {
        if n.x >= 0.0 { n } else { -n }
    } else {
        glam::DVec3::Z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        StepDatum, StepDimensionalLocation, StepDimensionalSize, StepGeometricTolerance,
        StepGeometricToleranceWithDatumReference,
        StepPropertyDefinitionRepr, StepReader,
    };
    use glam::DVec3;
    use rcad_modeling::make_box_brep;
    use std::io::Cursor;
    const HFSS_STEP: &str = include_str!("../../../assets/hfss.step");

    #[test]
    fn exports_full_box_and_reimports() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");
        // Default now uses gmsh_strict for GMSH/OCCT compatibility
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("exported STEP should parse");
        assert!(!reparsed.edges.is_empty());
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
        assert!(step.contains("CLOSED_SHELL"));
        // gmsh_strict mode does NOT include wireframe overlays
        assert!(!step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(!step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("SHAPE_DEFINITION_REPRESENTATION"));
    }

    #[test]
    fn exports_non_strict_includes_wireframe() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: false,  // Non-strict mode includes wireframe
                ..Default::default()
            },
        );

        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(step.contains("MANIFOLD_SOLID_BREP"));
        // Non-strict mode includes wireframe overlay for visualization
        assert!(step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(step.contains("GEOMETRIC_CURVE_SET"));
        assert!(step.contains("SHAPE_REPRESENTATION_RELATIONSHIP"));
        assert!(step.contains("SHAPE_DEFINITION_REPRESENTATION"));
    }

    #[test]
    fn exports_gmsh_strict_primary_brep_only_for_full_solid() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: true,
                ..Default::default()
            },
        );

        assert!(step.contains("ADVANCED_BREP_SHAPE_REPRESENTATION"));
        assert!(!step.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION"));
        assert!(!step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"));
        assert!(!step.contains("GEOMETRIC_CURVE_SET"));
    }

    #[test]
    fn gmsh_strict_output_has_correct_structure() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: true,
                ..Default::default()
            },
        );

        // GMSH strict should have MANIFOLD_SOLID_BREP
        assert!(step.contains("MANIFOLD_SOLID_BREP"), "should contain MANIFOLD_SOLID_BREP");
        assert!(step.contains("CLOSED_SHELL"), "should contain CLOSED_SHELL");
        assert!(step.contains("ADVANCED_FACE"), "should contain ADVANCED_FACE");

        // Should NOT have wireframe overlay entities
        assert!(!step.contains("GEOMETRIC_CURVE_SET"), "should not contain GEOMETRIC_CURVE_SET");
        assert!(!step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME"), "should not contain wireframe");
    }

    #[test]
    fn non_strict_mode_includes_geometry_entities() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: false,
                ..Default::default()
            },
        );

        // Non-strict should have both solid and wireframe
        assert!(step.contains("MANIFOLD_SOLID_BREP"), "should contain MANIFOLD_SOLID_BREP");
        assert!(step.contains("GEOMETRIC_CURVE_SET") || step.contains("SHELL_BASED_SURFACE_MODEL"),
                "should contain geometric entities");
    }

    #[test]
    fn exports_general_properties_when_provided() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let props = vec![
            StepGeneralProperty {
                name: "PartNumber".to_string(),
                description: Some("PN-001".to_string()),
            },
            StepGeneralProperty {
                name: "Revision".to_string(),
                description: Some("A".to_string()),
            },
        ];
        let step = StepWriter::write_string_with_properties(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &props,
            StepProtocol::Ap242,
        );

        assert!(step.contains("GENERAL_PROPERTY('PartNumber','PN-001',$)"));
        assert!(step.contains("GENERAL_PROPERTY('Revision','A',$)"));
        assert!(step.contains("PROPERTY_DEFINITION('PartNumber','PN-001',#"));
        assert!(step.contains("PROPERTY_DEFINITION('Revision','A',#"));
    }

    #[test]
    fn exports_ap242_metadata_entities_when_provided() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            property_definition_representations: vec![StepPropertyDefinitionRepr {
                entity_id: 0,
                property_definition_id: Some(10),
                representation_id: Some(20),
            }],
            dimensional_locations: vec![StepDimensionalLocation {
                entity_id: 0,
                name: Some("d_loc".into()),
                description: Some("desc".into()),
                from_entity_id: Some(30),
                to_entity_id: Some(31),
            }],
            dimensional_sizes: vec![StepDimensionalSize {
                entity_id: 0,
                name: Some("d_size".into()),
                description: None,
                shape_aspect_id: Some(40),
            }],
            geometric_tolerances: vec![StepGeometricTolerance {
                entity_id: 0,
                name: Some("flatness".into()),
                description: Some("gtol".into()),
                value_entity_id: Some(50),
                shape_aspect_id: Some(60),
            }],
            geometric_tolerances_with_datum_references: vec![
                StepGeometricToleranceWithDatumReference {
                    entity_id: 0,
                    name: Some("position".into()),
                    description: Some("gtol_datum".into()),
                    value_entity_id: Some(51),
                    shape_aspect_id: Some(61),
                    datum_system_id: Some(71),
                },
            ],
            datums: vec![StepDatum {
                entity_id: 0,
                name: Some("A".into()),
                description: Some("primary".into()),
                shape_aspect_id: Some(70),
            }],
            datum_systems: vec![StepDatumSystem {
                entity_id: 0,
                name: Some("A_SYS".into()),
                description: Some("primary_system".into()),
                datum_ids: vec![70],
            }],
            kinematic_pairs: vec![StepKinematicPair {
                entity_id: 0,
                entity_type: "REVOLUTE_PAIR".into(),
                name: Some("hinge".into()),
                description: Some("joint".into()),
                related_entity_ids: vec![81, 82, 83],
            }],
            // GDT extended fields
            dimensional_tolerances: vec![],
            tolerance_values: vec![],
            position_tolerances: vec![],
            orientation_tolerances: vec![],
            form_tolerances: vec![],
            runout_tolerances: vec![],
            profile_tolerances: vec![],
            datum_reference_frames: vec![],
            datum_targets: vec![],
            tolerance_zone_definitions_enhanced: vec![],
            // View and annotation fields
            views: vec![],
            cameras: vec![],
            view_volumes: vec![],
            notes: vec![],
            annotation_planes: vec![],
            annotation_occurrences: vec![],
            dimension_curves: vec![],
            terminator_symbols: vec![],
            datum_feature_callouts: vec![],
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );

        assert!(step.contains("PROPERTY_DEFINITION_REPRESENTATION(#10,#20)"));
        assert!(step.contains("DIMENSIONAL_LOCATION('d_loc','desc',#30,#31)"));
        assert!(step.contains("DIMENSIONAL_SIZE('d_size',$,#40)"));
        assert!(step.contains("GEOMETRIC_TOLERANCE('flatness','gtol',#50,#60)"));
        assert!(step.contains(
            "GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE('position','gtol_datum',#51,#61,#71)"
        ));
        assert!(step.contains("DATUM('A','primary',#70)"));
        assert!(step.contains("DATUM_SYSTEM('A_SYS','primary_system',(#70))"));
        assert!(step.contains("REVOLUTE_PAIR('hinge','joint',#81,#82,#83)"));
    }

    #[test]
    fn ap242_metadata_write_read_roundtrip_smoke() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            property_definition_representations: vec![StepPropertyDefinitionRepr {
                entity_id: 0,
                property_definition_id: Some(11),
                representation_id: Some(22),
            }],
            dimensional_locations: vec![StepDimensionalLocation {
                entity_id: 0,
                name: Some("L1".into()),
                description: Some("loc".into()),
                from_entity_id: Some(33),
                to_entity_id: Some(44),
            }],
            dimensional_sizes: vec![StepDimensionalSize {
                entity_id: 0,
                name: Some("S1".into()),
                description: Some("size".into()),
                shape_aspect_id: Some(55),
            }],
            geometric_tolerances: vec![StepGeometricTolerance {
                entity_id: 0,
                name: Some("parallelism".into()),
                description: Some("tol".into()),
                value_entity_id: Some(66),
                shape_aspect_id: Some(77),
            }],
            geometric_tolerances_with_datum_references: vec![
                StepGeometricToleranceWithDatumReference {
                    entity_id: 0,
                    name: Some("perpendicularity".into()),
                    description: Some("tol_datum".into()),
                    value_entity_id: Some(67),
                    shape_aspect_id: Some(78),
                    datum_system_id: Some(89),
                },
            ],
            datums: vec![StepDatum {
                entity_id: 0,
                name: Some("B".into()),
                description: Some("secondary".into()),
                shape_aspect_id: Some(88),
            }],
            datum_systems: vec![StepDatumSystem {
                entity_id: 0,
                name: Some("B_SYS".into()),
                description: Some("secondary_system".into()),
                datum_ids: vec![88],
            }],
            kinematic_pairs: vec![StepKinematicPair {
                entity_id: 0,
                entity_type: "PRISMATIC_PAIR".into(),
                name: Some("slider".into()),
                description: Some("guide".into()),
                related_entity_ids: vec![90, 91],
            }],
            // GDT extended fields
            dimensional_tolerances: vec![StepDimensionalTolerance {
                entity_id: 0,
                name: Some("diam_tol".into()),
                description: Some("diameter tolerance".into()),
                dimensional_characteristic_id: Some(100),
                upper_tolerance: Some(0.05),
                lower_tolerance: Some(-0.05),
                unit: Some("mm".into()),
            }],
            tolerance_values: vec![StepToleranceValue {
                entity_id: 0,
                name: Some("tol_val".into()),
                value: 0.025,
                unit: Some("mm".into()),
            }],
            position_tolerances: vec![StepPositionTolerance {
                entity_id: 0,
                name: Some("pos_tol".into()),
                description: Some("positional tolerance".into()),
                value_entity_id: Some(101),
                shape_aspect_id: Some(102),
                datum_system_id: Some(103),
                projected: false,
                projected_height: None,
            }],
            orientation_tolerances: vec![StepOrientationTolerance {
                entity_id: 0,
                name: Some("ang_tol".into()),
                description: Some("angularity".into()),
                value_entity_id: Some(104),
                shape_aspect_id: Some(105),
                datum_system_id: Some(106),
                orientation_type: OrientationToleranceType::Angularity,
            }],
            form_tolerances: vec![StepFormTolerance {
                entity_id: 0,
                name: Some("flat_tol".into()),
                description: Some("flatness".into()),
                value_entity_id: Some(107),
                shape_aspect_id: Some(108),
                form_type: FormToleranceType::Flatness,
            }],
            runout_tolerances: vec![StepRunoutTolerance {
                entity_id: 0,
                name: Some("cr_tol".into()),
                description: Some("circular runout".into()),
                value_entity_id: Some(109),
                shape_aspect_id: Some(110),
                datum_system_id: Some(111),
                runout_type: RunoutToleranceType::CircularRunout,
            }],
            profile_tolerances: vec![StepProfileTolerance {
                entity_id: 0,
                name: Some("lin_tol".into()),
                description: Some("profile of a line".into()),
                value_entity_id: Some(112),
                shape_aspect_id: Some(113),
                datum_system_id: None,
                profile_type: ProfileToleranceType::ProfileOfALine,
            }],
            datum_reference_frames: vec![StepDatumReferenceFrame {
                entity_id: 0,
                name: Some("DRF1".into()),
                description: Some("datum reference frame".into()),
                datum_system_ids: vec![88],
            }],
            datum_targets: vec![StepDatumTarget {
                entity_id: 0,
                name: Some("A1".into()),
                description: Some("datum target".into()),
                target_identifier: Some("A1".into()),
                datum_id: Some(88),
                target_type: DatumTargetType::Point,
                shape_aspect_id: Some(114),
            }],
            tolerance_zone_definitions_enhanced: vec![StepToleranceZoneDefinitionEnhanced {
                entity_id: 0,
                name: Some("cylindrical".into()),
                description: Some("symmetric".into()),
                tolerance_zone_id: Some(115),
                shape_aspect_id: Some(116),
                zone_shape: ToleranceZoneShape::Cylindrical,
                zone_position: ToleranceZonePosition::Symmetric,
                defining_shape_aspect_id: None,
            }],
            // View and annotation fields
            views: vec![],
            cameras: vec![],
            view_volumes: vec![],
            notes: vec![],
            annotation_planes: vec![],
            annotation_occurrences: vec![],
            dimension_curves: vec![],
            terminator_symbols: vec![],
            datum_feature_callouts: vec![],
        };

        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );

        let (_parsed_brep, doc_meta) =
            StepReader::parse_string_with_metadata(&step).expect("AP242 metadata STEP should parse");
        assert_eq!(doc_meta.property_definition_representations.len(), 1);
        assert_eq!(doc_meta.dimensional_locations.len(), 1);
        assert_eq!(doc_meta.dimensional_sizes.len(), 1);
        assert_eq!(doc_meta.geometric_tolerances.len(), 1);
        assert_eq!(doc_meta.geometric_tolerances_with_datum_references.len(), 1);
        assert_eq!(doc_meta.datums.len(), 1);
        assert_eq!(doc_meta.datum_systems.len(), 1);
        assert_eq!(doc_meta.kinematic_pairs.len(), 1);
        // GDT extended assertions
        assert_eq!(doc_meta.dimensional_tolerances.len(), 1);
        assert_eq!(doc_meta.tolerance_values.len(), 1);
        assert_eq!(doc_meta.position_tolerances.len(), 1);
        assert_eq!(doc_meta.orientation_tolerances.len(), 1);
        assert_eq!(doc_meta.form_tolerances.len(), 1);
        assert_eq!(doc_meta.runout_tolerances.len(), 1);
        assert_eq!(doc_meta.profile_tolerances.len(), 1);
        assert_eq!(doc_meta.datum_reference_frames.len(), 1);
        assert_eq!(doc_meta.datum_targets.len(), 1);
        assert_eq!(doc_meta.tolerance_zone_definitions_enhanced.len(), 1);
    }

    #[test]
    fn exports_selected_edges_without_faces() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[0, 1],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("edge-only export should parse");
        assert!(reparsed.solids.is_empty());
        assert_eq!(reparsed.edges.len(), 2);
    }

    #[test]
    fn stream_write_then_stream_read_roundtrip() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 2.0, 3.0)
            .expect("test box should be valid");

        let mut buf = Vec::<u8>::new();
        StepWriter::write_to(
            &mut buf,
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        )
        .expect("stream write should succeed");

        let reparsed =
            StepReader::parse_reader(Cursor::new(buf)).expect("stream read should parse");
        assert!(!reparsed.edges.is_empty());
        assert!(!reparsed.solids.is_empty());
    }

    #[test]
    fn exports_selected_faces_via_shell_based_surface_model() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[0],
                selected_edges: &[],
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("selected-face export should parse");
        assert!(!reparsed.solids.is_empty());
        assert!(step.contains("OPEN_SHELL"));
        assert!(step.contains("SHELL_BASED_SURFACE_MODEL"));
        assert!(step.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION"));
    }

    #[test]
    fn exports_analytic_surfaces_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        // All analytic surfaces should now be exported properly as ADVANCED_FACE
        // with their respective surface types, including seam faces on
        // spheres and cones.
        assert!(step.contains("SPHERICAL_SURFACE"));
        assert!(step.contains("CYLINDRICAL_SURFACE"));
        assert!(step.contains("TOROIDAL_SURFACE"));
        assert!(step.contains("CONICAL_SURFACE"));

        // Standalone 1D curves (GEOMETRIC_CURVE_SET) must also be exported
        // alongside the solid geometry.
        assert!(
            step.contains("GEOMETRIC_CURVE_SET"),
            "standalone wireframe edges should be exported"
        );
        assert!(
            step.contains("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION"),
            "wireframe shape representation should be present"
        );
    }

    #[test]
    fn round_trips_sphere_and_cone_surfaces() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        // Find the original cone half-angle and radius for comparison
        let mut orig_cone_angle = 0.0f64;
        let mut orig_cone_radius = 0.0f64;
        for surface in &brep.geom.surfaces {
            if let Surface3::Cone(c) = surface {
                orig_cone_angle = c.half_angle_rad;
                orig_cone_radius = c.radius;
            }
        }
        assert!(orig_cone_angle > 0.0, "should find a cone in hfss.step");

        // Use non-strict mode for complex geometry round-trip tests
        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                protocol: StepProtocol::Ap214,
                gmsh_strict: false,
                ..Default::default()
            },
        );

        let reparsed = StepReader::parse_string(&step).expect("re-exported STEP should parse");

        // Count faces with each surface type and verify cone parameters survive round-trip
        let mut sphere_count = 0usize;
        let mut cone_count = 0usize;
        for sid in reparsed.geom.face_surface.iter().flatten() {
            match reparsed.geom.surfaces.get(*sid) {
                Some(Surface3::Sphere(_)) => sphere_count += 1,
                Some(Surface3::Cone(c)) => {
                    cone_count += 1;
                    assert!(
                        (c.half_angle_rad - orig_cone_angle).abs() < 1e-6,
                        "cone half-angle drifted: original={} reparsed={}",
                        orig_cone_angle,
                        c.half_angle_rad,
                    );
                    assert!(
                        (c.radius - orig_cone_radius).abs() < 1e-6,
                        "cone radius drifted: original={} reparsed={}",
                        orig_cone_radius,
                        c.radius,
                    );
                }
                _ => {}
            }
        }
        assert!(
            sphere_count >= 1,
            "expected at least 1 sphere face after round-trip, got {}",
            sphere_count
        );
        assert!(
            cone_count >= 1,
            "expected at least 1 cone face after round-trip, got {}",
            cone_count
        );
    }

    #[test]
    fn exports_ellipsoid_surface_emits_semantic_tag() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Ellipsoid(rcad_kernel::EllipsoidalSurface {
            center: DVec3::new(0.5, 0.5, 0.5),
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            radius_x: 2.0,
            radius_y: 1.5,
            radius_z: 1.0,
        });

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                gmsh_strict: false,
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(step.contains("RCAD_ELLIPSOID"));

        let reparsed = StepReader::parse_string(&step).expect("ellipsoid fallback STEP should parse");
        let ellipsoid_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::Ellipsoid(_)))
            .count();
        assert!(
            ellipsoid_surfaces > 0,
            "expected at least one reparsed ellipsoid surface"
        );
    }

    #[test]
    fn exports_helicoid_surface_emits_semantic_tag() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Helicoid(rcad_kernel::HelicoidSurface {
            origin: DVec3::ZERO,
            axis: DVec3::Z,
            ref_dir: DVec3::X,
            pitch: 3.0,
        });

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                gmsh_strict: false,
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(step.contains("RCAD_HELICOID"));

        let reparsed = StepReader::parse_string(&step).expect("helicoid fallback STEP should parse");
        let helicoid_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::Helicoid(_)))
            .count();
        assert!(
            helicoid_surfaces > 0,
            "expected at least one reparsed helicoid surface"
        );
    }

    #[test]
    fn exports_coons_surface_via_bspline_fallback() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Coons(rcad_kernel::CoonsSurface {
            south: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::X,
            })),
            north: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 1.0, 1.0),
                direction: DVec3::X,
            })),
            west: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(0.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
            east: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::new(1.0, 0.0, 0.0),
                direction: DVec3::new(0.0, 1.0, 1.0),
            })),
        });

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                gmsh_strict: false,
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));
        assert!(!step.contains("RCAD_COONS"));

        let reparsed = StepReader::parse_string(&step).expect("Coons fallback STEP should parse");
        let bspline_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::BSpline(_)))
            .count();
        assert!(
            bspline_surfaces > 0,
            "expected at least one reparsed bspline surface from Coons fallback"
        );
    }

    #[test]
    fn exports_pipe_surface_via_bspline_fallback() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Pipe(rcad_kernel::PipeSurface {
            spine: Box::new(rcad_kernel::Curve3::Line(rcad_kernel::geom::Line3 {
                origin: DVec3::ZERO,
                direction: DVec3::Z,
            })),
            ref_dir: DVec3::X,
            radius: 1.25,
        });

        let step = StepWriter::write_string_with_options(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                gmsh_strict: false,
                ..Default::default()
            },
        );
        assert!(step.contains("B_SPLINE_SURFACE_WITH_KNOTS"));

        let reparsed = StepReader::parse_string(&step).expect("Pipe fallback STEP should parse");
        let bspline_surfaces = reparsed
            .geom
            .surfaces
            .iter()
            .filter(|surface| matches!(surface, Surface3::BSpline(_)))
            .count();
        assert!(
            bspline_surfaces > 0,
            "expected at least one reparsed bspline surface from Pipe fallback"
        );
    }

    // ── GDT Extended entity write tests ──────────────────────────────────────

    #[test]
    fn writes_dimensional_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            dimensional_tolerances: vec![StepDimensionalTolerance {
                entity_id: 0,
                name: Some("diam_tol".into()),
                description: Some("diameter tolerance".into()),
                dimensional_characteristic_id: Some(100),
                upper_tolerance: Some(0.05),
                lower_tolerance: Some(-0.05),
                unit: Some("mm".into()),
            }],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("DIMENSIONAL_TOLERANCE('diam_tol','diameter tolerance',#100,0.05,-0.05,'mm')"));
    }

    #[test]
    fn writes_tolerance_values() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            tolerance_values: vec![StepToleranceValue {
                entity_id: 0,
                name: Some("tol_value".into()),
                value: 0.025,
                unit: Some("mm".into()),
            }],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("MEASURE_REPRESENTATION_ITEM('tol_value',0.025,'mm')"));
    }

    #[test]
    fn writes_position_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            position_tolerances: vec![StepPositionTolerance {
                entity_id: 0,
                name: Some("pos_tol".into()),
                description: Some("positional tolerance".into()),
                value_entity_id: Some(20),
                shape_aspect_id: Some(30),
                datum_system_id: Some(40),
                projected: true,
                projected_height: Some(10.0),
            }],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("POSITION_TOLERANCE('pos_tol','positional tolerance',#20,#30,#40,.T.,10)"));
    }

    #[test]
    fn writes_orientation_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            orientation_tolerances: vec![
                StepOrientationTolerance {
                    entity_id: 0,
                    name: Some("ang_tol".into()),
                    description: Some("angularity".into()),
                    value_entity_id: Some(20),
                    shape_aspect_id: Some(30),
                    datum_system_id: Some(40),
                    orientation_type: OrientationToleranceType::Angularity,
                },
                StepOrientationTolerance {
                    entity_id: 0,
                    name: Some("perp_tol".into()),
                    description: Some("perpendicularity".into()),
                    value_entity_id: Some(21),
                    shape_aspect_id: Some(31),
                    datum_system_id: Some(41),
                    orientation_type: OrientationToleranceType::Perpendicularity,
                },
                StepOrientationTolerance {
                    entity_id: 0,
                    name: Some("para_tol".into()),
                    description: Some("parallelism".into()),
                    value_entity_id: Some(22),
                    shape_aspect_id: Some(32),
                    datum_system_id: Some(42),
                    orientation_type: OrientationToleranceType::Parallelism,
                },
            ],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("ANGULARITY_TOLERANCE('ang_tol','angularity',#20,#30,#40)"));
        assert!(step.contains("PERPENDICULARITY_TOLERANCE('perp_tol','perpendicularity',#21,#31,#41)"));
        assert!(step.contains("PARALLELISM_TOLERANCE('para_tol','parallelism',#22,#32,#42)"));
    }

    #[test]
    fn writes_form_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            form_tolerances: vec![
                StepFormTolerance {
                    entity_id: 0,
                    name: Some("flat_tol".into()),
                    description: Some("flatness".into()),
                    value_entity_id: Some(20),
                    shape_aspect_id: Some(30),
                    form_type: FormToleranceType::Flatness,
                },
                StepFormTolerance {
                    entity_id: 0,
                    name: Some("str_tol".into()),
                    description: Some("straightness".into()),
                    value_entity_id: Some(21),
                    shape_aspect_id: Some(31),
                    form_type: FormToleranceType::Straightness,
                },
                StepFormTolerance {
                    entity_id: 0,
                    name: Some("cir_tol".into()),
                    description: Some("circularity".into()),
                    value_entity_id: Some(22),
                    shape_aspect_id: Some(32),
                    form_type: FormToleranceType::Circularity,
                },
                StepFormTolerance {
                    entity_id: 0,
                    name: Some("cyl_tol".into()),
                    description: Some("cylindricity".into()),
                    value_entity_id: Some(23),
                    shape_aspect_id: Some(33),
                    form_type: FormToleranceType::Cylindricity,
                },
            ],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("FLATNESS_TOLERANCE('flat_tol','flatness',#20,#30)"));
        assert!(step.contains("STRAIGHTNESS_TOLERANCE('str_tol','straightness',#21,#31)"));
        assert!(step.contains("CIRCULARITY_TOLERANCE('cir_tol','circularity',#22,#32)"));
        assert!(step.contains("CYLINDRICITY_TOLERANCE('cyl_tol','cylindricity',#23,#33)"));
    }

    #[test]
    fn writes_runout_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            runout_tolerances: vec![
                StepRunoutTolerance {
                    entity_id: 0,
                    name: Some("cr_tol".into()),
                    description: Some("circular runout".into()),
                    value_entity_id: Some(20),
                    shape_aspect_id: Some(30),
                    datum_system_id: Some(40),
                    runout_type: RunoutToleranceType::CircularRunout,
                },
                StepRunoutTolerance {
                    entity_id: 0,
                    name: Some("tr_tol".into()),
                    description: Some("total runout".into()),
                    value_entity_id: Some(21),
                    shape_aspect_id: Some(31),
                    datum_system_id: Some(41),
                    runout_type: RunoutToleranceType::TotalRunout,
                },
            ],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("CIRCULAR_RUNOUT_TOLERANCE('cr_tol','circular runout',#20,#30,#40)"));
        assert!(step.contains("TOTAL_RUNOUT_TOLERANCE('tr_tol','total runout',#21,#31,#41)"));
    }

    #[test]
    fn writes_profile_tolerances() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            profile_tolerances: vec![
                StepProfileTolerance {
                    entity_id: 0,
                    name: Some("lin_tol".into()),
                    description: Some("profile of a line".into()),
                    value_entity_id: Some(20),
                    shape_aspect_id: Some(30),
                    datum_system_id: None,
                    profile_type: ProfileToleranceType::ProfileOfALine,
                },
                StepProfileTolerance {
                    entity_id: 0,
                    name: Some("surf_tol".into()),
                    description: Some("profile of a surface".into()),
                    value_entity_id: Some(21),
                    shape_aspect_id: Some(31),
                    datum_system_id: Some(41),
                    profile_type: ProfileToleranceType::ProfileOfASurface,
                },
            ],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("LINE_PROFILE_TOLERANCE('lin_tol','profile of a line',#20,#30,$)"));
        assert!(step.contains("SURFACE_PROFILE_TOLERANCE('surf_tol','profile of a surface',#21,#31,#41)"));
    }

    #[test]
    fn writes_datum_reference_frames() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            datum_reference_frames: vec![StepDatumReferenceFrame {
                entity_id: 0,
                name: Some("DRF1".into()),
                description: Some("primary datum reference frame".into()),
                datum_system_ids: vec![50, 51, 52],
            }],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("DATUM_REFERENCE_FRAME('DRF1','primary datum reference frame',(#50,#51,#52))"));
    }

    #[test]
    fn writes_datum_targets() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            datum_targets: vec![
                StepDatumTarget {
                    entity_id: 0,
                    name: Some("A1".into()),
                    description: Some("datum target point".into()),
                    target_identifier: Some("A1".into()),
                    datum_id: Some(80),
                    target_type: DatumTargetType::Point,
                    shape_aspect_id: Some(100),
                },
                StepDatumTarget {
                    entity_id: 0,
                    name: Some("B1".into()),
                    description: Some("datum target line".into()),
                    target_identifier: Some("B1".into()),
                    datum_id: Some(81),
                    target_type: DatumTargetType::Line,
                    shape_aspect_id: Some(101),
                },
                StepDatumTarget {
                    entity_id: 0,
                    name: Some("C1".into()),
                    description: Some("datum target area".into()),
                    target_identifier: Some("C1".into()),
                    datum_id: Some(82),
                    target_type: DatumTargetType::AreaCircle,
                    shape_aspect_id: Some(102),
                },
            ],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("DATUM_TARGET_POINT('A1','datum target point','A1',#80,#100)"));
        assert!(step.contains("DATUM_TARGET_LINE('B1','datum target line','B1',#81,#101)"));
        assert!(step.contains("DATUM_TARGET_CIRCLE('C1','datum target area','C1',#82,#102)"));
    }

    #[test]
    fn writes_enhanced_tolerance_zone_definitions() {
        let brep = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0)
            .expect("test box should be valid");
        let metadata = StepAp242Metadata {
            tolerance_zone_definitions_enhanced: vec![StepToleranceZoneDefinitionEnhanced {
                entity_id: 0,
                name: Some("cylindrical".into()),
                description: Some("symmetric".into()),
                tolerance_zone_id: Some(90),
                shape_aspect_id: Some(110),
                zone_shape: ToleranceZoneShape::Cylindrical,
                zone_position: ToleranceZonePosition::Symmetric,
                defining_shape_aspect_id: Some(120),
            }],
            ..Default::default()
        };
        let step = StepWriter::write_string_with_ap242_metadata(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &[],
            &metadata,
            StepProtocol::Ap242,
        );
        assert!(step.contains("TOLERANCE_ZONE_DEFINITION('cylindrical','symmetric',#90,#110,#120)"));
    }
}
