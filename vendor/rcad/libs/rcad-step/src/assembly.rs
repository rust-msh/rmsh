//! STEP assembly writer and reader.
//!
//! Writes a multi-component assembly where each component is a separate BRep
//! with an optional full affine transform and name.  Produces a STEP file with a
//! full NEXT_ASSEMBLY_USAGE_OCCURRENCE hierarchy so that importers (FreeCAD,
//! OCCT, etc.) can reconstruct the tree.
//!
//! # STEP assembly structure
//!
//! ```text
//! PRODUCT (assembly root)
//!   └─ PRODUCT_DEFINITION (assembly)
//!        └─ NEXT_ASSEMBLY_USAGE_OCCURRENCE → PRODUCT (component i)
//!                                               └─ PRODUCT_DEFINITION (component i)
//!                                                    └─ shape representation (geometry)
//! ```
//!
//! **Transform strategy**: transforms are baked into vertex/geometry coordinates
//! (`BRep::apply_transform`) rather than emitting `ITEM_DEFINED_TRANSFORMATION`
//! entities. This maximises compatibility with STEP readers that do not support
//! AP214 transformation entities.

use std::collections::{BTreeMap, HashMap, HashSet};

use glam::{DAffine3, DVec3};
use rcad_algorithms::{HealingOptions, HealingReport, analyze_and_heal};
use rcad_kernel::BRep;
use rcad_kernel::appearance::StepColor;
use serde::Serialize;

use crate::StepError;
use crate::writer::{ExportSelection, StepWriter};

// ─────────────────────────────────────────────────────────────────────────────
// AssemblyComponent
// ─────────────────────────────────────────────────────────────────────────────

/// A single component in an assembly.
///
/// The [`transform`][AssemblyComponent::transform] field is a full affine
/// transform (translation, rotation, uniform or non-uniform scale).  When
/// writing to STEP the transform is **baked** into vertex coordinates via
/// [`BRep::apply_transform`].
#[derive(Clone)]
pub struct AssemblyComponent {
    /// Human-readable part name.
    pub name: String,
    /// The geometry.
    pub brep: BRep,
    /// Full affine transform for this component (default: identity).
    pub transform: DAffine3,
    /// Optional RGB color for this component's faces.
    pub color: Option<rcad_kernel::appearance::Color>,
}

impl AssemblyComponent {
    /// Create a new component with an identity transform.
    pub fn new(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            brep,
            transform: DAffine3::IDENTITY,
            color: None,
        }
    }

    /// Set a full affine transform (replaces any previously set transform).
    pub fn with_transform(mut self, transform: DAffine3) -> Self {
        self.transform = transform;
        self
    }

    /// Convenience: set a pure translation transform.
    ///
    /// Equivalent to `with_transform(DAffine3::from_translation(t))`.
    pub fn with_translation(mut self, t: DVec3) -> Self {
        self.transform = DAffine3::from_translation(t);
        self
    }

    /// Set a color for this component.
    pub fn with_color(mut self, color: rcad_kernel::appearance::Color) -> Self {
        self.color = Some(color);
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AssemblyNode — tree-structured assembly
// ─────────────────────────────────────────────────────────────────────────────

/// A node in a hierarchical assembly tree.
///
/// - **Leaf** nodes carry a `BRep` geometry (and optional transform/color).
/// - **Branch** nodes carry child nodes and no geometry of their own.
///
/// Use [`write_assembly_tree`] / [`read_assembly_tree`] to round-trip the tree
/// through STEP.
#[derive(Clone)]
pub struct AssemblyNode {
    /// Human-readable name for this node.
    pub name: String,
    /// Geometry for leaf nodes; `None` for branch (sub-assembly) nodes.
    pub brep: Option<BRep>,
    /// Full affine transform applied before writing (baked into coordinates).
    pub transform: DAffine3,
    /// Optional color for leaf nodes.
    pub color: Option<rcad_kernel::appearance::Color>,
    /// Child nodes (empty for leaf nodes).
    pub children: Vec<AssemblyNode>,
}

impl AssemblyNode {
    /// Create a leaf node (part with geometry).
    pub fn leaf(name: impl Into<String>, brep: BRep) -> Self {
        Self {
            name: name.into(),
            brep: Some(brep),
            transform: DAffine3::IDENTITY,
            color: None,
            children: Vec::new(),
        }
    }

    /// Create a branch node (sub-assembly with children, no geometry).
    pub fn branch(name: impl Into<String>, children: Vec<AssemblyNode>) -> Self {
        Self {
            name: name.into(),
            brep: None,
            transform: DAffine3::IDENTITY,
            color: None,
            children,
        }
    }

    /// Set a transform (baked into coordinates on write).
    pub fn with_transform(mut self, t: DAffine3) -> Self {
        self.transform = t;
        self
    }

    /// Set a translation transform.
    pub fn with_translation(mut self, t: DVec3) -> Self {
        self.transform = DAffine3::from_translation(t);
        self
    }

    /// Set a color.
    pub fn with_color(mut self, color: rcad_kernel::appearance::Color) -> Self {
        self.color = Some(color);
        self
    }
}

/// Write a hierarchical assembly tree to a STEP string.
///
/// The tree is flattened into STEP `PRODUCT` / `PRODUCT_DEFINITION` /
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` entities.  Transforms are baked into
/// vertex coordinates.
pub fn write_assembly_tree(root_name: &str, root: &AssemblyNode) -> String {
    // Flatten the tree into a list of (parent_name, child_component) pairs
    // by recursively collecting all NAUO relationships.
    let mut records: Vec<String> = Vec::new();
    let mut next_id: u64 = 1;

    macro_rules! push {
        ($body:expr) => {{
            let id = next_id;
            next_id += 1;
            records.push(format!("#{}={};", id, $body));
            id
        }};
    }

    // Shared context entities (same as write_assembly_step)
    let app_ctx = push!("APPLICATION_CONTEXT('automotive_design')".to_string());
    push!(format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{})",
        app_ctx
    ));
    let prod_ctx = push!(format!(
        "PRODUCT_CONTEXT('part definition',#{},'mechanical')",
        app_ctx
    ));
    let def_ctx = push!(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{},'design')",
        app_ctx
    ));
    let len_unit = push!("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string());
    let rad_unit = push!("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
    let meas = push!(format!(
        "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})",
        rad_unit
    ));
    let dim_exp = push!("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)".to_string());
    let deg_unit = push!(format!(
        "( CONVERSION_BASED_UNIT('DEGREE',#{}) NAMED_UNIT(#{}) PLANE_ANGLE_UNIT() )",
        meas, dim_exp
    ));
    let sol_unit =
        push!("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string());
    let uncert = push!(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#{},'distance_accuracy_value','confusion accuracy')",
        len_unit
    ));
    let geom_ctx = push!(format!(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') )",
        uncert, len_unit, deg_unit, sol_unit
    ));

    // Recursively emit nodes; returns the PRODUCT_DEFINITION id for the node.
    fn emit_node(
        node: &AssemblyNode,
        records: &mut Vec<String>,
        next_id: &mut u64,
        prod_ctx: u64,
        def_ctx: u64,
        geom_ctx: u64,
        nauo_seq: &mut usize,
    ) -> u64 {
        // Apply transform to a cloned BRep if needed.
        let baked_brep: Option<BRep> = node.brep.as_ref().map(|b| {
            let mut b2 = b.clone();
            if node.transform != DAffine3::IDENTITY {
                b2.apply_transform(node.transform);
            }
            b2
        });

        // Emit geometry for leaf nodes.
        let shape_rep_id: Option<u64> = if let Some(brep) = &baked_brep {
            let colors = node.color.map(|c| StepColor::new().with_solid_color(c));
            let comp_step = if let Some(sc) = &colors {
                StepWriter::write_string_colored(brep, sc)
            } else {
                StepWriter::write_string(
                    brep,
                    ExportSelection { selected_faces: &[], selected_edges: &[] },
                )
            };
            let comp_records = extract_data_records(&comp_step);
            let id_offset = *next_id - 1;
            let comp_max_id = comp_records.iter().map(|(id, _)| *id).max().unwrap_or(1);
            let sr_id = comp_max_id + id_offset;
            for (orig_id, body) in &comp_records {
                let new_id = orig_id + id_offset;
                let renumbered_body = renumber_refs(body, id_offset);
                records.push(format!("#{}={};", new_id, renumbered_body));
            }
            *next_id = sr_id + 1;
            Some(sr_id)
        } else {
            None
        };

        // Emit PRODUCT / PRODUCT_DEFINITION for this node.
        let prod = push_record(records, next_id, format!(
            "PRODUCT('{}','{}','',( #{} ))",
            node.name, node.name, prod_ctx
        ));
        let formation = push_record(records, next_id, format!("PRODUCT_DEFINITION_FORMATION('','',#{})", prod));
        let pd = push_record(records, next_id, format!(
            "PRODUCT_DEFINITION('','',#{},#{})",
            formation, def_ctx
        ));
        let pds = push_record(records, next_id, format!("PRODUCT_DEFINITION_SHAPE('','',#{})", pd));

        if let Some(sr) = shape_rep_id {
            push_record(records, next_id, format!("SHAPE_DEFINITION_REPRESENTATION(#{},#{})", pds, sr));
        } else {
            // Branch node: emit an empty shape representation.
            let empty_sr = push_record(records, next_id, format!(
                "SHAPE_REPRESENTATION('{}',( #{} ),#{})",
                node.name, geom_ctx, geom_ctx
            ));
            push_record(records, next_id, format!("SHAPE_DEFINITION_REPRESENTATION(#{},#{})", pds, empty_sr));
        }

        // Recursively emit children and link via NAUO.
        for child in &node.children {
            let child_pd = emit_node(
                child, records, next_id, prod_ctx, def_ctx, geom_ctx, nauo_seq,
            );
            *nauo_seq += 1;
            let nauo = push_record(records, next_id, format!(
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','',#{},#{},$)",
                nauo_seq, child.name, pd, child_pd
            ));
            push_record(records, next_id, format!(
                "PRODUCT_DEFINITION_SHAPE('Acme','occurrence shape',#{})",
                nauo
            ));
        }

        pd
    }

    let mut nauo_seq = 0usize;

    // Emit the root node (which may itself be a branch).
    emit_node(
        root,
        &mut records,
        &mut next_id,
        prod_ctx,
        def_ctx,
        geom_ctx,
        &mut nauo_seq,
    );

    // Build output
    use std::fmt::Write as FmtWrite;
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    let _ = writeln!(
        out,
        "FILE_DESCRIPTION(('RCAD assembly tree: {}'),'2;1');",
        root_name
    );
    out.push_str(
        "FILE_NAME('rcad_assembly.step','2026-04-11T00:00:00',(''),(''),'RCAD','RCAD','');\n",
    );
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");
    for record in &records {
        out.push_str(record);
        out.push('\n');
    }
    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    out
}

/// Parse a STEP string into a hierarchical [`AssemblyNode`] tree.
///
/// The tree mirrors the `NEXT_ASSEMBLY_USAGE_OCCURRENCE` hierarchy in the file.
/// Leaf nodes carry isolated `BRep` geometry (via the same BFS approach as
/// [`read_assembly`]).  Branch nodes have `brep = None` and carry children.
///
/// Falls back to a single leaf node for plain single-part STEP files.
pub fn read_assembly_tree(step: &str) -> Result<AssemblyNode, StepError> {
    let entity_map = parse_entity_map(step);
    let reverse_map = build_reverse_map(&entity_map);

    // Build parent→children map from NAUO entities.
    // NAUO: ('seq','name','desc',#parent_pd,#child_pd,$)
    let mut children_of: HashMap<u64, Vec<(u64, String)>> = HashMap::new();
    let mut all_child_pds: HashSet<u64> = HashSet::new();

    for (_id, body) in &entity_map {
        if let Some(rest) = strip_entity_name(body, "NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            if let Some(args) = parse_args(rest) {
                let parent_pd = parse_ref(args.get(3).copied().unwrap_or(""));
                let child_pd = parse_ref(args.get(4).copied().unwrap_or(""));
                let rel_name = unquote(args.get(1).copied().unwrap_or(""));
                if parent_pd > 0 && child_pd > 0 {
                    children_of
                        .entry(parent_pd)
                        .or_default()
                        .push((child_pd, rel_name));
                    all_child_pds.insert(child_pd);
                }
            }
        }
    }

    // Find root PD(s): PDs that appear in NAUOs (as parent or child) but are
    // never a NAUO child.  This excludes PDs embedded in inlined geometry STEP
    // strings, which are not connected to any NAUO.
    let nauo_pds: HashSet<u64> = children_of
        .iter()
        .flat_map(|(&parent, children)| {
            std::iter::once(parent).chain(children.iter().map(|(child, _)| *child))
        })
        .collect();

    let root_pds: Vec<u64> = nauo_pds
        .iter()
        .copied()
        .filter(|id| !all_child_pds.contains(id))
        .collect();

    if root_pds.is_empty() || children_of.is_empty() {
        // Plain single-part STEP.
        let brep = crate::StepReader::parse_string(step)?;
        let name = find_root_product_name(&entity_map).unwrap_or_else(|| "part".to_string());
        return Ok(AssemblyNode::leaf(name, brep));
    }

    // Parse the whole STEP once as fallback for nodes whose geometry can't be isolated.
    let merged_brep = crate::StepReader::parse_string(step)?;

    // Recursively build the tree.
    fn build_node(
        pd_id: u64,
        entity_map: &HashMap<u64, String>,
        reverse_map: &HashMap<u64, Vec<u64>>,
        children_of: &HashMap<u64, Vec<(u64, String)>>,
        merged_brep: &BRep,
        step: &str,
    ) -> AssemblyNode {
        let name = resolve_product_name(pd_id, entity_map)
            .unwrap_or_else(|| format!("node_{}", pd_id));

        let child_entries = children_of.get(&pd_id);

        if let Some(entries) = child_entries {
            // Branch node: recurse into children.
            let children: Vec<AssemblyNode> = entries
                .iter()
                .map(|(child_pd, _)| {
                    build_node(
                        *child_pd,
                        entity_map,
                        reverse_map,
                        children_of,
                        merged_brep,
                        step,
                    )
                })
                .collect();
            AssemblyNode::branch(name, children)
        } else {
            // Leaf node: isolate geometry.
            let brep =
                if let Some(sr_id) = find_shape_rep_for_pd(pd_id, entity_map, reverse_map) {
                    let reachable = collect_reachable(sr_id, entity_map);
                    let comp_step = build_component_step(entity_map, &reachable);
                    crate::StepReader::parse_string(&comp_step)
                        .unwrap_or_else(|_| merged_brep.clone())
                } else {
                    merged_brep.clone()
                };
            AssemblyNode::leaf(name, brep)
        }
    }

    // If there are multiple roots (unusual), wrap them in a synthetic root.
    if root_pds.len() == 1 {
        Ok(build_node(
            root_pds[0],
            &entity_map,
            &reverse_map,
            &children_of,
            &merged_brep,
            step,
        ))
    } else {
        let children: Vec<AssemblyNode> = root_pds
            .iter()
            .map(|&pd| {
                build_node(
                    pd,
                    &entity_map,
                    &reverse_map,
                    &children_of,
                    &merged_brep,
                    step,
                )
            })
            .collect();
        Ok(AssemblyNode::branch("root".to_string(), children))
    }
}



/// Write a multi-component STEP assembly.
///
/// Each component is written as a separate `PRODUCT` + `PRODUCT_DEFINITION`
/// with its own geometry representation, linked into the root assembly via
/// `NEXT_ASSEMBLY_USAGE_OCCURRENCE`.
///
/// The [`AssemblyComponent::transform`] matrix is applied to vertex and geometry
/// coordinates before writing (baked into the STEP geometry, not stored as a
/// STEP transform entity).
pub fn write_assembly(assembly_name: &str, components: &[AssemblyComponent]) -> String {
    // Apply full DAffine3 transform into new BReps and collect colors.
    let prepared: Vec<(String, BRep, Option<rcad_kernel::appearance::Color>)> = components
        .iter()
        .map(|c| {
            let mut b = c.brep.clone();
            if c.transform != DAffine3::IDENTITY {
                b.apply_transform(c.transform);
            }
            (c.name.clone(), b, c.color)
        })
        .collect();

    write_assembly_step(assembly_name, &prepared)
}

fn push_record(records: &mut Vec<String>, next_id: &mut u64, body: String) -> u64 {
    let id = *next_id;
    *next_id += 1;
    records.push(format!("#{}={};", id, body));
    id
}

fn write_assembly_step(
    assembly_name: &str,
    components: &[(String, BRep, Option<rcad_kernel::appearance::Color>)],
) -> String {
    use std::fmt::Write as FmtWrite;

    let mut records: Vec<String> = Vec::new();
    let mut next_id: u64 = 1;

    macro_rules! push {
        ($body:expr) => {
            push_record(&mut records, &mut next_id, $body)
        };
    }

    // Shared context entities
    let app_ctx = push!("APPLICATION_CONTEXT('automotive_design')".to_string());
    push!(format!(
        "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#{})",
        app_ctx
    ));
    let prod_ctx = push!(format!(
        "PRODUCT_CONTEXT('part definition',#{},'mechanical')",
        app_ctx
    ));
    let def_ctx = push!(format!(
        "PRODUCT_DEFINITION_CONTEXT('part definition',#{},'design')",
        app_ctx
    ));

    // Measurement context
    let len_unit = push!("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string());
    let rad_unit = push!("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )".to_string());
    let meas = push!(format!(
        "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{})",
        rad_unit
    ));
    let dim_exp = push!("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)".to_string());
    let deg_unit = push!(format!(
        "( CONVERSION_BASED_UNIT('DEGREE',#{}) NAMED_UNIT(#{}) PLANE_ANGLE_UNIT() )",
        meas, dim_exp
    ));
    let sol_unit = push!("( NAMED_UNIT(*) SOLID_ANGLE_UNIT() SI_UNIT($,.STERADIAN.) )".to_string());
    let uncert = push!(format!(
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-6),#{},'distance_accuracy_value','confusion accuracy')",
        len_unit
    ));
    let geom_ctx = push!(format!(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{})) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{},#{},#{})) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') )",
        uncert, len_unit, deg_unit, sol_unit
    ));

    // Assembly root product
    let asm_product = push!(format!(
        "PRODUCT('{}','{}','',( #{} ))",
        assembly_name, assembly_name, prod_ctx
    ));
    let asm_formation = push!(format!(
        "PRODUCT_DEFINITION_FORMATION('','',#{})",
        asm_product
    ));
    let asm_definition = push!(format!(
        "PRODUCT_DEFINITION('','',#{},#{})",
        asm_formation, def_ctx
    ));
    let asm_shape = push!(format!(
        "PRODUCT_DEFINITION_SHAPE('','',#{})",
        asm_definition
    ));
    let asm_rep = push!(format!(
        "SHAPE_REPRESENTATION('{}',(#{}),#{})",
        assembly_name, geom_ctx, geom_ctx
    ));
    push!(format!(
        "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
        asm_shape, asm_rep
    ));

    // Component products + NAUO links
    for (i, (comp_name, brep, comp_color)) in components.iter().enumerate() {
        // Write the component geometry as a standalone STEP string,
        // then inline its DATA section records with offset IDs.
        let colors = comp_color.map(|c| StepColor::new().with_solid_color(c));
        let comp_step = if let Some(sc) = &colors {
            StepWriter::write_string_colored(brep, sc)
        } else {
            StepWriter::write_string(
                brep,
                ExportSelection {
                    selected_faces: &[],
                    selected_edges: &[],
                },
            )
        };

        // Extract DATA records from component STEP string and re-number them
        let comp_records = extract_data_records(&comp_step);
        let id_offset = next_id - 1;

        // Re-number and collect component records
        let renumbered: Vec<String> = comp_records
            .iter()
            .map(|(orig_id, body)| {
                let new_id = orig_id + id_offset;
                let renumbered_body = renumber_refs(body, id_offset);
                format!("#{}={};", new_id, renumbered_body)
            })
            .collect();

        // Find the highest ID in the component (= shape_representation or similar)
        let comp_max_id = comp_records.iter().map(|(id, _)| *id).max().unwrap_or(1);
        let shape_rep_id = comp_max_id + id_offset; // last record is typically the shape repr

        // Advance next_id past all component records
        for record in &renumbered {
            records.push(record.clone());
        }
        next_id = shape_rep_id + 1;

        // Component product/definition
        let comp_product = push!(format!(
            "PRODUCT('{}','{}','',( #{} ))",
            comp_name, comp_name, prod_ctx
        ));
        let comp_formation = push!(format!(
            "PRODUCT_DEFINITION_FORMATION('','',#{})",
            comp_product
        ));
        let comp_definition = push!(format!(
            "PRODUCT_DEFINITION('','',#{},#{})",
            comp_formation, def_ctx
        ));
        let comp_pds = push!(format!(
            "PRODUCT_DEFINITION_SHAPE('','',#{})",
            comp_definition
        ));
        push!(format!(
            "SHAPE_DEFINITION_REPRESENTATION(#{},#{})",
            comp_pds, shape_rep_id
        ));

        // NEXT_ASSEMBLY_USAGE_OCCURRENCE: link component to assembly
        let nauo = push!(format!(
            "NEXT_ASSEMBLY_USAGE_OCCURRENCE('{}','{}','',#{},#{},$)",
            i + 1,
            comp_name,
            asm_definition,
            comp_definition
        ));
        // PRODUCT_DEFINITION_SHAPE for the occurrence
        push!(format!(
            "PRODUCT_DEFINITION_SHAPE('Acme','occurrence shape',#{})",
            nauo
        ));
    }

    // Build output
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    let _ = writeln!(
        out,
        "FILE_DESCRIPTION(('RCAD assembly: {}'),'2;1');",
        assembly_name
    );
    out.push_str(
        "FILE_NAME('rcad_assembly.step','2026-04-11T00:00:00',(''),(''),'RCAD','RCAD','');\n",
    );
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");
    for record in &records {
        out.push_str(record);
        out.push('\n');
    }
    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Read
// ─────────────────────────────────────────────────────────────────────────────

/// Parse a STEP string that may contain a multi-component assembly and return
/// a flat list of [`AssemblyComponent`]s.
///
/// # Algorithm
///
/// 1. Scan the DATA section for `NEXT_ASSEMBLY_USAGE_OCCURRENCE` (NAUO) entities
///    to build a parent→child component map.
/// 2. For each NAUO child `PRODUCT_DEFINITION`, trace the chain:
///    `PRODUCT_DEFINITION` → `PRODUCT_DEFINITION_SHAPE` →
///    `SHAPE_DEFINITION_REPRESENTATION` → `SHAPE_REPRESENTATION`.
/// 3. BFS-collect all entity IDs reachable from that shape representation
///    (geometry, surfaces, curves, vertices, units, context, …).
/// 4. Build a self-contained sub-STEP string for each component and parse it
///    with [`crate::StepReader`] to obtain an isolated [`BRep`].
///
/// ## Fallback
///
/// If the file has no NAUO links (plain single-part STEP), the function returns
/// a single-element list containing the full parsed BRep with `transform =
/// IDENTITY`.
///
/// If a component's shape representation chain cannot be resolved (e.g. a
/// third-party file with non-standard structure), that component falls back to
/// the full merged BRep.
pub fn read_assembly(step: &str) -> Result<Vec<AssemblyComponent>, StepError> {
    let entity_map = parse_entity_map(step);
    let reverse_map = build_reverse_map(&entity_map);

    // Collect NAUO children: child PRODUCT_DEFINITION id + relation name
    let mut nauo_children: Vec<(u64, String)> = Vec::new();
    let mut has_nauo = false;

    for (_id, body) in &entity_map {
        if let Some(rest) = strip_entity_name(body, "NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            has_nauo = true;
            if let Some(args) = parse_args(rest) {
                let child_pd = parse_ref(args.get(4).copied().unwrap_or(""));
                let relation_name = unquote(args.get(1).copied().unwrap_or(""));
                if child_pd > 0 {
                    nauo_children.push((child_pd, relation_name));
                }
            }
        }
    }

    if !has_nauo || nauo_children.is_empty() {
        let brep = crate::StepReader::parse_string(step)?;
        let name = find_root_product_name(&entity_map).unwrap_or_else(|| "part".to_string());
        return Ok(vec![AssemblyComponent::new(name, brep)]);
    }

    // Parse the whole STEP once as a fallback for components whose geometry
    // cannot be isolated (e.g. third-party files with unusual structure).
    let merged_brep = crate::StepReader::parse_string(step)?;

    let mut components = Vec::new();
    for (pd_id, nauo_name) in &nauo_children {
        let name = resolve_product_name(*pd_id, &entity_map)
            .unwrap_or_else(|| nauo_name.clone());

        let brep = if let Some(sr_id) =
            find_shape_rep_for_pd(*pd_id, &entity_map, &reverse_map)
        {
            // BFS-collect all entities reachable from this component's shape
            // representation, then build a self-contained sub-STEP string.
            let reachable = collect_reachable(sr_id, &entity_map);
            let comp_step = build_component_step(&entity_map, &reachable);
            crate::StepReader::parse_string(&comp_step)
                .unwrap_or_else(|_| merged_brep.clone())
        } else {
            merged_brep.clone()
        };

        components.push(AssemblyComponent::new(name, brep));
    }

    Ok(components)
}

/// Parse STEP assembly and run healing pipeline on each component BRep.
///
/// Returns healed components plus one healing report per component in order.
pub fn read_assembly_with_healing(
    step: &str,
    options: HealingOptions,
) -> Result<(Vec<AssemblyComponent>, Vec<HealingReport>), StepError> {
    let mut components = read_assembly(step)?;
    let mut reports = Vec::with_capacity(components.len());

    for component in &mut components {
        let (healed, report) = analyze_and_heal(&component.brep, options);
        component.brep = healed;
        reports.push(report);
    }

    Ok((components, reports))
}

/// Stable JSON diagnostics payload for assembly import healing.
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyImportHealingJsonV1 {
    pub schema: &'static str,
    pub component_count: usize,
    pub clean_components: usize,
    pub failed_components: usize,
    pub issue_histogram: Vec<(String, usize)>,
}

/// Parse STEP assembly, run healing, and export a stable JSON diagnostics report.
pub fn read_assembly_with_healing_report_json(
    step: &str,
    options: HealingOptions,
) -> Result<(Vec<AssemblyComponent>, Vec<HealingReport>, String), StepError> {
    let (components, reports) = read_assembly_with_healing(step, options)?;

    let clean_components = reports.iter().filter(|r| r.final_result.is_valid()).count();
    let failed_components = reports.len().saturating_sub(clean_components);

    let mut issue_map: BTreeMap<String, usize> = BTreeMap::new();
    for report in &reports {
        for issue in &report.final_result.issues {
            *issue_map.entry(issue.to_string()).or_insert(0) += 1;
        }
    }

    let payload = AssemblyImportHealingJsonV1 {
        schema: "step.assembly.import.healing.v1",
        component_count: components.len(),
        clean_components,
        failed_components,
        issue_histogram: issue_map.into_iter().collect(),
    };

    let json = serde_json::to_string_pretty(&payload)
        .map_err(|e| StepError::InvalidFormat(format!("healing report JSON serialize failed: {e}")))?;

    Ok((components, reports, json))
}

/// Healing report for a tree node with geometry.
#[derive(Debug, Clone)]
pub struct AssemblyNodeHealingReport {
    /// Child-index path from tree root to this node.
    pub path: Vec<usize>,
    /// Node name snapshot for easier diagnostics.
    pub name: String,
    /// Healing analysis/repair report for this node's geometry.
    pub report: HealingReport,
}

/// Parse hierarchical STEP assembly and run healing on each geometric node.
///
/// Returns healed tree plus one report for each node that carries geometry
/// (`node.brep.is_some()`), in depth-first traversal order.
pub fn read_assembly_tree_with_healing(
    step: &str,
    options: HealingOptions,
) -> Result<(AssemblyNode, Vec<AssemblyNodeHealingReport>), StepError> {
    let mut root = read_assembly_tree(step)?;
    let mut reports = Vec::new();

    fn heal_node(
        node: &mut AssemblyNode,
        path: &mut Vec<usize>,
        options: HealingOptions,
        reports: &mut Vec<AssemblyNodeHealingReport>,
    ) {
        if let Some(brep) = node.brep.take() {
            let (healed, report) = analyze_and_heal(&brep, options);
            node.brep = Some(healed);
            reports.push(AssemblyNodeHealingReport {
                path: path.clone(),
                name: node.name.clone(),
                report,
            });
        }

        for (idx, child) in node.children.iter_mut().enumerate() {
            path.push(idx);
            heal_node(child, path, options, reports);
            path.pop();
        }
    }

    let mut path = Vec::new();
    heal_node(&mut root, &mut path, options, &mut reports);

    Ok((root, reports))
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal parsing helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Extract all `#N` reference IDs from a STEP entity body string.
fn extract_refs(body: &str) -> Vec<u64> {
    let mut refs = Vec::new();
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if i > start {
                if let Ok(n) = body[start..i].parse::<u64>() {
                    refs.push(n);
                }
            }
        } else {
            i += 1;
        }
    }
    refs
}

/// Build a reverse-reference map: referenced_id → list of entity IDs that
/// contain a `#referenced_id` in their body.
fn build_reverse_map(map: &HashMap<u64, String>) -> HashMap<u64, Vec<u64>> {
    let mut reverse: HashMap<u64, Vec<u64>> = HashMap::new();
    for (&id, body) in map {
        for ref_id in extract_refs(body) {
            reverse.entry(ref_id).or_default().push(id);
        }
    }
    reverse
}

/// BFS from `start_id`, following all `#N` references, and return the set of
/// all reachable entity IDs (including `start_id` itself).
fn collect_reachable(start_id: u64, map: &HashMap<u64, String>) -> HashSet<u64> {
    let mut visited: HashSet<u64> = HashSet::new();
    let mut queue = vec![start_id];
    while let Some(id) = queue.pop() {
        if !visited.insert(id) {
            continue;
        }
        if let Some(body) = map.get(&id) {
            for ref_id in extract_refs(body) {
                if !visited.contains(&ref_id) {
                    queue.push(ref_id);
                }
            }
        }
    }
    visited
}

/// Trace `PRODUCT_DEFINITION` → `PRODUCT_DEFINITION_SHAPE` →
/// `SHAPE_DEFINITION_REPRESENTATION` → shape representation ID.
///
/// Returns the shape representation entity ID, or `None` if the chain cannot
/// be resolved (e.g. the file uses a non-standard structure).
fn find_shape_rep_for_pd(
    pd_id: u64,
    map: &HashMap<u64, String>,
    reverse_map: &HashMap<u64, Vec<u64>>,
) -> Option<u64> {
    // Find PRODUCT_DEFINITION_SHAPE entities that reference pd_id
    let referencing = reverse_map.get(&pd_id)?;
    for &pds_id in referencing {
        let body = map.get(&pds_id)?;
        if !body.starts_with("PRODUCT_DEFINITION_SHAPE(") {
            continue;
        }
        // Find SHAPE_DEFINITION_REPRESENTATION entities that reference pds_id
        if let Some(referencing_pds) = reverse_map.get(&pds_id) {
            for &sdr_id in referencing_pds {
                if let Some(sdr_body) = map.get(&sdr_id) {
                    if let Some(args_str) =
                        sdr_body.strip_prefix("SHAPE_DEFINITION_REPRESENTATION(")
                    {
                        if let Some(args) = parse_args(args_str) {
                            let sr_id = parse_ref(args.get(1).copied().unwrap_or(""));
                            if sr_id > 0 {
                                return Some(sr_id);
                            }
                        }
                    }
                }
            }
        }
    }
    None
}

/// Build a self-contained STEP string containing only the entities in
/// `reachable`, using their original IDs from `map`.
fn build_component_step(map: &HashMap<u64, String>, reachable: &HashSet<u64>) -> String {
    let mut out = String::new();
    out.push_str("ISO-10303-21;\n");
    out.push_str("HEADER;\n");
    out.push_str("FILE_DESCRIPTION(('RCAD component'),'2;1');\n");
    out.push_str(
        "FILE_NAME('component.step','2026-04-11T00:00:00',(''),(''),'RCAD','RCAD','');\n",
    );
    out.push_str("FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\n");
    out.push_str("ENDSEC;\n");
    out.push_str("DATA;\n");

    let mut ids: Vec<u64> = reachable.iter().copied().collect();
    ids.sort_unstable();
    for id in ids {
        if let Some(body) = map.get(&id) {
            out.push_str(&format!("#{}={};\n", id, body));
        }
    }

    out.push_str("ENDSEC;\n");
    out.push_str("END-ISO-10303-21;\n");
    out
}

/// Build a map from entity ID → entity body (the part after `#id=` and before `;`).
fn parse_entity_map(step: &str) -> std::collections::HashMap<u64, String> {
    let mut map = std::collections::HashMap::new();
    let mut in_data = false;
    for line in step.lines() {
        let line = line.trim();
        if line == "DATA;" {
            in_data = true;
            continue;
        }
        if line == "ENDSEC;" {
            in_data = false;
            continue;
        }
        if !in_data {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('#') {
            if let Some(eq) = stripped.find('=') {
                let id_str = &stripped[..eq];
                let body = stripped[eq + 1..].trim_end_matches(';');
                if let Ok(id) = id_str.parse::<u64>() {
                    map.insert(id, body.to_string());
                }
            }
        }
    }
    map
}

/// If `body` starts with `ENTITY_NAME(`, return the argument string (after the
/// opening paren, before the matching closing paren).
fn strip_entity_name<'a>(body: &'a str, entity: &str) -> Option<&'a str> {
    let prefix = entity;
    if body.starts_with(prefix) {
        body[prefix.len()..].strip_prefix('(')
    } else if let Some(inner) = body.strip_prefix('(') {
        // compound entity: ( ENTITY_NAME() ... )
        // Find the sub-entity by scanning for the name inside compound parens
        let inner = inner.trim();
        if inner.starts_with(entity) {
            inner[entity.len()..].strip_prefix('(')
        } else {
            None
        }
    } else {
        None
    }
}

/// Split a STEP argument list (without outer parens) on `,`, respecting nested
/// parens and string literals.
fn parse_args(args_str: &str) -> Option<Vec<&str>> {
    let mut result = Vec::new();
    let mut depth = 0i32;
    let mut in_str = false;
    let mut start = 0;
    let bytes = args_str.as_bytes();
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'\'' if !in_str => in_str = true,
            b'\'' if in_str => in_str = false,
            b'(' if !in_str => depth += 1,
            b')' if !in_str => {
                if depth == 0 {
                    // closing paren of the whole arg list
                    result.push(args_str[start..i].trim());
                    return Some(result);
                }
                depth -= 1;
            }
            b',' if !in_str && depth == 0 => {
                result.push(args_str[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    if start < args_str.len() {
        result.push(args_str[start..].trim());
    }
    Some(result)
}

/// Parse `#N` → N, returning 0 on failure.
fn parse_ref(s: &str) -> u64 {
    s.trim().strip_prefix('#').and_then(|n| n.trim_end_matches(|c: char| !c.is_ascii_digit()).parse().ok()).unwrap_or(0)
}

/// Strip surrounding single-quotes from a STEP string literal.
fn unquote(s: &str) -> String {
    let s = s.trim();
    if s.starts_with('\'') && s.ends_with('\'') && s.len() >= 2 {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

/// Resolve PRODUCT_DEFINITION id → product name via PD → PD_FORMATION → PRODUCT.
fn resolve_product_name(
    pd_id: u64,
    map: &std::collections::HashMap<u64, String>,
) -> Option<String> {
    // PD body: PRODUCT_DEFINITION('','',#formation,#ctx)
    let pd_body = map.get(&pd_id)?;
    let pd_args = parse_args(pd_body.strip_prefix("PRODUCT_DEFINITION(")?.strip_suffix(')')?.as_ref())?;
    // Actually strip_entity_name handles compound; simpler approach:
    let formation_id = pd_args.get(2).map(|s| parse_ref(s))?;
    if formation_id == 0 {
        return None;
    }

    // Formation body: PRODUCT_DEFINITION_FORMATION('','',#product)
    let form_body = map.get(&formation_id)?;
    let form_args = parse_args(
        form_body
            .strip_prefix("PRODUCT_DEFINITION_FORMATION(")?.strip_suffix(')')?,
    )?;
    let prod_id = form_args.get(2).map(|s| parse_ref(s))?;
    if prod_id == 0 {
        return None;
    }

    // Product body: PRODUCT('id','name','desc',(#ctx))
    let prod_body = map.get(&prod_id)?;
    let prod_args = parse_args(prod_body.strip_prefix("PRODUCT(")?.strip_suffix(')')?.as_ref())?;
    // name is the second field (index 1)
    prod_args.get(1).map(|s| unquote(s))
}

/// Find the top-level PRODUCT name in a single-part STEP file (first PRODUCT entity).
fn find_root_product_name(map: &std::collections::HashMap<u64, String>) -> Option<String> {
    // Return the first PRODUCT entity's second field (name).
    let mut ids: Vec<u64> = map.keys().copied().collect();
    ids.sort();
    for id in ids {
        let body = &map[&id];
        if body.starts_with("PRODUCT(") {
            if let Some(args) = parse_args(body.strip_prefix("PRODUCT(")?.strip_suffix(')')?) {
                return args.get(1).map(|s| unquote(s));
            }
        }
    }
    None
}

// ─────────────────────────────────────────────────────────────────────────────
// Low-level helpers (retained from previous implementation)
// ─────────────────────────────────────────────────────────────────────────────

/// Extract (id, body) pairs from a STEP DATA section string.
fn extract_data_records(step: &str) -> Vec<(u64, String)> {
    let mut in_data = false;
    let mut result = Vec::new();
    for line in step.lines() {
        let line = line.trim();
        if line == "DATA;" {
            in_data = true;
            continue;
        }
        if line == "ENDSEC;" {
            in_data = false;
            continue;
        }
        if !in_data {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('#')
            && let Some(eq) = stripped.find('=')
        {
            let id_str = &stripped[..eq];
            let body = stripped[eq + 1..].trim_end_matches(';');
            if let Ok(id) = id_str.parse::<u64>() {
                result.push((id, body.to_string()));
            }
        }
    }
    result
}

/// Replace all `#N` references in a STEP entity body with `#(N + offset)`.
fn renumber_refs(body: &str, offset: u64) -> String {
    let mut out = String::with_capacity(body.len() + 16);
    let bytes = body.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if i > start {
                let num: u64 = body[start..i].parse().unwrap_or(0);
                out.push('#');
                out.push_str(&(num + offset).to_string());
            } else {
                out.push('#');
            }
        } else {
            out.push(bytes[i] as char);
            i += 1;
        }
    }
    out
}
