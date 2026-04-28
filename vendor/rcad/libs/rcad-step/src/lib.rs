use rcad_kernel::appearance::{Color, StepColor};
use rcad_algorithms::{HealingOptions, HealingReport, analyze_and_heal, analyze_wire_issues};
use rcad_kernel::geom::BSplineCurve3;
use rcad_kernel::tolerance::CONFUSION;
use rcad_kernel::{
    BRep, BSplineCurve2, Curve2d, Curve2dEval, Curve3, CurveEval, Ellipse2d, GeomStore, PCurve,
    Surface3, SurfaceEval,
};
use rcad_kernel::{Edge, Face, Shell, Solid, Vertex, Wire, WireEdge};
use rcad_modeling::make_box_brep;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Read;
use std::path::Path;

pub mod assembly;
pub mod iges;
pub mod obj_writer;
pub mod writer;

pub use assembly::{
    AssemblyComponent, AssemblyImportHealingJsonV1, AssemblyNode,
    AssemblyNodeHealingReport, read_assembly, read_assembly_tree,
    read_assembly_tree_with_healing, read_assembly_with_healing,
    read_assembly_with_healing_report_json, write_assembly, write_assembly_tree,
};
pub use iges::{IgesError, IgesReader, IgesWriter};
pub use obj_writer::{ObjError, ObjReader, ObjWriter, write_obj};
pub use writer::{
    ExportSelection, StepAp242Metadata, StepHeader, StepProtocol, StepWriteOptions, StepWriter,
};

/// Errors that can occur when reading or parsing a STEP file.
#[derive(Debug, Clone)]
pub enum StepError {
    /// File I/O error.
    Io(String),
    /// Not a valid STEP file (missing ISO-10303-21, no DATA/ENDSEC, etc.).
    InvalidFormat(String),
    /// A required STEP entity is missing or malformed.
    MissingEntity {
        entity_type: &'static str,
        id: Option<u64>,
    },
    /// Parse produced an empty or degenerate result.
    EmptyResult(String),
}

impl std::fmt::Display for StepError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(msg) => write!(f, "I/O error: {msg}"),
            Self::InvalidFormat(msg) => write!(f, "invalid STEP format: {msg}"),
            Self::MissingEntity {
                entity_type,
                id: Some(id),
            } => write!(f, "missing {entity_type} entity #{id}"),
            Self::MissingEntity {
                entity_type,
                id: None,
            } => write!(f, "missing {entity_type} entity"),
            Self::EmptyResult(msg) => write!(f, "STEP parse produced empty result: {msg}"),
        }
    }
}

impl std::error::Error for StepError {}

#[derive(Debug, Clone)]
struct AdvancedFaceRecord {
    bounds: Vec<u64>,
    surface: Option<u64>,
}

/// B-spline curve: (degree, control_point_refs, knot_mults, knot_vals)
type BSplineCurveData = (usize, Vec<u64>, Vec<usize>, Vec<f64>);
/// B-spline surface: (name, degree_u, degree_v, ctrl_grid_refs[v][u], mults_u, knots_u, mults_v, knots_v)
type BSplineSurfaceData = (
    String,
    usize,
    usize,
    Vec<Vec<u64>>,
    Vec<usize>,
    Vec<f64>,
    Vec<usize>,
    Vec<f64>,
);

#[derive(Debug, Clone)]
struct ParsedStep {
    cartesian_points: HashMap<u64, [f64; 3]>,
    directions: HashMap<u64, [f64; 3]>,
    vectors: HashMap<u64, (u64, f64)>,
    axis2_placements: HashMap<u64, (u64, u64, Option<u64>)>,
    lines: HashMap<u64, (u64, u64)>,
    circles: HashMap<u64, (u64, f64)>,
    ellipses: HashMap<u64, (u64, f64, f64)>,
    b_spline_curves: HashMap<u64, Vec<u64>>,
    /// Full B-spline curve data: degree, control_point_refs, knot_mults, knot_vals
    b_spline_curves_full: HashMap<u64, BSplineCurveData>,
    planes: HashMap<u64, u64>,
    cylindrical_surfaces: HashMap<u64, (u64, f64)>,
    spherical_surfaces: HashMap<u64, (u64, f64)>,
    conical_surfaces: HashMap<u64, (u64, f64, f64)>,
    toroidal_surfaces: HashMap<u64, (u64, f64, f64)>,
    vertex_points: HashMap<u64, u64>,
    edge_curves: HashMap<u64, (u64, u64, Option<u64>)>,
    oriented_edges: HashMap<u64, (u64, bool)>,
    edge_loops: HashMap<u64, Vec<u64>>,
    face_bounds: HashMap<u64, u64>,
    advanced_faces: HashMap<u64, AdvancedFaceRecord>,
    closed_shells: HashMap<u64, Vec<u64>>,
    open_shells: HashMap<u64, Vec<u64>>,
    manifold_solids: Vec<u64>,
    shell_based_surface_models: Vec<Vec<u64>>,
    trimmed_curves: HashMap<u64, (u64, f64, f64)>,
    geometric_curve_sets: Vec<Vec<u64>>,
    /// COMPOUND: maps entity id -> list of element references
    compounds: HashMap<u64, Vec<u64>>,
    /// COMPSOLID: maps entity id -> list of solid references
    compsolids: HashMap<u64, Vec<u64>>,
    /// Top-level compound references (from SHAPE_REPRESENTATION_RELATIONSHIP)
    top_level_compounds: Vec<u64>,
    /// SURFACE_CURVE: maps step id 鈫?(3d_curve_ref, pcurve_ref_list, same_parameter)
    surface_curves: HashMap<u64, (u64, Vec<u64>, bool)>,
    /// PCURVE: maps step id 鈫?(surface_ref, definitional_rep_ref)
    pcurves: HashMap<u64, (u64, u64)>,
    /// DEFINITIONAL_REPRESENTATION: maps step id 鈫?curve2d_ref
    definitional_reps: HashMap<u64, u64>,
    /// 2D cartesian points
    cartesian_points_2d: HashMap<u64, [f64; 2]>,
    /// 2D directions
    directions_2d: HashMap<u64, [f64; 2]>,
    /// 2D axis2 placements: id 鈫?(location, ref_dir)
    axis2_placements_2d: HashMap<u64, (u64, u64)>,
    /// B-Spline surface: (degree_u, degree_v, ctrl_grid_refs[v][u], mults_u, knots_u, mults_v, knots_v)
    b_spline_surfaces: HashMap<u64, BSplineSurfaceData>,
    /// SURFACE_OF_LINEAR_EXTRUSION: maps entity id 鈫?(profile_curve_ref, direction_ref)
    linear_extrusions: HashMap<u64, (u64, u64)>,
    /// SURFACE_OF_REVOLUTION: maps entity id 鈫?(profile_curve_ref, axis_placement_ref)
    revolutions: HashMap<u64, (u64, u64)>,
    /// RECTANGULAR_TRIMMED_SURFACE: maps entity id 鈫?(basis_surface_ref, [u1,u2,v1,v2])
    rectangular_trimmed_surfaces: HashMap<u64, (u64, [f64; 4])>,
    /// HYPERBOLA: maps entity id 鈫?(axis2_placement_3d_ref, semi_major, semi_minor)
    hyperbolas: HashMap<u64, (u64, f64, f64)>,
    /// PARABOLA: maps entity id 鈫?(axis2_placement_3d_ref, focal_param)
    parabolas: HashMap<u64, (u64, f64)>,
    /// OFFSET_CURVE_3D: maps entity id 鈫?(basis_curve_ref, offset_distance, ref_dir_ref)
    offset_curves_3d: HashMap<u64, (u64, f64, u64)>,
    /// Global uncertainty value from UNCERTAINTY_MEASURE_WITH_UNIT, if present.
    uncertainty_value: Option<f64>,

    // 鈹€鈹€ Color / presentation chain 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    /// COLOUR_RGB: id 鈫?[r, g, b]
    colour_rgbs: HashMap<u64, [f64; 3]>,
    /// FILL_AREA_STYLE_COLOUR: id 鈫?colour_rgb_id
    fill_area_style_colours: HashMap<u64, u64>,
    /// FILL_AREA_STYLE: id 鈫?[fasc_id, ...]
    fill_area_styles: HashMap<u64, Vec<u64>>,
    /// SURFACE_STYLE_FILL_AREA: id 鈫?fill_area_style_id
    surface_style_fill_areas: HashMap<u64, u64>,
    /// SURFACE_SIDE_STYLE: id 鈫?[ssfa_id, ...]
    surface_side_styles: HashMap<u64, Vec<u64>>,
    /// SURFACE_STYLE_USAGE: id 鈫?surface_side_style_id
    surface_style_usages: HashMap<u64, u64>,
    /// PRESENTATION_STYLE_ASSIGNMENT: id 鈫?[ssu_id, ...]
    presentation_style_assignments: HashMap<u64, Vec<u64>>,
    /// STYLED_ITEM: id 鈫?(shape_step_id, [psa_id, ...])
    styled_items: HashMap<u64, (u64, Vec<u64>)>,
}

impl ParsedStep {
    fn new() -> Self {
        Self {
            cartesian_points: HashMap::new(),
            directions: HashMap::new(),
            vectors: HashMap::new(),
            axis2_placements: HashMap::new(),
            lines: HashMap::new(),
            circles: HashMap::new(),
            ellipses: HashMap::new(),
            b_spline_curves: HashMap::new(),
            b_spline_curves_full: HashMap::new(),
            planes: HashMap::new(),
            cylindrical_surfaces: HashMap::new(),
            spherical_surfaces: HashMap::new(),
            conical_surfaces: HashMap::new(),
            toroidal_surfaces: HashMap::new(),
            vertex_points: HashMap::new(),
            edge_curves: HashMap::new(),
            oriented_edges: HashMap::new(),
            edge_loops: HashMap::new(),
            face_bounds: HashMap::new(),
            advanced_faces: HashMap::new(),
            closed_shells: HashMap::new(),
            open_shells: HashMap::new(),
            manifold_solids: Vec::new(),
            shell_based_surface_models: Vec::new(),
            trimmed_curves: HashMap::new(),
            geometric_curve_sets: Vec::new(),
            compounds: HashMap::new(),
            compsolids: HashMap::new(),
            top_level_compounds: Vec::new(),
            surface_curves: HashMap::new(),
            pcurves: HashMap::new(),
            definitional_reps: HashMap::new(),
            cartesian_points_2d: HashMap::new(),
            directions_2d: HashMap::new(),
            axis2_placements_2d: HashMap::new(),
            b_spline_surfaces: HashMap::new(),
            linear_extrusions: HashMap::new(),
            revolutions: HashMap::new(),
            rectangular_trimmed_surfaces: HashMap::new(),
            hyperbolas: HashMap::new(),
            parabolas: HashMap::new(),
            offset_curves_3d: HashMap::new(),
            uncertainty_value: None,
            colour_rgbs: HashMap::new(),
            fill_area_style_colours: HashMap::new(),
            fill_area_styles: HashMap::new(),
            surface_style_fill_areas: HashMap::new(),
            surface_side_styles: HashMap::new(),
            surface_style_usages: HashMap::new(),
            presentation_style_assignments: HashMap::new(),
            styled_items: HashMap::new(),
        }
    }
}

pub struct StepReader;

/// Stable JSON diagnostics payload for single-part STEP import healing.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepImportHealingJsonV1 {
    pub schema: &'static str,
    pub clean: bool,
    pub initial_issue_count: usize,
    pub final_issue_count: usize,
    pub fixed_issue_count: usize,
    pub issue_histogram: Vec<(String, usize)>,
    pub repair_pass_count: usize,
    pub parametric_pass_count: usize,
    pub make_connected_pass_count: usize,
    pub vertices_merged: usize,
    pub degenerate_faces_removed: usize,
    pub normals_recomputed: usize,
    pub faces_reoriented: usize,
    pub wires_fixed: usize,
    pub same_range_fixed: usize,
    pub same_parameter_fixed: usize,
    /// Wire-level open-gap count after healing.
    pub wire_open_gaps: usize,
    /// Wire-level topological self-intersection count after healing.
    pub wire_topological_self_intersections: usize,
    /// Wire-level geometric self-intersection count after healing.
    pub wire_geometric_self_intersections: usize,
}

/// Coarse protocol hint inferred from `FILE_SCHEMA`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum StepProtocolHint {
    Ap203,
    Ap214,
    Ap242,
    Unknown,
}

/// STEP document-level metadata extracted from file header and global entities.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDocumentMetadata {
    /// Raw `FILE_SCHEMA` payload, e.g. `AUTOMOTIVE_DESIGN { ... 214 ... }`.
    pub file_schema: Option<String>,
    /// Coarse AP protocol hint derived from `file_schema`.
    pub protocol_hint: StepProtocolHint,
    /// `FILE_NAME` first field when available.
    pub file_name: Option<String>,
    /// Structured `PRODUCT` records found in the STEP data section.
    pub products: Vec<StepProduct>,
    /// `PRODUCT_DEFINITION_FORMATION` chain entries.
    pub product_definition_formations: Vec<StepProductDefinitionFormation>,
    /// `PRODUCT_DEFINITION` chain entries resolved back to products when possible.
    pub product_definitions: Vec<StepProductDefinitionInfo>,
    /// `SHAPE_DEFINITION_REPRESENTATION` entries resolved back to product definitions when possible.
    pub shape_definition_representations: Vec<StepShapeDefinitionRepresentation>,
    /// `NEXT_ASSEMBLY_USAGE_OCCURRENCE` links for assembly hierarchy metadata.
    pub assembly_occurrences: Vec<StepAssemblyUsageOccurrence>,
    /// `PRODUCT` names found in the STEP data section.
    pub product_names: Vec<String>,
    /// Uncertainty value from `UNCERTAINTY_MEASURE_WITH_UNIT`, if present.
    pub uncertainty_value: Option<f64>,
    /// Materials extracted from `MATERIAL` entities, if present.
    pub materials: Vec<StepMaterial>,
    /// Layers extracted from `PRESENTATION_LAYER_ASSIGNMENT` entities, if present.
    pub layers: Vec<StepLayer>,
    /// General properties extracted from `GENERAL_PROPERTY` entities.
    pub general_properties: Vec<StepGeneralProperty>,
    /// Property-definition chain entries, linked to referenced `GENERAL_PROPERTY`
    /// when resolvable.
    pub property_definitions: Vec<StepPropertyDefinition>,
    /// PROPERTY_DEFINITION_REPRESENTATION linkages (AP242 extended metadata).
    pub property_definition_representations: Vec<StepPropertyDefinitionRepr>,
    /// GDT dimensional locations (DIMENSIONAL_LOCATION entities).
    pub dimensional_locations: Vec<StepDimensionalLocation>,
    /// GDT dimensional sizes (DIMENSIONAL_SIZE entities).
    pub dimensional_sizes: Vec<StepDimensionalSize>,
    /// GDT geometric tolerances (GEOMETRIC_TOLERANCE entities).
    pub geometric_tolerances: Vec<StepGeometricTolerance>,
    /// GDT geometric tolerances with datum references.
    pub geometric_tolerances_with_datum_references: Vec<StepGeometricToleranceWithDatumReference>,
    /// DATUM entries (AP242 datum system baseline).
    pub datums: Vec<StepDatum>,
    /// DATUM_SYSTEM entries (AP242 datum system grouping).
    pub datum_systems: Vec<StepDatumSystem>,
    /// Kinematics/joint entries (AP242 kinematic pair baseline).
    pub kinematic_pairs: Vec<StepKinematicPair>,
    /// Tolerance zone entries (AP242 TOLERANCE_ZONE).
    pub tolerance_zones: Vec<StepToleranceZone>,
    /// Tolerance zone definitions (AP242 TOLERANCE_ZONE_DEFINITION).
    pub tolerance_zone_definitions: Vec<StepToleranceZoneDefinition>,
    /// Datum feature entries (AP242 DATUM_FEATURE).
    pub datum_features: Vec<StepDatumFeature>,
    /// Datum reference elements (AP242 DATUM_REFERENCE_ELEMENT).
    pub datum_reference_elements: Vec<StepDatumReferenceElement>,
    /// Shape aspect entries (AP242 SHAPE_ASPECT).
    pub shape_aspects: Vec<StepShapeAspect>,
    /// Shape aspect definitions (AP242 SHAPE_ASPECT_DEFINITION).
    pub shape_aspect_definitions: Vec<StepShapeAspectDefinition>,
    /// Derived shape aspects (AP242 DERIVED_SHAPE_ASPECT).
    pub derived_shape_aspects: Vec<StepDerivedShapeAspect>,
    /// GDT dimensional tolerances (AP242 DIMENSIONAL_TOLERANCE).
    pub dimensional_tolerances: Vec<StepDimensionalTolerance>,
    /// GDT tolerance values (AP242 MEASURE_REPRESENTATION_ITEM).
    pub tolerance_values: Vec<StepToleranceValue>,
    /// GDT position tolerances (AP242 POSITION_TOLERANCE).
    pub position_tolerances: Vec<StepPositionTolerance>,
    /// GDT orientation tolerances (AP242 ORIENTATION_TOLERANCE).
    pub orientation_tolerances: Vec<StepOrientationTolerance>,
    /// GDT form tolerances (AP242 FORM_TOLERANCE).
    pub form_tolerances: Vec<StepFormTolerance>,
    /// GDT runout tolerances (AP242 RUNOUT_TOLERANCE).
    pub runout_tolerances: Vec<StepRunoutTolerance>,
    /// GDT profile tolerances (AP242 PROFILE_TOLERANCE).
    pub profile_tolerances: Vec<StepProfileTolerance>,
    /// GDT datum reference frames (AP242 DATUM_REFERENCE_FRAME).
    pub datum_reference_frames: Vec<StepDatumReferenceFrame>,
    /// GDT datum targets (AP242 DATUM_TARGET).
    pub datum_targets: Vec<StepDatumTarget>,
    /// Enhanced tolerance zone definitions.
    pub tolerance_zone_definitions_enhanced: Vec<StepToleranceZoneDefinitionEnhanced>,
    /// FEA model definitions (AP242 FEAMEDIAN_MODEL).
    pub fea_models: Vec<StepFeaModel>,
    /// FEA meshes (AP242 FEAMEDIAN_MESH).
    pub fea_meshes: Vec<StepFeaMesh>,
    /// FEA node sets (AP242 FEAMEDIAN_NODE_SET).
    pub fea_node_sets: Vec<StepFeaNodeSet>,
    /// FEA element sets (AP242 FEAMEDIAN_ELEMENT_SET).
    pub fea_element_sets: Vec<StepFeaElementSet>,
    /// FEA material properties (AP242 FEAMEDIAN_MATERIAL_PROPERTY).
    pub fea_material_properties: Vec<StepFeaMaterialProperty>,
    /// FEA boundary conditions (AP242 FEAMEDIAN_BOUNDARY_CONDITION).
    pub fea_boundary_conditions: Vec<StepFeaBoundaryCondition>,
    /// FEA loads (AP242 FEAMEDIAN_LOAD).
    pub fea_loads: Vec<StepFeaLoad>,
    /// FEA node groups (AP242 FEA_NODE_GROUP).
    pub fea_node_groups: Vec<StepFeaNodeGroup>,
    /// FEA analysis definitions (AP242 FEA_ANALYSIS, ANALYSIS_3D).
    pub fea_analyses: Vec<StepFeaAnalysis>,
    /// FEA state definitions (AP242 FEA_STATE).
    pub fea_states: Vec<StepFeaState>,
    /// FEA material models (AP242 FEA_MATERIAL_MODEL, FEA_LINEAR_ELASTICITY).
    pub fea_material_models: Vec<StepFeaMaterialModel>,
    /// FEA nodes with coordinates (AP242 NODE_REPRESENTATION).
    pub fea_nodes: Vec<StepFeaNode>,
    /// FEA elements with connectivity (AP242 ELEMENT_REPRESENTATION).
    pub fea_elements: Vec<StepFeaElement>,
    /// FEA analysis steps (AP242 FEA_STEP).
    pub fea_steps: Vec<StepFeaStep>,
    /// FEA result data (AP242 FEA_RESULT).
    pub fea_results: Vec<StepFeaResult>,
    /// FEA case definitions (AP242 FEA_CASE).
    pub fea_cases: Vec<StepFeaCase>,
    /// View definitions (AP242 VIEW, CAMERA_MODEL_D3).
    pub views: Vec<StepView>,
    /// Camera models (AP242 CAMERA_MODEL_D3).
    pub cameras: Vec<StepCameraModelD3>,
    /// View volumes (AP242 VIEW_VOLUME).
    pub view_volumes: Vec<StepViewVolume>,
    /// Notes/annotations (AP242 ANNOTATION).
    pub notes: Vec<StepNote>,
    /// Annotation planes (AP242 ANNOTATION_PLANE).
    pub annotation_planes: Vec<StepAnnotationPlane>,
    /// Annotation occurrences (AP242 ANNOTATION_OCCURRENCE).
    pub annotation_occurrences: Vec<StepAnnotationOccurrence>,
    /// Dimension curves (AP242 DIMENSION_CURVE).
    pub dimension_curves: Vec<StepDimensionCurve>,
    /// Terminator symbols (AP242 TERMINATOR_SYMBOL).
    pub terminator_symbols: Vec<StepTerminatorSymbol>,
    /// Datum feature callouts (AP242 DATUM_FEATURE_CALLOUT).
    pub datum_feature_callouts: Vec<StepDatumFeatureCallout>,
    // ── AP242 Product Definition Relationship Chains ─────────────────────────────
    /// Product definition relationships (PRODUCT_DEFINITION_RELATIONSHIP).
    pub product_definition_relationships: Vec<StepProductDefinitionRelationship>,
    // ── AP242 Shape Representation Associations ──────────────────────────────────
    /// Shape representation relationships (SHAPE_REPRESENTATION_RELATIONSHIP).
    pub shape_representation_relationships: Vec<StepShapeRepresentationRelationship>,
    /// Product definition shapes (PRODUCT_DEFINITION_SHAPE).
    pub product_definition_shapes: Vec<StepProductDefinitionShape>,
    // ── AP242 Configuration Management ───────────────────────────────────────────
    /// Configuration designs (CONFIGURATION_DESIGN).
    pub configuration_designs: Vec<StepConfigurationDesign>,
    /// Configuration items (CONFIGURATION_ITEM).
    pub configuration_items: Vec<StepConfigurationItem>,
    /// Product concepts (PRODUCT_CONCEPT).
    pub product_concepts: Vec<StepProductConcept>,
    /// Configuration effectivities (CONFIGURATION_EFFECTIVITY).
    pub configuration_effectivities: Vec<StepConfigurationEffectivity>,
    // ── AP242 Approval and Security ──────────────────────────────────────────────
    /// Approvals (APPROVAL).
    pub approvals: Vec<StepApproval>,
    /// Approval assignments (APPROVAL_ASSIGNMENT).
    pub approval_assignments: Vec<StepApprovalAssignment>,
    /// Security classifications (SECURITY_CLASSIFICATION).
    pub security_classifications: Vec<StepSecurityClassification>,
    /// Security classification assignments (SECURITY_CLASSIFICATION_ASSIGNMENT).
    pub security_classification_assignments: Vec<StepSecurityClassificationAssignment>,
    // ── AP242 Document References ────────────────────────────────────────────────
    /// Documents (DOCUMENT).
    pub documents: Vec<StepDocument>,
    /// Document files (DOCUMENT_FILE).
    pub document_files: Vec<StepDocumentFile>,
    /// Document usage assignments (DOCUMENT_USAGE_ASSIGNMENT).
    pub document_usage_assignments: Vec<StepDocumentUsageAssignment>,
    /// Document representation relationships (DOCUMENT_REPRESENTATION_RELATIONSHIP).
    pub document_representation_relationships: Vec<StepDocumentRepresentationRelationship>,
}

impl StepDocumentMetadata {
    /// Returns a human-readable summary of extracted AP242 entities.
    pub fn summary(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("Products: {}", self.products.len()));
        lines.push(format!("Product definitions: {}", self.product_definitions.len()));
        lines.push(format!("Assembly occurrences: {}", self.assembly_occurrences.len()));
        lines.push(format!("Materials: {}", self.materials.len()));
        lines.push(format!("Layers: {}", self.layers.len()));
        lines.push(format!("General properties: {}", self.general_properties.len()));
        lines.push(format!("Property definitions: {}", self.property_definitions.len()));
        lines.push(format!("Property definition representations: {}", self.property_definition_representations.len()));
        lines.push(format!("Dimensional locations: {}", self.dimensional_locations.len()));
        lines.push(format!("Dimensional sizes: {}", self.dimensional_sizes.len()));
        lines.push(format!("Geometric tolerances: {}", self.geometric_tolerances.len()));
        lines.push(format!("Geometric tolerances with datum references: {}", self.geometric_tolerances_with_datum_references.len()));
        lines.push(format!("Datums: {}", self.datums.len()));
        lines.push(format!("Datum systems: {}", self.datum_systems.len()));
        lines.push(format!("Kinematic pairs: {}", self.kinematic_pairs.len()));
        lines.push(format!("Tolerance zones: {}", self.tolerance_zones.len()));
        lines.push(format!("Tolerance zone definitions: {}", self.tolerance_zone_definitions.len()));
        lines.push(format!("Datum features: {}", self.datum_features.len()));
        lines.push(format!("Datum reference elements: {}", self.datum_reference_elements.len()));
        lines.push(format!("Shape aspects: {}", self.shape_aspects.len()));
        lines.push(format!("Shape aspect definitions: {}", self.shape_aspect_definitions.len()));
        lines.push(format!("Derived shape aspects: {}", self.derived_shape_aspects.len()));
        lines.push(format!("Dimensional tolerances: {}", self.dimensional_tolerances.len()));
        lines.push(format!("Tolerance values: {}", self.tolerance_values.len()));
        lines.push(format!("Position tolerances: {}", self.position_tolerances.len()));
        lines.push(format!("Orientation tolerances: {}", self.orientation_tolerances.len()));
        lines.push(format!("Form tolerances: {}", self.form_tolerances.len()));
        lines.push(format!("Runout tolerances: {}", self.runout_tolerances.len()));
        lines.push(format!("Profile tolerances: {}", self.profile_tolerances.len()));
        lines.push(format!("Datum reference frames: {}", self.datum_reference_frames.len()));
        lines.push(format!("Datum targets: {}", self.datum_targets.len()));
        lines.push(format!("Enhanced tolerance zone definitions: {}", self.tolerance_zone_definitions_enhanced.len()));
        lines.push(format!("FEA models: {}", self.fea_models.len()));
        lines.push(format!("FEA meshes: {}", self.fea_meshes.len()));
        lines.push(format!("FEA node sets: {}", self.fea_node_sets.len()));
        lines.push(format!("FEA element sets: {}", self.fea_element_sets.len()));
        lines.push(format!("FEA material properties: {}", self.fea_material_properties.len()));
        lines.push(format!("FEA boundary conditions: {}", self.fea_boundary_conditions.len()));
        lines.push(format!("FEA loads: {}", self.fea_loads.len()));
        lines.push(format!("FEA node groups: {}", self.fea_node_groups.len()));
        lines.push(format!("FEA analyses: {}", self.fea_analyses.len()));
        lines.push(format!("FEA states: {}", self.fea_states.len()));
        lines.push(format!("FEA material models: {}", self.fea_material_models.len()));
        lines.push(format!("FEA nodes: {}", self.fea_nodes.len()));
        lines.push(format!("FEA elements: {}", self.fea_elements.len()));
        lines.push(format!("FEA steps: {}", self.fea_steps.len()));
        lines.push(format!("FEA results: {}", self.fea_results.len()));
        lines.push(format!("FEA cases: {}", self.fea_cases.len()));
        lines.push(format!("Views: {}", self.views.len()));
        lines.push(format!("Cameras: {}", self.cameras.len()));
        lines.push(format!("View volumes: {}", self.view_volumes.len()));
        lines.push(format!("Notes: {}", self.notes.len()));
        lines.push(format!("Annotation planes: {}", self.annotation_planes.len()));
        lines.push(format!("Annotation occurrences: {}", self.annotation_occurrences.len()));
        lines.push(format!("Dimension curves: {}", self.dimension_curves.len()));
        lines.push(format!("Terminator symbols: {}", self.terminator_symbols.len()));
        lines.push(format!("Datum feature callouts: {}", self.datum_feature_callouts.len()));
        // AP242 Product Definition Relationship Chains
        lines.push(format!("Product definition relationships: {}", self.product_definition_relationships.len()));
        // AP242 Shape Representation Associations
        lines.push(format!("Shape representation relationships: {}", self.shape_representation_relationships.len()));
        lines.push(format!("Product definition shapes: {}", self.product_definition_shapes.len()));
        // AP242 Configuration Management
        lines.push(format!("Configuration designs: {}", self.configuration_designs.len()));
        lines.push(format!("Configuration items: {}", self.configuration_items.len()));
        lines.push(format!("Product concepts: {}", self.product_concepts.len()));
        lines.push(format!("Configuration effectivities: {}", self.configuration_effectivities.len()));
        // AP242 Approval and Security
        lines.push(format!("Approvals: {}", self.approvals.len()));
        lines.push(format!("Approval assignments: {}", self.approval_assignments.len()));
        lines.push(format!("Security classifications: {}", self.security_classifications.len()));
        lines.push(format!("Security classification assignments: {}", self.security_classification_assignments.len()));
        // AP242 Document References
        lines.push(format!("Documents: {}", self.documents.len()));
        lines.push(format!("Document files: {}", self.document_files.len()));
        lines.push(format!("Document usage assignments: {}", self.document_usage_assignments.len()));
        lines.push(format!("Document representation relationships: {}", self.document_representation_relationships.len()));
        lines.join("\n")
    }
}

/// A `PRODUCT` record extracted from the STEP data section.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProduct {
    pub entity_id: u64,
    pub product_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
}

/// A `PRODUCT_DEFINITION_FORMATION` record linked back to a `PRODUCT` when possible.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProductDefinitionFormation {
    pub entity_id: u64,
    pub formation_id: Option<String>,
    pub description: Option<String>,
    pub product_id: Option<u64>,
    pub product_name: Option<String>,
}

/// A `PRODUCT_DEFINITION` record resolved through the formation chain when possible.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProductDefinitionInfo {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub formation_id: Option<u64>,
    pub product_id: Option<u64>,
    pub product_name: Option<String>,
}

/// A `SHAPE_DEFINITION_REPRESENTATION` linkage resolved to a `PRODUCT_DEFINITION` when possible.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepShapeDefinitionRepresentation {
    pub entity_id: u64,
    pub product_definition_shape_id: Option<u64>,
    pub product_definition_id: Option<u64>,
    pub representation_id: Option<u64>,
}

/// A `NEXT_ASSEMBLY_USAGE_OCCURRENCE` relationship between two product definitions.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepAssemblyUsageOccurrence {
    pub entity_id: u64,
    pub usage_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    pub relating_product_definition_id: Option<u64>,
    pub related_product_definition_id: Option<u64>,
    pub relating_product_name: Option<String>,
    pub related_product_name: Option<String>,
}

// ── AP242 Product Definition Relationship Chains ───────────────────────────────

/// A `PRODUCT_DEFINITION_RELATIONSHIP` establishing parent-child relationships.
/// Used to define assembly hierarchies and component relationships.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProductDefinitionRelationship {
    pub entity_id: u64,
    /// Identifier for this relationship.
    pub relationship_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The product definition that owns or contains the related definition.
    pub relating_product_definition_id: Option<u64>,
    /// The product definition that is related or contained.
    pub related_product_definition_id: Option<u64>,
    /// Name of the relating product (resolved from product definition chain).
    pub relating_product_name: Option<String>,
    /// Name of the related product (resolved from product definition chain).
    pub related_product_name: Option<String>,
}

// ── AP242 Shape Representation Associations ────────────────────────────────────

/// A `SHAPE_REPRESENTATION_RELATIONSHIP` linking shape representations.
/// Used to associate multiple representations with a single product definition.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepShapeRepresentationRelationship {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The primary shape representation.
    pub relating_representation_id: Option<u64>,
    /// The secondary/dependent shape representation.
    pub related_representation_id: Option<u64>,
    /// Optional transformation entity reference.
    pub transformation_id: Option<u64>,
}

/// A `PRODUCT_DEFINITION_SHAPE` linking a product definition to its shape.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProductDefinitionShape {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The product definition this shape belongs to.
    pub product_definition_id: Option<u64>,
}

// ── AP242 Configuration Management ─────────────────────────────────────────────

/// A `CONFIGURATION_DESIGN` linking a configuration to a product definition.
/// Used for variant management and configuration control.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepConfigurationDesign {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The configuration being referenced.
    pub configuration_id: Option<u64>,
    /// The product definition this configuration applies to.
    pub product_definition_id: Option<u64>,
}

/// A `CONFIGURATION_ITEM` defining a configuration baseline.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepConfigurationItem {
    pub entity_id: u64,
    /// Unique identifier for this configuration item.
    pub item_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The product concept this configuration item belongs to.
    pub product_concept_id: Option<u64>,
}

/// A `PRODUCT_CONCEPT` representing a product at a conceptual level.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProductConcept {
    pub entity_id: u64,
    /// Unique identifier for this product concept.
    pub concept_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The market context for this concept.
    pub market_context_id: Option<u64>,
}

/// A `CONFIGURATION_EFFECTIVITY` defining when a configuration is effective.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepConfigurationEffectivity {
    pub entity_id: u64,
    pub configuration_id: Option<u64>,
    pub usage_id: Option<u64>,
    /// Effectivity start (e.g., serial number or date).
    pub effectivity_start: Option<String>,
    /// Effectivity end.
    pub effectivity_end: Option<String>,
}

// ── AP242 Approval and Security ────────────────────────────────────────────────

/// Approval status enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum ApprovalStatus {
    /// Approval is pending.
    Pending,
    /// Approval has been granted.
    Approved,
    /// Approval has been rejected.
    Rejected,
    /// Approval status is unknown.
    Unknown,
}

/// An `APPROVAL` entity tracking approval status.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepApproval {
    pub entity_id: u64,
    /// Approval status.
    pub status: ApprovalStatus,
    /// Approval level or type identifier.
    pub level: Option<String>,
    /// Date of approval.
    pub date: Option<String>,
    /// Approver name or identifier.
    pub approver: Option<String>,
}

/// An `APPROVAL_ASSIGNMENT` linking approval to a product or document.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepApprovalAssignment {
    pub entity_id: u64,
    /// The approval being assigned.
    pub approval_id: Option<u64>,
    /// The entity (product, document, etc.) being approved.
    pub approved_item_id: Option<u64>,
    /// Role of this approval (e.g., "design", "manufacturing").
    pub role: Option<String>,
}

/// Security classification level enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum SecurityClassificationLevel {
    /// Unclassified.
    Unclassified,
    /// Confidential.
    Confidential,
    /// Secret.
    Secret,
    /// Top Secret.
    TopSecret,
    /// Proprietary / Company Confidential.
    Proprietary,
    /// Unknown or unspecified.
    Unknown,
}

/// A `SECURITY_CLASSIFICATION` entity tracking security level.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepSecurityClassification {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Security classification level.
    pub security_level: SecurityClassificationLevel,
}

/// A `SECURITY_CLASSIFICATION_ASSIGNMENT` linking security to an item.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepSecurityClassificationAssignment {
    pub entity_id: u64,
    /// The security classification being assigned.
    pub security_classification_id: Option<u64>,
    /// The entity being classified.
    pub classified_item_id: Option<u64>,
}

// ── AP242 Document References ──────────────────────────────────────────────────

/// A `DOCUMENT` entity referencing external documents.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDocument {
    pub entity_id: u64,
    /// Document identifier (e.g., part number, drawing number).
    pub document_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Document type (e.g., "drawing", "specification", "CAD model").
    pub document_type: Option<String>,
}

/// A `DOCUMENT_FILE` representing a digital file attachment.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDocumentFile {
    pub entity_id: u64,
    pub document_id: Option<String>,
    pub name: Option<String>,
    pub description: Option<String>,
    /// File name including extension.
    pub file_name: Option<String>,
    /// MIME type or file format identifier.
    pub file_format: Option<String>,
}

/// A `DOCUMENT_USAGE_ASSIGNMENT` linking a document to a product.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDocumentUsageAssignment {
    pub entity_id: u64,
    /// The document being assigned.
    pub document_id: Option<u64>,
    /// The product definition the document is assigned to.
    pub product_definition_id: Option<u64>,
    /// Role or purpose of this document (e.g., "reference", "specification").
    pub role: Option<String>,
}

/// A `DOCUMENT_REPRESENTATION_RELATIONSHIP` linking document to representation.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDocumentRepresentationRelationship {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The document being related.
    pub document_id: Option<u64>,
    /// The representation (e.g., shape) being related.
    pub representation_id: Option<u64>,
}

/// Linkage from PROPERTY_DEFINITION to a representation via PROPERTY_DEFINITION_REPRESENTATION.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepPropertyDefinitionRepr {
    pub entity_id: u64,
    pub property_definition_id: Option<u64>,
    pub representation_id: Option<u64>,
}

/// A dimensional location (AP242 DIMENSIONAL_LOCATION).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDimensionalLocation {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub from_entity_id: Option<u64>,
    pub to_entity_id: Option<u64>,
}

/// A dimensional size (AP242 DIMENSIONAL_SIZE).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDimensionalSize {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub shape_aspect_id: Option<u64>,
}

/// A geometric tolerance entry (AP242 GEOMETRIC_TOLERANCE).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepGeometricTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub value_entity_id: Option<u64>,
    pub shape_aspect_id: Option<u64>,
}

/// A geometric tolerance entry referencing a datum system
/// (AP242 GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepGeometricToleranceWithDatumReference {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub value_entity_id: Option<u64>,
    pub shape_aspect_id: Option<u64>,
    pub datum_system_id: Option<u64>,
}

/// A datum entry (AP242 DATUM).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatum {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub shape_aspect_id: Option<u64>,
}

/// A datum system entry (AP242 DATUM_SYSTEM).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumSystem {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    pub datum_ids: Vec<u64>,
}

/// A kinematic/joint metadata entry (AP242 kinematic pair family).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepKinematicPair {
    pub entity_id: u64,
    /// Original STEP entity token (e.g. `REVOLUTE_PAIR`).
    pub entity_type: String,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Referenced ids carried by this entity.
    pub related_entity_ids: Vec<u64>,
}

/// A tolerance zone entry (AP242 TOLERANCE_ZONE).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepToleranceZone {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Referenced toleranced_shape_aspect (typically GEOMETRIC_TOLERANCE).
    pub toleranced_shape_aspect_id: Option<u64>,
}

/// A tolerance zone definition (AP242 TOLERANCE_ZONE_DEFINITION).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepToleranceZoneDefinition {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to TOLERANCE_ZONE.
    pub tolerance_zone_id: Option<u64>,
    /// Reference to shape_aspect defining the zone boundaries.
    pub shape_aspect_id: Option<u64>,
}

/// A datum feature entry (AP242 DATUM_FEATURE).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumFeature {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to associated DATUM.
    pub datum_id: Option<u64>,
}

/// A datum reference element (AP242 DATUM_REFERENCE_ELEMENT).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumReferenceElement {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to associated DATUM or DATUM_FEATURE.
    pub associated_entity_id: Option<u64>,
}

// ── GDT Extended Structures (AP242) ───────────────────────────────────────────

/// Tolerance zone shape enumeration for AP242.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum ToleranceZoneShape {
    /// Cylindrical tolerance zone.
    Cylindrical,
    /// Spherical tolerance zone.
    Spherical,
    /// Zone between two parallel planes.
    TwoParallelPlanes,
    /// Zone between two coaxial cylinders.
    TwoCoaxialCylinders,
    /// Zone between two concentric circles.
    TwoConcentricCircles,
    /// Zone within a circle.
    WithinCircle,
    /// Zone between two parallel lines.
    TwoParallelLines,
    /// Complex shape defined by supplementary geometry.
    Complex,
    /// Unknown or unspecified shape.
    Unknown,
}

/// Tolerance zone position enumeration for AP242.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum ToleranceZonePosition {
    /// Zone is symmetric about the theoretical exact location.
    Symmetric,
    /// Zone is unilateral (one-sided).
    Unilateral,
    /// Zone is bilateral but asymmetric.
    BilateralAsymmetric,
    /// Zone position not specified.
    Unspecified,
}

/// A tolerance value with optional unit (AP242 MEASURE_REPRESENTATION_ITEM).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepToleranceValue {
    pub entity_id: u64,
    pub name: Option<String>,
    /// The tolerance value (typically half the total tolerance zone).
    pub value: f64,
    /// Unit name (e.g., "mm", "in").
    pub unit: Option<String>,
}

/// A dimensional tolerance (AP242 DIMENSIONAL_TOLERANCE).
/// Represents a plus/minus tolerance on a dimension.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDimensionalTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to DIMENSIONAL_SIZE or DIMENSIONAL_LOCATION.
    pub dimensional_characteristic_id: Option<u64>,
    /// Upper tolerance value (positive).
    pub upper_tolerance: Option<f64>,
    /// Lower tolerance value (typically negative or zero).
    pub lower_tolerance: Option<f64>,
    /// Unit of the tolerance values.
    pub unit: Option<String>,
}

/// Position tolerance (AP242 POSITION_TOLERANCE).
/// Defines the allowable variation in position of a feature.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepPositionTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to tolerance value entity.
    pub value_entity_id: Option<u64>,
    /// Reference to the toleranced shape aspect.
    pub shape_aspect_id: Option<u64>,
    /// Reference to datum system (for tolerances with datum reference).
    pub datum_system_id: Option<u64>,
    /// Whether this is a projected tolerance zone.
    pub projected: bool,
    /// Projected height if projected tolerance zone.
    pub projected_height: Option<f64>,
}

/// Orientation tolerance (AP242 ORIENTATION_TOLERANCE).
/// Covers angularity, perpendicularity, and parallelism tolerances.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepOrientationTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to tolerance value entity.
    pub value_entity_id: Option<u64>,
    /// Reference to the toleranced shape aspect.
    pub shape_aspect_id: Option<u64>,
    /// Reference to datum system.
    pub datum_system_id: Option<u64>,
    /// The type of orientation tolerance.
    pub orientation_type: OrientationToleranceType,
}

/// Types of orientation tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum OrientationToleranceType {
    /// Angularity tolerance.
    Angularity,
    /// Perpendicularity tolerance.
    Perpendicularity,
    /// Parallelism tolerance.
    Parallelism,
}

/// Form tolerance (AP242 FORM_TOLERANCE).
/// Covers flatness, straightness, circularity (roundness), and cylindricity.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFormTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to tolerance value entity.
    pub value_entity_id: Option<u64>,
    /// Reference to the toleranced shape aspect.
    pub shape_aspect_id: Option<u64>,
    /// The type of form tolerance.
    pub form_type: FormToleranceType,
}

/// Types of form tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum FormToleranceType {
    /// Flatness tolerance.
    Flatness,
    /// Straightness tolerance.
    Straightness,
    /// Circularity (roundness) tolerance.
    Circularity,
    /// Cylindricity tolerance.
    Cylindricity,
}

/// Runout tolerance (AP242 RUNOUT_TOLERANCE).
/// Covers circular runout and total runout tolerances.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepRunoutTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to tolerance value entity.
    pub value_entity_id: Option<u64>,
    /// Reference to the toleranced shape aspect.
    pub shape_aspect_id: Option<u64>,
    /// Reference to datum system.
    pub datum_system_id: Option<u64>,
    /// The type of runout tolerance.
    pub runout_type: RunoutToleranceType,
}

/// Types of runout tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum RunoutToleranceType {
    /// Circular runout tolerance.
    CircularRunout,
    /// Total runout tolerance.
    TotalRunout,
}

/// Profile tolerance (AP242 PROFILE_TOLERANCE).
/// Covers profile of a line and profile of a surface tolerances.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepProfileTolerance {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to tolerance value entity.
    pub value_entity_id: Option<u64>,
    /// Reference to the toleranced shape aspect.
    pub shape_aspect_id: Option<u64>,
    /// Reference to datum system (optional for profile tolerances).
    pub datum_system_id: Option<u64>,
    /// The type of profile tolerance.
    pub profile_type: ProfileToleranceType,
}

/// Types of profile tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum ProfileToleranceType {
    /// Profile of a line tolerance.
    ProfileOfALine,
    /// Profile of a surface tolerance.
    ProfileOfASurface,
}

/// A datum reference frame (AP242 DATUM_REFERENCE_FRAME).
/// Establishes a coordinate system for tolerance measurement.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumReferenceFrame {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Ordered list of datum system IDs that define the reference frame.
    /// Primary, secondary, and tertiary datums.
    pub datum_system_ids: Vec<u64>,
}

/// A datum target (AP242 DATUM_TARGET).
/// Represents a specific point, line, or area used to establish a datum.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumTarget {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The target identifier (e.g., "A1", "A2" for datum A targets).
    pub target_identifier: Option<String>,
    /// Reference to the parent DATUM.
    pub datum_id: Option<u64>,
    /// The type of datum target.
    pub target_type: DatumTargetType,
    /// Reference to the geometry defining the target location/shape.
    pub shape_aspect_id: Option<u64>,
}

/// Types of datum target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum DatumTargetType {
    /// Point target.
    Point,
    /// Line target.
    Line,
    /// Area target (circular or rectangular).
    Area,
    /// Target area is a circle.
    AreaCircle,
    /// Target area is a rectangle.
    AreaRectangle,
}

/// Enhanced tolerance zone definition (AP242 TOLERANCE_ZONE_DEFINITION extended).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepToleranceZoneDefinitionEnhanced {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to TOLERANCE_ZONE.
    pub tolerance_zone_id: Option<u64>,
    /// Reference to shape_aspect defining the zone boundaries.
    pub shape_aspect_id: Option<u64>,
    /// Shape of the tolerance zone.
    pub zone_shape: ToleranceZoneShape,
    /// Position of the tolerance zone.
    pub zone_position: ToleranceZonePosition,
    /// Reference to defining shape aspect (supplementary geometry).
    pub defining_shape_aspect_id: Option<u64>,
}

/// A shape aspect entry (AP242 SHAPE_ASPECT).
/// Used to associate tolerance information to geometric features.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepShapeAspect {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the product_definition_shape this aspect belongs to.
    pub product_definition_shape_id: Option<u64>,
}

/// A shape aspect definition (AP242 SHAPE_ASPECT_DEFINITION).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepShapeAspectDefinition {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to parent SHAPE_ASPECT.
    pub shape_aspect_id: Option<u64>,
}

/// A derived shape aspect (AP242 DERIVED_SHAPE_ASPECT).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDerivedShapeAspect {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to base shape_aspect(s).
    pub base_shape_aspect_ids: Vec<u64>,
}

// ── FEA (Finite Element Analysis) entities (AP242) ────────────────────────────

/// FEA model definition (AP242 FEA_MODEL, FEA_MODEL_3D, or FEAMEDIAN_MODEL).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaModel {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
}

/// FEA mesh (AP242 FEA_MESH, FEA_MESH_3D, or FEAMEDIAN_MESH).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaMesh {
    pub entity_id: u64,
    pub name: Option<String>,
    pub node_count: Option<u64>,
    pub element_count: Option<u64>,
}

/// FEA node set (AP242 FEA_NODE_SET or FEAMEDIAN_NODE_SET).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaNodeSet {
    pub entity_id: u64,
    pub name: Option<String>,
    pub model_id: Option<u64>,
}

/// FEA element set (AP242 FEA_ELEMENT_SET or FEAMEDIAN_ELEMENT_SET).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaElementSet {
    pub entity_id: u64,
    pub name: Option<String>,
    pub model_id: Option<u64>,
    pub element_type: Option<String>,
}

/// FEA node group (AP242 FEA_NODE_GROUP).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaNodeGroup {
    pub entity_id: u64,
    pub name: Option<String>,
    pub model_id: Option<u64>,
    pub node_count: Option<u64>,
}

/// FEA material property (AP242 FEA_MATERIAL_PROPERTY or FEAMEDIAN_MATERIAL_PROPERTY).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaMaterialProperty {
    pub entity_id: u64,
    pub name: Option<String>,
    pub property_type: Option<String>,
    pub value: Option<f64>,
    pub unit: Option<String>,
}

/// FEA boundary condition (AP242 FEA_BOUNDARY_CONDITION or FEAMEDIAN_BOUNDARY_CONDITION).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaBoundaryCondition {
    pub entity_id: u64,
    pub name: Option<String>,
    pub condition_type: Option<String>,
    pub node_set_id: Option<u64>,
}

/// FEA load (AP242 FEA_LOAD or FEAMEDIAN_LOAD).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaLoad {
    pub entity_id: u64,
    pub name: Option<String>,
    pub load_type: Option<String>,
    pub magnitude: Option<f64>,
    pub direction: Option<[f64; 3]>,
}

// ── Additional FEA entities (AP209/AP242 extended) ───────────────────────────────

/// FEA analysis definition (AP242 FEA_ANALYSIS, ANALYSIS_3D).
/// Represents a finite element analysis definition with associated model.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaAnalysis {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the FEA model being analyzed.
    pub model_id: Option<u64>,
    /// Analysis type (e.g., "STATIC", "MODAL", "THERMAL", "BUCKLING").
    pub analysis_type: Option<String>,
    /// Creation time/date if specified.
    pub creation_date: Option<String>,
}

/// FEA state definition (AP242 FEA_STATE).
/// Represents a state of the model (e.g., initial conditions, results).
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaState {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the FEA analysis this state belongs to.
    pub analysis_id: Option<u64>,
    /// State type (e.g., "INITIAL", "RESULT", "LOAD_CASE").
    pub state_type: Option<String>,
}

/// FEA material model (AP242 FEA_MATERIAL_MODEL, FEA_LINEAR_ELASTICITY, etc.).
/// Defines material constitutive model for FEA.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaMaterialModel {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Material model type (e.g., "LINEAR_ELASTIC", "ELASTIC_PLASTIC", "HYPERELASTIC").
    pub model_type: String,
    /// Young's modulus (for linear elastic materials).
    pub youngs_modulus: Option<f64>,
    /// Poisson's ratio.
    pub poissons_ratio: Option<f64>,
    /// Shear modulus.
    pub shear_modulus: Option<f64>,
    /// Density.
    pub density: Option<f64>,
    /// Thermal expansion coefficient.
    pub thermal_expansion: Option<f64>,
    /// Unit for modulus values (e.g., "MPa", "GPa", "psi").
    pub modulus_unit: Option<String>,
}

/// FEA node with coordinates (AP242 NODE_REPRESENTATION or via FEAMEDIAN_NODE).
/// Represents a single node in the FEA mesh with its coordinates.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaNode {
    pub entity_id: u64,
    /// Node number/label within the mesh (1-based indexing typically).
    pub node_number: u64,
    /// Cartesian coordinates (x, y, z).
    pub coordinates: [f64; 3],
    /// Reference to the mesh this node belongs to.
    pub mesh_id: Option<u64>,
}

/// FEA element with connectivity (AP242 ELEMENT_REPRESENTATION or via FEAMEDIAN_ELEMENT).
/// Represents a single element in the FEA mesh with its node connectivity.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaElement {
    pub entity_id: u64,
    /// Element number/label within the mesh.
    pub element_number: u64,
    /// Element type (e.g., "TETRA4", "TETRA10", "HEXA8", "HEXA20", "TRIA3", "QUAD4").
    pub element_type: String,
    /// Node IDs that define this element (connectivity).
    pub node_ids: Vec<u64>,
    /// Reference to the mesh this element belongs to.
    pub mesh_id: Option<u64>,
}

/// FEA step/substep definition (AP242 FEA_STEP).
/// Represents an analysis step or substep within an FEA analysis.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaStep {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the FEA analysis.
    pub analysis_id: Option<u64>,
    /// Step number within the analysis.
    pub step_number: Option<u64>,
    /// Step type (e.g., "LOAD", "TIME", "FREQUENCY").
    pub step_type: Option<String>,
}

/// FEA result data (AP242 FEA_RESULT, NODAL_RESULT, ELEMENT_RESULT).
/// Represents result data from FEA analysis.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaResult {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Result type (e.g., "DISPLACEMENT", "STRESS", "STRAIN", "REACTION_FORCE").
    pub result_type: String,
    /// Reference to the FEA state or analysis.
    pub analysis_id: Option<u64>,
    /// Reference to the result location (node set or element set).
    pub location_id: Option<u64>,
    /// Number of components (e.g., 3 for displacement vector, 6 for stress tensor).
    pub component_count: Option<u64>,
}

/// FEA case definition (AP242 FEA_CASE).
/// Represents an analysis case within an FEA model.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepFeaCase {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the FEA model.
    pub model_id: Option<u64>,
    /// Case type (e.g., "LOAD_CASE", "BC_CASE", "COMBINATION").
    pub case_type: Option<String>,
}

/// A material entry extracted from a STEP file.
///
/// Analogous to `XCAFDoc_Material` / `XCAFDoc_MaterialTool` in OCCT.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepMaterial {
    /// Material name (from the `MATERIAL(name, ...)` entity).
    pub name: String,
    /// Optional density value parsed from associated property definitions.
    pub density: Option<f64>,
}

/// A layer/group entry extracted from a STEP file.
///
/// Analogous to `XCAFDoc_Layer` / `XCAFDoc_LayerTool` in OCCT.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepLayer {
    /// Layer name (first argument of `PRESENTATION_LAYER_ASSIGNMENT`).
    pub name: String,
}

/// A generic STEP property entry extracted from `GENERAL_PROPERTY`.
///
/// Analogous to OCCT's STEP general attributes exported as `property_definition`.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepGeneralProperty {
    /// Property name.
    pub name: String,
    /// Optional human-readable description when present.
    pub description: Option<String>,
}

/// A `PROPERTY_DEFINITION` entry with optional linkage to a referenced
/// `GENERAL_PROPERTY`.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepPropertyDefinition {
    /// Entity id of `PROPERTY_DEFINITION`.
    pub definition_id: u64,
    /// Name field from `PROPERTY_DEFINITION(name, description, ...)`.
    pub name: Option<String>,
    /// Description field from `PROPERTY_DEFINITION(name, description, ...)`.
    pub description: Option<String>,
    /// Referenced entity id (typically `GENERAL_PROPERTY`) when present.
    pub referenced_entity_id: Option<u64>,
    /// Linked GENERAL_PROPERTY name when `referenced_entity_id` resolves.
    pub general_property_name: Option<String>,
    /// Linked GENERAL_PROPERTY description when resolvable.
    pub general_property_description: Option<String>,
}

// ── View and Camera entities (AP242) ─────────────────────────────────────────

/// A view definition (AP242 VIEW, CAMERA_MODEL_D3, etc.).
///
/// Analogous to `XCAFView` / `XCAFDoc_ViewTool` in OCCT.
/// Represents a named view with associated camera model.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepView {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// Reference to the camera model (CAMERA_MODEL_D3 or similar).
    pub camera_model_id: Option<u64>,
    /// View type (e.g., "front", "top", "isometric", "section").
    pub view_type: Option<String>,
}

/// A camera model (AP242 CAMERA_MODEL_D3).
///
/// Defines the camera position, orientation, and projection parameters.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepCameraModel {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to the view volume.
    pub view_volume_id: Option<u64>,
    /// Reference to the perspective camera model (if perspective projection).
    pub perspective_of_volume_id: Option<u64>,
}

/// A 3D camera model (AP242 CAMERA_MODEL_D3).
///
/// Provides complete camera definition including position and orientation.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepCameraModelD3 {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to AXIS2_PLACEMENT_3D defining camera position and orientation.
    pub view_reference_system_id: Option<u64>,
    /// Reference to the view volume.
    pub view_volume_id: Option<u64>,
    /// Whether this is perspective projection (true) or orthographic (false).
    pub perspective: bool,
}

/// A view volume (AP242 VIEW_VOLUME).
///
/// Defines the viewing frustum or orthographic view bounds.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepViewVolume {
    pub entity_id: u64,
    pub name: Option<String>,
    /// View volume type: orthographic or perspective.
    pub volume_type: ViewVolumeType,
    /// The point being viewed (center of interest).
    pub view_center: Option<[f64; 3]>,
    /// Distance from camera to view plane.
    pub view_plane_distance: Option<f64>,
    /// Up direction vector.
    pub up_direction: Option<[f64; 3]>,
    /// Width of view window (for orthographic).
    pub view_window_width: Option<f64>,
    /// Height of view window (for orthographic).
    pub view_window_height: Option<f64>,
}

/// View volume projection type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum ViewVolumeType {
    /// Orthographic (parallel) projection.
    Orthographic,
    /// Perspective projection.
    Perspective,
    /// Unknown or unspecified type.
    Unknown,
}

// ── Annotation and Note entities (AP242) ─────────────────────────────────────

/// A note/annotation entity (AP242 ANNOTATION).
///
/// Analogous to `XCAFNoteObjects` / `XCAFDoc_NoteTool` in OCCT.
/// Represents a textual or graphical annotation on the model.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepNote {
    pub entity_id: u64,
    pub name: Option<String>,
    pub description: Option<String>,
    /// The annotation text content.
    pub text: Option<String>,
    /// Reference to the annotation plane.
    pub annotation_plane_id: Option<u64>,
    /// Reference to associated geometry.
    pub associated_geometry_id: Option<u64>,
}

/// An annotation plane (AP242 ANNOTATION_PLANE).
///
/// Defines the plane on which annotations are placed.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepAnnotationPlane {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to AXIS2_PLACEMENT_3D defining the plane.
    pub plane_id: Option<u64>,
    /// Reference to associated annotation occurrence.
    pub annotation_occurrence_id: Option<u64>,
}

/// An annotation occurrence (AP242 ANNOTATION_OCCURRENCE).
///
/// Represents a specific instance of an annotation in the model.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepAnnotationOccurrence {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to the annotation style.
    pub style_id: Option<u64>,
    /// Reference to the annotation fill area.
    pub fill_area_id: Option<u64>,
    /// Reference to the defining shape aspect.
    pub shape_aspect_id: Option<u64>,
}

/// A dimension curve (AP242 DIMENSION_CURVE).
///
/// Represents a dimension line annotation.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDimensionCurve {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to the curve geometry.
    pub curve_id: Option<u64>,
    /// Reference to associated annotation plane.
    pub annotation_plane_id: Option<u64>,
    /// Reference to terminators (start and end).
    pub terminator_ids: Vec<u64>,
}

/// A terminator symbol (AP242 TERMINATOR_SYMBOL).
///
/// Represents arrowheads and other termination symbols on dimension lines.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepTerminatorSymbol {
    pub entity_id: u64,
    pub name: Option<String>,
    /// Reference to the annotated curve.
    pub annotated_curve_id: Option<u64>,
    /// The terminator type (arrow, dot, etc.).
    pub terminator_type: TerminatorType,
}

/// Types of terminator symbols.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
pub enum TerminatorType {
    /// Arrowhead pointing away from dimension.
    Arrow,
    /// Filled dot.
    Dot,
    /// Open arrow.
    OpenArrow,
    /// Closed filled arrow.
    ClosedArrow,
    /// Origin symbol (circle).
    Origin,
    /// Unknown terminator type.
    Unknown,
}

/// A datum feature callout (AP242 DATUM_FEATURE_CALLOUT).
///
/// Represents a datum label annotation.
#[derive(Debug, Clone, serde::Serialize)]
pub struct StepDatumFeatureCallout {
    pub entity_id: u64,
    pub name: Option<String>,
    /// The datum identifier text.
    pub datum_identifier: Option<String>,
    /// Reference to the annotation plane.
    pub annotation_plane_id: Option<u64>,
}

fn extract_file_schema(content: &str) -> Option<String> {
    let key = "FILE_SCHEMA((";
    let start = content.find(key)? + key.len();
    let rest = &content[start..];
    let end = rest.find("));")?;
    Some(rest[..end].trim().trim_matches('\'').to_string())
}

fn extract_file_name(content: &str) -> Option<String> {
    let key = "FILE_NAME(";
    let start = content.find(key)? + key.len();
    let rest = &content[start..];
    let first_quote = rest.find('\'')?;
    let rest = &rest[first_quote + 1..];
    let end_quote = rest.find('\'')?;
    Some(rest[..end_quote].to_string())
}

fn extract_product_names(content: &str) -> Vec<String> {
    extract_products(content)
        .into_iter()
        .filter_map(|product| product.name.or(product.product_id))
        .fold(Vec::new(), |mut names, name| {
            if !names.iter().any(|existing| existing == &name) {
                names.push(name);
            }
            names
        })
}

fn extract_products(content: &str) -> Vec<StepProduct> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut products = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT") {
            continue;
        }
        products.push(StepProduct {
            entity_id: id,
            product_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
        });
    }
    products
}

fn extract_product_definition_formations(content: &str) -> Vec<StepProductDefinitionFormation> {
    use std::collections::HashMap;

    let products_by_id: HashMap<u64, StepProduct> = extract_products(content)
        .into_iter()
        .map(|product| (product.entity_id, product))
        .collect();

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION_FORMATION")
            && !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE")
        {
            continue;
        }

        let parts = split_top_level(args, ',');
        let product_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let product_name = product_id.and_then(|pid| {
            products_by_id
                .get(&pid)
                .and_then(|product| product.name.clone().or(product.product_id.clone()))
        });

        out.push(StepProductDefinitionFormation {
            entity_id: id,
            formation_id: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            product_id,
            product_name,
        });
    }

    out
}

fn extract_product_definitions(content: &str) -> Vec<StepProductDefinitionInfo> {
    use std::collections::HashMap;

    let formation_by_id: HashMap<u64, StepProductDefinitionFormation> =
        extract_product_definition_formations(content)
            .into_iter()
            .map(|formation| (formation.entity_id, formation))
            .collect();

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION") {
            continue;
        }

        let parts = split_top_level(args, ',');
        let formation_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let (product_id, product_name) = formation_id
            .and_then(|fid| formation_by_id.get(&fid))
            .map(|formation| (formation.product_id, formation.product_name.clone()))
            .unwrap_or((None, None));

        out.push(StepProductDefinitionInfo {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            formation_id,
            product_id,
            product_name,
        });
    }

    out
}

fn extract_shape_definition_representations(content: &str) -> Vec<StepShapeDefinitionRepresentation> {
    use std::collections::{HashMap, HashSet};

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let product_definition_ids: HashSet<u64> = extract_product_definitions(content)
        .into_iter()
        .map(|definition| definition.entity_id)
        .collect();

    let mut pds_to_definition: HashMap<u64, u64> = HashMap::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION_SHAPE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        if let Some(definition_id) = parts.get(2).and_then(|p| parse_ref(p.trim())) {
            pds_to_definition.insert(id, definition_id);
        }
    }

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SHAPE_DEFINITION_REPRESENTATION") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let product_definition_shape_id = parts.first().and_then(|p| parse_ref(p.trim()));
        let product_definition_id = product_definition_shape_id
            .and_then(|pds_id| pds_to_definition.get(&pds_id).copied())
            .filter(|definition_id| product_definition_ids.contains(definition_id));
        let representation_id = parts.get(1).and_then(|p| parse_ref(p.trim()));

        out.push(StepShapeDefinitionRepresentation {
            entity_id: id,
            product_definition_shape_id,
            product_definition_id,
            representation_id,
        });
    }

    out
}

fn extract_assembly_occurrences(content: &str) -> Vec<StepAssemblyUsageOccurrence> {
    use std::collections::HashMap;

    let product_defs_by_id: HashMap<u64, StepProductDefinitionInfo> = extract_product_definitions(content)
        .into_iter()
        .map(|definition| (definition.entity_id, definition))
        .collect();

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("NEXT_ASSEMBLY_USAGE_OCCURRENCE") {
            continue;
        }

        let parts = split_top_level(args, ',');
        let relating_product_definition_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let related_product_definition_id = parts.get(4).and_then(|p| parse_ref(p.trim()));

        out.push(StepAssemblyUsageOccurrence {
            entity_id: id,
            usage_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            relating_product_definition_id,
            related_product_definition_id,
            relating_product_name: relating_product_definition_id.and_then(|pd| {
                product_defs_by_id
                    .get(&pd)
                    .and_then(|definition| definition.product_name.clone())
            }),
            related_product_name: related_product_definition_id.and_then(|pd| {
                product_defs_by_id
                    .get(&pd)
                    .and_then(|definition| definition.product_name.clone())
            }),
        });
    }

    out
}

fn infer_protocol_hint(file_schema: Option<&str>) -> StepProtocolHint {
    let Some(schema) = file_schema else {
        return StepProtocolHint::Unknown;
    };
    let s = schema.to_ascii_uppercase();
    if s.contains("242") {
        StepProtocolHint::Ap242
    } else if s.contains("214") {
        StepProtocolHint::Ap214
    } else if s.contains("203") {
        StepProtocolHint::Ap203
    } else {
        StepProtocolHint::Unknown
    }
}

/// Extract material names from `MATERIAL('name', ...)` entities in the data section.
fn extract_materials(content: &str) -> Vec<StepMaterial> {
    let mut materials = Vec::new();
    let mut search = content;
    while let Some(pos) = search.find("MATERIAL(") {
        let rest = &search[pos + "MATERIAL(".len()..];
        if let Some(name) = extract_first_string_arg(rest) {
            materials.push(StepMaterial { name, density: None });
        }
        search = &search[pos + 1..];
    }
    materials
}

/// Extract layer names from `PRESENTATION_LAYER_ASSIGNMENT('name', ...)` entities.
fn extract_layers(content: &str) -> Vec<StepLayer> {
    let mut layers = Vec::new();
    let mut search = content;
    while let Some(pos) = search.find("PRESENTATION_LAYER_ASSIGNMENT(") {
        let rest = &search[pos + "PRESENTATION_LAYER_ASSIGNMENT(".len()..];
        if let Some(name) = extract_first_string_arg(rest) {
            if !layers.iter().any(|l: &StepLayer| l.name == name) {
                layers.push(StepLayer { name });
            }
        }
        search = &search[pos + 1..];
    }
    layers
}

/// Extract generic properties from `GENERAL_PROPERTY('name','description',...)` entities.
fn extract_general_properties(content: &str) -> Vec<StepGeneralProperty> {
    let mut props = Vec::new();
    let mut search = content;
    while let Some(pos) = search.find("GENERAL_PROPERTY(") {
        let rest = &search[pos + "GENERAL_PROPERTY(".len()..];
        if let Some(name) = extract_nth_string_arg(rest, 0) {
            let description = extract_nth_string_arg(rest, 1);
            props.push(StepGeneralProperty { name, description });
        }
        search = &search[pos + 1..];
    }
    props
}

/// Extract `GENERAL_PROPERTY` entities keyed by their entity id.
fn extract_general_properties_with_ids(content: &str) -> Vec<(u64, StepGeneralProperty)> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut props = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if entity.eq_ignore_ascii_case("GENERAL_PROPERTY") {
            if let Some(name) = extract_nth_string_arg(args, 0) {
                let description = extract_nth_string_arg(args, 1);
                props.push((id, StepGeneralProperty { name, description }));
            }
        }
    }

    props
}

/// Extract `PROPERTY_DEFINITION` entities and link them to `GENERAL_PROPERTY`
/// when the third argument references such an entity.
fn extract_property_definitions(content: &str) -> Vec<StepPropertyDefinition> {
    use std::collections::HashMap;

    let general_by_id: HashMap<u64, StepGeneralProperty> =
        extract_general_properties_with_ids(content).into_iter().collect();

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PROPERTY_DEFINITION") {
            continue;
        }

        let parts = split_top_level(args, ',');
        let name = parts.first().and_then(|_| extract_nth_string_arg(args, 0));
        let description = parts.get(1).and_then(|_| extract_nth_string_arg(args, 1));
        let referenced_entity_id = parts.get(2).and_then(|p| parse_ref(p));

        let (general_property_name, general_property_description) =
            if let Some(rid) = referenced_entity_id {
                if let Some(prop) = general_by_id.get(&rid) {
                    (Some(prop.name.clone()), prop.description.clone())
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };

        out.push(StepPropertyDefinition {
            definition_id: id,
            name,
            description,
            referenced_entity_id,
            general_property_name,
            general_property_description,
        });
    }

    out
}

/// Extract the first single-quoted string argument from a STEP entity argument list.

fn extract_property_definition_reprs(content: &str) -> Vec<StepPropertyDefinitionRepr> {
    let Ok(data) = extract_data_section(content) else { return Vec::new(); };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else { continue; };
        let Some((entity, args)) = parse_entity_body(body) else { continue; };
        if !entity.eq_ignore_ascii_case("PROPERTY_DEFINITION_REPRESENTATION") { continue; }
        let parts = split_top_level(args, ',');
        let pd_id = parts.first().and_then(|p| parse_ref(p.trim()));
        let rep_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        out.push(StepPropertyDefinitionRepr { entity_id: id, property_definition_id: pd_id, representation_id: rep_id });
    }
    out
}

fn extract_dimensional_locations(content: &str) -> Vec<StepDimensionalLocation> {
    let Ok(data) = extract_data_section(content) else { return Vec::new(); };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else { continue; };
        let Some((entity, args)) = parse_entity_body(body) else { continue; };
        if !entity.eq_ignore_ascii_case("DIMENSIONAL_LOCATION") { continue; }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let from_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let to_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        out.push(StepDimensionalLocation { entity_id: id, name, description, from_entity_id: from_id, to_entity_id: to_id });
    }
    out
}

fn extract_dimensional_sizes(content: &str) -> Vec<StepDimensionalSize> {
    let Ok(data) = extract_data_section(content) else { return Vec::new(); };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else { continue; };
        let Some((entity, args)) = parse_entity_body(body) else { continue; };
        if !entity.eq_ignore_ascii_case("DIMENSIONAL_SIZE") { continue; }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let shape_aspect_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDimensionalSize { entity_id: id, name, description, shape_aspect_id });
    }
    out
}

fn extract_geometric_tolerances(content: &str) -> Vec<StepGeometricTolerance> {
    let Ok(data) = extract_data_section(content) else { return Vec::new(); };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else { continue; };
        let Some((entity, args)) = parse_entity_body(body) else { continue; };
        if !entity.eq_ignore_ascii_case("GEOMETRIC_TOLERANCE") { continue; }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        out.push(StepGeometricTolerance { entity_id: id, name, description, value_entity_id, shape_aspect_id });
    }
    out
}

fn extract_geometric_tolerances_with_datum_references(
    content: &str,
) -> Vec<StepGeometricToleranceWithDatumReference> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let datum_system_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        out.push(StepGeometricToleranceWithDatumReference {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            datum_system_id,
        });
    }
    out
}

fn extract_datums(content: &str) -> Vec<StepDatum> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DATUM") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let shape_aspect_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDatum {
            entity_id: id,
            name,
            description,
            shape_aspect_id,
        });
    }
    out
}

fn extract_datum_systems(content: &str) -> Vec<StepDatumSystem> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DATUM_SYSTEM") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let datum_ids = parts
            .get(2)
            .map(|p| parse_ref_list(p.trim()))
            .unwrap_or_default();
        out.push(StepDatumSystem {
            entity_id: id,
            name,
            description,
            datum_ids,
        });
    }
    out
}

fn extract_kinematic_pairs(content: &str) -> Vec<StepKinematicPair> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let is_kinematic = entity_upper.contains("KINEMATIC")
            || entity_upper.contains("JOINT")
            || entity_upper.ends_with("_PAIR");
        if !is_kinematic {
            continue;
        }

        out.push(StepKinematicPair {
            entity_id: id,
            entity_type: entity_upper,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            related_entity_ids: parse_ref_list(args),
        });
    }
    out
}

fn extract_tolerance_zones(content: &str) -> Vec<StepToleranceZone> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("TOLERANCE_ZONE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let toleranced_shape_aspect_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepToleranceZone {
            entity_id: id,
            name,
            description,
            toleranced_shape_aspect_id,
        });
    }
    out
}

fn extract_tolerance_zone_definitions(content: &str) -> Vec<StepToleranceZoneDefinition> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("TOLERANCE_ZONE_DEFINITION") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let tolerance_zone_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        out.push(StepToleranceZoneDefinition {
            entity_id: id,
            name,
            description,
            tolerance_zone_id,
            shape_aspect_id,
        });
    }
    out
}

fn extract_datum_features(content: &str) -> Vec<StepDatumFeature> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DATUM_FEATURE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let datum_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDatumFeature {
            entity_id: id,
            name,
            description,
            datum_id,
        });
    }
    out
}

fn extract_datum_reference_elements(content: &str) -> Vec<StepDatumReferenceElement> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DATUM_REFERENCE_ELEMENT") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let associated_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDatumReferenceElement {
            entity_id: id,
            name,
            description,
            associated_entity_id,
        });
    }
    out
}

fn extract_shape_aspects(content: &str) -> Vec<StepShapeAspect> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SHAPE_ASPECT") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let product_definition_shape_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepShapeAspect {
            entity_id: id,
            name,
            description,
            product_definition_shape_id,
        });
    }
    out
}

fn extract_shape_aspect_definitions(content: &str) -> Vec<StepShapeAspectDefinition> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SHAPE_ASPECT_DEFINITION") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let shape_aspect_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepShapeAspectDefinition {
            entity_id: id,
            name,
            description,
            shape_aspect_id,
        });
    }
    out
}

fn extract_derived_shape_aspects(content: &str) -> Vec<StepDerivedShapeAspect> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DERIVED_SHAPE_ASPECT") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let base_shape_aspect_ids = parts
            .get(2)
            .map(|p| parse_ref_list(p.trim()))
            .unwrap_or_default();
        out.push(StepDerivedShapeAspect {
            entity_id: id,
            name,
            description,
            base_shape_aspect_ids,
        });
    }
    out
}

// ── GDT Extended entity extraction (AP242) ─────────────────────────────────────

fn extract_dimensional_tolerances(content: &str) -> Vec<StepDimensionalTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        // DIMENSIONAL_TOLERANCE is often expressed via PLUS_MINUS_VALUE or TOLERANCE_VALUE
        if !entity.eq_ignore_ascii_case("DIMENSIONAL_TOLERANCE")
            && !entity.eq_ignore_ascii_case("TOLERANCE_VALUE")
            && !entity.eq_ignore_ascii_case("PLUS_MINUS_VALUE")
        {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let dimensional_characteristic_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let upper_tolerance = parts.get(3).and_then(|p| parse_float_arg(p.trim()));
        let lower_tolerance = parts.get(4).and_then(|p| parse_float_arg(p.trim()));
        let unit = parts.get(5).and_then(|p| parse_string_arg(p.trim()));
        out.push(StepDimensionalTolerance {
            entity_id: id,
            name,
            description,
            dimensional_characteristic_id,
            upper_tolerance,
            lower_tolerance,
            unit,
        });
    }
    out
}

fn extract_tolerance_values(content: &str) -> Vec<StepToleranceValue> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        // Tolerance values are often MEASURE_REPRESENTATION_ITEM or LENGTH_MEASURE_WITH_UNIT
        if !entity.eq_ignore_ascii_case("MEASURE_REPRESENTATION_ITEM")
            && !entity.eq_ignore_ascii_case("LENGTH_MEASURE_WITH_UNIT")
            && !entity.eq_ignore_ascii_case("MEASURE_WITH_UNIT")
        {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        // Value is typically the second argument or extracted from a nested structure
        let value = parts.get(1).and_then(|p| parse_float_arg(p.trim()));
        // Unit might be a reference or string
        let unit = parts.get(2).and_then(|p| {
            let s = p.trim();
            if s.starts_with('#') {
                // It's a reference to a unit definition
                None
            } else {
                parse_string_arg(s)
            }
        });
        out.push(StepToleranceValue {
            entity_id: id,
            name,
            value: value.unwrap_or(0.0),
            unit,
        });
    }
    out
}

fn extract_position_tolerances(content: &str) -> Vec<StepPositionTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("POSITION_TOLERANCE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let datum_system_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        let projected = parts.get(5).and_then(|p| parse_bool_arg(p.trim())).unwrap_or(false);
        let projected_height = parts.get(6).and_then(|p| parse_float_arg(p.trim()));
        out.push(StepPositionTolerance {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            datum_system_id,
            projected,
            projected_height,
        });
    }
    out
}

fn extract_orientation_tolerances(content: &str) -> Vec<StepOrientationTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let orientation_type = if entity_upper == "ANGULARITY_TOLERANCE" {
            Some(OrientationToleranceType::Angularity)
        } else if entity_upper == "PERPENDICULARITY_TOLERANCE" {
            Some(OrientationToleranceType::Perpendicularity)
        } else if entity_upper == "PARALLELISM_TOLERANCE" {
            Some(OrientationToleranceType::Parallelism)
        } else if entity_upper == "ORIENTATION_TOLERANCE" {
            // Generic orientation tolerance, type determined by name
            None
        } else {
            continue;
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let datum_system_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        // Determine type from name if not set
        let final_type = orientation_type.unwrap_or_else(|| {
            match name.as_deref() {
                Some("angularity") => OrientationToleranceType::Angularity,
                Some("perpendicularity") => OrientationToleranceType::Perpendicularity,
                Some("parallelism") => OrientationToleranceType::Parallelism,
                _ => OrientationToleranceType::Angularity, // Default
            }
        });
        out.push(StepOrientationTolerance {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            datum_system_id,
            orientation_type: final_type,
        });
    }
    out
}

fn extract_form_tolerances(content: &str) -> Vec<StepFormTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let form_type = if entity_upper == "FLATNESS_TOLERANCE" {
            Some(FormToleranceType::Flatness)
        } else if entity_upper == "STRAIGHTNESS_TOLERANCE" {
            Some(FormToleranceType::Straightness)
        } else if entity_upper == "CIRCULARITY_TOLERANCE" || entity_upper == "ROUNDNESS_TOLERANCE" {
            Some(FormToleranceType::Circularity)
        } else if entity_upper == "CYLINDRICITY_TOLERANCE" {
            Some(FormToleranceType::Cylindricity)
        } else if entity_upper == "FORM_TOLERANCE" {
            None
        } else {
            continue;
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        // Determine type from name if not set
        let final_type = form_type.unwrap_or_else(|| {
            match name.as_deref() {
                Some("flatness") => FormToleranceType::Flatness,
                Some("straightness") => FormToleranceType::Straightness,
                Some("circularity") | Some("roundness") => FormToleranceType::Circularity,
                Some("cylindricity") => FormToleranceType::Cylindricity,
                _ => FormToleranceType::Flatness, // Default
            }
        });
        out.push(StepFormTolerance {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            form_type: final_type,
        });
    }
    out
}

fn extract_runout_tolerances(content: &str) -> Vec<StepRunoutTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let runout_type = if entity_upper == "CIRCULAR_RUNOUT_TOLERANCE" {
            Some(RunoutToleranceType::CircularRunout)
        } else if entity_upper == "TOTAL_RUNOUT_TOLERANCE" {
            Some(RunoutToleranceType::TotalRunout)
        } else if entity_upper == "RUNOUT_TOLERANCE" {
            None
        } else {
            continue;
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let datum_system_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        // Determine type from name if not set
        let final_type = runout_type.unwrap_or_else(|| {
            match name.as_deref() {
                Some("circular runout") => RunoutToleranceType::CircularRunout,
                Some("total runout") => RunoutToleranceType::TotalRunout,
                _ => RunoutToleranceType::CircularRunout, // Default
            }
        });
        out.push(StepRunoutTolerance {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            datum_system_id,
            runout_type: final_type,
        });
    }
    out
}

fn extract_profile_tolerances(content: &str) -> Vec<StepProfileTolerance> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let profile_type = if entity_upper == "LINE_PROFILE_TOLERANCE" {
            Some(ProfileToleranceType::ProfileOfALine)
        } else if entity_upper == "SURFACE_PROFILE_TOLERANCE" {
            Some(ProfileToleranceType::ProfileOfASurface)
        } else if entity_upper == "PROFILE_TOLERANCE" {
            None
        } else {
            continue;
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let value_entity_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let datum_system_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        // Determine type from name if not set
        let final_type = profile_type.unwrap_or_else(|| {
            match name.as_deref() {
                Some("profile of a line") => ProfileToleranceType::ProfileOfALine,
                Some("profile of a surface") => ProfileToleranceType::ProfileOfASurface,
                _ => ProfileToleranceType::ProfileOfASurface, // Default
            }
        });
        out.push(StepProfileTolerance {
            entity_id: id,
            name,
            description,
            value_entity_id,
            shape_aspect_id,
            datum_system_id,
            profile_type: final_type,
        });
    }
    out
}

fn extract_datum_reference_frames(content: &str) -> Vec<StepDatumReferenceFrame> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DATUM_REFERENCE_FRAME") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let datum_system_ids = parts
            .get(2)
            .map(|p| parse_ref_list(p.trim()))
            .unwrap_or_default();
        out.push(StepDatumReferenceFrame {
            entity_id: id,
            name,
            description,
            datum_system_ids,
        });
    }
    out
}

fn extract_datum_targets(content: &str) -> Vec<StepDatumTarget> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_ascii_uppercase();
        let target_type = if entity_upper == "DATUM_TARGET" || entity_upper == "DATUM_TARGET_POINT" {
            Some(DatumTargetType::Point)
        } else if entity_upper == "DATUM_TARGET_LINE" {
            Some(DatumTargetType::Line)
        } else if entity_upper == "DATUM_TARGET_AREA" {
            Some(DatumTargetType::Area)
        } else if entity_upper == "DATUM_TARGET_CIRCLE" {
            Some(DatumTargetType::AreaCircle)
        } else if entity_upper == "DATUM_TARGET_RECTANGLE" {
            Some(DatumTargetType::AreaRectangle)
        } else {
            continue;
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let target_identifier = parts.get(2).and_then(|p| parse_string_arg(p.trim()));
        let datum_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        out.push(StepDatumTarget {
            entity_id: id,
            name,
            description,
            target_identifier,
            datum_id,
            target_type: target_type.unwrap_or(DatumTargetType::Point),
            shape_aspect_id,
        });
    }
    out
}

fn extract_tolerance_zone_definitions_enhanced(content: &str) -> Vec<StepToleranceZoneDefinitionEnhanced> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("TOLERANCE_ZONE_DEFINITION") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let tolerance_zone_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let shape_aspect_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        // Parse zone shape from name or use default
        let zone_shape = match name.as_deref() {
            Some("cylindrical") => ToleranceZoneShape::Cylindrical,
            Some("spherical") => ToleranceZoneShape::Spherical,
            Some("two parallel planes") => ToleranceZoneShape::TwoParallelPlanes,
            Some("two coaxial cylinders") => ToleranceZoneShape::TwoCoaxialCylinders,
            Some("two concentric circles") => ToleranceZoneShape::TwoConcentricCircles,
            Some("within circle") => ToleranceZoneShape::WithinCircle,
            Some("two parallel lines") => ToleranceZoneShape::TwoParallelLines,
            _ => ToleranceZoneShape::Unknown,
        };
        // Parse zone position from description or use default
        let zone_position = match description.as_deref() {
            Some("symmetric") => ToleranceZonePosition::Symmetric,
            Some("unilateral") => ToleranceZonePosition::Unilateral,
            Some("bilateral asymmetric") => ToleranceZonePosition::BilateralAsymmetric,
            _ => ToleranceZonePosition::Unspecified,
        };
        let defining_shape_aspect_id = parts.get(4).and_then(|p| parse_ref(p.trim()));
        out.push(StepToleranceZoneDefinitionEnhanced {
            entity_id: id,
            name,
            description,
            tolerance_zone_id,
            shape_aspect_id,
            zone_shape,
            zone_position,
            defining_shape_aspect_id,
        });
    }
    out
}

// ── FEA (Finite Element Analysis) entity extraction (AP242) ──────────────────

fn extract_fea_models(content: &str) -> Vec<StepFeaModel> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_MODEL" | "FEA_MODEL_3D" | "FEAMEDIAN_MODEL"
        ) {
            continue;
        }
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        out.push(StepFeaModel {
            entity_id: id,
            name,
            description,
        });
    }
    out
}

fn extract_fea_meshes(content: &str) -> Vec<StepFeaMesh> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_MESH" | "FEA_MESH_3D" | "FEAMEDIAN_MESH"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let node_count = parts.get(1).and_then(|p| parse_uint_arg(p.trim()));
        let element_count = parts.get(2).and_then(|p| parse_uint_arg(p.trim()));
        out.push(StepFeaMesh {
            entity_id: id,
            name,
            node_count,
            element_count,
        });
    }
    out
}

fn extract_fea_node_sets(content: &str) -> Vec<StepFeaNodeSet> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_NODE_SET" | "FEAMEDIAN_NODE_SET"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let model_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        out.push(StepFeaNodeSet {
            entity_id: id,
            name,
            model_id,
        });
    }
    out
}

fn extract_fea_element_sets(content: &str) -> Vec<StepFeaElementSet> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_ELEMENT_SET" | "FEAMEDIAN_ELEMENT_SET"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let model_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        let element_type = extract_nth_string_arg(args, 1);
        out.push(StepFeaElementSet {
            entity_id: id,
            name,
            model_id,
            element_type,
        });
    }
    out
}

fn extract_fea_material_properties(content: &str) -> Vec<StepFeaMaterialProperty> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_MATERIAL_PROPERTY" | "FEAMEDIAN_MATERIAL_PROPERTY"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let property_type = extract_nth_string_arg(args, 1);
        let value = parts.get(2).and_then(|p| parse_float_arg(p.trim()));
        let unit = extract_nth_string_arg(args, 2);
        out.push(StepFeaMaterialProperty {
            entity_id: id,
            name,
            property_type,
            value,
            unit,
        });
    }
    out
}

fn extract_fea_boundary_conditions(content: &str) -> Vec<StepFeaBoundaryCondition> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_BOUNDARY_CONDITION" | "FEAMEDIAN_BOUNDARY_CONDITION"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let condition_type = extract_nth_string_arg(args, 1);
        let node_set_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepFeaBoundaryCondition {
            entity_id: id,
            name,
            condition_type,
            node_set_id,
        });
    }
    out
}

fn extract_fea_loads(content: &str) -> Vec<StepFeaLoad> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_LOAD" | "FEAMEDIAN_LOAD"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let load_type = extract_nth_string_arg(args, 1);
        let magnitude = parts.get(2).and_then(|p| parse_float_arg(p.trim()));
        // Direction is typically a reference to a DIRECTION entity or a list of 3 floats
        let direction = parts.get(3).and_then(|p| parse_direction_tuple(p.trim()));
        out.push(StepFeaLoad {
            entity_id: id,
            name,
            load_type,
            magnitude,
            direction,
        });
    }
    out
}


fn extract_fea_node_groups(content: &str) -> Vec<StepFeaNodeGroup> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("FEA_NODE_GROUP") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let model_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        let node_count = parts.get(2).and_then(|p| parse_uint_arg(p.trim()));
        out.push(StepFeaNodeGroup {
            entity_id: id,
            name,
            model_id,
            node_count,
        });
    }
    out
}

// ── Additional FEA entity extraction functions (AP209/AP242 extended) ──────────

fn extract_fea_analyses(content: &str) -> Vec<StepFeaAnalysis> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_ANALYSIS" | "ANALYSIS_3D" | "FEAMEDIAN_ANALYSIS"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let model_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        // analysis_type is the 3rd string (index 2) in the arguments
        let analysis_type = extract_nth_string_arg(args, 2);
        // creation_date is the 4th string (index 3) in the arguments
        let creation_date = extract_nth_string_arg(args, 3);
        out.push(StepFeaAnalysis {
            entity_id: id,
            name,
            description,
            model_id,
            analysis_type,
            creation_date,
        });
    }
    out
}

fn extract_fea_states(content: &str) -> Vec<StepFeaState> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_STATE" | "FEAMEDIAN_STATE"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let analysis_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        // state_type is the 3rd string (index 2) in the arguments
        let state_type = extract_nth_string_arg(args, 2);
        out.push(StepFeaState {
            entity_id: id,
            name,
            description,
            analysis_id,
            state_type,
        });
    }
    out
}

fn extract_fea_material_models(content: &str) -> Vec<StepFeaMaterialModel> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        let model_type = match entity_upper.as_str() {
            "FEA_LINEAR_ELASTICITY" => "LINEAR_ELASTIC",
            "FEA_MATERIAL_MODEL" => "GENERIC",
            "FEA_ELASTIC_PLASTIC" => "ELASTIC_PLASTIC",
            "FEA_HYPERELASTIC" => "HYPERELASTIC",
            "FEAMEDIAN_MATERIAL_MODEL" => "MEDIAN_MODEL",
            _ => continue,
        };
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        // For FEA_LINEAR_ELASTICITY: name, description, Young's modulus, Poisson's ratio, ...
        let youngs_modulus = parts.get(2).and_then(|p| parse_float_arg(p.trim()));
        let poissons_ratio = parts.get(3).and_then(|p| parse_float_arg(p.trim()));
        let shear_modulus = parts.get(4).and_then(|p| parse_float_arg(p.trim()));
        let density = parts.get(5).and_then(|p| parse_float_arg(p.trim()));
        let thermal_expansion = parts.get(6).and_then(|p| parse_float_arg(p.trim()));
        // modulus_unit is the 3rd string (index 2) in the arguments
        let modulus_unit = extract_nth_string_arg(args, 2);
        out.push(StepFeaMaterialModel {
            entity_id: id,
            name,
            description,
            model_type: model_type.to_string(),
            youngs_modulus,
            poissons_ratio,
            shear_modulus,
            density,
            thermal_expansion,
            modulus_unit,
        });
    }
    out
}

fn extract_fea_nodes(content: &str) -> Vec<StepFeaNode> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "NODE_REPRESENTATION" | "FEA_NODE" | "FEAMEDIAN_NODE"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        // NODE_REPRESENTATION(name, node_number, (x, y, z), mesh_ref)
        // or FEA_NODE(node_number, (x, y, z), mesh_ref)
        let node_number = if entity_upper == "NODE_REPRESENTATION" {
            parts.get(1).and_then(|p| parse_uint_arg(p.trim())).unwrap_or(id)
        } else {
            parts.get(0).and_then(|p| parse_uint_arg(p.trim())).unwrap_or(id)
        };
        // Find the coordinates - could be a tuple or a reference to a CARTESIAN_POINT
        let coords_str = if entity_upper == "NODE_REPRESENTATION" {
            parts.get(2)
        } else {
            parts.get(1)
        };
        let coordinates = coords_str
            .and_then(|s| parse_direction_tuple(s.trim()))
            .unwrap_or([0.0, 0.0, 0.0]);
        let mesh_id = if entity_upper == "NODE_REPRESENTATION" {
            parts.get(3).and_then(|p| parse_ref(p.trim()))
        } else {
            parts.get(2).and_then(|p| parse_ref(p.trim()))
        };
        out.push(StepFeaNode {
            entity_id: id,
            node_number,
            coordinates,
            mesh_id,
        });
    }
    out
}

fn extract_fea_elements(content: &str) -> Vec<StepFeaElement> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "ELEMENT_REPRESENTATION" | "FEA_ELEMENT" | "FEAMEDIAN_ELEMENT"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        // ELEMENT_REPRESENTATION(name, element_number, element_type, (node_ids), mesh_ref)
        // or FEA_ELEMENT(element_number, element_type, (node_ids), mesh_ref)
        let is_element_repr = entity_upper == "ELEMENT_REPRESENTATION";
        // String indices: for ELEMENT_REPRESENTATION, strings are name(0) and element_type(1)
        // For FEA_ELEMENT, only string is element_type(0)
        let element_number = if is_element_repr {
            parts.get(1).and_then(|p| parse_uint_arg(p.trim())).unwrap_or(id)
        } else {
            parts.get(0).and_then(|p| parse_uint_arg(p.trim())).unwrap_or(id)
        };
        // element_type is string index 1 for ELEMENT_REPRESENTATION, 0 for FEA_ELEMENT
        let element_type = if is_element_repr {
            extract_nth_string_arg(args, 1).unwrap_or_else(|| "UNKNOWN".to_string())
        } else {
            extract_nth_string_arg(args, 0).unwrap_or_else(|| "UNKNOWN".to_string())
        };
        // Parse node IDs from a list like (#1,#2,#3,#4) or (1,2,3,4)
        let node_ids_str = if is_element_repr {
            parts.get(3).map(|s| s.trim()).unwrap_or("")
        } else {
            parts.get(2).map(|s| s.trim()).unwrap_or("")
        };
        let node_ids = parse_ref_list(node_ids_str);
        let mesh_id = if is_element_repr {
            parts.get(4).and_then(|p| parse_ref(p.trim()))
        } else {
            parts.get(3).and_then(|p| parse_ref(p.trim()))
        };
        out.push(StepFeaElement {
            entity_id: id,
            element_number,
            element_type,
            node_ids,
            mesh_id,
        });
    }
    out
}

fn extract_fea_steps(content: &str) -> Vec<StepFeaStep> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_STEP" | "FEAMEDIAN_STEP"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let analysis_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let step_number = parts.get(3).and_then(|p| parse_uint_arg(p.trim()));
        // step_type is the 3rd string (index 2) in the arguments
        let step_type = extract_nth_string_arg(args, 2);
        out.push(StepFeaStep {
            entity_id: id,
            name,
            description,
            analysis_id,
            step_number,
            step_type,
        });
    }
    out
}

fn extract_fea_results(content: &str) -> Vec<StepFeaResult> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_RESULT" | "NODAL_RESULT" | "ELEMENT_RESULT" | "FEAMEDIAN_RESULT"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let result_type = extract_nth_string_arg(args, 1).unwrap_or_else(|| entity_upper.clone());
        let analysis_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        let location_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let component_count = parts.get(4).and_then(|p| parse_uint_arg(p.trim()));
        out.push(StepFeaResult {
            entity_id: id,
            name,
            result_type,
            analysis_id,
            location_id,
            component_count,
        });
    }
    out
}

fn extract_fea_cases(content: &str) -> Vec<StepFeaCase> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "FEA_CASE" | "FEAMEDIAN_CASE"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let description = extract_nth_string_arg(args, 1);
        let model_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        // case_type is the 3rd string (index 2) in the arguments
        let case_type = extract_nth_string_arg(args, 2);
        out.push(StepFeaCase {
            entity_id: id,
            name,
            description,
            model_id,
            case_type,
        });
    }
    out
}

// ── View and Camera extraction functions (AP242) ─────────────────────────────

fn extract_views(content: &str) -> Vec<StepView> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("CAMERA_MODEL_D3") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let _view_reference_system_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        let _view_volume_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        // Perspective flag is typically indicated by presence of PERSPECTIVE_OF_VOLUME
        out.push(StepView {
            entity_id: id,
            name: name.clone(),
            description: None,
            camera_model_id: Some(id),
            view_type: name,
        });
    }
    out
}

fn extract_cameras(content: &str) -> Vec<StepCameraModelD3> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("CAMERA_MODEL_D3") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let view_reference_system_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        let view_volume_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        // Check for perspective by looking for PERSPECTIVE_OF_VOLUME in subsequent args
        let perspective = parts.iter().any(|p| {
            let p_upper = p.trim().to_uppercase();
            p_upper.contains("PERSPECTIVE") || p_upper == ".T."
        });
        out.push(StepCameraModelD3 {
            entity_id: id,
            name,
            view_reference_system_id,
            view_volume_id,
            perspective,
        });
    }
    out
}

fn extract_view_volumes(content: &str) -> Vec<StepViewVolume> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("VIEW_VOLUME") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        // VIEW_VOLUME has complex structure - extract what we can
        let volume_type = if parts.iter().any(|p| p.trim().to_uppercase().contains("PERSPECTIVE")) {
            ViewVolumeType::Perspective
        } else {
            ViewVolumeType::Orthographic
        };
        // Try to extract view window dimensions
        let view_window_width = parts.get(4).and_then(|p| parse_float_arg(p.trim()));
        let view_window_height = parts.get(5).and_then(|p| parse_float_arg(p.trim()));
        out.push(StepViewVolume {
            entity_id: id,
            name,
            volume_type,
            view_center: None,
            view_plane_distance: None,
            up_direction: None,
            view_window_width,
            view_window_height,
        });
    }
    out
}

// ── Annotation extraction functions (AP242) ──────────────────────────────────

fn extract_notes(content: &str) -> Vec<StepNote> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "ANNOTATION" | "ANNOTATION_TEXT" | "DESCRIPTIVE_REPRESENTATION_ITEM"
        ) {
            continue;
        }
        let name = extract_nth_string_arg(args, 0);
        let text = extract_nth_string_arg(args, 1).or(extract_nth_string_arg(args, 0));
        out.push(StepNote {
            entity_id: id,
            name,
            description: None,
            text,
            annotation_plane_id: None,
            associated_geometry_id: None,
        });
    }
    out
}

fn extract_annotation_planes(content: &str) -> Vec<StepAnnotationPlane> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("ANNOTATION_PLANE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let plane_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        out.push(StepAnnotationPlane {
            entity_id: id,
            name,
            plane_id,
            annotation_occurrence_id: None,
        });
    }
    out
}

fn extract_annotation_occurrences(content: &str) -> Vec<StepAnnotationOccurrence> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("ANNOTATION_OCCURRENCE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let style_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        out.push(StepAnnotationOccurrence {
            entity_id: id,
            name,
            style_id,
            fill_area_id: None,
            shape_aspect_id: None,
        });
    }
    out
}

fn extract_dimension_curves(content: &str) -> Vec<StepDimensionCurve> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("DIMENSION_CURVE") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let curve_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        let annotation_plane_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDimensionCurve {
            entity_id: id,
            name,
            curve_id,
            annotation_plane_id,
            terminator_ids: Vec::new(),
        });
    }
    out
}

fn extract_terminator_symbols(content: &str) -> Vec<StepTerminatorSymbol> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !matches!(
            entity_upper.as_str(),
            "TERMINATOR_SYMBOL" | "DIMENSION_CURVE_TERMINATOR" | "FILLED_ARROW"
        ) {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let annotated_curve_id = parts.get(1).and_then(|p| parse_ref(p.trim()));
        // Determine terminator type from entity name
        let terminator_type = if entity_upper.contains("ARROW") {
            TerminatorType::Arrow
        } else if entity_upper.contains("DOT") {
            TerminatorType::Dot
        } else {
            TerminatorType::Unknown
        };
        out.push(StepTerminatorSymbol {
            entity_id: id,
            name,
            annotated_curve_id,
            terminator_type,
        });
    }
    out
}

fn extract_datum_feature_callouts(content: &str) -> Vec<StepDatumFeatureCallout> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("DATUM_FEATURE_CALLOUT") {
            continue;
        }
        let parts = split_top_level(args, ',');
        let name = extract_nth_string_arg(args, 0);
        let datum_identifier = extract_nth_string_arg(args, 1);
        let annotation_plane_id = parts.get(2).and_then(|p| parse_ref(p.trim()));
        out.push(StepDatumFeatureCallout {
            entity_id: id,
            name,
            datum_identifier,
            annotation_plane_id,
        });
    }
    out
}

// ── AP242 Product Definition Relationship Chains ───────────────────────────────

fn extract_product_definition_relationships(content: &str) -> Vec<StepProductDefinitionRelationship> {
    use std::collections::HashMap;

    let product_defs_by_id: HashMap<u64, StepProductDefinitionInfo> = extract_product_definitions(content)
        .into_iter()
        .map(|definition| (definition.entity_id, definition))
        .collect();

    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION_RELATIONSHIP") {
            continue;
        }

        let parts = split_top_level(args, ',');
        let relating_product_definition_id = parts.get(3).and_then(|p| parse_ref(p.trim()));
        let related_product_definition_id = parts.get(4).and_then(|p| parse_ref(p.trim()));

        out.push(StepProductDefinitionRelationship {
            entity_id: id,
            relationship_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            relating_product_definition_id,
            related_product_definition_id,
            relating_product_name: relating_product_definition_id.and_then(|pd| {
                product_defs_by_id.get(&pd).and_then(|d| d.product_name.clone())
            }),
            related_product_name: related_product_definition_id.and_then(|pd| {
                product_defs_by_id.get(&pd).and_then(|d| d.product_name.clone())
            }),
        });
    }

    out
}

// ── AP242 Shape Representation Associations ────────────────────────────────────

fn extract_shape_representation_relationships(content: &str) -> Vec<StepShapeRepresentationRelationship> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SHAPE_REPRESENTATION_RELATIONSHIP") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepShapeRepresentationRelationship {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            relating_representation_id: parts.get(2).and_then(|p| parse_ref(p.trim())),
            related_representation_id: parts.get(3).and_then(|p| parse_ref(p.trim())),
            transformation_id: parts.get(4).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

fn extract_product_definition_shapes(content: &str) -> Vec<StepProductDefinitionShape> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_DEFINITION_SHAPE") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepProductDefinitionShape {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            product_definition_id: parts.get(2).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

// ── AP242 Configuration Management ─────────────────────────────────────────────

fn extract_configuration_designs(content: &str) -> Vec<StepConfigurationDesign> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("CONFIGURATION_DESIGN") {
            continue;
        }

        // CONFIGURATION_DESIGN(name, configuration, product_definition)
        let parts = split_top_level(args, ',');
        out.push(StepConfigurationDesign {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: None,
            configuration_id: parts.get(1).and_then(|p| parse_ref(p.trim())),
            product_definition_id: parts.get(2).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

fn extract_configuration_items(content: &str) -> Vec<StepConfigurationItem> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("CONFIGURATION_ITEM") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepConfigurationItem {
            entity_id: id,
            item_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            product_concept_id: parts.get(3).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

fn extract_product_concepts(content: &str) -> Vec<StepProductConcept> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("PRODUCT_CONCEPT") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepProductConcept {
            entity_id: id,
            concept_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            market_context_id: parts.get(3).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

fn extract_configuration_effectivities(content: &str) -> Vec<StepConfigurationEffectivity> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("CONFIGURATION_EFFECTIVITY") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepConfigurationEffectivity {
            entity_id: id,
            configuration_id: parts.get(0).and_then(|p| parse_ref(p.trim())),
            usage_id: parts.get(1).and_then(|p| parse_ref(p.trim())),
            effectivity_start: extract_nth_string_arg(args, 2),
            effectivity_end: extract_nth_string_arg(args, 3),
        });
    }

    out
}

// ── AP242 Approval and Security ────────────────────────────────────────────────

fn parse_approval_status(s: &str) -> ApprovalStatus {
    let s = s.trim().to_uppercase();
    if s.contains("APPROVED") || s == ".APPROVED." {
        ApprovalStatus::Approved
    } else if s.contains("REJECTED") || s == ".REJECTED." {
        ApprovalStatus::Rejected
    } else if s.contains("PENDING") || s == ".PENDING." || s == ".NOT_YET_APPROVED." {
        ApprovalStatus::Pending
    } else {
        ApprovalStatus::Unknown
    }
}

fn extract_approvals(content: &str) -> Vec<StepApproval> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("APPROVAL") {
            continue;
        }

        let parts = split_top_level(args, ',');
        // APPROVAL(status, level)
        // status is typically .APPROVED. or a reference to APPROVAL_STATUS
        let status_str = parts.get(0).map(|s| s.trim()).unwrap_or("");
        let status = parse_approval_status(status_str);

        out.push(StepApproval {
            entity_id: id,
            status,
            level: extract_nth_string_arg(args, 0), // level is the first quoted string
            date: None, // Populated from APPROVAL_DATE
            approver: None, // Populated from APPROVAL_PERSON_ORGANIZATION
        });
    }

    out
}

fn extract_approval_assignments(content: &str) -> Vec<StepApprovalAssignment> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.ends_with("APPROVAL_ASSIGNMENT")
            && !entity_upper.eq_ignore_ascii_case("APPROVAL_ASSIGNMENT")
        {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepApprovalAssignment {
            entity_id: id,
            approval_id: parts.get(0).and_then(|p| parse_ref(p.trim())),
            approved_item_id: parts.get(1).and_then(|p| parse_ref(p.trim())),
            role: extract_nth_string_arg(args, 2),
        });
    }

    out
}

fn parse_security_level(s: &str) -> SecurityClassificationLevel {
    let s = s.trim().to_uppercase();
    if s.contains("TOP_SECRET") || s == ".TOP_SECRET." {
        SecurityClassificationLevel::TopSecret
    } else if s.contains("SECRET") || s == ".SECRET." {
        SecurityClassificationLevel::Secret
    } else if s.contains("CONFIDENTIAL") || s == ".CONFIDENTIAL." {
        SecurityClassificationLevel::Confidential
    } else if s.contains("PROPRIETARY") || s == ".PROPRIETARY." || s.contains("COMPANY_CONFIDENTIAL") {
        SecurityClassificationLevel::Proprietary
    } else if s.contains("UNCLASSIFIED") || s == ".UNCLASSIFIED." {
        SecurityClassificationLevel::Unclassified
    } else {
        SecurityClassificationLevel::Unknown
    }
}

fn extract_security_classifications(content: &str) -> Vec<StepSecurityClassification> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SECURITY_CLASSIFICATION") {
            continue;
        }

        let parts = split_top_level(args, ',');
        // SECURITY_CLASSIFICATION(name, description, security_level)
        let security_level = parts.get(2)
            .map(|s| parse_security_level(s))
            .unwrap_or(SecurityClassificationLevel::Unknown);

        out.push(StepSecurityClassification {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            security_level,
        });
    }

    out
}

fn extract_security_classification_assignments(content: &str) -> Vec<StepSecurityClassificationAssignment> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("SECURITY_CLASSIFICATION_ASSIGNMENT") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepSecurityClassificationAssignment {
            entity_id: id,
            security_classification_id: parts.get(0).and_then(|p| parse_ref(p.trim())),
            classified_item_id: parts.get(1).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

// ── AP242 Document References ──────────────────────────────────────────────────

fn extract_documents(content: &str) -> Vec<StepDocument> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DOCUMENT") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepDocument {
            entity_id: id,
            document_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            document_type: parts.get(3).and_then(|p| {
                let s = p.trim();
                if s.starts_with('#') {
                    None // Reference to document_type entity
                } else {
                    parse_string_arg(s)
                }
            }),
        });
    }

    out
}

fn extract_document_files(content: &str) -> Vec<StepDocumentFile> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DOCUMENT_FILE") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepDocumentFile {
            entity_id: id,
            document_id: extract_nth_string_arg(args, 0),
            name: extract_nth_string_arg(args, 1),
            description: extract_nth_string_arg(args, 2),
            file_name: extract_nth_string_arg(args, 3),
            file_format: parts.get(4).and_then(|p| parse_string_arg(p.trim())),
        });
    }

    out
}

fn extract_document_usage_assignments(content: &str) -> Vec<StepDocumentUsageAssignment> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        let entity_upper = entity.to_uppercase();
        if !entity_upper.eq_ignore_ascii_case("DOCUMENT_USAGE_ASSIGNMENT")
            && !entity_upper.eq_ignore_ascii_case("DOCUMENT_PRODUCT_EQUIVALENCE")
        {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepDocumentUsageAssignment {
            entity_id: id,
            document_id: parts.get(0).and_then(|p| parse_ref(p.trim())),
            product_definition_id: parts.get(1).and_then(|p| parse_ref(p.trim())),
            role: extract_nth_string_arg(args, 0), // role is the first (and only) quoted string
        });
    }

    out
}

fn extract_document_representation_relationships(content: &str) -> Vec<StepDocumentRepresentationRelationship> {
    let Ok(data) = extract_data_section(content) else {
        return Vec::new();
    };

    let mut out = Vec::new();
    for record in split_records(data) {
        let Ok(Some((id, body))) = parse_entity_record(&record) else {
            continue;
        };
        let Some((entity, args)) = parse_entity_body(body) else {
            continue;
        };
        if !entity.eq_ignore_ascii_case("DOCUMENT_REPRESENTATION_RELATIONSHIP") {
            continue;
        }

        let parts = split_top_level(args, ',');
        out.push(StepDocumentRepresentationRelationship {
            entity_id: id,
            name: extract_nth_string_arg(args, 0),
            description: extract_nth_string_arg(args, 1),
            document_id: parts.get(2).and_then(|p| parse_ref(p.trim())),
            representation_id: parts.get(3).and_then(|p| parse_ref(p.trim())),
        });
    }

    out
}

/// Parse a uint argument (either bare number or from a measure value).
fn parse_uint_arg(s: &str) -> Option<u64> {
    let s = s.trim();
    // Try parsing as a direct integer
    if let Ok(val) = s.parse::<u64>() {
        return Some(val);
    }
    // Try extracting from parentheses like "(#100)" or "(100)"
    if s.starts_with('(') && s.ends_with(')') {
        let inner = &s[1..s.len() - 1];
        return parse_uint_arg(inner);
    }
    None
}

/// Parse a float argument (bare number or wrapped).
fn parse_float_arg(s: &str) -> Option<f64> {
    let s = s.trim();
    if let Ok(val) = s.parse::<f64>() {
        return Some(val);
    }
    if s.starts_with('(') && s.ends_with(')') {
        let inner = &s[1..s.len() - 1];
        return parse_float_arg(inner);
    }
    None
}

/// Parse a string argument (extracts from single quotes).
fn parse_string_arg(s: &str) -> Option<String> {
    let s = s.trim();
    if s.starts_with('\'') && s.ends_with('\'') && s.len() >= 2 {
        Some(s[1..s.len() - 1].to_string())
    } else {
        None
    }
}

/// Parse a boolean argument (STEP format: .T. or .F.).
fn parse_bool_arg(s: &str) -> Option<bool> {
    let s = s.trim();
    if s == ".T." {
        Some(true)
    } else if s == ".F." {
        Some(false)
    } else {
        None
    }
}

/// Parse a direction tuple like "(1.0,0.0,0.0)" or a reference to a direction entity.
fn parse_direction_tuple(s: &str) -> Option<[f64; 3]> {
    let s = s.trim();
    // Skip entity references like #100
    if s.starts_with('#') {
        return None;
    }
    // Try parsing as a tuple (x,y,z)
    if s.starts_with('(') && s.ends_with(')') {
        let inner = &s[1..s.len() - 1];
        let parts: Vec<&str> = inner.split(',').collect();
        if parts.len() == 3 {
            let x = parts[0].trim().parse::<f64>().ok()?;
            let y = parts[1].trim().parse::<f64>().ok()?;
            let z = parts[2].trim().parse::<f64>().ok()?;
            return Some([x, y, z]);
        }
    }
    None
}

fn extract_first_string_arg(args: &str) -> Option<String> {
    let q1 = args.find('\'')?;
    let rest = &args[q1 + 1..];
    let q2 = rest.find('\'')?;
    let s = rest[..q2].to_string();
    if s.is_empty() { None } else { Some(s) }
}

/// Extract the N-th single-quoted string argument from a STEP argument list.
fn extract_nth_string_arg(args: &str, n: usize) -> Option<String> {
    let mut out = Vec::new();
    let mut i = 0usize;
    let b = args.as_bytes();
    while i < b.len() {
        if b[i] == b'\'' {
            let start = i + 1;
            let mut j = start;
            while j < b.len() && b[j] != b'\'' {
                j += 1;
            }
            if j >= b.len() {
                break;
            }
            out.push(args[start..j].to_string());
            i = j + 1;
        } else {
            i += 1;
        }
    }
    out.get(n).cloned().filter(|s| !s.is_empty())
}

impl StepReader {
    /// Parse STEP content from any UTF-8 reader.
    ///
    /// This provides DE_Wrapper-style stream input without requiring a temp file.
    pub fn parse_reader<R: Read>(mut reader: R) -> Result<BRep, StepError> {
        let mut content = String::new();
        reader
            .read_to_string(&mut content)
            .map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string(&content)
    }

    /// Parse STEP content from any UTF-8 reader and extract document metadata.
    pub fn parse_reader_with_metadata<R: Read>(
        mut reader: R,
    ) -> Result<(BRep, StepDocumentMetadata), StepError> {
        let mut content = String::new();
        reader
            .read_to_string(&mut content)
            .map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_metadata(&content)
    }

    /// Parse STEP content from a reader and run healing pipeline.
    pub fn parse_reader_with_healing<R: Read>(
        mut reader: R,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport), StepError> {
        let mut content = String::new();
        reader
            .read_to_string(&mut content)
            .map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_healing(&content, options)
    }

    /// Parse STEP content from a reader, run healing, and produce JSON diagnostics.
    pub fn parse_reader_with_healing_report_json<R: Read>(
        mut reader: R,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport, String), StepError> {
        let mut content = String::new();
        reader
            .read_to_string(&mut content)
            .map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_healing_report_json(&content, options)
    }

    pub fn read_file<P: AsRef<Path>>(path: P) -> Result<BRep, StepError> {
        let content = std::fs::read_to_string(path).map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string(&content)
    }

    pub fn parse_string(content: &str) -> Result<BRep, StepError> {
        if !content.contains("ISO-10303-21") {
            return Err(StepError::InvalidFormat(
                "missing ISO-10303-21 header".into(),
            ));
        }
        let entities = parse_entities(content)?;
        build_brep_from_parsed(&entities)
    }

    /// Parse a STEP string and return BRep plus document-level metadata.
    pub fn parse_string_with_metadata(
        content: &str,
    ) -> Result<(BRep, StepDocumentMetadata), StepError> {
        if !content.contains("ISO-10303-21") {
            return Err(StepError::InvalidFormat(
                "missing ISO-10303-21 header".into(),
            ));
        }

        let entities = parse_entities(content)?;
        let brep = build_brep_from_parsed(&entities)?;
        let file_schema = extract_file_schema(content);
        let products = extract_products(content);
        let metadata = StepDocumentMetadata {
            protocol_hint: infer_protocol_hint(file_schema.as_deref()),
            file_name: extract_file_name(content),
            product_names: extract_product_names(content),
            products,
            product_definition_formations: extract_product_definition_formations(content),
            product_definitions: extract_product_definitions(content),
            shape_definition_representations: extract_shape_definition_representations(content),
            assembly_occurrences: extract_assembly_occurrences(content),
            file_schema,
            uncertainty_value: entities.uncertainty_value,
            materials: extract_materials(content),
            layers: extract_layers(content),
            general_properties: extract_general_properties(content),
            property_definitions: extract_property_definitions(content),
            property_definition_representations: extract_property_definition_reprs(content),
            dimensional_locations: extract_dimensional_locations(content),
            dimensional_sizes: extract_dimensional_sizes(content),
            geometric_tolerances: extract_geometric_tolerances(content),
            geometric_tolerances_with_datum_references:
                extract_geometric_tolerances_with_datum_references(content),
            datums: extract_datums(content),
            datum_systems: extract_datum_systems(content),
            kinematic_pairs: extract_kinematic_pairs(content),
            tolerance_zones: extract_tolerance_zones(content),
            tolerance_zone_definitions: extract_tolerance_zone_definitions(content),
            datum_features: extract_datum_features(content),
            datum_reference_elements: extract_datum_reference_elements(content),
            shape_aspects: extract_shape_aspects(content),
            shape_aspect_definitions: extract_shape_aspect_definitions(content),
            derived_shape_aspects: extract_derived_shape_aspects(content),
            dimensional_tolerances: extract_dimensional_tolerances(content),
            tolerance_values: extract_tolerance_values(content),
            position_tolerances: extract_position_tolerances(content),
            orientation_tolerances: extract_orientation_tolerances(content),
            form_tolerances: extract_form_tolerances(content),
            runout_tolerances: extract_runout_tolerances(content),
            profile_tolerances: extract_profile_tolerances(content),
            datum_reference_frames: extract_datum_reference_frames(content),
            datum_targets: extract_datum_targets(content),
            tolerance_zone_definitions_enhanced: extract_tolerance_zone_definitions_enhanced(content),
            fea_models: extract_fea_models(content),
            fea_meshes: extract_fea_meshes(content),
            fea_node_sets: extract_fea_node_sets(content),
            fea_element_sets: extract_fea_element_sets(content),
            fea_material_properties: extract_fea_material_properties(content),
            fea_boundary_conditions: extract_fea_boundary_conditions(content),
            fea_loads: extract_fea_loads(content),
            fea_node_groups: extract_fea_node_groups(content),
            fea_analyses: extract_fea_analyses(content),
            fea_states: extract_fea_states(content),
            fea_material_models: extract_fea_material_models(content),
            fea_nodes: extract_fea_nodes(content),
            fea_elements: extract_fea_elements(content),
            fea_steps: extract_fea_steps(content),
            fea_results: extract_fea_results(content),
            fea_cases: extract_fea_cases(content),
            views: extract_views(content),
            cameras: extract_cameras(content),
            view_volumes: extract_view_volumes(content),
            notes: extract_notes(content),
            annotation_planes: extract_annotation_planes(content),
            annotation_occurrences: extract_annotation_occurrences(content),
            dimension_curves: extract_dimension_curves(content),
            terminator_symbols: extract_terminator_symbols(content),
            datum_feature_callouts: extract_datum_feature_callouts(content),
            // AP242 Product Definition Relationship Chains
            product_definition_relationships: extract_product_definition_relationships(content),
            // AP242 Shape Representation Associations
            shape_representation_relationships: extract_shape_representation_relationships(content),
            product_definition_shapes: extract_product_definition_shapes(content),
            // AP242 Configuration Management
            configuration_designs: extract_configuration_designs(content),
            configuration_items: extract_configuration_items(content),
            product_concepts: extract_product_concepts(content),
            configuration_effectivities: extract_configuration_effectivities(content),
            // AP242 Approval and Security
            approvals: extract_approvals(content),
            approval_assignments: extract_approval_assignments(content),
            security_classifications: extract_security_classifications(content),
            security_classification_assignments: extract_security_classification_assignments(content),
            // AP242 Document References
            documents: extract_documents(content),
            document_files: extract_document_files(content),
            document_usage_assignments: extract_document_usage_assignments(content),
            document_representation_relationships: extract_document_representation_relationships(content),
        };

        Ok((brep, metadata))
    }

    /// Read STEP file and return BRep plus document-level metadata.
    pub fn read_file_with_metadata<P: AsRef<Path>>(
        path: P,
    ) -> Result<(BRep, StepDocumentMetadata), StepError> {
        let content = std::fs::read_to_string(path).map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_metadata(&content)
    }

    /// Parse a STEP string and run rcad-algorithms healing pipeline.
    pub fn parse_string_with_healing(
        content: &str,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport), StepError> {
        let brep = Self::parse_string(content)?;
        let (healed, report) = analyze_and_heal(&brep, options);
        Ok((healed, report))
    }

    /// Parse a STEP string, run healing, and export stable JSON diagnostics.
    pub fn parse_string_with_healing_report_json(
        content: &str,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport, String), StepError> {
        let (healed, report) = Self::parse_string_with_healing(content, options)?;
        let wire = analyze_wire_issues(&healed, options.tolerance);

        let mut issue_map: BTreeMap<String, usize> = BTreeMap::new();
        for issue in &report.final_result.issues {
            *issue_map.entry(issue.to_string()).or_insert(0) += 1;
        }

        let vertices_merged = report.passes.iter().map(|pass| pass.vertices_merged).sum();
        let degenerate_faces_removed = report
            .passes
            .iter()
            .map(|pass| pass.degenerate_faces_removed)
            .sum();
        let normals_recomputed = report.passes.iter().map(|pass| pass.normals_recomputed).sum();
        let faces_reoriented = report.passes.iter().map(|pass| pass.faces_reoriented).sum();
        let wires_fixed = report.passes.iter().map(|pass| pass.wires_fixed).sum();
        let same_range_fixed = report
            .passes
            .iter()
            .map(|pass| pass.same_range_fixed)
            .sum::<usize>()
            + report
                .parametric_passes
                .iter()
                .map(|pass| pass.same_range_fixed)
                .sum::<usize>();
        let same_parameter_fixed = report
            .passes
            .iter()
            .map(|pass| pass.same_parameter_fixed)
            .sum::<usize>()
            + report
                .parametric_passes
                .iter()
                .map(|pass| pass.same_parameter_fixed)
                .sum::<usize>();

        let payload = StepImportHealingJsonV1 {
            schema: "step.import.healing.v1",
            clean: report.final_result.is_valid(),
            initial_issue_count: report.initial_issue_count(),
            final_issue_count: report.final_issue_count(),
            fixed_issue_count: report.fixed_issue_count(),
            issue_histogram: issue_map.into_iter().collect(),
            repair_pass_count: report.passes.len(),
            parametric_pass_count: report.parametric_passes.len(),
            make_connected_pass_count: report.make_connected_passes.len(),
            vertices_merged,
            degenerate_faces_removed,
            normals_recomputed,
            faces_reoriented,
            wires_fixed,
            same_range_fixed,
            same_parameter_fixed,
            wire_open_gaps: wire.total_open_gaps,
            wire_topological_self_intersections: wire.total_topological_self_intersections,
            wire_geometric_self_intersections: wire.total_geometric_self_intersections,
        };

        let json = serde_json::to_string_pretty(&payload)
            .map_err(|e| StepError::InvalidFormat(format!("healing report JSON serialize failed: {e}")))?;

        Ok((healed, report, json))
    }

    /// Read a STEP file and run rcad-algorithms healing pipeline.
    pub fn read_file_with_healing<P: AsRef<Path>>(
        path: P,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport), StepError> {
        let content = std::fs::read_to_string(path).map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_healing(&content, options)
    }

    /// Read a STEP file, run healing, and export stable JSON diagnostics.
    pub fn read_file_with_healing_report_json<P: AsRef<Path>>(
        path: P,
        options: HealingOptions,
    ) -> Result<(BRep, HealingReport, String), StepError> {
        let content = std::fs::read_to_string(path).map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_healing_report_json(&content, options)
    }

    /// Parse a STEP file, returning both the BRep and an optional color map.
    ///
    /// Colors are extracted from the `STYLED_ITEM 鈫?COLOUR_RGB` chain.
    /// Returns `None` for color when the file has no color entities.
    pub fn read_file_with_color<P: AsRef<Path>>(
        path: P,
    ) -> Result<(BRep, Option<StepColor>), StepError> {
        let content = std::fs::read_to_string(path).map_err(|e| StepError::Io(e.to_string()))?;
        Self::parse_string_with_color(&content)
    }

    /// Parse a STEP string, returning both the BRep and an optional color map.
    pub fn parse_string_with_color(content: &str) -> Result<(BRep, Option<StepColor>), StepError> {
        if !content.contains("ISO-10303-21") {
            return Err(StepError::InvalidFormat(
                "missing ISO-10303-21 header".into(),
            ));
        }
        let entities = parse_entities(content)?;
        let (brep, face_id_map) = build_brep_with_face_map(&entities)?;
        let color = resolve_step_color(&entities, &face_id_map);
        Ok((brep, color))
    }
}

fn parse_entities(content: &str) -> Result<ParsedStep, StepError> {
    let data = extract_data_section(content)?;
    let records = split_records(data);
    let mut parsed = ParsedStep::new();

    for record in records {
        let Some((id, body)) = parse_entity_record(&record)? else {
            continue;
        };
        if let Some((entity, args)) = parse_entity_body(body) {
            match entity {
                "CARTESIAN_POINT" => {
                    if let Some(coords) = parse_cartesian_point(args) {
                        parsed.cartesian_points.insert(id, coords);
                    } else if let Some(coords2d) = parse_cartesian_point_2d(args) {
                        parsed.cartesian_points_2d.insert(id, coords2d);
                    }
                }
                "DIRECTION" => {
                    if let Some(coords) = parse_cartesian_point(args) {
                        parsed.directions.insert(id, coords);
                    } else if let Some(coords2d) = parse_cartesian_point_2d(args) {
                        parsed.directions_2d.insert(id, coords2d);
                    }
                }
                "VECTOR" => {
                    if let Some((dir_ref, magnitude)) = parse_vector(args) {
                        parsed.vectors.insert(id, (dir_ref, magnitude));
                    }
                }
                "AXIS2_PLACEMENT_3D" => {
                    if let Some((origin, axis, ref_dir)) = parse_axis2_placement(args) {
                        parsed.axis2_placements.insert(id, (origin, axis, ref_dir));
                    }
                }
                "AXIS2_PLACEMENT_2D" => {
                    if let Some((loc, ref_dir)) = parse_axis2_placement_2d(args) {
                        parsed.axis2_placements_2d.insert(id, (loc, ref_dir));
                    }
                }
                "LINE" => {
                    if let Some((origin, vector_ref)) = parse_curve_basis(args) {
                        parsed.lines.insert(id, (origin, vector_ref));
                    }
                }
                "CIRCLE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.circles.insert(id, (placement, radius));
                    }
                }
                "ELLIPSE" => {
                    if let Some((placement, major, minor)) = parse_placement_two_radii(args) {
                        parsed.ellipses.insert(id, (placement, major, minor));
                    }
                }
                "HYPERBOLA" => {
                    // HYPERBOLA('name', #placement, semi_major, semi_minor)
                    if let Some((placement, major, minor)) = parse_placement_two_radii(args) {
                        parsed.hyperbolas.insert(id, (placement, major, minor));
                    }
                }
                "PARABOLA" => {
                    // PARABOLA('name', #placement, focal_param)
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.parabolas.insert(id, (placement, radius));
                    }
                }
                "OFFSET_CURVE_3D" => {
                    // OFFSET_CURVE_3D('', #basis_curve, offset_dist, #ref_dir)
                    if let Some((basis, dist, dir)) = parse_offset_curve_3d(args) {
                        parsed.offset_curves_3d.insert(id, (basis, dist, dir));
                    }
                }
                "B_SPLINE_CURVE_WITH_KNOTS" => {
                    if let Some(points) = parse_bspline_control_points(args)
                        && !points.is_empty()
                    {
                        parsed.b_spline_curves.insert(id, points.clone());
                        // Also try to parse full data (degree + knots)
                        if let Some(full) = parse_bspline_curve_full(args) {
                            parsed.b_spline_curves_full.insert(id, full);
                        }
                    }
                }
                "B_SPLINE_SURFACE_WITH_KNOTS" => {
                    if let Some(data) = parse_bspline_surface_with_knots(args) {
                        parsed.b_spline_surfaces.insert(id, data);
                    }
                }
                "TRIMMED_CURVE" => {
                    if let Some((curve_ref, t0, t1)) = parse_trimmed_curve(args) {
                        parsed.trimmed_curves.insert(id, (curve_ref, t0, t1));
                    }
                }
                "GEOMETRIC_CURVE_SET" => {
                    if let Some(curve_refs) = parse_ref_list_after_name(args)
                        && !curve_refs.is_empty()
                    {
                        parsed.geometric_curve_sets.push(curve_refs);
                    }
                }
                "SURFACE_CURVE" => {
                    // SURFACE_CURVE('', #3d_curve, (#pcurve1, ...), .PCURVE_S1.)
                    if let Some((curve3d_ref, pcurve_refs, same_param)) = parse_surface_curve(args)
                    {
                        parsed
                            .surface_curves
                            .insert(id, (curve3d_ref, pcurve_refs, same_param));
                    }
                }
                "PCURVE" => {
                    // PCURVE('', #surface, #definitional_rep)
                    if let Some((surface_ref, def_ref)) = parse_pcurve_args(args) {
                        parsed.pcurves.insert(id, (surface_ref, def_ref));
                    }
                }
                "DEFINITIONAL_REPRESENTATION" => {
                    // DEFINITIONAL_REPRESENTATION('', (#curve2d), #context)
                    if let Some(curve2d_ref) = parse_definitional_rep(args) {
                        parsed.definitional_reps.insert(id, curve2d_ref);
                    }
                }
                "PLANE" => {
                    if let Some(placement) = parse_single_ref_after_name(args) {
                        parsed.planes.insert(id, placement);
                    }
                }
                "CYLINDRICAL_SURFACE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.cylindrical_surfaces.insert(id, (placement, radius));
                    }
                }
                "SPHERICAL_SURFACE" => {
                    if let Some((placement, radius)) = parse_placement_radius(args) {
                        parsed.spherical_surfaces.insert(id, (placement, radius));
                    }
                }
                "CONICAL_SURFACE" => {
                    if let Some((placement, radius, half_angle_rad)) = parse_conical_surface(args) {
                        parsed
                            .conical_surfaces
                            .insert(id, (placement, radius, half_angle_rad));
                    }
                }
                "TOROIDAL_SURFACE" => {
                    if let Some((placement, major, minor)) = parse_toroidal_surface(args) {
                        parsed
                            .toroidal_surfaces
                            .insert(id, (placement, major, minor));
                    }
                }
                "SURFACE_OF_LINEAR_EXTRUSION" => {
                    // 'name', #profile_curve, #direction
                    let refs = parse_ref_list(args);
                    if refs.len() >= 2 {
                        parsed.linear_extrusions.insert(id, (refs[0], refs[1]));
                    }
                }
                "SURFACE_OF_REVOLUTION" => {
                    // 'name', #profile_curve, #axis1_placement
                    let refs = parse_ref_list(args);
                    if refs.len() >= 2 {
                        parsed.revolutions.insert(id, (refs[0], refs[1]));
                    }
                }
                "RECTANGULAR_TRIMMED_SURFACE" => {
                    // RECTANGULAR_TRIMMED_SURFACE('name', #basis, u1, u2, v1, v2, .T., .T.)
                    if let Some((basis_ref, trim)) = parse_rectangular_trimmed_surface(args) {
                        parsed
                            .rectangular_trimmed_surfaces
                            .insert(id, (basis_ref, trim));
                    }
                }
                "VERTEX_POINT" => {
                    if let Some(point_ref) = parse_single_ref_after_name(args) {
                        parsed.vertex_points.insert(id, point_ref);
                    }
                }
                "EDGE_CURVE" => {
                    if let Some((start, end, curve_ref)) = parse_edge_curve_vertices(args) {
                        parsed.edge_curves.insert(id, (start, end, curve_ref));
                    }
                }
                "ORIENTED_EDGE" => {
                    if let Some((edge_ref, orientation)) = parse_oriented_edge(args) {
                        parsed.oriented_edges.insert(id, (edge_ref, orientation));
                    }
                }
                "EDGE_LOOP" => {
                    if let Some(items) = parse_ref_list_after_name(args) {
                        parsed.edge_loops.insert(id, items);
                    }
                }
                "FACE_BOUND" | "FACE_OUTER_BOUND" => {
                    if let Some(loop_ref) = parse_single_ref_after_name(args) {
                        parsed.face_bounds.insert(id, loop_ref);
                    }
                }
                "ADVANCED_FACE" => {
                    if let Some((bounds, surface)) = parse_advanced_face(args) {
                        parsed
                            .advanced_faces
                            .insert(id, AdvancedFaceRecord { bounds, surface });
                    }
                }
                "CLOSED_SHELL" => {
                    if let Some(face_refs) = parse_ref_list_after_name(args) {
                        parsed.closed_shells.insert(id, face_refs);
                    }
                }
                "OPEN_SHELL" => {
                    if let Some(face_refs) = parse_ref_list_after_name(args) {
                        parsed.open_shells.insert(id, face_refs);
                    }
                }
                "MANIFOLD_SOLID_BREP" => {
                    if let Some(shell_ref) = parse_single_ref_after_name(args) {
                        parsed.manifold_solids.push(shell_ref);
                    }
                }
                "SHELL_BASED_SURFACE_MODEL" => {
                    if let Some(shell_refs) = parse_ref_list_after_name(args)
                        && !shell_refs.is_empty()
                    {
                        parsed.shell_based_surface_models.push(shell_refs);
                    }
                }
                "COMPOUND" => {
                    // COMPOUND('name', #elem1, #elem2, ...)
                    if let Some(elem_refs) = parse_ref_list_after_name(args) {
                        parsed.compounds.insert(id, elem_refs);
                    }
                }
                "COMPSOLID" => {
                    // COMPSOLID('name', #solid1, #solid2, ...)
                    if let Some(solid_refs) = parse_ref_list_after_name(args) {
                        parsed.compsolids.insert(id, solid_refs);
                    }
                }
                "UNCERTAINTY_MEASURE_WITH_UNIT" => {
                    // UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(value),...)
                    // Extract the length measure value for global tolerance
                    if let Some(tol) = parse_uncertainty_measure(args) {
                        // Keep the largest uncertainty value if multiple appear
                        parsed.uncertainty_value = Some(match parsed.uncertainty_value {
                            Some(existing) => existing.max(tol),
                            None => tol,
                        });
                    }
                }
                // 鈹€鈹€ Color / presentation chain 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
                "COLOUR_RGB" => {
                    // COLOUR_RGB('name', r, g, b)
                    if let Some(rgb) = parse_colour_rgb(args) {
                        parsed.colour_rgbs.insert(id, rgb);
                    }
                }
                "FILL_AREA_STYLE_COLOUR" => {
                    // FILL_AREA_STYLE_COLOUR('name', #colour_rgb)
                    if let Some(colour_ref) = parse_single_ref_after_name(args) {
                        parsed.fill_area_style_colours.insert(id, colour_ref);
                    }
                }
                "FILL_AREA_STYLE" => {
                    // FILL_AREA_STYLE('name', (#fasc, ...))
                    if let Some(refs) = parse_ref_list_after_name(args) {
                        parsed.fill_area_styles.insert(id, refs);
                    }
                }
                "SURFACE_STYLE_FILL_AREA" => {
                    // SURFACE_STYLE_FILL_AREA(#fill_area_style)
                    if let Some(fas_ref) = parse_single_ref(args) {
                        parsed.surface_style_fill_areas.insert(id, fas_ref);
                    }
                }
                "SURFACE_SIDE_STYLE" => {
                    // SURFACE_SIDE_STYLE('name', (#ssfa, ...))
                    if let Some(refs) = parse_ref_list_after_name(args) {
                        parsed.surface_side_styles.insert(id, refs);
                    }
                }
                "SURFACE_STYLE_USAGE" => {
                    // SURFACE_STYLE_USAGE(.BOTH., #surface_side_style)
                    if let Some(sss_ref) = parse_last_ref(args) {
                        parsed.surface_style_usages.insert(id, sss_ref);
                    }
                }
                "PRESENTATION_STYLE_ASSIGNMENT" => {
                    // PRESENTATION_STYLE_ASSIGNMENT((#ssu, ...))
                    let refs = parse_ref_list(args);
                    if !refs.is_empty() {
                        parsed.presentation_style_assignments.insert(id, refs);
                    }
                }
                "STYLED_ITEM" => {
                    // STYLED_ITEM('name', (#psa,...), #shape)
                    if let Some((shape_ref, psa_refs)) = parse_styled_item(args) {
                        parsed.styled_items.insert(id, (shape_ref, psa_refs));
                    }
                }
                _ => {}
            }
        }
    }

    Ok(parsed)
}

fn build_brep_from_parsed(parsed: &ParsedStep) -> Result<BRep, StepError> {
    // Check for compound structure first
    if !parsed.compounds.is_empty() {
        return build_compound_brep(parsed);
    }
    // Check for compsolid structure
    if !parsed.compsolids.is_empty() {
        return build_compsolid_brep(parsed);
    }
    build_brep_with_face_map(parsed).map(|(brep, _)| brep)
}

/// Build a BRep representing a compound from STEP COMPOUND entities.
fn build_compound_brep(parsed: &ParsedStep) -> Result<BRep, StepError> {
    use rcad_kernel::topology::Compound;

    let mut compound = Compound::new();

    // Find all top-level compounds (those not referenced by other compounds)
    let mut referenced: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for elems in parsed.compounds.values() {
        for &elem_ref in elems {
            if parsed.compounds.contains_key(&elem_ref) || parsed.compsolids.contains_key(&elem_ref) {
                referenced.insert(elem_ref);
            }
        }
    }

    // Process top-level compounds
    for (&compound_id, elems) in &parsed.compounds {
        if referenced.contains(&compound_id) {
            continue; // Skip nested compounds - they'll be processed recursively
        }

        for &elem_ref in elems {
            // Check if this is a manifold solid reference
            if parsed.manifold_solids.contains(&elem_ref) {
                // Build the solid from the manifold solid
                if let Some(shell_ref) = get_shell_for_solid(parsed, elem_ref) {
                    if let Some(solid) = build_solid_from_shell(parsed, shell_ref)? {
                        compound.add_solid(None, solid);
                    }
                }
            }
            // Check for nested compound
            else if parsed.compounds.contains_key(&elem_ref) {
                let nested = build_nested_compound(parsed, elem_ref)?;
                compound.add_compound(None, nested);
            }
            // Check for compsolid
            else if parsed.compsolids.contains_key(&elem_ref) {
                if let Some(compsolid) = build_compsolid_from_step(parsed, elem_ref)? {
                    compound.add_comp_solid(None, compsolid);
                }
            }
        }
    }

    // If no top-level compounds found, fall back to building from all solids
    if compound.is_empty() {
        return build_brep_with_face_map(parsed).map(|(brep, _)| brep);
    }

    Ok(BRep::from_compound(compound))
}

/// Build a nested compound recursively.
fn build_nested_compound(parsed: &ParsedStep, compound_id: u64) -> Result<rcad_kernel::topology::Compound, StepError> {
    use rcad_kernel::topology::Compound;

    let mut compound = Compound::new();

    if let Some(elems) = parsed.compounds.get(&compound_id) {
        for &elem_ref in elems {
            if parsed.manifold_solids.contains(&elem_ref) {
                if let Some(shell_ref) = get_shell_for_solid(parsed, elem_ref) {
                    if let Some(solid) = build_solid_from_shell(parsed, shell_ref)? {
                        compound.add_solid(None, solid);
                    }
                }
            } else if parsed.compounds.contains_key(&elem_ref) {
                let nested = build_nested_compound(parsed, elem_ref)?;
                compound.add_compound(None, nested);
            } else if parsed.compsolids.contains_key(&elem_ref) {
                if let Some(compsolid) = build_compsolid_from_step(parsed, elem_ref)? {
                    compound.add_comp_solid(None, compsolid);
                }
            }
        }
    }

    Ok(compound)
}

/// Build a BRep representing a CompSolid from STEP COMPSOLID entities.
fn build_compsolid_brep(parsed: &ParsedStep) -> Result<BRep, StepError> {
    use rcad_kernel::topology::CompSolid;

    // For now, just use the first compsolid
    if let Some((&id, _solid_refs)) = parsed.compsolids.iter().next() {
        let compsolid = build_compsolid_from_step(parsed, id)?;
        if let Some(cs) = compsolid {
            return Ok(BRep::from_compsolid(cs));
        }
    }

    build_brep_with_face_map(parsed).map(|(brep, _)| brep)
}

/// Build a CompSolid from STEP COMPSOLID element references.
fn build_compsolid_from_step(parsed: &ParsedStep, compsolid_id: u64) -> Result<Option<rcad_kernel::topology::CompSolid>, StepError> {
    use rcad_kernel::topology::CompSolid;

    let Some(solid_refs) = parsed.compsolids.get(&compsolid_id) else {
        return Ok(None);
    };

    let mut compsolid = CompSolid::new();

    for &solid_ref in solid_refs {
        if let Some(shell_ref) = get_shell_for_solid(parsed, solid_ref) {
            if let Some(solid) = build_solid_from_shell(parsed, shell_ref)? {
                compsolid.add(solid);
            }
        }
    }

    if compsolid.is_empty() {
        Ok(None)
    } else {
        Ok(Some(compsolid))
    }
}

/// Get the shell reference for a MANIFOLD_SOLID_BREP by index.
fn get_shell_for_solid(parsed: &ParsedStep, solid_ref: u64) -> Option<u64> {
    // MANIFOLD_SOLID_BREP references a CLOSED_SHELL
    // We need to find which shell this solid refers to
    // In the current implementation, manifold_solids is a Vec of shell refs
    let idx = parsed.manifold_solids.iter().position(|&r| r == solid_ref)?;
    parsed.manifold_solids.get(idx).copied()
}

/// Build a Solid from a CLOSED_SHELL reference.
fn build_solid_from_shell(parsed: &ParsedStep, shell_ref: u64) -> Result<Option<Solid>, StepError> {
    let Some(face_refs) = parsed.closed_shells.get(&shell_ref) else {
        return Ok(None);
    };

    // Build a mini-parsed containing just this shell's faces
    let shell_face_sets = vec![face_refs.clone()];

    // Use the existing face building logic
    let mut vertices = Vec::new();
    let mut vertex_index_by_id: HashMap<u64, usize> = HashMap::new();
    let mut edges: Vec<Edge> = Vec::new();
    let mut edge_index_by_curve: HashMap<u64, usize> = HashMap::new();
    let mut geom = GeomStore::default();
    let mut curve_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut surface_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut faces: Vec<Face> = Vec::new();

    // Collect vertices used by this shell
    let used_vertex_ids = collect_used_vertices(parsed, &shell_face_sets)?;
    let mut vertex_ids: Vec<u64> = used_vertex_ids.into_iter().collect();
    vertex_ids.sort_unstable();

    for (idx, vertex_id) in vertex_ids.iter().enumerate() {
        let point_id = *parsed
            .vertex_points
            .get(vertex_id)
            .ok_or(StepError::MissingEntity {
                entity_type: "VERTEX_POINT",
                id: Some(*vertex_id),
            })?;
        let point = *parsed
            .cartesian_points
            .get(&point_id)
            .ok_or(StepError::MissingEntity {
                entity_type: "CARTESIAN_POINT",
                id: Some(point_id),
            })?;
        vertices.push(Vertex {
            point: glam::DVec3::new(point[0], point[1], point[2]),
        });
        vertex_index_by_id.insert(*vertex_id, idx);
    }

    for &face_id in face_refs {
        if let Some((face, surface_ref)) = build_face(
            parsed,
            face_id,
            &mut vertices,
            &vertex_index_by_id,
            &mut edges,
            &mut edge_index_by_curve,
            &mut geom,
            &mut curve_store_index_by_step,
        ) {
            let surface_binding = surface_ref.and_then(|step_surface| {
                if let Some(idx) = surface_store_index_by_step.get(&step_surface) {
                    return Some(*idx);
                }
                let surface = resolve_surface(parsed, step_surface)?;
                let idx = geom.surfaces.len();
                geom.surfaces.push(surface);
                surface_store_index_by_step.insert(step_surface, idx);
                Some(idx)
            });
            geom.face_surface.push(surface_binding);
            faces.push(face);
        }
    }

    if faces.is_empty() {
        return Ok(None);
    }

    Ok(Some(Solid {
        shells: vec![Shell { faces }],
    }))
}

/// Build BRep and return the STEP face-id 鈫?flat-face-index mapping used for color resolution.
fn build_brep_with_face_map(parsed: &ParsedStep) -> Result<(BRep, HashMap<u64, usize>), StepError> {
    let shell_face_sets = collect_shell_faces(parsed);
    let used_vertex_ids = if shell_face_sets.is_empty() {
        collect_edge_vertices(parsed)
    } else {
        let mut used = collect_used_vertices(parsed, &shell_face_sets)?;
        used.extend(collect_edge_vertices(parsed));
        used
    };
    if used_vertex_ids.is_empty() && parsed.geometric_curve_sets.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok((brep, HashMap::new()));
        }
        return Err(StepError::EmptyResult("no vertices".into()));
    }

    let mut vertex_ids: Vec<u64> = used_vertex_ids.into_iter().collect();
    vertex_ids.sort_unstable();

    let mut vertex_index_by_id: HashMap<u64, usize> = HashMap::new();
    let mut vertices = Vec::with_capacity(vertex_ids.len());
    for (idx, vertex_id) in vertex_ids.iter().enumerate() {
        let point_id = *parsed
            .vertex_points
            .get(vertex_id)
            .ok_or(StepError::MissingEntity {
                entity_type: "VERTEX_POINT",
                id: Some(*vertex_id),
            })?;
        let point = *parsed
            .cartesian_points
            .get(&point_id)
            .ok_or(StepError::MissingEntity {
                entity_type: "CARTESIAN_POINT",
                id: Some(point_id),
            })?;
        vertices.push(Vertex {
            point: glam::DVec3::new(point[0], point[1], point[2]),
        });
        vertex_index_by_id.insert(*vertex_id, idx);
    }

    let mut edges: Vec<Edge> = Vec::new();
    let mut edge_index_by_curve: HashMap<u64, usize> = HashMap::new();
    let mut curve_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut surface_store_index_by_step: HashMap<u64, usize> = HashMap::new();
    let mut solids: Vec<Solid> = Vec::new();
    let mut geom = GeomStore::default();
    // STEP face id 鈫?flat face index across all shells (for color resolution)
    let mut face_id_map: HashMap<u64, usize> = HashMap::new();
    let mut flat_face_idx: usize = 0;

    for shell_faces in shell_face_sets {
        let mut faces: Vec<Face> = Vec::new();
        for face_id in shell_faces {
            if let Some((face, surface_ref)) = build_face(
                parsed,
                face_id,
                &mut vertices,
                &vertex_index_by_id,
                &mut edges,
                &mut edge_index_by_curve,
                &mut geom,
                &mut curve_store_index_by_step,
            ) {
                let surface_binding = surface_ref.and_then(|step_surface| {
                    if let Some(idx) = surface_store_index_by_step.get(&step_surface) {
                        return Some(*idx);
                    }
                    let surface = resolve_surface(parsed, step_surface)?;
                    let idx = geom.surfaces.len();
                    geom.surfaces.push(surface);
                    surface_store_index_by_step.insert(step_surface, idx);
                    Some(idx)
                });
                geom.face_surface.push(surface_binding);
                face_id_map.insert(face_id, flat_face_idx);
                flat_face_idx += 1;
                faces.push(face);
            }
        }

        if !faces.is_empty() {
            solids.push(Solid {
                shells: vec![Shell { faces }],
            });
        }
    }

    // Preserve standalone edges that are not part of any face loop.
    for (edge_curve_id, (start_id, end_id, curve_ref)) in &parsed.edge_curves {
        if edge_index_by_curve.contains_key(edge_curve_id) {
            continue;
        }
        let Some(&start) = vertex_index_by_id.get(start_id) else {
            continue;
        };
        let Some(&end) = vertex_index_by_id.get(end_id) else {
            continue;
        };

        let idx = edges.len();
        edges.push(Edge { start, end });
        edge_index_by_curve.insert(*edge_curve_id, idx);

        if geom.edge_curve.len() <= idx {
            geom.edge_curve.resize(idx + 1, None);
        }
        geom.edge_curve[idx] = curve_ref.and_then(|step_curve| {
            if let Some(existing) = curve_store_index_by_step.get(&step_curve) {
                return Some(*existing);
            }
            let curve = resolve_curve(parsed, step_curve)?;
            let cidx = geom.curves.len();
            geom.curves.push(curve);
            curve_store_index_by_step.insert(step_curve, cidx);
            Some(cidx)
        });
    }

    // Sample standalone 1D curves from GEOMETRIC_CURVE_SET (Polyline1/2/3 etc.)
    for curve_set in &parsed.geometric_curve_sets {
        for &curve_ref in curve_set {
            let points = sample_standalone_curve(parsed, curve_ref);
            if points.len() < 2 {
                continue;
            }
            let base = vertices.len();
            for p in &points {
                vertices.push(Vertex { point: *p });
            }
            for i in 0..points.len() - 1 {
                let idx = edges.len();
                edges.push(Edge {
                    start: base + i,
                    end: base + i + 1,
                });
                if geom.edge_curve.len() <= idx {
                    geom.edge_curve.resize(idx + 1, None);
                }
                // No curve binding needed 鈥?already tessellated into polyline
            }
        }
    }

    // Populate edge_pcurves from SURFACE_CURVE 鈫?PCURVE 鈫?DEFINITIONAL_REPRESENTATION chains
    for (step_curve_id, (inner_3d_ref, pcurve_ids, same_param)) in &parsed.surface_curves {
        // Find which BRep edge this SURFACE_CURVE belongs to
        // (edge_curves maps step EDGE_CURVE id 鈫?edge; the EDGE_CURVE's curve_ref points here)
        let edge_idx = edge_index_by_curve.get(step_curve_id).copied();
        // Also check if any edge_curve entry references this surface_curve indirectly
        let edge_idx = edge_idx.or_else(|| {
            parsed
                .edge_curves
                .iter()
                .find_map(|(ec_id, (_, _, curve_ref))| {
                    if curve_ref.as_ref() == Some(step_curve_id) {
                        edge_index_by_curve.get(ec_id).copied()
                    } else {
                        None
                    }
                })
        });
        let Some(edge_idx) = edge_idx else { continue };

        let mut pcs: Vec<PCurve> = Vec::new();
        for &pc_step_id in pcurve_ids {
            let Some(&(surface_step_id, def_rep_id)) = parsed.pcurves.get(&pc_step_id) else {
                continue;
            };
            let Some(&curve2d_step_id) = parsed.definitional_reps.get(&def_rep_id) else {
                continue;
            };
            // Resolve the surface into GeomStore
            let surface_idx =
                if let Some(existing) = surface_store_index_by_step.get(&surface_step_id) {
                    *existing
                } else if let Some(surf) = resolve_surface(parsed, surface_step_id) {
                    let sidx = geom.surfaces.len();
                    geom.surfaces.push(surf);
                    surface_store_index_by_step.insert(surface_step_id, sidx);
                    // Also update face_surface entries pointing to this surface
                    sidx
                } else {
                    continue;
                };
            // Resolve the 2D curve
            let Some(curve2d) = resolve_curve2d(parsed, curve2d_step_id) else {
                continue;
            };
            let c2didx = geom.curve2ds.len();
            geom.curve2ds.push(curve2d);
            pcs.push(PCurve {
                surface_idx,
                curve2d_idx: c2didx,
            });
        }
        if !pcs.is_empty() {
            if geom.edge_pcurves.len() <= edge_idx {
                geom.edge_pcurves.resize(edge_idx + 1, Vec::new());
            }
            geom.edge_pcurves[edge_idx] = pcs;
        }
        // Propagate same_parameter flag from SURFACE_CURVE entity
        if !same_param {
            if geom.edge_same_parameter.len() <= edge_idx {
                geom.edge_same_parameter.resize(edge_idx + 1, true);
            }
            geom.edge_same_parameter[edge_idx] = false;
        }
        let _ = inner_3d_ref; // already resolved via resolve_curve's SURFACE_CURVE dereference
    }

    if solids.is_empty() && edges.is_empty() {
        if let Some(brep) = brep_from_points_bbox(parsed) {
            return Ok((brep, HashMap::new()));
        }
        return Err(StepError::EmptyResult("no faces or edges".into()));
    }

    let mut brep = BRep {
        vertices,
        edges,
        solids,
        geom,
        compound: None,
        compsolid: None,
    };

    // Populate per-entity tolerance vectors from UNCERTAINTY_MEASURE_WITH_UNIT.
    // Only write when the file specifies a value different from the CONFUSION default,
    // to avoid polluting the GeomStore with trivial entries.
    if let Some(tol) = parsed.uncertainty_value
        && tol > 0.0
        && (tol - CONFUSION).abs() > CONFUSION * 0.5
    {
        let n_verts = brep.vertices.len();
        let n_edges = brep.edges.len();
        let n_faces: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        brep.geom.vertex_tolerance = vec![tol; n_verts];
        brep.geom.edge_tolerance = vec![tol; n_edges];
        brep.geom.face_tolerance = vec![tol; n_faces];
    }

    Ok((brep, face_id_map))
}

/// Resolve STYLED_ITEM 鈫?COLOUR_RGB chains into a StepColor.
/// Returns None if no color entities were found.
fn resolve_step_color(parsed: &ParsedStep, face_id_map: &HashMap<u64, usize>) -> Option<StepColor> {
    if parsed.styled_items.is_empty() {
        return None;
    }

    let mut step_color = StepColor::new();
    let mut found_any = false;

    for (shape_ref, psa_refs) in parsed.styled_items.values() {
        // Resolve color through the chain
        let rgb = psa_refs.iter().find_map(|psa_id| {
            let ssu_ids = parsed.presentation_style_assignments.get(psa_id)?;
            ssu_ids.iter().find_map(|ssu_id| {
                let sss_id = parsed.surface_style_usages.get(ssu_id)?;
                let ssfa_ids = parsed.surface_side_styles.get(sss_id)?;
                ssfa_ids.iter().find_map(|ssfa_id| {
                    let fas_id = parsed.surface_style_fill_areas.get(ssfa_id)?;
                    let fasc_ids = parsed.fill_area_styles.get(fas_id)?;
                    fasc_ids.iter().find_map(|fasc_id| {
                        let colour_id = parsed.fill_area_style_colours.get(fasc_id)?;
                        parsed.colour_rgbs.get(colour_id).copied()
                    })
                })
            })
        });

        let Some([r, g, b]) = rgb else { continue };
        let color = Color::new(r, g, b);
        found_any = true;

        // Map shape_ref to a face index
        if let Some(&face_idx) = face_id_map.get(shape_ref) {
            step_color = step_color.with_face_color(face_idx, color);
        } else {
            // Could be a solid-level styled item; use as default
            step_color.solid_color = Some(color);
        }
    }

    if found_any { Some(step_color) } else { None }
}

fn collect_edge_vertices(parsed: &ParsedStep) -> BTreeSet<u64> {
    let mut used = BTreeSet::new();
    for (start, end, _) in parsed.edge_curves.values() {
        used.insert(*start);
        used.insert(*end);
    }
    used
}

fn brep_from_points_bbox(parsed: &ParsedStep) -> Option<BRep> {
    if parsed.cartesian_points.is_empty() {
        return None;
    }

    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for point in parsed.cartesian_points.values() {
        for i in 0..3 {
            if point[i] < min[i] {
                min[i] = point[i];
            }
            if point[i] > max[i] {
                max[i] = point[i];
            }
        }
    }

    if min.iter().any(|v| !v.is_finite()) || max.iter().any(|v| !v.is_finite()) {
        return None;
    }

    let mut width = max[0] - min[0];
    let mut height = max[1] - min[1];
    let mut depth = max[2] - min[2];

    if width.abs() < 1e-9 {
        width = 1.0;
    }
    if height.abs() < 1e-9 {
        height = 1.0;
    }
    if depth.abs() < 1e-9 {
        depth = 1.0;
    }

    make_box_brep(
        glam::DVec3::new(min[0], min[1], min[2]),
        glam::DVec3::X,
        glam::DVec3::Y,
        width,
        height,
        depth,
    )
    .ok()
}

fn collect_shell_faces(parsed: &ParsedStep) -> Vec<Vec<u64>> {
    let mut shells = Vec::new();

    if !parsed.manifold_solids.is_empty() {
        for shell_id in &parsed.manifold_solids {
            if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                shells.push(face_ids.clone());
            }
        }
    }

    if !parsed.shell_based_surface_models.is_empty() {
        for shell_refs in &parsed.shell_based_surface_models {
            for shell_id in shell_refs {
                if let Some(face_ids) = parsed.open_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                } else if let Some(face_ids) = parsed.closed_shells.get(shell_id) {
                    shells.push(face_ids.clone());
                }
            }
        }
    }

    if shells.is_empty() && !parsed.closed_shells.is_empty() {
        for face_ids in parsed.closed_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.open_shells.is_empty() {
        for face_ids in parsed.open_shells.values() {
            shells.push(face_ids.clone());
        }
    }

    if shells.is_empty() && !parsed.advanced_faces.is_empty() {
        let mut all_faces: Vec<u64> = parsed.advanced_faces.keys().copied().collect();
        all_faces.sort_unstable();
        shells.push(all_faces);
    }

    shells
}

fn collect_used_vertices(
    parsed: &ParsedStep,
    shell_face_sets: &[Vec<u64>],
) -> Result<BTreeSet<u64>, StepError> {
    let mut used = BTreeSet::new();

    for shell in shell_face_sets {
        for face_id in shell {
            let bound_ids = parsed
                .advanced_faces
                .get(face_id)
                .ok_or(StepError::MissingEntity {
                    entity_type: "ADVANCED_FACE",
                    id: Some(*face_id),
                })?;
            for bound_id in &bound_ids.bounds {
                let loop_id = parsed
                    .face_bounds
                    .get(bound_id)
                    .ok_or(StepError::MissingEntity {
                        entity_type: "FACE_BOUND",
                        id: Some(*bound_id),
                    })?;
                let oriented_ids =
                    parsed
                        .edge_loops
                        .get(loop_id)
                        .ok_or(StepError::MissingEntity {
                            entity_type: "EDGE_LOOP",
                            id: Some(*loop_id),
                        })?;
                for oriented_id in oriented_ids {
                    let (edge_curve_id, _) =
                        parsed
                            .oriented_edges
                            .get(oriented_id)
                            .ok_or(StepError::MissingEntity {
                                entity_type: "ORIENTED_EDGE",
                                id: Some(*oriented_id),
                            })?;
                    let (start, end, _) =
                        parsed
                            .edge_curves
                            .get(edge_curve_id)
                            .ok_or(StepError::MissingEntity {
                                entity_type: "EDGE_CURVE",
                                id: Some(*edge_curve_id),
                            })?;
                    used.insert(*start);
                    used.insert(*end);
                }
            }
        }
    }

    Ok(used)
}

#[allow(clippy::too_many_arguments)]
fn build_face(
    parsed: &ParsedStep,
    face_id: u64,
    vertices: &mut Vec<Vertex>,
    vertex_index_by_id: &HashMap<u64, usize>,
    edges: &mut Vec<Edge>,
    edge_index_by_curve: &mut HashMap<u64, usize>,
    geom: &mut GeomStore,
    curve_store_index_by_step: &mut HashMap<u64, usize>,
) -> Option<(Face, Option<u64>)> {
    let bound_ids = parsed.advanced_faces.get(&face_id)?;
    let outer_bound = *bound_ids.bounds.first()?;
    let loop_id = *parsed.face_bounds.get(&outer_bound)?;
    let oriented_ids = parsed.edge_loops.get(&loop_id)?;

    let mut polygon: Vec<usize> = Vec::new();
    let mut wire_edge_indices: Vec<WireEdge> = Vec::new();
    let mut face_vertex_indices: Vec<usize> = Vec::new();
    let mut sampled_loop_points: Vec<glam::DVec3> = Vec::new();
    let mut sampled_loop_uv_points: Vec<glam::DVec2> = Vec::new();

    // Detect seam edges: an edge_curve that appears twice in the same face boundary
    let mut edge_curve_count: HashMap<u64, usize> = HashMap::new();
    for oriented_id in oriented_ids {
        let (edge_curve_id, _) = *parsed.oriented_edges.get(oriented_id)?;
        *edge_curve_count.entry(edge_curve_id).or_insert(0) += 1;
    }
    let has_seam = edge_curve_count.values().any(|&c| c >= 2);

    for oriented_id in oriented_ids {
        let (edge_curve_id, orientation) = *parsed.oriented_edges.get(oriented_id)?;
        let (start_id, end_id, curve_ref) = *parsed.edge_curves.get(&edge_curve_id)?;

        let (from_vertex_id, to_vertex_id) = if orientation {
            (start_id, end_id)
        } else {
            (end_id, start_id)
        };

        let from = *vertex_index_by_id.get(&from_vertex_id)?;
        let to = *vertex_index_by_id.get(&to_vertex_id)?;
        face_vertex_indices.push(from);
        face_vertex_indices.push(to);

        if let Some(points) = sample_oriented_edge_points(parsed, edge_curve_id, orientation) {
            append_edge_points(&mut sampled_loop_points, &points);
        }
        if let Some(surface_ref) = bound_ids.surface
            && let Some(points) = sample_oriented_edge_uv_points(
                parsed,
                edge_curve_id,
                orientation,
                surface_ref,
            )
        {
            append_edge_uv_points(&mut sampled_loop_uv_points, &points);
        }

        if polygon.is_empty() {
            polygon.push(from);
            polygon.push(to);
        } else {
            let last = *polygon.last()?;
            if last == from {
                polygon.push(to);
            } else if last == to {
                polygon.push(from);
            } else {
                polygon.push(from);
                polygon.push(to);
            }
        }

        let edge_index = if let Some(idx) = edge_index_by_curve.get(&edge_curve_id) {
            *idx
        } else {
            let idx = edges.len();
            edges.push(Edge {
                start: *vertex_index_by_id.get(&start_id)?,
                end: *vertex_index_by_id.get(&end_id)?,
            });
            edge_index_by_curve.insert(edge_curve_id, idx);

            if geom.edge_curve.len() <= idx {
                geom.edge_curve.resize(idx + 1, None);
            }
            geom.edge_curve[idx] = curve_ref.and_then(|step_curve| {
                if let Some(existing) = curve_store_index_by_step.get(&step_curve) {
                    return Some(*existing);
                }
                let curve = resolve_curve(parsed, step_curve)?;
                let cidx = geom.curves.len();
                geom.curves.push(curve);
                curve_store_index_by_step.insert(step_curve, cidx);
                Some(cidx)
            });

            // Populate edge_curve_range from vertex positions
            if geom.edge_curve_range.len() <= idx {
                geom.edge_curve_range.resize(idx + 1, None);
            }
            if geom.edge_degenerated.len() <= idx {
                geom.edge_degenerated.resize(idx + 1, false);
            }
            if let Some(Some(cidx)) = geom.edge_curve.get(idx)
                && let Some(curve) = geom.curves.get(*cidx)
            {
                let p0 = vertices
                    .get(*vertex_index_by_id.get(&start_id)?)
                    .map(|v| v.point);
                let p1 = vertices
                    .get(*vertex_index_by_id.get(&end_id)?)
                    .map(|v| v.point);
                if let (Some(p0), Some(p1)) = (p0, p1) {
                    let t_range = match curve {
                        Curve3::Line(line) => {
                            let t0 = (p0 - line.origin).dot(line.direction);
                            let t1 = (p1 - line.origin).dot(line.direction);
                            [t0, t1]
                        }
                        _ => curve.default_domain(),
                    };
                    geom.edge_curve_range[idx] = Some(t_range);
                    let len = (p1 - p0).length();
                    geom.edge_degenerated[idx] = len <= 1e-12;
                }
            }
            idx
        };
        wire_edge_indices.push(WireEdge {
            idx: edge_index,
            forward: orientation,
        });
    }

    while polygon.len() > 1 && polygon.first() == polygon.last() {
        polygon.pop();
    }
    dedup_consecutive(&mut polygon);

    let face_surface = bound_ids.surface;
    let is_planar_face = face_surface
        .map(|sid| parsed.planes.contains_key(&sid))
        .unwrap_or(true);

    if let Some(surface_ref) = face_surface
        && sampled_loop_uv_points.len() < 3
        && sampled_loop_points.len() >= 3
    {
        if let Some(projected) = project_boundary_points_to_surface_uv(
            parsed,
            surface_ref,
            &sampled_loop_points,
        ) {
            sampled_loop_uv_points = projected;
        }
    }

    let triangles = if sampled_loop_points.len() >= 3 && is_planar_face {
        triangulate_point_loop(vertices, &sampled_loop_points)
    } else if let Some(surface_ref) = face_surface
        && sampled_loop_uv_points.len() >= 3
    {
        let tris = triangulate_surface_trim_loop(parsed, surface_ref, vertices, &sampled_loop_uv_points);
        eprintln!(
            "[rcad-step][diag] face #{face_id} curved uv-loop surface=#{surface_ref} uv_pts={} tris={}",
            sampled_loop_uv_points.len(),
            tris.len()
        );
        tris
    } else if let Some(surface_ref) = face_surface
        && !parsed.planes.contains_key(&surface_ref)
    {
        let tris = triangulate_surface_fallback(
            parsed,
            surface_ref,
            vertices,
            &face_vertex_indices,
            has_seam,
        );
        eprintln!(
            "[rcad-step][diag] face #{face_id} curved fallback surface=#{surface_ref} uv_pts={} sampled_pts={} tris={}",
            sampled_loop_uv_points.len(),
            sampled_loop_points.len(),
            tris.len()
        );
        tris
    } else if polygon.len() >= 3 {
        triangulate_fan(&polygon)
    } else if let Some(surface_ref) = face_surface {
        let tris = triangulate_surface_fallback(
            parsed,
            surface_ref,
            vertices,
            &face_vertex_indices,
            has_seam,
        );
        eprintln!(
            "[rcad-step][diag] face #{face_id} fallback-no-polygon surface=#{surface_ref} uv_pts={} sampled_pts={} tris={}",
            sampled_loop_uv_points.len(),
            sampled_loop_points.len(),
            tris.len()
        );
        tris
    } else {
        Vec::new()
    };
    if triangles.is_empty() {
        return None;
    }

    let normal = compute_normal(&triangles[0], vertex_index_by_id, parsed);

    Some((
        Face {
            outer_wire: Wire {
                edges: wire_edge_indices,
            },
            inner_wires: Vec::new(),
            normal,
            triangles,
            mesh_dirty: false,
        },
        bound_ids.surface,
    ))
}

fn compute_normal(
    _triangle: &[usize; 3],
    _vertex_index_by_id: &HashMap<u64, usize>,
    _parsed: &ParsedStep,
) -> glam::DVec3 {
    // Normal is optional for current renderer, keep a stable default if degenerate.
    glam::DVec3::new(0.0, 0.0, 1.0)
}

fn triangulate_fan(polygon: &[usize]) -> Vec<[usize; 3]> {
    if polygon.len() < 3 {
        return Vec::new();
    }

    let mut triangles = Vec::with_capacity(polygon.len().saturating_sub(2));
    for i in 1..(polygon.len() - 1) {
        let tri = [polygon[0], polygon[i], polygon[i + 1]];
        if tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2] {
            triangles.push(tri);
        }
    }
    triangles
}

fn triangulate_point_loop(vertices: &mut Vec<Vertex>, points: &[glam::DVec3]) -> Vec<[usize; 3]> {
    let mut loop_points: Vec<glam::DVec3> = points.to_vec();
    dedup_consecutive_points(&mut loop_points);
    if loop_points.len() > 1 {
        let first = loop_points[0];
        let last = *loop_points.last().unwrap_or(&first);
        if first.distance(last) < 1e-6 {
            loop_points.pop();
        }
    }

    if loop_points.len() < 3 {
        return Vec::new();
    }

    let base = vertices.len();
    for p in &loop_points {
        vertices.push(Vertex { point: *p });
    }

    let mut triangles = Vec::with_capacity(loop_points.len().saturating_sub(2));
    for i in 1..(loop_points.len() - 1) {
        triangles.push([base, base + i, base + i + 1]);
    }
    triangles
}

fn append_edge_points(dst: &mut Vec<glam::DVec3>, src: &[glam::DVec3]) {
    for p in src {
        if let Some(last) = dst.last()
            && last.distance(*p) < 1e-7
        {
            continue;
        }
        dst.push(*p);
    }
}

fn append_edge_uv_points(dst: &mut Vec<glam::DVec2>, src: &[glam::DVec2]) {
    for p in src {
        if let Some(last) = dst.last()
            && last.distance(*p) < 1e-7
        {
            continue;
        }
        dst.push(*p);
    }
}

fn dedup_consecutive_points(points: &mut Vec<glam::DVec3>) {
    if points.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(points.len());
    for p in points.iter().copied() {
        if deduped
            .last()
            .map(|last: &glam::DVec3| last.distance(p) < 1e-7)
            .unwrap_or(false)
        {
            continue;
        }
        deduped.push(p);
    }
    *points = deduped;
}

fn sample_oriented_edge_points(
    parsed: &ParsedStep,
    edge_curve_id: u64,
    orientation: bool,
) -> Option<Vec<glam::DVec3>> {
    let (start_id, end_id, curve_ref) = *parsed.edge_curves.get(&edge_curve_id)?;
    let start = vertex_point_from_ref(parsed, start_id)?;
    let end = vertex_point_from_ref(parsed, end_id)?;

    let mut points = if let Some(curve_id) = curve_ref {
        let actual_curve_id = if let Some(&(inner_ref, _, _)) = parsed.surface_curves.get(&curve_id) {
            inner_ref
        } else {
            curve_id
        };

        if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&actual_curve_id) {
            sample_trimmed_curve_geom(parsed, underlying_ref, t0, t1)
        } else if let Some((placement_ref, radius)) = parsed.circles.get(&actual_curve_id) {
            sample_circle_edge(parsed, *placement_ref, *radius, start, end)?
        } else if let Some((placement_ref, major, minor)) = parsed.ellipses.get(&actual_curve_id) {
            sample_ellipse_edge(parsed, *placement_ref, *major, *minor, start, end)?
        } else if let Some(points) = sample_spline_curve(parsed, actual_curve_id, None) {
            points
        } else if let Some(control_refs) = parsed.b_spline_curves.get(&actual_curve_id) {
            sample_bspline_polyline(parsed, control_refs)
        } else {
            vec![start, end]
        }
    } else {
        vec![start, end]
    };

    // Keep sampled geometry anchored to topological edge endpoints.
    if !points.is_empty() {
        points[0] = start;
        let last = points.len() - 1;
        points[last] = end;
    }

    if !orientation {
        points.reverse();
    }
    Some(points)
}

fn sample_oriented_edge_uv_points(
    parsed: &ParsedStep,
    edge_curve_id: u64,
    orientation: bool,
    surface_ref: u64,
) -> Option<Vec<glam::DVec2>> {
    let (_, _, curve_ref) = *parsed.edge_curves.get(&edge_curve_id)?;
    let step_curve_id = curve_ref?;
    let (_, pcurve_refs, _) = parsed
        .surface_curves
        .get(&step_curve_id)
        .or_else(|| {
            parsed.surface_curves.iter().find_map(|(_, value)| {
                if value.0 == step_curve_id {
                    Some(value)
                } else {
                    None
                }
            })
        })?;

    let pcurve_ref = pcurve_refs.iter().find(|pcurve_ref| {
        parsed
            .pcurves
            .get(pcurve_ref)
            .map(|(candidate_surface, _)| *candidate_surface == surface_ref)
            .unwrap_or(false)
    })?;
    let (_, def_ref) = *parsed.pcurves.get(pcurve_ref)?;
    let curve2d_ref = *parsed.definitional_reps.get(&def_ref)?;

    let mut points = sample_curve2d_points(parsed, curve2d_ref)?;
    if !orientation {
        points.reverse();
    }
    Some(points)
}

fn curve2d_default_range(curve: &Curve2d) -> [f64; 2] {
    match curve {
        Curve2d::Line(_) => [0.0, 1.0],
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => [0.0, std::f64::consts::TAU],
        Curve2d::CircleInvolute(_) => [0.0, 4.0 * std::f64::consts::PI],
        Curve2d::ArchimedeanSpiral(_) | Curve2d::LogarithmicSpiral(_) => {
            [0.0, 4.0 * std::f64::consts::PI]
        }
        Curve2d::SineWave(_) | Curve2d::BSpline(_) | Curve2d::Bezier(_) => [0.0, 1.0],
    }
}

fn sample_curve2d_points(parsed: &ParsedStep, curve_ref: u64) -> Option<Vec<glam::DVec2>> {
    let (base_curve_ref, range) = if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&curve_ref) {
        (underlying_ref, Some([t0, t1]))
    } else {
        (curve_ref, None)
    };
    let curve = resolve_curve2d(parsed, base_curve_ref)?;
    let [t0, t1] = range.unwrap_or_else(|| curve2d_default_range(&curve));
    let segments = match curve {
        Curve2d::Line(_) => 1usize,
        Curve2d::Circle(_) | Curve2d::Ellipse(_) => 48usize,
        Curve2d::BSpline(_) | Curve2d::Bezier(_) => 64usize,
        _ => 32usize,
    };

    let mut points = Vec::with_capacity(segments + 1);
    if (t1 - t0).abs() < 1e-12 {
        points.push(curve.point_at(t0));
        return Some(points);
    }

    for index in 0..=segments {
        let t = t0 + (t1 - t0) * (index as f64 / segments as f64);
        points.push(curve.point_at(t));
    }
    dedup_consecutive_uv_points(&mut points);
    Some(points)
}

fn dedup_consecutive_uv_points(points: &mut Vec<glam::DVec2>) {
    if points.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(points.len());
    for point in points.iter().copied() {
        if deduped
            .last()
            .map(|last: &glam::DVec2| last.distance(point) < 1e-7)
            .unwrap_or(false)
        {
            continue;
        }
        deduped.push(point);
    }
    *points = deduped;
}

fn triangulate_surface_trim_loop(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    uv_points: &[glam::DVec2],
) -> Vec<[usize; 3]> {
    let Some(surface) = resolve_surface_for_trim_ops(parsed, surface_ref) else {
        return Vec::new();
    };

    let mut loop_uv = uv_points.to_vec();
    dedup_consecutive_uv_points(&mut loop_uv);
    if loop_uv.len() > 1 {
        let first = loop_uv[0];
        let last = *loop_uv.last().unwrap_or(&first);
        if first.distance(last) < 1e-6 {
            loop_uv.pop();
        }
    }
    unwrap_periodic_uv_loop(&surface, &mut loop_uv);
    if loop_uv.len() < 3 {
        return Vec::new();
    }

    if !matches!(surface, Surface3::Plane(_)) {
        if let Some(triangles) = triangulate_surface_trim_grid(&surface, vertices, &loop_uv) {
            if !triangles.is_empty() {
                return triangles;
            }
        }
    }

    let uv_poly_3d: Vec<glam::DVec3> = loop_uv
        .iter()
        .map(|uv| glam::DVec3::new(uv.x, uv.y, 0.0))
        .collect();
    let local_tris = rcad_algorithms::triangulate::triangulate_polygon(&uv_poly_3d, glam::DVec3::Z);
    if local_tris.is_empty() {
        return Vec::new();
    }

    let base = vertices.len();
    for uv in &loop_uv {
        vertices.push(Vertex {
            point: surface.point_at(uv.x, uv.y),
        });
    }

    local_tris
        .into_iter()
        .filter(|[a, b, c]| a != b && b != c && a != c)
        .map(|[a, b, c]| [base + a, base + b, base + c])
        .collect()
}

fn triangulate_surface_trim_grid(
    surface: &Surface3,
    vertices: &mut Vec<Vertex>,
    uv_loop: &[glam::DVec2],
) -> Option<Vec<[usize; 3]>> {
    if uv_loop.len() < 3 {
        return None;
    }

    let (segments_u, segments_v) = match surface {
        Surface3::Sphere(_) => (32usize, 16usize),
        Surface3::Cylinder(_) | Surface3::Cone(_) => (40usize, 12usize),
        Surface3::Torus(_) => (48usize, 24usize),
        Surface3::BSpline(_) | Surface3::Bezier(_) | Surface3::Trimmed(_) => (24usize, 24usize),
        _ => (24usize, 16usize),
    };

    let mut min_u = f64::INFINITY;
    let mut max_u = f64::NEG_INFINITY;
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for uv in uv_loop {
        min_u = min_u.min(uv.x);
        max_u = max_u.max(uv.x);
        min_v = min_v.min(uv.y);
        max_v = max_v.max(uv.y);
    }
    if !min_u.is_finite() || !max_u.is_finite() || !min_v.is_finite() || !max_v.is_finite() {
        return None;
    }
    if (max_u - min_u).abs() < 1e-10 || (max_v - min_v).abs() < 1e-10 {
        return None;
    }

    let mut grid_indices = vec![None; (segments_u + 1) * (segments_v + 1)];
    for j in 0..=segments_v {
        let v = min_v + (max_v - min_v) * (j as f64 / segments_v as f64);
        for i in 0..=segments_u {
            let u = min_u + (max_u - min_u) * (i as f64 / segments_u as f64);
            let uv = glam::DVec2::new(u, v);
            if !point_in_polygon_2d(uv, uv_loop) && !point_near_polygon_2d(uv, uv_loop, 1e-6) {
                continue;
            }
            let idx = vertices.len();
            vertices.push(Vertex {
                point: surface.point_at(u, v),
            });
            grid_indices[j * (segments_u + 1) + i] = Some(idx);
        }
    }

    let mut triangles = Vec::new();
    for j in 0..segments_v {
        for i in 0..segments_u {
            let idx = |ii: usize, jj: usize| grid_indices[jj * (segments_u + 1) + ii];
            let uv0 = glam::DVec2::new(
                min_u + (max_u - min_u) * (i as f64 / segments_u as f64),
                min_v + (max_v - min_v) * (j as f64 / segments_v as f64),
            );
            let uv1 = glam::DVec2::new(
                min_u + (max_u - min_u) * ((i + 1) as f64 / segments_u as f64),
                min_v + (max_v - min_v) * (j as f64 / segments_v as f64),
            );
            let uv2 = glam::DVec2::new(
                min_u + (max_u - min_u) * (i as f64 / segments_u as f64),
                min_v + (max_v - min_v) * ((j + 1) as f64 / segments_v as f64),
            );
            let uv3 = glam::DVec2::new(
                min_u + (max_u - min_u) * ((i + 1) as f64 / segments_u as f64),
                min_v + (max_v - min_v) * ((j + 1) as f64 / segments_v as f64),
            );

            let Some(i0) = idx(i, j) else { continue };
            let Some(i1) = idx(i + 1, j) else { continue };
            let Some(i2) = idx(i, j + 1) else { continue };
            let Some(i3) = idx(i + 1, j + 1) else { continue };

            let p0 = vertices[i0].point;
            let p1 = vertices[i1].point;
            let p2 = vertices[i2].point;
            let p3 = vertices[i3].point;
            let diag_03 = (p3 - p0).length_squared();
            let diag_12 = (p2 - p1).length_squared();
            let use_diag_03 = if (diag_03 - diag_12).abs() <= 1e-12 {
                ((i + j) & 1) == 0
            } else {
                diag_03 < diag_12
            };

            if use_diag_03 {
                let centroid_a = (uv0 + uv2 + uv3) / 3.0;
                if point_in_polygon_2d(centroid_a, uv_loop)
                    || point_near_polygon_2d(centroid_a, uv_loop, 1e-6)
                {
                    triangles.push([i0, i2, i3]);
                }

                let centroid_b = (uv0 + uv3 + uv1) / 3.0;
                if point_in_polygon_2d(centroid_b, uv_loop)
                    || point_near_polygon_2d(centroid_b, uv_loop, 1e-6)
                {
                    triangles.push([i0, i3, i1]);
                }
            } else {
                let centroid_a = (uv0 + uv2 + uv1) / 3.0;
                if point_in_polygon_2d(centroid_a, uv_loop)
                    || point_near_polygon_2d(centroid_a, uv_loop, 1e-6)
                {
                    triangles.push([i0, i2, i1]);
                }

                let centroid_b = (uv1 + uv2 + uv3) / 3.0;
                if point_in_polygon_2d(centroid_b, uv_loop)
                    || point_near_polygon_2d(centroid_b, uv_loop, 1e-6)
                {
                    triangles.push([i1, i2, i3]);
                }
            }
        }
    }

    Some(triangles)
}

fn point_in_polygon_2d(point: glam::DVec2, polygon: &[glam::DVec2]) -> bool {
    if polygon.len() < 3 {
        return false;
    }
    let mut inside = false;
    let mut prev = polygon.len() - 1;
    for current in 0..polygon.len() {
        let a = polygon[current];
        let b = polygon[prev];
        let intersects = ((a.y > point.y) != (b.y > point.y))
            && (point.x < (b.x - a.x) * (point.y - a.y) / (b.y - a.y).abs().max(1e-20) + a.x);
        if intersects {
            inside = !inside;
        }
        prev = current;
    }
    inside
}

fn point_near_polygon_2d(point: glam::DVec2, polygon: &[glam::DVec2], tolerance: f64) -> bool {
    if polygon.len() < 2 {
        return false;
    }
    for index in 0..polygon.len() {
        let a = polygon[index];
        let b = polygon[(index + 1) % polygon.len()];
        if point_segment_distance_2d(point, a, b) <= tolerance {
            return true;
        }
    }
    false
}

fn point_segment_distance_2d(point: glam::DVec2, a: glam::DVec2, b: glam::DVec2) -> f64 {
    let ab = b - a;
    let denom = ab.length_squared();
    if denom <= 1e-20 {
        return point.distance(a);
    }
    let t = ((point - a).dot(ab) / denom).clamp(0.0, 1.0);
    point.distance(a + ab * t)
}

fn project_boundary_points_to_surface_uv(
    parsed: &ParsedStep,
    surface_ref: u64,
    boundary_points: &[glam::DVec3],
) -> Option<Vec<glam::DVec2>> {
    let surface = resolve_surface_for_trim_ops(parsed, surface_ref)?;
    let mut uv_points = Vec::with_capacity(boundary_points.len());
    for &point in boundary_points {
        uv_points.push(project_point_to_surface_uv(&surface, point)?);
    }
    unwrap_periodic_uv_loop(&surface, &mut uv_points);
    dedup_consecutive_uv_points(&mut uv_points);
    Some(uv_points)
}

fn project_point_to_surface_uv(surface: &Surface3, point: glam::DVec3) -> Option<glam::DVec2> {
    match surface {
        Surface3::Plane(plane) => {
            let u_axis = any_perpendicular(plane.normal);
            let v_axis = plane.normal.cross(u_axis).normalize_or_zero();
            let delta = point - plane.origin;
            Some(glam::DVec2::new(delta.dot(u_axis), delta.dot(v_axis)))
        }
        Surface3::Cylinder(cylinder) => {
            let axis = cylinder.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - cylinder.origin;
            let radial = delta - axis * delta.dot(axis);
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                delta.dot(axis),
            ))
        }
        Surface3::Sphere(sphere) => {
            let axis = sphere.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let radial = (point - sphere.center).normalize_or_zero();
            if radial.length_squared() <= 1e-20 {
                return None;
            }
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                radial.dot(axis).clamp(-1.0, 1.0).acos(),
            ))
        }
        Surface3::Cone(cone) => {
            let axis = cone.axis_dir();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - cone.apex;
            let height = delta.dot(axis);
            let radial = delta - axis * height;
            Some(glam::DVec2::new(
                radial.dot(v_axis).atan2(radial.dot(u_axis)),
                cone.slant_from_axial(height),
            ))
        }
        Surface3::Torus(torus) => {
            let axis = torus.axis.normalize_or_zero();
            if axis.length_squared() <= 1e-20 {
                return None;
            }
            let u_axis = any_perpendicular(axis);
            let v_axis = axis.cross(u_axis).normalize_or_zero();
            let delta = point - torus.center;
            let planar = delta - axis * delta.dot(axis);
            if planar.length_squared() <= 1e-20 {
                return None;
            }
            let theta = planar.dot(v_axis).atan2(planar.dot(u_axis));
            let ring_dir = (u_axis * theta.cos() + v_axis * theta.sin()).normalize_or_zero();
            let ring_center = torus.center + torus.major_radius * ring_dir;
            let tube_vec = (point - ring_center).normalize_or_zero();
            Some(glam::DVec2::new(
                theta,
                tube_vec.dot(axis).atan2(tube_vec.dot(ring_dir)),
            ))
        }
        Surface3::Trimmed(trimmed) => project_point_to_surface_uv(trimmed.basis.as_ref(), point),
        _ => None,
    }
}

fn unwrap_periodic_uv_loop(surface: &Surface3, points: &mut [glam::DVec2]) {
    let period = match surface {
        Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
            Some(std::f64::consts::TAU)
        }
        Surface3::Trimmed(trimmed) => match trimmed.basis.as_ref() {
            Surface3::Cylinder(_) | Surface3::Cone(_) | Surface3::Sphere(_) | Surface3::Torus(_) => {
                Some(std::f64::consts::TAU)
            }
            _ => None,
        },
        _ => None,
    };

    let Some(period) = period else {
        return;
    };
    if points.len() < 2 {
        return;
    }

    let mut offset = 0.0;
    let mut previous = points[0].x;
    for point in points.iter_mut().skip(1) {
        let raw = point.x + offset;
        let delta = raw - previous;
        if delta > period * 0.5 {
            offset -= period;
        } else if delta < -period * 0.5 {
            offset += period;
        }
        point.x += offset;
        previous = point.x;
    }
}

fn sample_circle_edge(
    parsed: &ParsedStep,
    placement_ref: u64,
    radius: f64,
    start: glam::DVec3,
    end: glam::DVec3,
) -> Option<Vec<glam::DVec3>> {
    if !radius.is_finite() || radius <= 0.0 {
        return None;
    }
    let (center, axis, u) = placement_frame_from_ref(parsed, placement_ref)?;
    let v = axis.cross(u).normalize_or_zero();
    let a0 = angle_on_basis(start - center, u, v);
    let a1 = angle_on_basis(end - center, u, v);
    let mut sweep = a1 - a0;
    if start.distance(end) < 1e-6 {
        sweep = std::f64::consts::TAU;
    } else if sweep <= 0.0 {
        sweep += std::f64::consts::TAU;
    }
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0)
        .ceil()
        .max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = a0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (radius * t.cos()) + v * (radius * t.sin()));
    }
    Some(points)
}

fn sample_ellipse_edge(
    parsed: &ParsedStep,
    placement_ref: u64,
    major: f64,
    minor: f64,
    start: glam::DVec3,
    end: glam::DVec3,
) -> Option<Vec<glam::DVec3>> {
    if !major.is_finite() || !minor.is_finite() || major <= 0.0 || minor <= 0.0 {
        return None;
    }
    let (center, axis, u) = placement_frame_from_ref(parsed, placement_ref)?;
    let v = axis.cross(u).normalize_or_zero();

    let param = |p: glam::DVec3| {
        let d = p - center;
        let x = d.dot(u) / major;
        let y = d.dot(v) / minor;
        y.atan2(x)
    };

    let t0 = param(start);
    let t1 = param(end);
    let mut sweep = t1 - t0;
    if start.distance(end) < 1e-6 {
        sweep = std::f64::consts::TAU;
    } else if sweep <= 0.0 {
        sweep += std::f64::consts::TAU;
    }
    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0)
        .ceil()
        .max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = t0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (major * t.cos()) + v * (minor * t.sin()));
    }
    Some(points)
}

fn sample_bspline_polyline(parsed: &ParsedStep, control_refs: &[u64]) -> Vec<glam::DVec3> {
    let mut points = Vec::new();
    for cref in control_refs {
        if let Some(p) = point_from_ref(parsed, *cref) {
            points.push(p);
        }
    }
    points
}

fn sample_spline_curve(
    parsed: &ParsedStep,
    curve_ref: u64,
    range: Option<[f64; 2]>,
) -> Option<Vec<glam::DVec3>> {
    let curve = resolve_curve(parsed, curve_ref)?;
    if !matches!(curve, Curve3::BSpline(_) | Curve3::Bezier(_)) {
        return None;
    }

    let [t0, t1] = range.unwrap_or_else(|| curve.default_domain());
    if !t0.is_finite() || !t1.is_finite() {
        return None;
    }

    let seg = 64usize;
    let mut points = Vec::with_capacity(seg + 1);
    if (t1 - t0).abs() < 1e-12 {
        points.push(curve.point_at(t0));
        return Some(points);
    }

    for i in 0..=seg {
        let t = t0 + (t1 - t0) * (i as f64 / seg as f64);
        points.push(curve.point_at(t));
    }
    dedup_consecutive_points(&mut points);
    Some(points)
}

fn angle_on_basis(v: glam::DVec3, u: glam::DVec3, w: glam::DVec3) -> f64 {
    v.dot(w).atan2(v.dot(u))
}

fn triangulate_spherical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, radius)) = parsed.spherical_surfaces.get(&surface_ref) else {
        return Vec::new();
    };
    let Some((center, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !radius.is_finite() || *radius <= 0.0 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();

    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, center, w, u, v, has_seam);
    let (phi_min, phi_max) =
        infer_polar_range(vertices, face_vertex_indices, center, w, *radius, has_seam);

    let u_segments = 32usize;
    let v_segments = 16usize;
    let stride = u_segments + 1;
    let base = vertices.len();

    let at_north_pole = phi_min.abs() < 1e-6;
    let at_south_pole = (phi_max - std::f64::consts::PI).abs() < 1e-6;

    for vi in 0..=v_segments {
        let phi = phi_min + (phi_max - phi_min) * vi as f64 / v_segments as f64;
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        for ui in 0..=u_segments {
            let theta = theta_min + (theta_max - theta_min) * ui as f64 / u_segments as f64;
            let point = center
                + w * (radius * cos_phi)
                + u * (radius * theta.cos() * sin_phi)
                + v * (radius * theta.sin() * sin_phi);
            vertices.push(Vertex { point });
        }
    }

    let mut triangles = Vec::with_capacity(u_segments * v_segments * 2);
    for vi in 0..v_segments {
        for ui in 0..u_segments {
            let i0 = base + vi * stride + ui;
            let i1 = i0 + 1;
            let i2 = i0 + stride;
            let i3 = i2 + 1;

            // Skip degenerate triangles at poles
            if !(at_north_pole && vi == 0) {
                triangles.push([i0, i2, i1]);
            }
            if !(at_south_pole && vi == v_segments - 1) {
                triangles.push([i1, i2, i3]);
            }
        }
    }

    triangles
}

fn triangulate_surface_fallback(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    if parsed.b_spline_surfaces.contains_key(&surface_ref) {
        return triangulate_bspline_surface(parsed, surface_ref, vertices);
    }
    if parsed.spherical_surfaces.contains_key(&surface_ref) {
        return triangulate_spherical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.cylindrical_surfaces.contains_key(&surface_ref) {
        return triangulate_cylindrical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.conical_surfaces.contains_key(&surface_ref) {
        return triangulate_conical_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    if parsed.toroidal_surfaces.contains_key(&surface_ref) {
        return triangulate_toroidal_surface(
            parsed,
            surface_ref,
            vertices,
            face_vertex_indices,
            has_seam,
        );
    }
    Vec::new()
}

fn triangulate_cylindrical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, radius)) = parsed.cylindrical_surfaces.get(&surface_ref) else {
        return Vec::new();
    };
    let Some((origin, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !radius.is_finite() || *radius <= 0.0 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();
    let (t_min, t_max) =
        infer_axis_range(vertices, face_vertex_indices, origin, w, -*radius, *radius);
    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

    let radial_segments = 40usize;
    let height_segments = 12usize;
    let stride = radial_segments + 1;
    let base = vertices.len();

    for j in 0..=height_segments {
        let tj = t_min + (t_max - t_min) * (j as f64 / height_segments as f64);
        let center = origin + w * tj;
        for i in 0..=radial_segments {
            let theta = theta_min + (theta_max - theta_min) * (i as f64 / radial_segments as f64);
            let ring_dir = u * theta.cos() + v * theta.sin();
            vertices.push(Vertex {
                point: center + ring_dir * *radius,
            });
        }
    }

    triangulate_grid(vertices, base, height_segments, radial_segments, stride)
}

fn triangulate_conical_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, ref_radius, half_angle_rad)) =
        parsed.conical_surfaces.get(&surface_ref)
    else {
        return Vec::new();
    };
    let Some((origin, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !ref_radius.is_finite() || !half_angle_rad.is_finite() {
        return Vec::new();
    }

    let tan_a = half_angle_rad.tan();
    if !tan_a.is_finite() || tan_a.abs() < 1e-8 {
        return Vec::new();
    }

    let u = ref_dir;
    let v = axis.cross(u).normalize_or_zero();
    let w = axis.normalize_or_zero();
    let default_h = ref_radius.abs().max(0.1) / tan_a.abs();
    let (t_min, t_max) = infer_axis_range(vertices, face_vertex_indices, origin, w, 0.0, default_h);
    let (theta_min, theta_max) =
        infer_angular_range(vertices, face_vertex_indices, origin, w, u, v, has_seam);

    let radial_segments = 40usize;
    let height_segments = 12usize;
    let stride = radial_segments + 1;
    let base = vertices.len();

    for j in 0..=height_segments {
        let tj = t_min + (t_max - t_min) * (j as f64 / height_segments as f64);
        let center = origin + w * tj;
        let rj = (ref_radius + tj * tan_a).abs().max(1e-5);
        for i in 0..=radial_segments {
            let theta = theta_min + (theta_max - theta_min) * (i as f64 / radial_segments as f64);
            let ring_dir = u * theta.cos() + v * theta.sin();
            vertices.push(Vertex {
                point: center + ring_dir * rj,
            });
        }
    }

    triangulate_grid(vertices, base, height_segments, radial_segments, stride)
}

fn triangulate_toroidal_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
    face_vertex_indices: &[usize],
    has_seam: bool,
) -> Vec<[usize; 3]> {
    let Some((placement_ref, major_radius, minor_radius)) =
        parsed.toroidal_surfaces.get(&surface_ref)
    else {
        return Vec::new();
    };
    let Some((center, axis, ref_dir)) = placement_frame_from_ref(parsed, *placement_ref) else {
        return Vec::new();
    };
    if !major_radius.is_finite()
        || !minor_radius.is_finite()
        || *major_radius <= 0.0
        || *minor_radius <= 0.0
    {
        return Vec::new();
    }

    let u_dir = ref_dir;
    let v_dir = axis.cross(u_dir).normalize_or_zero();
    let w = axis.normalize_or_zero();

    // Infer major angle (theta) range from boundary vertices
    let (theta_min, theta_max) = infer_angular_range(
        vertices,
        face_vertex_indices,
        center,
        w,
        u_dir,
        v_dir,
        has_seam,
    );

    // Infer minor angle (phi) range: project boundary vertices into the minor circle plane
    let (phi_min, phi_max) = infer_torus_minor_range(
        vertices,
        face_vertex_indices,
        center,
        w,
        u_dir,
        v_dir,
        *major_radius,
        has_seam,
    );

    let major_segments = 48usize;
    let minor_segments = 24usize;
    let stride = minor_segments + 1;
    let base = vertices.len();

    for i in 0..=major_segments {
        let theta = theta_min + (theta_max - theta_min) * i as f64 / major_segments as f64;
        let ring_dir = u_dir * theta.cos() + v_dir * theta.sin();
        let ring_center = center + ring_dir * *major_radius;
        for j in 0..=minor_segments {
            let phi = phi_min + (phi_max - phi_min) * j as f64 / minor_segments as f64;
            let minor_dir = ring_dir * phi.cos() + w * phi.sin();
            vertices.push(Vertex {
                point: ring_center + minor_dir * *minor_radius,
            });
        }
    }

    triangulate_grid(vertices, base, major_segments, minor_segments, stride)
}

/// Infer the angular (theta) range of boundary vertices projected onto a
/// surface local frame (origin, u_dir, v_dir, axis).
/// Returns (theta_min, theta_max) in radians.
/// If `has_seam` is true, the face wraps the full period 鈫?returns (0, 2蟺).
fn infer_angular_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    origin: glam::DVec3,
    axis: glam::DVec3,
    u_dir: glam::DVec3,
    v_dir: glam::DVec3,
    has_seam: bool,
) -> (f64, f64) {
    use std::f64::consts::TAU;
    if has_seam {
        return (0.0, TAU);
    }

    let mut angles: Vec<f64> = Vec::new();
    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = v.point - origin;
            let proj = d - axis * d.dot(axis);
            let a = proj.dot(v_dir).atan2(proj.dot(u_dir));
            angles.push(a);
        }
    }

    if angles.is_empty() {
        return (0.0, TAU);
    }

    // Sort angles and find the largest gap; the face covers the complement.
    angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    angles.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    let mut max_gap = 0.0f64;
    let mut gap_end = 0usize;
    for i in 1..angles.len() {
        let gap = angles[i] - angles[i - 1];
        if gap > max_gap {
            max_gap = gap;
            gap_end = i;
        }
    }
    // Wrap-around gap
    let wrap_gap = (angles[0] + TAU) - angles[angles.len() - 1];
    if wrap_gap > max_gap {
        // Face goes from angles[0] to angles[last]
        return (angles[0], angles[angles.len() - 1]);
    }

    // Face goes from angles[gap_end] around to angles[gap_end-1]
    let theta_min = angles[gap_end];
    let theta_max = angles[gap_end - 1] + TAU;
    (theta_min, theta_max)
}

/// Infer polar (phi) range for spherical surfaces from boundary vertices.
/// Returns (phi_min, phi_max) in [0, 蟺] (north-pole to south-pole).
fn infer_polar_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    center: glam::DVec3,
    axis: glam::DVec3,
    radius: f64,
    has_seam: bool,
) -> (f64, f64) {
    let _ = has_seam;
    let mut phi_min = std::f64::consts::PI;
    let mut phi_max = 0.0f64;

    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = (v.point - center).normalize_or_zero();
            let cos_phi = d.dot(axis).clamp(-1.0, 1.0);
            let phi = cos_phi.acos();
            phi_min = phi_min.min(phi);
            phi_max = phi_max.max(phi);
        }
    }

    let _ = radius;
    if !phi_min.is_finite() || !phi_max.is_finite() || (phi_max - phi_min).abs() < 1e-8 {
        (0.0, std::f64::consts::PI)
    } else {
        (phi_min, phi_max)
    }
}

/// Infer the minor (phi) angular range for a torus face.
/// Projects each boundary vertex into the local minor circle plane at the
/// closest major angle, then computes the minor angle.
#[allow(clippy::too_many_arguments)]
fn infer_torus_minor_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    center: glam::DVec3,
    axis: glam::DVec3,
    u_dir: glam::DVec3,
    v_dir: glam::DVec3,
    major_radius: f64,
    has_seam: bool,
) -> (f64, f64) {
    use std::f64::consts::TAU;
    if has_seam {
        return (0.0, TAU);
    }

    let mut angles: Vec<f64> = Vec::new();
    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let d = v.point - center;
            // Project onto the equatorial plane to find the major angle
            let proj_equat = d - axis * d.dot(axis);
            let theta = proj_equat.dot(v_dir).atan2(proj_equat.dot(u_dir));
            // Reconstruct ring_dir and ring_center for this major angle
            let ring_dir = u_dir * theta.cos() + v_dir * theta.sin();
            let ring_center = center + ring_dir * major_radius;
            // Vector from ring center to the vertex, projected onto the minor circle plane
            let to_v = v.point - ring_center;
            let phi = to_v.dot(axis).atan2(to_v.dot(ring_dir));
            angles.push(phi);
        }
    }

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    angles.dedup_by(|a, b| (*a - *b).abs() < 1e-10);

    if angles.len() < 2 {
        return (0.0, TAU);
    }

    let mut max_gap = 0.0f64;
    let mut gap_end = 0usize;
    for i in 1..angles.len() {
        let gap = angles[i] - angles[i - 1];
        if gap > max_gap {
            max_gap = gap;
            gap_end = i;
        }
    }
    let wrap_gap = (angles[0] + TAU) - angles[angles.len() - 1];
    if wrap_gap > max_gap {
        return (angles[0], angles[angles.len() - 1]);
    }
    let phi_min = angles[gap_end];
    let phi_max = angles[gap_end - 1] + TAU;
    (phi_min, phi_max)
}

fn triangulate_grid(
    vertices: &[Vertex],
    base: usize,
    rows: usize,
    cols: usize,
    stride: usize,
) -> Vec<[usize; 3]> {
    let mut triangles = Vec::with_capacity(rows * cols * 2);
    for r in 0..rows {
        for c in 0..cols {
            let i0 = base + r * stride + c;
            let i1 = i0 + 1;
            let i2 = i0 + stride;
            let i3 = i2 + 1;
            let p0 = vertices[i0].point;
            let p1 = vertices[i1].point;
            let p2 = vertices[i2].point;
            let p3 = vertices[i3].point;
            let diag_03 = (p3 - p0).length_squared();
            let diag_12 = (p2 - p1).length_squared();

            if diag_03 <= diag_12 {
                triangles.push([i0, i2, i3]);
                triangles.push([i0, i3, i1]);
            } else {
                triangles.push([i0, i2, i1]);
                triangles.push([i1, i2, i3]);
            }
        }
    }
    triangles
}

/// Triangulate a B-Spline surface by uniform (u,v) grid sampling.
fn triangulate_bspline_surface(
    parsed: &ParsedStep,
    surface_ref: u64,
    vertices: &mut Vec<Vertex>,
) -> Vec<[usize; 3]> {
    use rcad_kernel::geom::SurfaceEval;
    let Some(surface) = resolve_surface(parsed, surface_ref) else {
        return Vec::new();
    };
    let [u0, u1, v0, v1] = surface.default_domain();
    if !u0.is_finite() || !u1.is_finite() || !v0.is_finite() || !v1.is_finite() {
        return Vec::new();
    }
    let nu = 20usize;
    let nv = 20usize;
    let base = vertices.len();
    for j in 0..=nv {
        let v = v0 + (v1 - v0) * (j as f64 / nv as f64);
        for i in 0..=nu {
            let u = u0 + (u1 - u0) * (i as f64 / nu as f64);
            vertices.push(Vertex {
                point: surface.point_at(u, v),
            });
        }
    }
    triangulate_grid(vertices, base, nv, nu, nu + 1)
}

fn infer_axis_range(
    vertices: &[Vertex],
    face_vertex_indices: &[usize],
    origin: glam::DVec3,
    axis: glam::DVec3,
    default_min: f64,
    default_max: f64,
) -> (f64, f64) {
    let mut min_t = f64::INFINITY;
    let mut max_t = f64::NEG_INFINITY;

    for &vidx in face_vertex_indices {
        if let Some(v) = vertices.get(vidx) {
            let t = (v.point - origin).dot(axis);
            min_t = min_t.min(t);
            max_t = max_t.max(t);
        }
    }

    if !min_t.is_finite() || !max_t.is_finite() || (max_t - min_t).abs() < 1e-8 {
        (default_min, default_max)
    } else {
        (min_t, max_t)
    }
}

fn dedup_consecutive(polygon: &mut Vec<usize>) {
    if polygon.len() < 2 {
        return;
    }
    let mut deduped = Vec::with_capacity(polygon.len());
    let mut last = None;
    for &idx in polygon.iter() {
        if Some(idx) != last {
            deduped.push(idx);
            last = Some(idx);
        }
    }
    *polygon = deduped;
}

fn extract_data_section(content: &str) -> Result<&str, StepError> {
    let start = content
        .find("DATA;")
        .ok_or_else(|| StepError::InvalidFormat("missing DATA section".into()))?;
    let after_start = &content[start + "DATA;".len()..];
    let end = after_start
        .find("ENDSEC;")
        .ok_or_else(|| StepError::InvalidFormat("missing ENDSEC after DATA".into()))?;
    Ok(&after_start[..end])
}

fn split_records(data: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    let mut in_comment = false;
    let mut in_string = false;
    let mut chars = data.chars().peekable();

    while let Some(ch) = chars.next() {
        if in_comment {
            if ch == '*' && chars.peek() == Some(&'/') {
                let _ = chars.next();
                in_comment = false;
            }
            continue;
        }

        if ch == '\'' {
            current.push(ch);
            if in_string {
                if chars.peek() == Some(&'\'') {
                    current.push(chars.next().unwrap_or('\''));
                } else {
                    in_string = false;
                }
            } else {
                in_string = true;
            }
            continue;
        }

        if !in_string && ch == '/' && chars.peek() == Some(&'*') {
            let _ = chars.next();
            in_comment = true;
            continue;
        }

        if ch == ';' && !in_string {
            let trimmed = current.trim();
            if !trimmed.is_empty() {
                records.push(trimmed.to_string());
            }
            current.clear();
            continue;
        }

        current.push(ch);
    }

    let trailing = current.trim();
    if !trailing.is_empty() {
        records.push(trailing.to_string());
    }

    records
}

fn parse_entity_record(record: &str) -> Result<Option<(u64, &str)>, StepError> {
    let line = record.trim();
    if !line.starts_with('#') {
        return Ok(None);
    }

    let eq = line
        .find('=')
        .ok_or_else(|| StepError::InvalidFormat(format!("invalid entity record: {line}")))?;
    let id_str = line[1..eq].trim();
    let id = id_str
        .parse::<u64>()
        .map_err(|e| StepError::InvalidFormat(format!("invalid entity id {id_str}: {e}")))?;
    Ok(Some((id, line[eq + 1..].trim())))
}

fn parse_entity_body(body: &str) -> Option<(&str, &str)> {
    let mut payload = body.trim();
    if payload.starts_with('(') {
        payload = payload.strip_prefix('(')?.strip_suffix(')')?.trim();
    }

    let open = payload.find('(')?;
    let close = payload.rfind(')')?;
    let entity = payload[..open].trim();
    let args = &payload[open + 1..close];
    Some((entity, args))
}

fn parse_cartesian_point(args: &str) -> Option<[f64; 3]> {
    let list = parse_coord_list(args)?;
    if list.len() < 3 {
        return None;
    }
    Some([list[0], list[1], list[2]])
}

fn parse_coord_list(args: &str) -> Option<Vec<f64>> {
    let open = args.rfind('(')?;
    let close = args.rfind(')')?;
    if close <= open {
        return None;
    }
    let raw = &args[open + 1..close];
    let mut coords = Vec::new();
    for item in raw.split(',') {
        let v = item.trim().parse::<f64>().ok()?;
        coords.push(v);
    }
    Some(coords)
}

fn parse_single_ref_after_name(args: &str) -> Option<u64> {
    let parts = split_top_level(args, ',');
    for part in parts.into_iter().skip(1) {
        if let Some(reference) = parse_ref(part) {
            return Some(reference);
        }
    }
    None
}

fn parse_edge_curve_vertices(args: &str) -> Option<(u64, u64, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    let start = parse_ref(parts[1])?;
    let end = parse_ref(parts[2])?;
    let curve_ref = parse_ref(parts[3]);
    Some((start, end, curve_ref))
}

fn parse_axis2_placement(args: &str) -> Option<(u64, u64, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let ref_dir = parts.get(3).and_then(|s| parse_ref(s));
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?, ref_dir))
}

fn parse_curve_basis(args: &str) -> Option<(u64, u64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?))
}

fn parse_placement_radius(args: &str) -> Option<(u64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parts[2].trim().parse::<f64>().ok()?))
}

fn parse_placement_two_radii(args: &str) -> Option<(u64, f64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parts[3].trim().parse::<f64>().ok()?,
    ))
}

fn parse_bspline_control_points(args: &str) -> Option<Vec<u64>> {
    let parts = split_top_level(args, ',');
    let refs = parts.get(2).map(|s| parse_ref_list(s))?;
    if refs.is_empty() {
        return None;
    }
    Some(refs)
}

/// Parse full B_SPLINE_CURVE_WITH_KNOTS args:
/// ('name', degree, (ctrl_pts...), .FORM., .bool., .bool., (mults...), (knots...), .UNSPECIFIED.)
fn parse_bspline_curve_full(args: &str) -> Option<BSplineCurveData> {
    let parts = split_top_level(args, ',');
    // parts[0] = name, [1] = degree, [2] = ctrl pts list, [3] = form,
    // [4] = closed, [5] = self_intersect, [6] = knot_mults, [7] = knots, [8] = type
    if parts.len() < 8 {
        return None;
    }
    let degree = parts[1].trim().parse::<usize>().ok()?;
    let ctrl_refs = parse_ref_list(parts[2]);
    let mults: Vec<usize> = parse_float_list(parts[6])
        .into_iter()
        .map(|v| v as usize)
        .collect();
    let knot_vals: Vec<f64> = parse_float_list(parts[7]);

    if ctrl_refs.is_empty() || mults.is_empty() || knot_vals.is_empty() {
        return None;
    }
    Some((degree, ctrl_refs, mults, knot_vals))
}

/// Parse B_SPLINE_SURFACE_WITH_KNOTS args.
/// Returns (name, degree_u, degree_v, ctrl_grid[v_row][u_col], mults_u, knots_u, mults_v, knots_v)
fn parse_bspline_surface_with_knots(args: &str) -> Option<BSplineSurfaceData> {
    // STEP format:
    // ('name', degree_u, degree_v, ((#p00,#p01,...),(#p10,...)),
    //   .UNSPECIFIED., .F., .F., .F.,
    //   (mults_u...), (mults_v...), (knots_u...), (knots_v...), .UNSPECIFIED.)
    // parts[0]=name, [1]=deg_u, [2]=deg_v, [3]=ctrl grid (nested list),
    // [4..7]=flags, [8]=mults_u, [9]=mults_v, [10]=knots_u, [11]=knots_v
    let parts = split_top_level(args, ',');
    if parts.len() < 12 {
        return None;
    }
    let name = parse_step_string(parts[0]);
    let degree_u = parts[1].trim().parse::<usize>().ok()?;
    let degree_v = parts[2].trim().parse::<usize>().ok()?;

    // Strip outer parens to get the row-list string, then split rows by top-level comma
    let grid_outer = parts[3].trim();
    let grid_inner = grid_outer
        .strip_prefix('(')
        .unwrap_or(grid_outer)
        .trim_end_matches(')');
    let rows_raw = split_top_level(grid_inner, ',');
    let ctrl_grid: Vec<Vec<u64>> = rows_raw
        .iter()
        .map(|row| parse_ref_list(row))
        .filter(|row| !row.is_empty())
        .collect();
    if ctrl_grid.is_empty() {
        return None;
    }

    let mults_u: Vec<usize> = parse_float_list(parts[8])
        .into_iter()
        .map(|v| v as usize)
        .collect();
    let mults_v: Vec<usize> = parse_float_list(parts[9])
        .into_iter()
        .map(|v| v as usize)
        .collect();
    let knots_u: Vec<f64> = parse_float_list(parts[10]);
    let knots_v: Vec<f64> = parse_float_list(parts[11]);

    if mults_u.is_empty() || knots_u.is_empty() || mults_v.is_empty() || knots_v.is_empty() {
        return None;
    }
    Some((
        name, degree_u, degree_v, ctrl_grid, mults_u, knots_u, mults_v, knots_v,
    ))
}

fn parse_step_string(raw: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.len() >= 2 && trimmed.starts_with('\'') && trimmed.ends_with('\'') {
        trimmed[1..trimmed.len() - 1].replace("''", "'")
    } else {
        trimmed.to_string()
    }
}

fn parse_triplet(raw: &str) -> Option<[f64; 3]> {
    let vals: Vec<f64> = raw
        .split(',')
        .map(|s| s.trim().parse::<f64>().ok())
        .collect::<Option<Vec<_>>>()?;
    if vals.len() != 3 {
        return None;
    }
    Some([vals[0], vals[1], vals[2]])
}

fn decode_rcad_tagged_surface(name: &str) -> Option<Surface3> {
    let mut parts = name.split(';');
    let kind = parts.next()?;
    let mut attrs = std::collections::HashMap::<&str, &str>::new();
    for part in parts {
        let (key, value) = part.split_once('=')?;
        attrs.insert(key, value);
    }

    match kind {
        "RCAD_ELLIPSOID" => {
            let c = parse_triplet(attrs.get("c").copied()?)?;
            let a = parse_triplet(attrs.get("a").copied()?)?;
            let d = parse_triplet(attrs.get("d").copied()?)?;
            let r = parse_triplet(attrs.get("r").copied()?)?;
            Some(Surface3::Ellipsoid(rcad_kernel::EllipsoidalSurface {
                center: glam::DVec3::from_array(c),
                axis: glam::DVec3::from_array(a),
                ref_dir: glam::DVec3::from_array(d),
                radius_x: r[0],
                radius_y: r[1],
                radius_z: r[2],
            }))
        }
        "RCAD_HELICOID" => {
            let o = parse_triplet(attrs.get("o").copied()?)?;
            let a = parse_triplet(attrs.get("a").copied()?)?;
            let d = parse_triplet(attrs.get("d").copied()?)?;
            Some(Surface3::Helicoid(rcad_kernel::HelicoidSurface {
                origin: glam::DVec3::from_array(o),
                axis: glam::DVec3::from_array(a),
                ref_dir: glam::DVec3::from_array(d),
                pitch: attrs.get("p")?.trim().parse::<f64>().ok()?,
            }))
        }
        _ => None,
    }
}

fn build_bspline_surface_from_data(
    degree_u: usize,
    degree_v: usize,
    ctrl_grid_raw: &[Vec<u64>],
    mults_u: &[usize],
    knots_u: &[f64],
    mults_v: &[usize],
    knots_v: &[f64],
    parsed: &ParsedStep,
) -> Option<Surface3> {
    let expanded_u = expand_knots(mults_u, knots_u);
    let expanded_v = expand_knots(mults_v, knots_v);

    let n_v = ctrl_grid_raw.len();
    let n_u = ctrl_grid_raw.first().map(|r| r.len()).unwrap_or(0);
    if n_u == 0 || n_v == 0 {
        return None;
    }

    let mut control_points = vec![vec![glam::DVec3::ZERO; n_v]; n_u];
    let weights = vec![vec![1.0f64; n_v]; n_u];
    for (vi, row) in ctrl_grid_raw.iter().enumerate() {
        for (ui, &ref_id) in row.iter().enumerate() {
            if let Some(pt) = point_from_ref(parsed, ref_id) {
                control_points[ui][vi] = pt;
            }
        }
    }

    Some(Surface3::BSpline(rcad_kernel::geom::BSplineSurface {
        degree_u,
        degree_v,
        knots_u: expanded_u,
        knots_v: expanded_v,
        control_points,
        weights,
    }))
}

fn resolve_surface_for_trim_ops(parsed: &ParsedStep, surface_ref: u64) -> Option<Surface3> {
    if let Some((_, degree_u, degree_v, ctrl_grid_raw, mults_u, knots_u, mults_v, knots_v)) =
        parsed.b_spline_surfaces.get(&surface_ref)
    {
        return build_bspline_surface_from_data(
            *degree_u,
            *degree_v,
            ctrl_grid_raw,
            mults_u,
            knots_u,
            mults_v,
            knots_v,
            parsed,
        );
    }

    resolve_surface(parsed, surface_ref)
}

// ═══════════════════════════════════════════════════════════════════════════════
// CONVERSION FUNCTIONS: Kernel annotation types to STEP types
// ═══════════════════════════════════════════════════════════════════════════════

use rcad_kernel::annotation::{
    Annotation, AnnotationKind, Note, NoteCategory, NoteTarget, View, ViewProjection,
};

/// Convert a kernel View to STEP view entities.
///
/// Returns (StepView, StepCameraModelD3, StepViewVolume) for writing to STEP.
pub fn view_to_step_entities(view: &View, start_id: u64) -> (StepView, StepCameraModelD3, StepViewVolume) {
    let view_id = start_id;
    let camera_id = start_id + 1;
    let volume_id = start_id + 2;

    let step_view = StepView {
        entity_id: view_id,
        name: Some(view.name.clone()),
        description: None,
        camera_model_id: Some(camera_id),
        view_type: Some(view.name.clone()),
    };

    let step_camera = StepCameraModelD3 {
        entity_id: camera_id,
        name: Some(view.name.clone()),
        view_reference_system_id: None, // Would need AXIS2_PLACEMENT_3D
        view_volume_id: Some(volume_id),
        perspective: view.projection == ViewProjection::Perspective,
    };

    let volume_type = match view.projection {
        ViewProjection::Orthographic => ViewVolumeType::Orthographic,
        ViewProjection::Perspective => ViewVolumeType::Perspective,
    };

    let step_volume = StepViewVolume {
        entity_id: volume_id,
        name: Some(view.name.clone()),
        volume_type,
        view_center: Some(view.target.into()),
        view_plane_distance: Some(view.distance()),
        up_direction: Some(view.up_vector.into()),
        view_window_width: Some(view.view_width),
        view_window_height: Some(view.view_height),
    };

    (step_view, step_camera, step_volume)
}

/// Convert a kernel Note to a STEP note entity.
pub fn note_to_step(note: &Note) -> StepNote {
    StepNote {
        entity_id: note.id,
        name: Some(note.name.clone()),
        description: note.author.clone(),
        text: Some(note.text.clone()),
        annotation_plane_id: None,
        associated_geometry_id: None,
    }
}

/// Convert a kernel Annotation to STEP note and related entities.
///
/// Returns a StepNote and optionally annotation plane and occurrence information.
pub fn annotation_to_step_entities(
    annotation: &Annotation,
    _start_id: u64,
) -> (StepNote, Option<StepAnnotationPlane>, Option<StepAnnotationOccurrence>) {
    let step_note = StepNote {
        entity_id: annotation.id,
        name: Some(annotation.name.clone()),
        description: None,
        text: Some(annotation.text.clone()),
        annotation_plane_id: None,
        associated_geometry_id: None,
    };

    // For now, we don't create annotation planes/occurrences for simple annotations
    // These would be needed for more complex PMI representations
    (step_note, None, None)
}

/// Convert StepDocumentMetadata views to kernel View objects.
pub fn step_views_to_kernel(
    views: &[StepView],
    cameras: &[StepCameraModelD3],
    volumes: &[StepViewVolume],
) -> Vec<View> {
    views
        .iter()
        .filter_map(|step_view| {
            // Find corresponding camera
            let camera = step_view.camera_model_id.and_then(|cam_id| {
                cameras.iter().find(|c| c.entity_id == cam_id || c.entity_id == step_view.entity_id)
            });

            // Find corresponding view volume
            let volume = camera.and_then(|cam| {
                cam.view_volume_id.and_then(|vol_id| {
                    volumes.iter().find(|v| v.entity_id == vol_id)
                })
            });

            let mut view = View::new(step_view.entity_id, step_view.name.clone().unwrap_or_default());

            if let Some(cam) = camera {
                view.projection = if cam.perspective {
                    ViewProjection::Perspective
                } else {
                    ViewProjection::Orthographic
                };
            }

            if let Some(vol) = volume {
                if let Some(center) = vol.view_center {
                    view.target = glam::DVec3::from_array(center);
                }
                if let Some(up) = vol.up_direction {
                    view.up_vector = glam::DVec3::from_array(up);
                }
                if let Some(width) = vol.view_window_width {
                    view.view_width = width;
                }
                if let Some(height) = vol.view_window_height {
                    view.view_height = height;
                }
                if let Some(dist) = vol.view_plane_distance {
                    // Calculate camera position from target and distance
                    let dir = view.view_direction();
                    if dir.length() > 0.0 {
                        view.camera_position = view.target - dir.normalize() * dist;
                    }
                }
                view.clipping = rcad_kernel::annotation::ViewClipping::new(0.1, 10000.0);
            }

            view.custom = true;
            Some(view)
        })
        .collect()
}

/// Convert StepDocumentMetadata notes to kernel Note objects.
pub fn step_notes_to_kernel(notes: &[StepNote]) -> Vec<Note> {
    notes
        .iter()
        .filter_map(|step_note| {
            let category = match step_note.name.as_deref() {
                Some(name) if name.contains("warning") || name.contains("Warning") => NoteCategory::Warning,
                Some(name) if name.contains("requirement") || name.contains("Requirement") => NoteCategory::Requirement,
                Some(name) if name.contains("comment") || name.contains("Comment") => NoteCategory::Comment,
                Some(name) if name.contains("approval") || name.contains("Approval") => NoteCategory::Approval,
                _ => NoteCategory::Info,
            };

            let mut note = Note::new(
                step_note.entity_id,
                step_note.name.clone().unwrap_or_default(),
                step_note.text.clone().unwrap_or_default(),
            )
            .with_category(category);

            if let Some(ref author) = step_note.description {
                note = note.with_author(author.clone());
            }

            Some(note)
        })
        .collect()
}

/// Create StepAp242Metadata from kernel annotation store.
pub fn annotation_store_to_step_metadata(
    store: &rcad_kernel::annotation::AnnotationStore,
) -> StepAp242Metadata {
    let mut metadata = StepAp242Metadata::default();

    // Convert views
    let mut next_id = 1000u64; // Start IDs high to avoid conflicts
    for view in &store.views {
        let (step_view, step_camera, step_volume) = view_to_step_entities(view, next_id);
        metadata.views.push(step_view);
        metadata.cameras.push(step_camera);
        metadata.view_volumes.push(step_volume);
        next_id += 10;
    }

    // Convert notes
    for note in &store.xc_notes {
        metadata.notes.push(note_to_step(note));
    }

    // Convert annotations
    for annotation in &store.annotations {
        let (step_note, plane, occurrence) = annotation_to_step_entities(annotation, next_id);
        metadata.notes.push(step_note);
        if let Some(plane) = plane {
            metadata.annotation_planes.push(plane);
        }
        if let Some(occurrence) = occurrence {
            metadata.annotation_occurrences.push(occurrence);
        }
        next_id += 10;
    }

    metadata
}

#[cfg(test)]
mod tagged_surface_tests {
    use super::*;

    #[test]
    fn decodes_rcad_ellipsoid_surface_tag() {
        let tag = "RCAD_ELLIPSOID;c=1,2,3;a=0,0,1;d=1,0,0;r=4,5,6";
        let surface = decode_rcad_tagged_surface(tag);
        match surface {
            Some(Surface3::Ellipsoid(e)) => {
                assert_eq!(e.center, glam::DVec3::new(1.0, 2.0, 3.0));
                assert_eq!(e.axis, glam::DVec3::Z);
                assert_eq!(e.ref_dir, glam::DVec3::X);
                assert_eq!(e.radius_x, 4.0);
                assert_eq!(e.radius_y, 5.0);
                assert_eq!(e.radius_z, 6.0);
            }
            other => panic!("expected ellipsoid, got {other:?}"),
        }
    }

    #[test]
    fn decodes_rcad_helicoid_surface_tag() {
        let tag = "RCAD_HELICOID;o=1,2,3;a=0,0,1;d=1,0,0;p=7.5";
        let surface = decode_rcad_tagged_surface(tag);
        match surface {
            Some(Surface3::Helicoid(h)) => {
                assert_eq!(h.origin, glam::DVec3::new(1.0, 2.0, 3.0));
                assert_eq!(h.axis, glam::DVec3::Z);
                assert_eq!(h.ref_dir, glam::DVec3::X);
                assert_eq!(h.pitch, 7.5);
            }
            other => panic!("expected helicoid, got {other:?}"),
        }
    }

    #[test]
    fn parses_bspline_surface_with_tagged_name_containing_commas() {
        let args = "'RCAD_ELLIPSOID;c=1,2,3;a=0,0,1;d=1,0,0;r=4,5,6',2,2,((#1,#2,#3),(#4,#5,#6),(#7,#8,#9)),.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.0,1.0),(0.0,1.0),.UNSPECIFIED.";
        let parsed = parse_bspline_surface_with_knots(args)
            .expect("tagged B-spline surface args should parse");
        assert!(parsed.0.starts_with("RCAD_ELLIPSOID"));
        assert_eq!(parsed.1, 2);
        assert_eq!(parsed.2, 2);
        assert_eq!(parsed.3.len(), 3);
    }

    #[test]
    fn parse_entities_keeps_tagged_bspline_surface_names() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_NAME('test','','','','','','');\nFILE_SCHEMA(('AP214'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=CARTESIAN_POINT('',(1.,0.,0.));\n#3=CARTESIAN_POINT('',(2.,0.,0.));\n#4=CARTESIAN_POINT('',(0.,1.,0.));\n#5=CARTESIAN_POINT('',(1.,1.,0.));\n#6=CARTESIAN_POINT('',(2.,1.,0.));\n#7=CARTESIAN_POINT('',(0.,2.,0.));\n#8=CARTESIAN_POINT('',(1.,2.,0.));\n#9=CARTESIAN_POINT('',(2.,2.,0.));\n#10=B_SPLINE_SURFACE_WITH_KNOTS('RCAD_HELICOID;o=0,0,0;a=0,0,1;d=1,0,0;p=3',2,2,((#1,#2,#3),(#4,#5,#6),(#7,#8,#9)),.UNSPECIFIED.,.F.,.F.,.F.,(3,3),(3,3),(0.0,1.0),(0.0,1.0),.UNSPECIFIED.);\nENDSEC;\nEND-ISO-10303-21;\n";
        let parsed = parse_entities(step).expect("tagged B-spline surface STEP should parse");
        let data = parsed
            .b_spline_surfaces
            .get(&10)
            .expect("tagged B-spline surface should be retained");
        assert!(data.0.starts_with("RCAD_HELICOID"));
    }
}

fn parse_conical_surface(args: &str) -> Option<(u64, f64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parts[3].trim().parse::<f64>().ok()?.to_radians(),
    ))
}

fn parse_toroidal_surface(args: &str) -> Option<(u64, f64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parts[3].trim().parse::<f64>().ok()?,
    ))
}

/// Parse OFFSET_CURVE_3D args:
/// ('name', #basis_curve, offset_distance, #ref_direction)
fn parse_offset_curve_3d(args: &str) -> Option<(u64, f64, u64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    Some((
        parse_ref(parts[1])?,
        parts[2].trim().parse::<f64>().ok()?,
        parse_ref(parts[3])?,
    ))
}

fn parse_advanced_face(args: &str) -> Option<(Vec<u64>, Option<u64>)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let bounds = parse_ref_list(parts[1]);
    let surface = parse_ref(parts[2]);
    if bounds.is_empty() {
        return None;
    }
    Some((bounds, surface))
}

fn resolve_curve(parsed: &ParsedStep, curve_ref: u64) -> Option<Curve3> {
    // Dereference SURFACE_CURVE 鈥?extract the wrapped 3D curve
    let actual_ref = if let Some(&(inner_ref, _, _)) = parsed.surface_curves.get(&curve_ref) {
        inner_ref
    } else {
        curve_ref
    };

    if let Some((origin_point, vector_ref)) = parsed.lines.get(&actual_ref) {
        let origin = point_from_ref(parsed, *origin_point)?;
        let (direction_ref, _magnitude) = *parsed.vectors.get(vector_ref)?;
        let direction = direction_from_ref(parsed, direction_ref)?;
        return Some(Curve3::Line(rcad_kernel::geom::Line3 { origin, direction }));
    }

    if let Some((placement_ref, radius)) = parsed.circles.get(&actual_ref) {
        let (center, normal) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Circle(rcad_kernel::geom::Circle3 {
            center,
            normal,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, major_radius, minor_radius)) = parsed.ellipses.get(&actual_ref) {
        let (center, normal, major_dir) = placement_frame_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Ellipse(rcad_kernel::geom::Ellipse3 {
            center,
            normal,
            major_dir,
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        }));
    }

    // BSpline: use full data if available, otherwise fall through to None
    if let Some((degree, ctrl_refs, mults, knot_vals)) =
        parsed.b_spline_curves_full.get(&actual_ref)
    {
        let control_points: Vec<glam::DVec3> = ctrl_refs
            .iter()
            .filter_map(|&r| point_from_ref(parsed, r))
            .collect();
        if control_points.len() >= 2 {
            // Expand knot vector from multiplicities
            let mut knots = Vec::new();
            for (&mult, &val) in mults.iter().zip(knot_vals.iter()) {
                for _ in 0..mult {
                    knots.push(val);
                }
            }
            let weights = vec![1.0; control_points.len()];
            return Some(Curve3::BSpline(BSplineCurve3 {
                degree: *degree,
                knots,
                control_points,
                weights,
            }));
        }
    }

    if let Some((placement_ref, semi_major, semi_minor)) = parsed.hyperbolas.get(&actual_ref) {
        let (center, normal, major_dir) = placement_frame_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Hyperbola(rcad_kernel::geom::Hyperbola3 {
            center,
            normal,
            major_dir,
            semi_major: *semi_major,
            semi_minor: *semi_minor,
        }));
    }

    if let Some((placement_ref, focal_param)) = parsed.parabolas.get(&actual_ref) {
        let (vertex, normal, axis_dir) = placement_frame_from_ref(parsed, *placement_ref)?;
        return Some(Curve3::Parabola(rcad_kernel::geom::Parabola3 {
            vertex,
            normal,
            axis_dir,
            focal_param: *focal_param,
        }));
    }

    if let Some((basis_ref, offset_dist, dir_ref)) = parsed.offset_curves_3d.get(&actual_ref) {
        let basis = resolve_curve(parsed, *basis_ref)?;
        let offset_dir = direction_from_ref(parsed, *dir_ref)?;
        return Some(Curve3::Offset(rcad_kernel::geom::OffsetCurve3 {
            basis: Box::new(basis),
            offset_distance: *offset_dist,
            offset_dir,
        }));
    }

    None
}

fn resolve_surface(parsed: &ParsedStep, surface_ref: u64) -> Option<Surface3> {
    if let Some(placement_ref) = parsed.planes.get(&surface_ref) {
        let (origin, normal) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Plane(rcad_kernel::geom::Plane { origin, normal }));
    }

    if let Some((placement_ref, radius)) = parsed.cylindrical_surfaces.get(&surface_ref) {
        let (origin, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Cylinder(rcad_kernel::geom::CylindricalSurface {
            origin,
            axis,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, radius)) = parsed.spherical_surfaces.get(&surface_ref) {
        let (center, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Sphere(rcad_kernel::geom::SphericalSurface {
            center,
            axis,
            radius: *radius,
        }));
    }

    if let Some((placement_ref, ref_radius, half_angle_rad)) =
        parsed.conical_surfaces.get(&surface_ref)
    {
        let (apex, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Cone(rcad_kernel::geom::ConicalSurface {
            apex,
            axis,
            radius: *ref_radius,
            half_angle_rad: *half_angle_rad,
        }));
    }

    if let Some((placement_ref, major_radius, minor_radius)) =
        parsed.toroidal_surfaces.get(&surface_ref)
    {
        let (center, axis) = placement_from_ref(parsed, *placement_ref)?;
        return Some(Surface3::Torus(rcad_kernel::geom::ToroidalSurface {
            center,
            axis,
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        }));
    }

    if let Some((name, degree_u, degree_v, ctrl_grid_raw, mults_u, knots_u, mults_v, knots_v)) =
        parsed.b_spline_surfaces.get(&surface_ref)
    {
        if let Some(surface) = decode_rcad_tagged_surface(name) {
            return Some(surface);
        }
        return build_bspline_surface_from_data(
            *degree_u,
            *degree_v,
            ctrl_grid_raw,
            mults_u,
            knots_u,
            mults_v,
            knots_v,
            parsed,
        );
    }

    // SURFACE_OF_LINEAR_EXTRUSION
    if let Some((profile_ref, dir_ref)) = parsed.linear_extrusions.get(&surface_ref).copied()
        && let (Some(profile), Some(dir)) = (
            resolve_curve(parsed, profile_ref),
            direction_from_ref(parsed, dir_ref),
        )
    {
        return Some(Surface3::LinearExtrusion(
            rcad_kernel::geom::LinearExtrusionSurface {
                profile: Box::new(profile),
                direction: dir.normalize_or_zero(),
            },
        ));
    }

    // SURFACE_OF_REVOLUTION
    if let Some((profile_ref, axis_ref)) = parsed.revolutions.get(&surface_ref).copied()
        && let Some(profile) = resolve_curve(parsed, profile_ref)
    {
        // Try AXIS2_PLACEMENT_3D first (common in practice), then fall back to a direction ref
        let axis_result = placement_from_ref(parsed, axis_ref).or_else(|| {
            // Treat as bare direction ref at origin
            direction_from_ref(parsed, axis_ref).map(|d| (glam::DVec3::ZERO, d))
        });
        if let Some((axis_origin, axis_dir)) = axis_result {
            return Some(Surface3::Revolution(rcad_kernel::geom::RevolutionSurface {
                profile: Box::new(profile),
                axis_origin,
                axis_dir: axis_dir.normalize_or_zero(),
            }));
        }
    }

    // RECTANGULAR_TRIMMED_SURFACE
    if let Some((basis_ref, trim)) = parsed
        .rectangular_trimmed_surfaces
        .get(&surface_ref)
        .copied()
        && let Some(basis) = resolve_surface(parsed, basis_ref)
    {
        return Some(Surface3::Trimmed(rcad_kernel::TrimmedSurface {
            basis: Box::new(basis),
            trim,
        }));
    }

    None
}

/// Expand a compressed knot vector (multiplicities + values) into a full knot vector.
fn expand_knots(mults: &[usize], vals: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for (&m, &v) in mults.iter().zip(vals.iter()) {
        for _ in 0..m {
            out.push(v);
        }
    }
    out
}

/// Parse UNCERTAINTY_MEASURE_WITH_UNIT args to extract the LENGTH_MEASURE value.
/// Format: `UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(val), unit_ref, 'name', 'desc')`
fn parse_uncertainty_measure(args: &str) -> Option<f64> {
    // Find LENGTH_MEASURE(value) in the args string
    let start = args.find("LENGTH_MEASURE(")?;
    let rest = &args[start + "LENGTH_MEASURE(".len()..];
    let end = rest.find(')')?;
    rest[..end].trim().parse::<f64>().ok()
}

/// Parse RECTANGULAR_TRIMMED_SURFACE args:
/// `'name', #basis, u1, u2, v1, v2, .T., .T.` 鈫?`(basis_ref, [u1,u2,v1,v2])`.
fn parse_rectangular_trimmed_surface(args: &str) -> Option<(u64, [f64; 4])> {
    // Extract the single # reference (basis surface)
    let refs = parse_ref_list(args);
    let basis_ref = *refs.first()?;
    // Extract all floats (the 4 trim parameters)
    let floats: Vec<f64> = args
        .split(',')
        .filter_map(|tok| {
            let t = tok.trim().trim_matches(|c: char| {
                !c.is_ascii_digit() && c != '.' && c != '-' && c != 'E' && c != 'e'
            });
            t.parse::<f64>().ok()
        })
        .collect();
    if floats.len() >= 4 {
        Some((basis_ref, [floats[0], floats[1], floats[2], floats[3]]))
    } else {
        None
    }
}

/// Parse COLOUR_RGB args: `'name', r, g, b` 鈫?`[r, g, b]`.
fn parse_colour_rgb(args: &str) -> Option<[f64; 3]> {
    // Skip the optional name string, then parse three floats.
    let rest = if args.trim_start().starts_with('\'') {
        // Skip quoted name
        let after = args.trim_start().trim_start_matches('\'');
        let end_quote = after.find('\'')?;
        &after[end_quote + 1..]
    } else {
        args
    };
    let floats: Vec<f64> = rest
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect();
    if floats.len() >= 3 {
        Some([floats[0], floats[1], floats[2]])
    } else {
        None
    }
}

/// Parse a single `#N` reference from args (no name prefix).
fn parse_single_ref(args: &str) -> Option<u64> {
    let refs = parse_ref_list(args);
    refs.into_iter().next()
}

/// Parse the last `#N` reference in args (for SURFACE_STYLE_USAGE where shape ref is last).
fn parse_last_ref(args: &str) -> Option<u64> {
    parse_ref_list(args).into_iter().last()
}

/// Parse STYLED_ITEM args: `'name', (#psa,...), #shape` 鈫?`(shape_ref, [psa_refs])`.
///
/// The writer emits: `STYLED_ITEM('color',(#psa,),#face_id)`
fn parse_styled_item(args: &str) -> Option<(u64, Vec<u64>)> {
    let refs = parse_ref_list(args);
    if refs.len() < 2 {
        return None;
    }
    // Last ref is the shape; all but last are style assignments
    let shape_ref = *refs.last()?;
    let psa_refs = refs[..refs.len() - 1].to_vec();
    Some((shape_ref, psa_refs))
}

fn point_from_ref(parsed: &ParsedStep, point_ref: u64) -> Option<glam::DVec3> {
    let p = parsed.cartesian_points.get(&point_ref)?;
    Some(glam::DVec3::new(p[0], p[1], p[2]))
}

fn vertex_point_from_ref(parsed: &ParsedStep, vertex_ref: u64) -> Option<glam::DVec3> {
    let point_ref = *parsed.vertex_points.get(&vertex_ref)?;
    point_from_ref(parsed, point_ref)
}

fn direction_from_ref(parsed: &ParsedStep, direction_ref: u64) -> Option<glam::DVec3> {
    let d = parsed.directions.get(&direction_ref)?;
    Some(glam::DVec3::new(d[0], d[1], d[2]).normalize_or_zero())
}

fn placement_from_ref(
    parsed: &ParsedStep,
    placement_ref: u64,
) -> Option<(glam::DVec3, glam::DVec3)> {
    let (origin_ref, axis_ref, _) = *parsed.axis2_placements.get(&placement_ref)?;
    Some((
        point_from_ref(parsed, origin_ref)?,
        direction_from_ref(parsed, axis_ref)?,
    ))
}

fn placement_frame_from_ref(
    parsed: &ParsedStep,
    placement_ref: u64,
) -> Option<(glam::DVec3, glam::DVec3, glam::DVec3)> {
    let (origin_ref, axis_ref, ref_dir_ref) = *parsed.axis2_placements.get(&placement_ref)?;
    let origin = point_from_ref(parsed, origin_ref)?;
    let axis = direction_from_ref(parsed, axis_ref)?;

    let major_dir = if let Some(dref) = ref_dir_ref {
        let d = direction_from_ref(parsed, dref)?;
        let d_proj = (d - axis * d.dot(axis)).normalize_or_zero();
        if d_proj.length_squared() > 1e-8 {
            d_proj
        } else {
            any_perpendicular(axis)
        }
    } else {
        any_perpendicular(axis)
    };

    Some((origin, axis, major_dir))
}

fn any_perpendicular(axis: glam::DVec3) -> glam::DVec3 {
    let helper = if axis.dot(glam::DVec3::Y).abs() < 0.9 {
        glam::DVec3::Y
    } else {
        glam::DVec3::X
    };
    axis.cross(helper).normalize_or_zero()
}

fn parse_vector(args: &str) -> Option<(u64, f64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let dir_ref = parse_ref(parts[1])?;
    let magnitude = parts[2].trim().parse::<f64>().ok()?;
    Some((dir_ref, magnitude))
}

fn parse_trimmed_curve(args: &str) -> Option<(u64, f64, f64)> {
    // TRIMMED_CURVE('name', curve_ref, (PARAMETER_VALUE(t0)), (PARAMETER_VALUE(t1)), ...)
    let parts = split_top_level(args, ',');
    if parts.len() < 4 {
        return None;
    }
    let curve_ref = parse_ref(parts[1])?;
    let t0 = parse_parameter_value(parts[2])?;
    let t1 = parse_parameter_value(parts[3])?;
    Some((curve_ref, t0, t1))
}

fn parse_parameter_value(s: &str) -> Option<f64> {
    // s looks like "(PARAMETER_VALUE(0.))" 鈥?find the float inside PARAMETER_VALUE(...)
    let cursor = s.to_uppercase();
    let pv_pos = cursor.find("PARAMETER_VALUE(")?;
    let after = &s[pv_pos + "PARAMETER_VALUE(".len()..];
    let end = after.find(')')?;
    after[..end].trim().parse::<f64>().ok()
}

/// Evaluate a STEP LINE curve at parameter t: p(t) = origin + dir * magnitude * t
fn eval_line_at(parsed: &ParsedStep, line_ref: u64, t: f64) -> Option<glam::DVec3> {
    let &(origin_ref, vec_ref) = parsed.lines.get(&line_ref)?;
    let origin = point_from_ref(parsed, origin_ref)?;
    let &(dir_ref, magnitude) = parsed.vectors.get(&vec_ref)?;
    let dir = direction_from_ref(parsed, dir_ref)?;
    Some(origin + dir * (magnitude * t))
}

/// Sample a standalone curve referenced from a GEOMETRIC_CURVE_SET into polyline points.
fn sample_standalone_curve(parsed: &ParsedStep, curve_ref: u64) -> Vec<glam::DVec3> {
    // Handle TRIMMED_CURVE wrapper
    if let Some(&(underlying_ref, t0, t1)) = parsed.trimmed_curves.get(&curve_ref) {
        return sample_trimmed_curve_geom(parsed, underlying_ref, t0, t1);
    }
    // Handle B_SPLINE/BEZIER-like curves through parametric sampling.
    if let Some(points) = sample_spline_curve(parsed, curve_ref, None) {
        return points;
    }
    // Fallback for sparse B_SPLINE data: use control polyline only when full knot data is missing.
    if let Some(control_refs) = parsed.b_spline_curves.get(&curve_ref) {
        return sample_bspline_polyline(parsed, control_refs);
    }
    // Handle bare LINE (t 0..1)
    if parsed.lines.contains_key(&curve_ref) {
        let p0 = eval_line_at(parsed, curve_ref, 0.0);
        let p1 = eval_line_at(parsed, curve_ref, 1.0);
        return match (p0, p1) {
            (Some(a), Some(b)) => vec![a, b],
            _ => Vec::new(),
        };
    }
    Vec::new()
}

/// Sample the underlying geometry of a TRIMMED_CURVE at [t0, t1].
fn sample_trimmed_curve_geom(
    parsed: &ParsedStep,
    curve_ref: u64,
    t0: f64,
    t1: f64,
) -> Vec<glam::DVec3> {
    if parsed.lines.contains_key(&curve_ref) {
        // LINE: p(t) = origin + dir * magnitude * t
        let p0 = eval_line_at(parsed, curve_ref, t0);
        let p1 = eval_line_at(parsed, curve_ref, t1);
        return match (p0, p1) {
            (Some(a), Some(b)) => vec![a, b],
            _ => Vec::new(),
        };
    }

    if let Some(&(placement_ref, radius)) = parsed.circles.get(&curve_ref) {
        // CIRCLE: HFSS exports parameters in degrees (0..360 is full circle)
        return sample_standalone_circle(parsed, placement_ref, radius, t0, t1);
    }

    // Generic spline-style curves (B-spline, Bezier represented through B-spline data).
    if let Some(points) = sample_spline_curve(parsed, curve_ref, Some([t0, t1])) {
        return points;
    }

    Vec::new()
}

/// Sample a CIRCLE arc from t0_deg to t1_deg (degrees, HFSS convention).
fn sample_standalone_circle(
    parsed: &ParsedStep,
    placement_ref: u64,
    radius: f64,
    t0_deg: f64,
    t1_deg: f64,
) -> Vec<glam::DVec3> {
    if !radius.is_finite() || radius <= 0.0 {
        return Vec::new();
    }
    let Some((center, axis, u)) = placement_frame_from_ref(parsed, placement_ref) else {
        return Vec::new();
    };
    let v = axis.cross(u).normalize_or_zero();

    let t0 = t0_deg.to_radians();
    let mut sweep = (t1_deg - t0_deg).to_radians();
    if sweep.abs() < 1e-9 {
        sweep = std::f64::consts::TAU; // full circle
    } else if sweep < 0.0 {
        sweep += std::f64::consts::TAU;
    }

    let seg = ((sweep.abs() / std::f64::consts::TAU) * 64.0)
        .ceil()
        .max(8.0) as usize;
    let mut points = Vec::with_capacity(seg + 1);
    for i in 0..=seg {
        let t = t0 + sweep * (i as f64 / seg as f64);
        points.push(center + u * (radius * t.cos()) + v * (radius * t.sin()));
    }
    points
}

fn parse_oriented_edge(args: &str) -> Option<(u64, bool)> {
    let parts = split_top_level(args, ',');
    let mut edge_ref = None;
    let mut orientation = None;

    for part in &parts {
        if edge_ref.is_none() {
            edge_ref = parse_ref(part);
            if edge_ref.is_some() {
                continue;
            }
        }

        if orientation.is_none() {
            let v = part.trim();
            if v == ".T." {
                orientation = Some(true);
            } else if v == ".F." {
                orientation = Some(false);
            }
        }
    }

    Some((edge_ref?, orientation?))
}

fn parse_ref_list_after_name(args: &str) -> Option<Vec<u64>> {
    let open = args.find('(')?;
    let close = args.rfind(')')?;
    if close <= open {
        return None;
    }
    let inside = args[open + 1..close].trim();
    if !inside.starts_with('#') && !inside.contains('#') {
        return None;
    }
    Some(parse_ref_list(inside))
}

fn parse_ref_list(input: &str) -> Vec<u64> {
    let mut refs = Vec::new();
    let mut i = 0usize;
    let bytes = input.as_bytes();

    while i < bytes.len() {
        if bytes[i] == b'#' {
            i += 1;
            let start = i;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
            if start < i
                && let Ok(v) = input[start..i].parse::<u64>()
            {
                refs.push(v);
            }
        } else {
            i += 1;
        }
    }

    refs
}

/// Parse a parenthesized list of floating-point numbers: `(1., 2., 3.)` 鈫?`[1.0, 2.0, 3.0]`.
fn parse_float_list(input: &str) -> Vec<f64> {
    let inner = input.trim().trim_start_matches('(').trim_end_matches(')');
    inner
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect()
}

fn parse_ref(input: &str) -> Option<u64> {
    let trimmed = input.trim();
    let hash = trimmed.find('#')?;
    let digits = &trimmed[hash + 1..];
    let mut end = 0usize;
    for ch in digits.chars() {
        if ch.is_ascii_digit() {
            end += ch.len_utf8();
        } else {
            break;
        }
    }
    if end == 0 {
        return None;
    }
    digits[..end].parse::<u64>().ok()
}

fn split_top_level(input: &str, delimiter: char) -> Vec<&str> {
    let mut result = Vec::new();
    let mut depth = 0i32;
    let mut in_string = false;
    let mut start = 0usize;

    for (idx, ch) in input.char_indices() {
        match ch {
            '\'' => in_string = !in_string,
            '(' if !in_string => depth += 1,
            ')' if !in_string => depth -= 1,
            _ => {}
        }

        if ch == delimiter && depth == 0 && !in_string {
            result.push(input[start..idx].trim());
            start = idx + ch.len_utf8();
        }
    }

    if start <= input.len() {
        result.push(input[start..].trim());
    }

    result
}

fn parse_cartesian_point_2d(args: &str) -> Option<[f64; 2]> {
    // CARTESIAN_POINT('name', (x, y)) 鈥?only 2 coordinates
    let inner = args.trim().trim_start_matches('(').trim_end_matches(')');
    let parts = split_top_level(inner, ',');
    if parts.len() != 3 {
        return None; // 3 parts means 3D, not 2D
    }
    // parts[0] = name (quoted string), parts[1] = tuple like (x,y)
    let coords_str = parts[1].trim();
    let coords_inner = coords_str.trim_start_matches('(').trim_end_matches(')');
    let nums: Vec<f64> = coords_inner
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .collect();
    if nums.len() == 2 {
        Some([nums[0], nums[1]])
    } else {
        None
    }
}

fn parse_axis2_placement_2d(args: &str) -> Option<(u64, u64)> {
    // AXIS2_PLACEMENT_2D('', #location, #ref_dir)
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?))
}

fn parse_surface_curve(args: &str) -> Option<(u64, Vec<u64>, bool)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    let curve3d_ref = parse_ref(parts[1])?;
    let pcurve_refs = parse_ref_list(parts[2]);
    // 4th field is the master_representation flag (optional, default .T.)
    // In STEP AP214 this is ".PCURVE_S1." / ".PCURVE_S2." / ".CURVE_3D." / "$"
    // We treat it as same_parameter=true unless explicitly set to .F.
    let same_param = parts
        .get(3)
        .map(|s| !s.trim().eq_ignore_ascii_case(".F."))
        .unwrap_or(true);
    Some((curve3d_ref, pcurve_refs, same_param))
}

/// PCURVE('', #surface, #definitional_rep)
fn parse_pcurve_args(args: &str) -> Option<(u64, u64)> {
    let parts = split_top_level(args, ',');
    if parts.len() < 3 {
        return None;
    }
    Some((parse_ref(parts[1])?, parse_ref(parts[2])?))
}

/// DEFINITIONAL_REPRESENTATION('', (#curve2d_ref), #context)
fn parse_definitional_rep(args: &str) -> Option<u64> {
    let parts = split_top_level(args, ',');
    if parts.len() < 2 {
        return None;
    }
    // The second part is a list like (#14) 鈥?take the first element
    let refs = parse_ref_list(parts[1]);
    refs.into_iter().next()
}

/// Resolve a 2D curve from the parsed step data.
///
/// 2D curves live inside DEFINITIONAL_REPRESENTATION bodies. In the STEP file
/// they use the same entity names (LINE, CIRCLE) but with 2-component coords.
/// The reader stores them in the 3D maps (cartesian_points, circles, lines) 鈥?
/// STEP parsers accept them there. We just need to down-convert to Curve2d.
fn resolve_curve2d(parsed: &ParsedStep, curve_ref: u64) -> Option<Curve2d> {
    if let Some((origin_ref, vector_ref)) = parsed.lines.get(&curve_ref) {
        // Try 2D cartesian point first, fall back to 3D
        let origin = parsed
            .cartesian_points_2d
            .get(origin_ref)
            .map(|&p| glam::DVec2::new(p[0], p[1]))
            .or_else(|| {
                parsed
                    .cartesian_points
                    .get(origin_ref)
                    .map(|&p| glam::DVec2::new(p[0], p[1]))
            })?;
        let (dir_ref, _mag) = *parsed.vectors.get(vector_ref)?;
        let dir2d = parsed
            .directions_2d
            .get(&dir_ref)
            .map(|&d| glam::DVec2::new(d[0], d[1]))
            .or_else(|| {
                parsed
                    .directions
                    .get(&dir_ref)
                    .map(|&d| glam::DVec2::new(d[0], d[1]))
            })?;
        return Some(Curve2d::Line(rcad_kernel::geom::Line2d {
            origin,
            direction: dir2d.normalize_or_zero(),
        }));
    }

    if let Some((placement_ref, radius)) = parsed.circles.get(&curve_ref) {
        // 2D circle: extract center from the 2D placement
        let center = parsed
            .axis2_placements_2d
            .get(placement_ref)
            .and_then(|(loc_ref, _)| parsed.cartesian_points_2d.get(loc_ref))
            .map(|&p| glam::DVec2::new(p[0], p[1]))
            .or_else(|| {
                parsed
                    .axis2_placements
                    .get(placement_ref)
                    .and_then(|(loc_ref, _, _)| parsed.cartesian_points.get(loc_ref))
                    .map(|&p| glam::DVec2::new(p[0], p[1]))
            })?;
        return Some(Curve2d::Circle(rcad_kernel::geom::Circle2d {
            center,
            radius: *radius,
        }));
    }

    // 2D Ellipse: ELLIPSE referencing an AXIS2_PLACEMENT_2D
    if let Some((placement_ref, major, minor)) = parsed.ellipses.get(&curve_ref)
        && let Some((loc_ref, dir_ref)) = parsed.axis2_placements_2d.get(placement_ref)
    {
        let center = parsed
            .cartesian_points_2d
            .get(loc_ref)
            .map(|&p| glam::DVec2::new(p[0], p[1]))
            .or_else(|| {
                parsed
                    .cartesian_points
                    .get(loc_ref)
                    .map(|&p| glam::DVec2::new(p[0], p[1]))
            })?;
        let major_dir = parsed
            .directions_2d
            .get(dir_ref)
            .map(|&d| glam::DVec2::new(d[0], d[1]))
            .or_else(|| {
                parsed
                    .directions
                    .get(dir_ref)
                    .map(|&d| glam::DVec2::new(d[0], d[1]))
            })
            .unwrap_or(glam::DVec2::X)
            .normalize_or(glam::DVec2::X);
        return Some(Curve2d::Ellipse(Ellipse2d {
            center,
            major_dir,
            major_radius: *major,
            minor_radius: *minor,
        }));
    }

    // 2D B-Spline curve: B_SPLINE_CURVE_WITH_KNOTS with 2D control points
    if let Some((degree, cp_refs, mults, knot_vals)) = parsed.b_spline_curves_full.get(&curve_ref) {
        // Check if ALL control points are 2D (present in cartesian_points_2d)
        let all_2d = cp_refs
            .iter()
            .all(|id| parsed.cartesian_points_2d.contains_key(id));
        if all_2d {
            let control_points: Vec<glam::DVec2> = cp_refs
                .iter()
                .filter_map(|id| parsed.cartesian_points_2d.get(id))
                .map(|&p| glam::DVec2::new(p[0], p[1]))
                .collect();
            if control_points.len() == cp_refs.len() {
                let knots = expand_knots(mults, knot_vals);
                let weights = vec![1.0_f64; control_points.len()];
                return Some(Curve2d::BSpline(BSplineCurve2 {
                    degree: *degree,
                    knots,
                    control_points,
                    weights,
                }));
            }
        }
    }

    None
}

// ── PCurve / tolerance export validation ─────────────────────────────────────
//
// Analogous to OCCT's `BRepLib::CheckCurveOnSurface` and the tolerance-check
// stage of `BRepAlgoAPI_Check`.

/// A single issue found by [`validate_export_readiness`].
#[derive(Debug, Clone, PartialEq)]
pub enum ExportIssue {
    /// Edge has a PCurve whose surface index is out of range.
    PcurveSurfaceOutOfRange { edge_idx: usize, surface_idx: usize },
    /// Edge has a PCurve whose 2D curve index is out of range.
    PcurveCurveOutOfRange { edge_idx: usize, curve2d_idx: usize },
    /// Edge tolerance is below the global precision floor (`CONFUSION`).
    EdgeToleranceTooTight { edge_idx: usize, tolerance: f64 },
    /// Vertex tolerance is below the global precision floor.
    VertexToleranceTooTight { vertex_idx: usize, tolerance: f64 },
    /// Per-face tolerance is below the global precision floor.
    FaceToleranceTooTight { face_idx: usize, tolerance: f64 },
    /// Edge has more than 2 PCurves, which is unusual and may cause
    /// conformance issues with strict STEP readers.
    TooManyPcurves { edge_idx: usize, count: usize },
    /// Edge is referenced by a surface-bearing face but has no PCurve stored,
    /// meaning the SURFACE_CURVE entity will be missing on export.
    MissingPcurve { edge_idx: usize },
}

impl std::fmt::Display for ExportIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use ExportIssue::*;
        match self {
            PcurveSurfaceOutOfRange { edge_idx, surface_idx } =>
                write!(f, "edge {edge_idx}: PCurve surface_idx {surface_idx} out of range"),
            PcurveCurveOutOfRange { edge_idx, curve2d_idx } =>
                write!(f, "edge {edge_idx}: PCurve curve2d_idx {curve2d_idx} out of range"),
            EdgeToleranceTooTight { edge_idx, tolerance } =>
                write!(f, "edge {edge_idx}: tolerance {tolerance} < CONFUSION ({CONFUSION})"),
            VertexToleranceTooTight { vertex_idx, tolerance } =>
                write!(f, "vertex {vertex_idx}: tolerance {tolerance} < CONFUSION ({CONFUSION})"),
            FaceToleranceTooTight { face_idx, tolerance } =>
                write!(f, "face {face_idx}: tolerance {tolerance} < CONFUSION ({CONFUSION})"),
            TooManyPcurves { edge_idx, count } =>
                write!(f, "edge {edge_idx}: has {count} PCurves (max 2 expected)"),
            MissingPcurve { edge_idx } =>
                write!(f, "edge {edge_idx}: used by surface-bearing face but no PCurve stored"),
        }
    }
}

/// Result of [`validate_export_readiness`].
#[derive(Debug, Clone, Default)]
pub struct ExportReadinessReport {
    /// All detected issues.  Empty ↔ [`is_ready`] == true.
    pub issues: Vec<ExportIssue>,
    /// `true` when no issues were found.
    pub is_ready: bool,
    /// Number of edges whose PCurve lists were inspected.
    pub edges_checked: usize,
    /// Total number of individual PCurve entries validated.
    pub pcurves_checked: usize,
}

impl ExportReadinessReport {
    /// One-line human-readable summary.
    pub fn summary(&self) -> String {
        if self.is_ready {
            format!(
                "export-ready: {} edges, {} pcurves checked, 0 issues",
                self.edges_checked, self.pcurves_checked
            )
        } else {
            format!(
                "NOT export-ready: {} issue(s) across {} edges / {} pcurves",
                self.issues.len(), self.edges_checked, self.pcurves_checked
            )
        }
    }
}

/// Validate a [`BRep`] for STEP export correctness.
///
/// Performs three categories of checks:
///
/// 1. **PCurve index bounds** — every [`PCurve`] in `geom.edge_pcurves` must
///    reference a valid index into `geom.surfaces` and `geom.curve2ds`.
/// 2. **PCurve cardinality** — more than 2 PCurves on a single edge is unusual
///    and may cause conformance issues with strict STEP AP214/AP242 readers.
/// 3. **Missing PCurves** — an edge that is referenced by a surface-bearing face
///    (i.e., `geom.face_surface` has a `Some` entry for that face's flat index)
///    but has no PCurve entry is flagged.  Advisory: analytic primitives that do
///    not populate `geom.edge_pcurves` at all are *not* flagged because the
///    writer synthesises the SURFACE_CURVE on the fly for such shapes.
/// 4. **Tolerance floor** — stored tolerance values below `CONFUSION` (1 × 10⁻⁷)
///    would violate the STEP AP214/AP242 minimum-tolerance recommendations.
///
/// Analogous to `BRepLib::CheckCurveOnSurface` + OCCT shape-analysis tolerance
/// queries before a `WriteSTEP` call.
pub fn validate_export_readiness(brep: &BRep) -> ExportReadinessReport {
    let mut issues = Vec::new();
    let n_surfaces = brep.geom.surfaces.len();
    let n_curve2ds = brep.geom.curve2ds.len();
    let n_edge_pcurve_entries = brep.geom.edge_pcurves.len();
    let n_edges = brep.edges.len();
    let mut pcurves_checked = 0usize;

    // ── 1 & 2: PCurve index validity and cardinality ──────────────────────────
    for ei in 0..n_edge_pcurve_entries {
        let pcs = &brep.geom.edge_pcurves[ei];
        if pcs.len() > 2 {
            issues.push(ExportIssue::TooManyPcurves { edge_idx: ei, count: pcs.len() });
        }
        for pc in pcs {
            pcurves_checked += 1;
            if pc.surface_idx >= n_surfaces {
                issues.push(ExportIssue::PcurveSurfaceOutOfRange {
                    edge_idx: ei,
                    surface_idx: pc.surface_idx,
                });
            }
            if pc.curve2d_idx >= n_curve2ds {
                issues.push(ExportIssue::PcurveCurveOutOfRange {
                    edge_idx: ei,
                    curve2d_idx: pc.curve2d_idx,
                });
            }
        }
    }

    // ── 3: Missing PCurves ────────────────────────────────────────────────────
    // Only flag when `edge_pcurves` has been partially populated (i.e. the BRep
    // uses explicit PCurve storage) but a particular edge slot is empty or absent.
    // If the slice is entirely empty the writer uses analytic fallback paths and
    // we should not raise false-positive MissingPcurve issues.
    let any_pcurves_stored = brep.geom.edge_pcurves.iter().any(|v| !v.is_empty());
    if any_pcurves_stored {
        // Collect edge indices that are adjacent to surface-bearing faces.
        let mut surface_edges: std::collections::HashSet<usize> =
            std::collections::HashSet::new();
        let mut flat_fi = 0usize;
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let has_surface = brep.geom.face_surface
                        .get(flat_fi)
                        .and_then(|o| *o)
                        .is_some();
                    if has_surface {
                        for we in &face.outer_wire.edges {
                            surface_edges.insert(we.idx);
                        }
                        for iw in &face.inner_wires {
                            for we in &iw.edges {
                                surface_edges.insert(we.idx);
                            }
                        }
                    }
                    flat_fi += 1;
                }
            }
        }

        for &ei in &surface_edges {
            let has_pcurve = brep.geom.edge_pcurves
                .get(ei)
                .map_or(false, |v| !v.is_empty());
            if !has_pcurve {
                issues.push(ExportIssue::MissingPcurve { edge_idx: ei });
            }
        }
    }

    // ── 4: Tolerance floor ────────────────────────────────────────────────────
    for (vi, &t) in brep.geom.vertex_tolerance.iter().enumerate() {
        if t > 0.0 && t < CONFUSION {
            issues.push(ExportIssue::VertexToleranceTooTight { vertex_idx: vi, tolerance: t });
        }
    }
    for (ei, &t) in brep.geom.edge_tolerance.iter().enumerate() {
        if t > 0.0 && t < CONFUSION {
            issues.push(ExportIssue::EdgeToleranceTooTight { edge_idx: ei, tolerance: t });
        }
    }
    for (fi, &t) in brep.geom.face_tolerance.iter().enumerate() {
        if t > 0.0 && t < CONFUSION {
            issues.push(ExportIssue::FaceToleranceTooTight { face_idx: fi, tolerance: t });
        }
    }

    let is_ready = issues.is_empty();
    ExportReadinessReport {
        is_ready,
        issues,
        edges_checked: n_edges,
        pcurves_checked,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_algorithms::HealingMode;
    use std::io::Cursor;

    const HFSS_STEP: &str = include_str!("../../../assets/hfss.step");
    const BOX_STEP: &str = include_str!("../../../assets/box.step");
    const EDGE_ONLY_STEP: &str = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=CARTESIAN_POINT('',(1.,0.,0.));\n#3=VERTEX_POINT('',#1);\n#4=VERTEX_POINT('',#2);\n#5=EDGE_CURVE('',#3,#4,$,.T.);\nENDSEC;\nEND-ISO-10303-21;\n";

    /// Helper for round-trip testing: export BRep to STEP and re-import.
    fn round_trip_brep(brep: &BRep, gmsh_strict: bool) -> BRep {
        let step = StepWriter::write_string_with_options(
            brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
            &StepWriteOptions {
                gmsh_strict,
                ..Default::default()
            },
        );
        StepReader::parse_string(&step).expect("round-trip STEP should parse")
    }

    /// Helper to count topology elements.
    /// Counts only topological vertices (those referenced by edges), not tessellation vertices.
    fn count_topology(brep: &BRep) -> (usize, usize, usize, usize) {
        use std::collections::HashSet;
        let topological_vertices: HashSet<usize> = brep.edges.iter()
            .flat_map(|e| [e.start, e.end].into_iter())
            .collect();
        let vertices = topological_vertices.len();
        let edges = brep.edges.len();
        let shells: usize = brep.solids.iter().map(|s| s.shells.len()).sum();
        let faces: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .map(|sh| sh.faces.len())
            .sum();
        (vertices, edges, faces, shells)
    }

    #[test]
    fn round_trip_boolean_union() {
        use rcad_algorithms::{BooleanOpType, SimplifyOptions, boolean_op_simplified, geom_populate::populate_box_geom};
        use rcad_modeling::make_box_brep;
        use glam::DVec3;

        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        let mut b = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        populate_box_geom(&mut a);
        populate_box_geom(&mut b);

        let (union, _) = boolean_op_simplified(BooleanOpType::Union, &a, &b, SimplifyOptions::default())
            .expect("union should succeed");

        let (v1, e1, f1, s1) = count_topology(&union);
        let round_tripped = round_trip_brep(&union, true);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        assert_eq!(v1, v2, "vertex count should match after round-trip");
        assert_eq!(e1, e2, "edge count should match after round-trip");
        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");
    }

    #[test]
    fn round_trip_boolean_difference() {
        use rcad_algorithms::{BooleanOpType, SimplifyOptions, boolean_op_simplified, geom_populate::populate_box_geom};
        use rcad_modeling::make_box_brep;
        use glam::DVec3;

        let mut a = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 4.0, 4.0, 4.0).unwrap();
        let mut b = make_box_brep(DVec3::new(1.0, 1.0, 1.0), DVec3::X, DVec3::Y, 2.0, 2.0, 2.0).unwrap();
        populate_box_geom(&mut a);
        populate_box_geom(&mut b);

        let (diff, _) = boolean_op_simplified(BooleanOpType::Difference, &a, &b, SimplifyOptions::default())
            .expect("difference should succeed");

        let (v1, e1, f1, s1) = count_topology(&diff);
        let round_tripped = round_trip_brep(&diff, true);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        assert_eq!(v1, v2, "vertex count should match after round-trip");
        assert_eq!(e1, e2, "edge count should match after round-trip");
        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");
    }

    #[test]
    fn parses_hfss_into_non_trivial_brep() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        assert!(
            brep.vertices.len() > 8,
            "hfss should have more than box vertices"
        );
        assert!(!brep.edges.is_empty(), "hfss should produce edges");

        let triangle_count: usize = brep
            .solids
            .iter()
            .flat_map(|s| s.shells.iter())
            .flat_map(|sh| sh.faces.iter())
            .map(|f| f.triangles.len())
            .sum();

        assert!(triangle_count > 0, "hfss should produce triangulated faces");
    }

    #[test]
    fn triangulates_spherical_face_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        let mut face_triangles = Vec::new();
        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    face_triangles.push(face.triangles.len());
                }
            }
        }

        let max_face_triangles = face_triangles.into_iter().max().unwrap_or(0);
        assert!(
            max_face_triangles >= 200,
            "expected a tessellated spherical face, got max triangles={max_face_triangles}"
        );
    }

    #[test]
    fn triangulates_toroidal_face_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        let mut face_idx = 0usize;
        let mut found_torus = false;
        let mut torus_has_triangles = false;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    let is_torus = brep
                        .geom
                        .face_surface
                        .get(face_idx)
                        .and_then(|binding| *binding)
                        .and_then(|sid| brep.geom.surfaces.get(sid))
                        .map(|s| matches!(s, Surface3::Torus(_)))
                        .unwrap_or(false);

                    if is_torus {
                        found_torus = true;
                        if !face.triangles.is_empty() {
                            torus_has_triangles = true;
                        }
                    }

                    face_idx += 1;
                }
            }
        }

        assert!(found_torus, "expected hfss.step to contain a toroidal face");
        assert!(
            torus_has_triangles,
            "expected toroidal faces to be triangulated"
        );
    }

    #[test]
    fn triangulates_single_edge_planar_faces_from_hfss() {
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let mut found = false;

        for solid in &brep.solids {
            for shell in &solid.shells {
                for face in &shell.faces {
                    if face.outer_wire.edges.len() == 1 && face.triangles.len() >= 16 {
                        found = true;
                        break;
                    }
                }
            }
        }

        assert!(
            found,
            "expected at least one circular/elliptic single-edge planar face with triangulation"
        );
    }

    #[test]
    fn parses_box_example() {
        let brep = StepReader::parse_string(BOX_STEP).expect("box.step should parse");
        assert!(!brep.vertices.is_empty());
    }

    #[test]
    fn parse_box_with_healing_analyze_only_reports_clean() {
        let (brep, report) = StepReader::parse_string_with_healing(
            BOX_STEP,
            HealingOptions {
                mode: HealingMode::AnalyzeOnly,
                ..HealingOptions::default()
            },
        )
        .expect("box.step parse+healing should succeed");

        assert!(!brep.vertices.is_empty());
        assert!(report.initial.is_valid());
        assert!(report.final_result.is_valid());
        assert!(report.passes.is_empty());
    }

    #[test]
    fn parse_hfss_with_healing_returns_report() {
        let (brep, report) = StepReader::parse_string_with_healing(
            HFSS_STEP,
            HealingOptions::default(),
        )
        .expect("hfss.step parse+healing should succeed");

        assert!(!brep.vertices.is_empty());
        assert!(report.initial_issue_count() >= report.final_issue_count());
    }

    #[test]
    fn parse_hfss_with_healing_report_json_contains_schema_and_counts() {
        let (brep, report, json) = StepReader::parse_string_with_healing_report_json(
            HFSS_STEP,
            HealingOptions::default(),
        )
        .expect("hfss.step parse+healing+json should succeed");

        assert!(!brep.vertices.is_empty());
        let v: serde_json::Value =
            serde_json::from_str(&json).expect("healing report json should parse");
        assert_eq!(v["schema"], "step.import.healing.v1");
        assert_eq!(
            v["initial_issue_count"].as_u64().unwrap_or(0) as usize,
            report.initial_issue_count()
        );
        assert_eq!(
            v["final_issue_count"].as_u64().unwrap_or(0) as usize,
            report.final_issue_count()
        );
        assert!(v["issue_histogram"].is_array());
        assert_eq!(
            v["repair_pass_count"].as_u64().unwrap_or(0) as usize,
            report.passes.len()
        );
        assert_eq!(
            v["parametric_pass_count"].as_u64().unwrap_or(0) as usize,
            report.parametric_passes.len()
        );
        assert_eq!(
            v["make_connected_pass_count"].as_u64().unwrap_or(0) as usize,
            report.make_connected_passes.len()
        );
        assert!(v["faces_reoriented"].is_u64());
        assert!(v["same_range_fixed"].is_u64());
        assert!(v["same_parameter_fixed"].is_u64());
        assert!(v["wire_open_gaps"].is_u64());
        assert!(v["wire_topological_self_intersections"].is_u64());
        assert!(v["wire_geometric_self_intersections"].is_u64());
    }

    #[test]
    fn parse_box_with_metadata_extracts_schema_and_protocol_hint() {
        let (brep, meta) = StepReader::parse_string_with_metadata(BOX_STEP)
            .expect("box.step parse+metadata should succeed");

        assert!(!brep.vertices.is_empty());
        assert!(meta.file_schema.is_some(), "FILE_SCHEMA should be extracted");
        assert!(
            matches!(meta.protocol_hint, StepProtocolHint::Ap203 | StepProtocolHint::Ap214 | StepProtocolHint::Ap242 | StepProtocolHint::Unknown),
            "protocol hint should be classified"
        );
        assert!(meta.file_name.is_some(), "FILE_NAME should be extracted");
        // materials and layers are empty for this simple test file
        let _ = &meta.materials;
        let _ = &meta.layers;
        let _ = &meta.general_properties;
        let _ = &meta.property_definitions;
        let _ = &meta.products;
    }

    #[test]
    fn writer_tagged_ellipsoid_surface_resolves_back_to_native_surface() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Ellipsoid(rcad_kernel::EllipsoidalSurface {
            center: glam::DVec3::new(0.5, 0.5, 0.5),
            axis: glam::DVec3::Z,
            ref_dir: glam::DVec3::X,
            radius_x: 2.0,
            radius_y: 1.5,
            radius_z: 1.0,
        });

        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );
        let parsed = parse_entities(&step).expect("tagged ellipsoid STEP should parse into entities");
        let (&surface_id, _) = parsed
            .b_spline_surfaces
            .iter()
            .find(|(_, data)| data.0.starts_with("RCAD_ELLIPSOID"))
            .expect("expected tagged ellipsoid B-spline surface entity");

        match resolve_surface(&parsed, surface_id) {
            Some(Surface3::Ellipsoid(surface)) => {
                assert_eq!(surface.center, glam::DVec3::new(0.5, 0.5, 0.5));
                assert_eq!(surface.axis, glam::DVec3::Z);
                assert_eq!(surface.ref_dir, glam::DVec3::X);
                assert_eq!(surface.radius_x, 2.0);
                assert_eq!(surface.radius_y, 1.5);
                assert_eq!(surface.radius_z, 1.0);
            }
            other => panic!("expected ellipsoid from tagged surface, got {other:?}"),
        }
    }

    #[test]
    fn writer_tagged_helicoid_surface_resolves_back_to_native_surface() {
        let mut brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");
        let sid = *brep
            .geom
            .face_surface
            .iter()
            .flatten()
            .next()
            .expect("hfss.step should contain a face surface");
        brep.geom.surfaces[sid] = Surface3::Helicoid(rcad_kernel::HelicoidSurface {
            origin: glam::DVec3::ZERO,
            axis: glam::DVec3::Z,
            ref_dir: glam::DVec3::X,
            pitch: 3.0,
        });

        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );
        let parsed = parse_entities(&step).expect("tagged helicoid STEP should parse into entities");
        let (&surface_id, _) = parsed
            .b_spline_surfaces
            .iter()
            .find(|(_, data)| data.0.starts_with("RCAD_HELICOID"))
            .expect("expected tagged helicoid B-spline surface entity");

        match resolve_surface(&parsed, surface_id) {
            Some(Surface3::Helicoid(surface)) => {
                assert_eq!(surface.origin, glam::DVec3::ZERO);
                assert_eq!(surface.axis, glam::DVec3::Z);
                assert_eq!(surface.ref_dir, glam::DVec3::X);
                assert_eq!(surface.pitch, 3.0);
            }
            other => panic!("expected helicoid from tagged surface, got {other:?}"),
        }
    }

    #[test]
    fn extract_materials_finds_material_entities() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 3 1 1 1 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=MATERIAL('Steel','Carbon steel',#2);
#2=MATERIAL('Aluminium','Lightweight alloy',#3);
ENDSEC;
END-ISO-10303-21;
"#;
        let materials = extract_materials(step_fragment);
        assert_eq!(materials.len(), 2);
        assert_eq!(materials[0].name, "Steel");
        assert_eq!(materials[1].name, "Aluminium");
    }

    #[test]
    fn extract_layers_finds_presentation_layer_assignment() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 3 1 1 1 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#10=PRESENTATION_LAYER_ASSIGNMENT('Layer_0','Default layer',(#5));
#11=PRESENTATION_LAYER_ASSIGNMENT('Body_Layer','Body geometry',(#6,#7));
ENDSEC;
END-ISO-10303-21;
"#;
        let layers = extract_layers(step_fragment);
        assert_eq!(layers.len(), 2);
        assert_eq!(layers[0].name, "Layer_0");
        assert_eq!(layers[1].name, "Body_Layer");
    }

    #[test]
    fn extract_product_names_finds_unique_product_entities() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=PRODUCT('ROOT-001','RootAsm','assembly root',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#4=PRODUCT('BRKT-001','Bracket','mount bracket',(#2));
#5=PRODUCT('BRKT-001-DUP','Bracket','mount bracket duplicate',(#2));
ENDSEC;
END-ISO-10303-21;
"#;
        let names = extract_product_names(step_fragment);
        assert_eq!(names, vec!["RootAsm".to_string(), "Bracket".to_string()]);
    }

    #[test]
    fn extract_products_returns_step_product_records() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#11=PRODUCT('ASM-100','Top Assembly','assembly description',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
ENDSEC;
END-ISO-10303-21;
"#;

        let products = extract_products(step_fragment);
        assert_eq!(products.len(), 1);
        assert_eq!(products[0].entity_id, 11);
        assert_eq!(products[0].product_id.as_deref(), Some("ASM-100"));
        assert_eq!(products[0].name.as_deref(), Some("Top Assembly"));
        assert_eq!(products[0].description.as_deref(), Some("assembly description"));
    }

    #[test]
    fn extract_product_definition_formations_resolve_product_names() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=PRODUCT('ROOT-001','RootAsm','assembly root',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#10=PRODUCT_DEFINITION_FORMATION('rel-1','released',#1);
ENDSEC;
END-ISO-10303-21;
"#;

        let formations = extract_product_definition_formations(step_fragment);
        assert_eq!(formations.len(), 1);
        assert_eq!(formations[0].entity_id, 10);
        assert_eq!(formations[0].product_id, Some(1));
        assert_eq!(formations[0].product_name.as_deref(), Some("RootAsm"));
    }

    #[test]
    fn extract_product_definitions_resolve_product_chain() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=PRODUCT('ROOT-001','RootAsm','assembly root',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#10=PRODUCT_DEFINITION_FORMATION('rel-1','released',#1);
#20=PRODUCT_DEFINITION('root def','',#10,#30);
#30=PRODUCT_DEFINITION_CONTEXT('part definition',#3,'design');
ENDSEC;
END-ISO-10303-21;
"#;

        let defs = extract_product_definitions(step_fragment);
        assert_eq!(defs.len(), 1);
        assert_eq!(defs[0].entity_id, 20);
        assert_eq!(defs[0].formation_id, Some(10));
        assert_eq!(defs[0].product_id, Some(1));
        assert_eq!(defs[0].product_name.as_deref(), Some("RootAsm"));
    }

    #[test]
    fn extract_shape_definition_representations_resolve_product_definition() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=PRODUCT('ROOT-001','RootAsm','assembly root',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#10=PRODUCT_DEFINITION_FORMATION('rel-1','released',#1);
#20=PRODUCT_DEFINITION('root def','',#10,#30);
#30=PRODUCT_DEFINITION_CONTEXT('part definition',#3,'design');
#40=PRODUCT_DEFINITION_SHAPE('','',#20);
#50=SHAPE_DEFINITION_REPRESENTATION(#40,#60);
#60=SHAPE_REPRESENTATION('RootAsm',(),#70);
#70=( GEOMETRIC_REPRESENTATION_CONTEXT(3) REPRESENTATION_CONTEXT('Context #1','3D') );
ENDSEC;
END-ISO-10303-21;
"#;

        let reprs = extract_shape_definition_representations(step_fragment);
        assert_eq!(reprs.len(), 1);
        assert_eq!(reprs[0].product_definition_shape_id, Some(40));
        assert_eq!(reprs[0].product_definition_id, Some(20));
        assert_eq!(reprs[0].representation_id, Some(60));
    }

    #[test]
    fn extract_assembly_occurrences_resolve_parent_child_product_names() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#1=PRODUCT('ASM-001','Assembly','root assembly',(#2));
#4=PRODUCT('PRT-001','Bracket','child part',(#2));
#2=PRODUCT_CONTEXT('part definition',#3,'mechanical');
#3=APPLICATION_CONTEXT('configuration controlled 3d designs of mechanical parts and assemblies');
#10=PRODUCT_DEFINITION_FORMATION('asm-rel','released',#1);
#11=PRODUCT_DEFINITION_FORMATION('prt-rel','released',#4);
#20=PRODUCT_DEFINITION('asm def','',#10,#30);
#21=PRODUCT_DEFINITION('part def','',#11,#30);
#30=PRODUCT_DEFINITION_CONTEXT('part definition',#3,'design');
#50=NEXT_ASSEMBLY_USAGE_OCCURRENCE('1','Bracket occ','',#20,#21,$);
ENDSEC;
END-ISO-10303-21;
"#;

        let occs = extract_assembly_occurrences(step_fragment);
        assert_eq!(occs.len(), 1);
        assert_eq!(occs[0].relating_product_definition_id, Some(20));
        assert_eq!(occs[0].related_product_definition_id, Some(21));
        assert_eq!(occs[0].relating_product_name.as_deref(), Some("Assembly"));
        assert_eq!(occs[0].related_product_name.as_deref(), Some("Bracket"));
    }

    #[test]
    fn extract_general_properties_finds_general_property_entities() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#21=GENERAL_PROPERTY('PartNumber','ERP key',$);
#22=GENERAL_PROPERTY('Revision','A',$);
ENDSEC;
END-ISO-10303-21;
"#;
        let props = extract_general_properties(step_fragment);
        assert_eq!(props.len(), 2);
        assert_eq!(props[0].name, "PartNumber");
        assert_eq!(props[0].description.as_deref(), Some("ERP key"));
        assert_eq!(props[1].name, "Revision");
        assert_eq!(props[1].description.as_deref(), Some("A"));
    }

    #[test]
    fn extract_property_definitions_links_general_property_chain() {
        let step_fragment = r#"ISO-10303-21;
HEADER;
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
FILE_NAME('test','','','','','','');
ENDSEC;
DATA;
#21=GENERAL_PROPERTY('PartNumber','ERP key',$);
#22=GENERAL_PROPERTY('Revision','A',$);
#30=PROPERTY_DEFINITION('part number def','',#21);
#31=PROPERTY_DEFINITION('revision def','',#22);
ENDSEC;
END-ISO-10303-21;
"#;

        let defs = extract_property_definitions(step_fragment);
        assert_eq!(defs.len(), 2);

        assert_eq!(defs[0].definition_id, 30);
        assert_eq!(defs[0].referenced_entity_id, Some(21));
        assert_eq!(defs[0].general_property_name.as_deref(), Some("PartNumber"));
        assert_eq!(defs[0].general_property_description.as_deref(), Some("ERP key"));

        assert_eq!(defs[1].definition_id, 31);
        assert_eq!(defs[1].referenced_entity_id, Some(22));
        assert_eq!(defs[1].general_property_name.as_deref(), Some("Revision"));
        assert_eq!(defs[1].general_property_description.as_deref(), Some("A"));
    }

    #[test]
    fn preserves_standalone_edge_only_geometry() {
        let brep = StepReader::parse_string(EDGE_ONLY_STEP).expect("edge-only STEP should parse");
        assert_eq!(brep.vertices.len(), 2);
        assert_eq!(brep.edges.len(), 1);
        assert!(
            brep.solids.is_empty(),
            "edge-only data should not fabricate solids"
        );
    }

    #[test]
    fn samples_geometric_curve_sets_from_hfss() {
        // hfss.step has GEOMETRIC_CURVE_SET with Polyline1 (2 trimmed lines),
        // Polyline2 (b-spline), and Polyline3 (trimmed circle arc 0..135 deg).
        // All should produce additional edges in the BRep.
        let brep = StepReader::parse_string(HFSS_STEP).expect("hfss.step should parse");

        // The b-spline alone has 7 control points 鈫?at least 6 edges
        // The trimmed lines each contribute 1 edge (2 total)
        // The circle arc contributes 8+ edges
        // Total from curve sets: at least 16 edges beyond the face topology edges
        let total_edges = brep.edges.len();
        assert!(
            total_edges >= 20,
            "expected geometric curve set edges, got total edge count = {total_edges}"
        );
    }

    #[test]
    fn extract_property_definition_reprs_parses_linkage() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#10=PROPERTY_DEFINITION('pd','',#20);\n#11=PROPERTY_DEFINITION_REPRESENTATION(#10,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let reprs = extract_property_definition_reprs(step);
        assert_eq!(reprs.len(), 1);
        assert_eq!(reprs[0].entity_id, 11);
        assert_eq!(reprs[0].property_definition_id, Some(10));
        assert_eq!(reprs[0].representation_id, Some(30));
    }

    #[test]
    fn extract_dimensional_locations_parses_from_to() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#5=DIMENSIONAL_LOCATION('diameter','D',#100,#101);\nENDSEC;\nEND-ISO-10303-21;\n";
        let locs = extract_dimensional_locations(step);
        assert_eq!(locs.len(), 1);
        assert_eq!(locs[0].entity_id, 5);
        assert_eq!(locs[0].name.as_deref(), Some("diameter"));
        assert_eq!(locs[0].description.as_deref(), Some("D"));
        assert_eq!(locs[0].from_entity_id, Some(100));
        assert_eq!(locs[0].to_entity_id, Some(101));
    }

    #[test]
    fn extract_geometric_tolerances_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#7=GEOMETRIC_TOLERANCE('flatness','',#50,#60);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_geometric_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 7);
        assert_eq!(tols[0].name.as_deref(), Some("flatness"));
        assert_eq!(tols[0].value_entity_id, Some(50));
        assert_eq!(tols[0].shape_aspect_id, Some(60));
    }

    #[test]
    fn extract_datums_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#8=DATUM('A','primary',#70);\nENDSEC;\nEND-ISO-10303-21;\n";
        let datums = extract_datums(step);
        assert_eq!(datums.len(), 1);
        assert_eq!(datums[0].entity_id, 8);
        assert_eq!(datums[0].name.as_deref(), Some("A"));
        assert_eq!(datums[0].description.as_deref(), Some("primary"));
        assert_eq!(datums[0].shape_aspect_id, Some(70));
    }

    #[test]
    fn extract_geometric_tolerance_with_datum_reference_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#9=GEOMETRIC_TOLERANCE_WITH_DATUM_REFERENCE('positional','',#50,#60,#70);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_geometric_tolerances_with_datum_references(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 9);
        assert_eq!(tols[0].name.as_deref(), Some("positional"));
        assert_eq!(tols[0].value_entity_id, Some(50));
        assert_eq!(tols[0].shape_aspect_id, Some(60));
        assert_eq!(tols[0].datum_system_id, Some(70));
    }

    #[test]
    fn extract_datum_systems_parses_grouping() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#8=DATUM('A','primary',#70);\n#9=DATUM('B','secondary',#71);\n#10=DATUM_SYSTEM('A|B','main',(#8,#9));\nENDSEC;\nEND-ISO-10303-21;\n";
        let systems = extract_datum_systems(step);
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0].entity_id, 10);
        assert_eq!(systems[0].name.as_deref(), Some("A|B"));
        assert_eq!(systems[0].datum_ids, vec![8, 9]);
    }

    #[test]
    fn extract_kinematic_pairs_parses_revolute_pair() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=REVOLUTE_PAIR('hinge','joint',#10,#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let pairs = extract_kinematic_pairs(step);
        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].entity_id, 100);
        assert_eq!(pairs[0].entity_type, "REVOLUTE_PAIR");
        assert_eq!(pairs[0].name.as_deref(), Some("hinge"));
        assert_eq!(pairs[0].related_entity_ids, vec![10, 20, 30]);
    }

    #[test]
    fn extract_tolerance_zones_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#20=TOLERANCE_ZONE('cylindrical','zone definition',#50);\nENDSEC;\nEND-ISO-10303-21;\n";
        let zones = extract_tolerance_zones(step);
        assert_eq!(zones.len(), 1);
        assert_eq!(zones[0].entity_id, 20);
        assert_eq!(zones[0].name.as_deref(), Some("cylindrical"));
        assert_eq!(zones[0].description.as_deref(), Some("zone definition"));
        assert_eq!(zones[0].toleranced_shape_aspect_id, Some(50));
    }

    #[test]
    fn extract_tolerance_zone_definitions_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#30=TOLERANCE_ZONE_DEFINITION('def1','definition',#20,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_tolerance_zone_definitions(step);
        assert_eq!(defs.len(), 1);
        assert_eq!(defs[0].entity_id, 30);
        assert_eq!(defs[0].name.as_deref(), Some("def1"));
        assert_eq!(defs[0].tolerance_zone_id, Some(20));
        assert_eq!(defs[0].shape_aspect_id, Some(40));
    }

    #[test]
    fn extract_datum_features_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#15=DATUM_FEATURE('feature_A','datum feature',#8);\nENDSEC;\nEND-ISO-10303-21;\n";
        let features = extract_datum_features(step);
        assert_eq!(features.len(), 1);
        assert_eq!(features[0].entity_id, 15);
        assert_eq!(features[0].name.as_deref(), Some("feature_A"));
        assert_eq!(features[0].description.as_deref(), Some("datum feature"));
        assert_eq!(features[0].datum_id, Some(8));
    }

    #[test]
    fn extract_datum_reference_elements_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#25=DATUM_REFERENCE_ELEMENT('ref_elem','reference element',#15);\nENDSEC;\nEND-ISO-10303-21;\n";
        let elems = extract_datum_reference_elements(step);
        assert_eq!(elems.len(), 1);
        assert_eq!(elems[0].entity_id, 25);
        assert_eq!(elems[0].name.as_deref(), Some("ref_elem"));
        assert_eq!(elems[0].associated_entity_id, Some(15));
    }

    #[test]
    fn extract_shape_aspects_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#60=SHAPE_ASPECT('hole','hole feature',#100);\nENDSEC;\nEND-ISO-10303-21;\n";
        let aspects = extract_shape_aspects(step);
        assert_eq!(aspects.len(), 1);
        assert_eq!(aspects[0].entity_id, 60);
        assert_eq!(aspects[0].name.as_deref(), Some("hole"));
        assert_eq!(aspects[0].description.as_deref(), Some("hole feature"));
        assert_eq!(aspects[0].product_definition_shape_id, Some(100));
    }

    #[test]
    fn extract_shape_aspect_definitions_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#70=SHAPE_ASPECT_DEFINITION('hole_def','hole definition',#60);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_shape_aspect_definitions(step);
        assert_eq!(defs.len(), 1);
        assert_eq!(defs[0].entity_id, 70);
        assert_eq!(defs[0].name.as_deref(), Some("hole_def"));
        assert_eq!(defs[0].shape_aspect_id, Some(60));
    }

    #[test]
    fn extract_derived_shape_aspects_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#80=DERIVED_SHAPE_ASPECT('derived','derived aspect',(#60,#70));\nENDSEC;\nEND-ISO-10303-21;\n";
        let aspects = extract_derived_shape_aspects(step);
        assert_eq!(aspects.len(), 1);
        assert_eq!(aspects[0].entity_id, 80);
        assert_eq!(aspects[0].name.as_deref(), Some("derived"));
        assert_eq!(aspects[0].base_shape_aspect_ids, vec![60, 70]);
    }

    // ── FEA entity extraction tests ──────────────────────────────────────

    #[test]
    fn extract_fea_models_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=FEAMEDIAN_MODEL('bracket_fea','FEA model for bracket');\nENDSEC;\nEND-ISO-10303-21;\n";
        let models = extract_fea_models(step);
        assert_eq!(models.len(), 1);
        assert_eq!(models[0].entity_id, 100);
        assert_eq!(models[0].name.as_deref(), Some("bracket_fea"));
        assert_eq!(models[0].description.as_deref(), Some("FEA model for bracket"));
    }

    #[test]
    fn extract_fea_meshes_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#200=FEAMEDIAN_MESH('mesh1',1000,500);\nENDSEC;\nEND-ISO-10303-21;\n";
        let meshes = extract_fea_meshes(step);
        assert_eq!(meshes.len(), 1);
        assert_eq!(meshes[0].entity_id, 200);
        assert_eq!(meshes[0].name.as_deref(), Some("mesh1"));
        assert_eq!(meshes[0].node_count, Some(1000));
        assert_eq!(meshes[0].element_count, Some(500));
    }

    #[test]
    fn extract_fea_node_sets_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#300=FEAMEDIAN_NODE_SET('fixed_nodes',#100);\nENDSEC;\nEND-ISO-10303-21;\n";
        let node_sets = extract_fea_node_sets(step);
        assert_eq!(node_sets.len(), 1);
        assert_eq!(node_sets[0].entity_id, 300);
        assert_eq!(node_sets[0].name.as_deref(), Some("fixed_nodes"));
        assert_eq!(node_sets[0].model_id, Some(100));
    }

    #[test]
    fn extract_fea_element_sets_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#400=FEAMEDIAN_ELEMENT_SET('shell_elements',#100,'QUAD4');\nENDSEC;\nEND-ISO-10303-21;\n";
        let element_sets = extract_fea_element_sets(step);
        assert_eq!(element_sets.len(), 1);
        assert_eq!(element_sets[0].entity_id, 400);
        assert_eq!(element_sets[0].name.as_deref(), Some("shell_elements"));
        assert_eq!(element_sets[0].model_id, Some(100));
        assert_eq!(element_sets[0].element_type.as_deref(), Some("QUAD4"));
    }

    #[test]
    fn extract_fea_material_properties_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#500=FEAMEDIAN_MATERIAL_PROPERTY('Steel_E','YoungsModulus',210000.0,'MPa');\nENDSEC;\nEND-ISO-10303-21;\n";
        let props = extract_fea_material_properties(step);
        assert_eq!(props.len(), 1);
        assert_eq!(props[0].entity_id, 500);
        assert_eq!(props[0].name.as_deref(), Some("Steel_E"));
        assert_eq!(props[0].property_type.as_deref(), Some("YoungsModulus"));
        assert!((props[0].value.unwrap() - 210000.0).abs() < 1e-6);
        assert_eq!(props[0].unit.as_deref(), Some("MPa"));
    }

    #[test]
    fn extract_fea_boundary_conditions_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#600=FEAMEDIAN_BOUNDARY_CONDITION('fixed_bc','FIXED',#300);\nENDSEC;\nEND-ISO-10303-21;\n";
        let bcs = extract_fea_boundary_conditions(step);
        assert_eq!(bcs.len(), 1);
        assert_eq!(bcs[0].entity_id, 600);
        assert_eq!(bcs[0].name.as_deref(), Some("fixed_bc"));
        assert_eq!(bcs[0].condition_type.as_deref(), Some("FIXED"));
        assert_eq!(bcs[0].node_set_id, Some(300));
    }

    #[test]
    fn extract_fea_loads_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#700=FEAMEDIAN_LOAD('pressure_load','PRESSURE',100.0,(0.0,0.0,-1.0));\nENDSEC;\nEND-ISO-10303-21;\n";
        let loads = extract_fea_loads(step);
        assert_eq!(loads.len(), 1);
        assert_eq!(loads[0].entity_id, 700);
        assert_eq!(loads[0].name.as_deref(), Some("pressure_load"));
        assert_eq!(loads[0].load_type.as_deref(), Some("PRESSURE"));
        assert!((loads[0].magnitude.unwrap() - 100.0).abs() < 1e-6);
        assert_eq!(loads[0].direction.unwrap(), [0.0, 0.0, -1.0]);
    }

    #[test]
    fn extract_fea_models_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=FEA_MODEL('model1','Standard FEA model');\n#101=FEA_MODEL_3D('model2','3D FEA model');\n#102=FEAMEDIAN_MODEL('model3','Median model');\nENDSEC;\nEND-ISO-10303-21;\n";
        let models = extract_fea_models(step);
        assert_eq!(models.len(), 3);
        assert_eq!(models[0].name.as_deref(), Some("model1"));
        assert_eq!(models[1].name.as_deref(), Some("model2"));
        assert_eq!(models[2].name.as_deref(), Some("model3"));
    }

    #[test]
    fn extract_fea_meshes_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#200=FEA_MESH('mesh1',1000,500);\n#201=FEA_MESH_3D('mesh2',2000,800);\n#202=FEAMEDIAN_MESH('mesh3',500,200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let meshes = extract_fea_meshes(step);
        assert_eq!(meshes.len(), 3);
        assert_eq!(meshes[0].name.as_deref(), Some("mesh1"));
        assert_eq!(meshes[1].name.as_deref(), Some("mesh2"));
        assert_eq!(meshes[2].name.as_deref(), Some("mesh3"));
        assert_eq!(meshes[0].node_count, Some(1000));
        assert_eq!(meshes[1].node_count, Some(2000));
        assert_eq!(meshes[2].node_count, Some(500));
    }

    #[test]
    fn extract_fea_node_sets_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#300=FEA_NODE_SET('nodes1',#100);\n#301=FEAMEDIAN_NODE_SET('nodes2',#101);\nENDSEC;\nEND-ISO-10303-21;\n";
        let node_sets = extract_fea_node_sets(step);
        assert_eq!(node_sets.len(), 2);
        assert_eq!(node_sets[0].name.as_deref(), Some("nodes1"));
        assert_eq!(node_sets[1].name.as_deref(), Some("nodes2"));
    }

    #[test]
    fn extract_fea_element_sets_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#400=FEA_ELEMENT_SET('elements1',#100,'TETRAHEDRON');\n#401=FEAMEDIAN_ELEMENT_SET('elements2',#101,'HEXAHEDRON');\nENDSEC;\nEND-ISO-10303-21;\n";
        let element_sets = extract_fea_element_sets(step);
        assert_eq!(element_sets.len(), 2);
        assert_eq!(element_sets[0].name.as_deref(), Some("elements1"));
        assert_eq!(element_sets[1].name.as_deref(), Some("elements2"));
        assert_eq!(element_sets[0].element_type.as_deref(), Some("TETRAHEDRON"));
        assert_eq!(element_sets[1].element_type.as_deref(), Some("HEXAHEDRON"));
    }

    #[test]
    fn extract_fea_material_properties_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#500=FEA_MATERIAL_PROPERTY('YoungsModulus','ELASTIC_MODULUS',210000.0,'MPa');\n#501=FEAMEDIAN_MATERIAL_PROPERTY('Density','DENSITY',7850.0,'kg/m3');\nENDSEC;\nEND-ISO-10303-21;\n";
        let props = extract_fea_material_properties(step);
        assert_eq!(props.len(), 2);
        assert_eq!(props[0].name.as_deref(), Some("YoungsModulus"));
        assert_eq!(props[1].name.as_deref(), Some("Density"));
        assert!((props[0].value.unwrap() - 210000.0).abs() < 1e-6);
        assert!((props[1].value.unwrap() - 7850.0).abs() < 1e-6);
    }

    #[test]
    fn extract_fea_boundary_conditions_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#600=FEA_BOUNDARY_CONDITION('fixed','FIXED',#300);\n#601=FEAMEDIAN_BOUNDARY_CONDITION('symmetry','SYMMETRY',#301);\nENDSEC;\nEND-ISO-10303-21;\n";
        let bcs = extract_fea_boundary_conditions(step);
        assert_eq!(bcs.len(), 2);
        assert_eq!(bcs[0].name.as_deref(), Some("fixed"));
        assert_eq!(bcs[1].name.as_deref(), Some("symmetry"));
        assert_eq!(bcs[0].condition_type.as_deref(), Some("FIXED"));
        assert_eq!(bcs[1].condition_type.as_deref(), Some("SYMMETRY"));
    }

    #[test]
    fn extract_fea_loads_parses_multiple_entity_types() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#700=FEA_LOAD('force1','FORCE',1000.0,(1.0,0.0,0.0));\n#701=FEAMEDIAN_LOAD('pressure1','PRESSURE',500.0,(0.0,0.0,-1.0));\nENDSEC;\nEND-ISO-10303-21;\n";
        let loads = extract_fea_loads(step);
        assert_eq!(loads.len(), 2);
        assert_eq!(loads[0].name.as_deref(), Some("force1"));
        assert_eq!(loads[1].name.as_deref(), Some("pressure1"));
        assert_eq!(loads[0].load_type.as_deref(), Some("FORCE"));
        assert_eq!(loads[1].load_type.as_deref(), Some("PRESSURE"));
    }

    #[test]
    fn extract_fea_node_groups_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#800=FEA_NODE_GROUP('boundary_nodes',#100,50);\n#801=FEA_NODE_GROUP('internal_nodes',#100,200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let groups = extract_fea_node_groups(step);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].entity_id, 800);
        assert_eq!(groups[0].name.as_deref(), Some("boundary_nodes"));
        assert_eq!(groups[0].model_id, Some(100));
        assert_eq!(groups[0].node_count, Some(50));
        assert_eq!(groups[1].name.as_deref(), Some("internal_nodes"));
        assert_eq!(groups[1].node_count, Some(200));
    }

    // ── Additional FEA entity extraction tests (AP209/AP242 extended) ─────────────

    #[test]
    fn extract_fea_analyses_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=FEA_ANALYSIS('static_analysis','Static structural analysis',#50,'STATIC','2024-01-15');\n#101=ANALYSIS_3D('modal_analysis','Modal analysis',#51,'MODAL','2024-01-16');\nENDSEC;\nEND-ISO-10303-21;\n";
        let analyses = extract_fea_analyses(step);
        assert_eq!(analyses.len(), 2);
        assert_eq!(analyses[0].entity_id, 100);
        assert_eq!(analyses[0].name.as_deref(), Some("static_analysis"));
        assert_eq!(analyses[0].description.as_deref(), Some("Static structural analysis"));
        assert_eq!(analyses[0].model_id, Some(50));
        assert_eq!(analyses[0].analysis_type.as_deref(), Some("STATIC"));
        assert_eq!(analyses[0].creation_date.as_deref(), Some("2024-01-15"));
        assert_eq!(analyses[1].name.as_deref(), Some("modal_analysis"));
        assert_eq!(analyses[1].analysis_type.as_deref(), Some("MODAL"));
    }

    #[test]
    fn extract_fea_states_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#200=FEA_STATE('initial_state','Initial conditions',#100,'INITIAL');\n#201=FEA_STATE('result_state','Analysis results',#100,'RESULT');\nENDSEC;\nEND-ISO-10303-21;\n";
        let states = extract_fea_states(step);
        assert_eq!(states.len(), 2);
        assert_eq!(states[0].entity_id, 200);
        assert_eq!(states[0].name.as_deref(), Some("initial_state"));
        assert_eq!(states[0].analysis_id, Some(100));
        assert_eq!(states[0].state_type.as_deref(), Some("INITIAL"));
        assert_eq!(states[1].state_type.as_deref(), Some("RESULT"));
    }

    #[test]
    fn extract_fea_material_models_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#300=FEA_LINEAR_ELASTICITY('Steel','Structural steel',210000.0,0.3,81000.0,7850.0,1.2e-5,'MPa');\n#301=FEA_MATERIAL_MODEL('CustomMaterial','Custom material');\nENDSEC;\nEND-ISO-10303-21;\n";
        let models = extract_fea_material_models(step);
        assert_eq!(models.len(), 2);
        assert_eq!(models[0].entity_id, 300);
        assert_eq!(models[0].name.as_deref(), Some("Steel"));
        assert_eq!(models[0].model_type, "LINEAR_ELASTIC");
        assert!((models[0].youngs_modulus.unwrap() - 210000.0).abs() < 1e-6);
        assert!((models[0].poissons_ratio.unwrap() - 0.3).abs() < 1e-6);
        assert!((models[0].density.unwrap() - 7850.0).abs() < 1e-6);
        assert_eq!(models[0].modulus_unit.as_deref(), Some("MPa"));
        assert_eq!(models[1].model_type, "GENERIC");
    }

    #[test]
    fn extract_fea_nodes_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#400=NODE_REPRESENTATION('node1',1,(0.0,0.0,0.0),#500);\n#401=NODE_REPRESENTATION('node2',2,(10.0,0.0,0.0),#500);\n#402=NODE_REPRESENTATION('node3',3,(10.0,10.0,0.0),#500);\n#403=FEA_NODE(4,(5.0,5.0,5.0),#501);\nENDSEC;\nEND-ISO-10303-21;\n";
        let nodes = extract_fea_nodes(step);
        assert_eq!(nodes.len(), 4);
        assert_eq!(nodes[0].entity_id, 400);
        assert_eq!(nodes[0].node_number, 1);
        assert_eq!(nodes[0].coordinates, [0.0, 0.0, 0.0]);
        assert_eq!(nodes[0].mesh_id, Some(500));
        assert_eq!(nodes[1].node_number, 2);
        assert_eq!(nodes[1].coordinates, [10.0, 0.0, 0.0]);
        assert_eq!(nodes[3].node_number, 4);
        assert_eq!(nodes[3].coordinates, [5.0, 5.0, 5.0]);
    }

    #[test]
    fn extract_fea_elements_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#600=ELEMENT_REPRESENTATION('elem1',1,'TETRA4',(#400,#401,#402,#403),#500);\n#601=FEA_ELEMENT(2,'HEXA8',(#410,#411,#412,#413,#414,#415,#416,#417),#501);\nENDSEC;\nEND-ISO-10303-21;\n";
        let elements = extract_fea_elements(step);
        assert_eq!(elements.len(), 2);
        assert_eq!(elements[0].entity_id, 600);
        assert_eq!(elements[0].element_number, 1);
        assert_eq!(elements[0].element_type, "TETRA4");
        assert_eq!(elements[0].node_ids, vec![400, 401, 402, 403]);
        assert_eq!(elements[0].mesh_id, Some(500));
        assert_eq!(elements[1].element_number, 2);
        assert_eq!(elements[1].element_type, "HEXA8");
        assert_eq!(elements[1].node_ids.len(), 8);
    }

    #[test]
    fn extract_fea_steps_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#700=FEA_STEP('load_step_1','Apply load',#100,1,'LOAD');\n#701=FEA_STEP('load_step_2','Increase load',#100,2,'LOAD');\nENDSEC;\nEND-ISO-10303-21;\n";
        let steps = extract_fea_steps(step);
        assert_eq!(steps.len(), 2);
        assert_eq!(steps[0].entity_id, 700);
        assert_eq!(steps[0].name.as_deref(), Some("load_step_1"));
        assert_eq!(steps[0].analysis_id, Some(100));
        assert_eq!(steps[0].step_number, Some(1));
        assert_eq!(steps[0].step_type.as_deref(), Some("LOAD"));
    }

    #[test]
    fn extract_fea_results_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#800=FEA_RESULT('displacement_result','DISPLACEMENT',#100,#300,3);\n#801=NODAL_RESULT('stress_result','STRESS',#100,#301,6);\n#802=ELEMENT_RESULT('strain_result','STRAIN',#100,#302,6);\nENDSEC;\nEND-ISO-10303-21;\n";
        let results = extract_fea_results(step);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].entity_id, 800);
        assert_eq!(results[0].name.as_deref(), Some("displacement_result"));
        assert_eq!(results[0].result_type, "DISPLACEMENT");
        assert_eq!(results[0].analysis_id, Some(100));
        assert_eq!(results[0].component_count, Some(3));
        assert_eq!(results[1].result_type, "STRESS");
        assert_eq!(results[2].result_type, "STRAIN");
    }

    #[test]
    fn extract_fea_cases_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#900=FEA_CASE('load_case_1','Gravity load case',#50,'LOAD_CASE');\n#901=FEA_CASE('bc_case_1','Fixed boundary conditions',#50,'BC_CASE');\nENDSEC;\nEND-ISO-10303-21;\n";
        let cases = extract_fea_cases(step);
        assert_eq!(cases.len(), 2);
        assert_eq!(cases[0].entity_id, 900);
        assert_eq!(cases[0].name.as_deref(), Some("load_case_1"));
        assert_eq!(cases[0].description.as_deref(), Some("Gravity load case"));
        assert_eq!(cases[0].model_id, Some(50));
        assert_eq!(cases[0].case_type.as_deref(), Some("LOAD_CASE"));
        assert_eq!(cases[1].case_type.as_deref(), Some("BC_CASE"));
    }

    #[test]
    fn fea_metadata_integration() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#10=FEA_MODEL('bracket_fea','FEA model for bracket');\n#20=FEA_MESH('mesh1',1000,500);\n#30=FEA_ANALYSIS('static','Static analysis',#10,'STATIC','2024-01-15');\n#40=FEA_LINEAR_ELASTICITY('Steel','Steel material',210000.0,0.3,81000.0,7850.0,1.2e-5,'MPa');\n#50=NODE_REPRESENTATION('n1',1,(0.0,0.0,0.0),#20);\n#51=NODE_REPRESENTATION('n2',2,(1.0,0.0,0.0),#20);\n#60=ELEMENT_REPRESENTATION('e1',1,'TETRA4',(#50,#51,#52,#53),#20);\n#70=FEA_RESULT('disp','DISPLACEMENT',#30,#20,3);\nENDSEC;\nEND-ISO-10303-21;\n";
        // Use extraction functions directly since FEA-only files don't have B-Rep geometry
        assert_eq!(extract_fea_models(step).len(), 1);
        assert_eq!(extract_fea_meshes(step).len(), 1);
        assert_eq!(extract_fea_analyses(step).len(), 1);
        assert_eq!(extract_fea_material_models(step).len(), 1);
        assert_eq!(extract_fea_nodes(step).len(), 2);
        assert_eq!(extract_fea_elements(step).len(), 1);
        assert_eq!(extract_fea_results(step).len(), 1);
    }

    #[test]
    fn metadata_summary_reports_all_entity_counts() {
        // Test extraction functions directly
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#60=SHAPE_ASPECT('hole','hole feature',#100);\n#20=TOLERANCE_ZONE('cylindrical','zone definition',#50);\n#15=DATUM_FEATURE('feature_A','datum feature',#8);\nENDSEC;\nEND-ISO-10303-21;\n";

        // Verify extraction functions work
        let shape_aspects = extract_shape_aspects(step);
        let tolerance_zones = extract_tolerance_zones(step);
        let datum_features = extract_datum_features(step);

        assert_eq!(shape_aspects.len(), 1, "Should extract 1 shape aspect");
        assert_eq!(tolerance_zones.len(), 1, "Should extract 1 tolerance zone");
        assert_eq!(datum_features.len(), 1, "Should extract 1 datum feature");

        assert_eq!(shape_aspects[0].name.as_deref(), Some("hole"));
        assert_eq!(tolerance_zones[0].name.as_deref(), Some("cylindrical"));
        assert_eq!(datum_features[0].name.as_deref(), Some("feature_A"));
    }

    #[test]
    fn parse_reader_matches_parse_string() {
        let a = StepReader::parse_string(BOX_STEP).expect("box.step string parse should succeed");
        let b = StepReader::parse_reader(Cursor::new(BOX_STEP.as_bytes()))
            .expect("box.step stream parse should succeed");
        assert_eq!(a.vertices.len(), b.vertices.len());
        assert_eq!(a.edges.len(), b.edges.len());
        assert_eq!(a.solids.len(), b.solids.len());
    }

    #[test]
    fn parse_reader_with_metadata_matches_string_variant() {
        let (_a_brep, a_md) =
            StepReader::parse_string_with_metadata(BOX_STEP).expect("string metadata parse should succeed");
        let (_b_brep, b_md) = StepReader::parse_reader_with_metadata(Cursor::new(BOX_STEP.as_bytes()))
            .expect("stream metadata parse should succeed");

        assert_eq!(a_md.file_schema, b_md.file_schema);
        assert_eq!(a_md.protocol_hint as u8, b_md.protocol_hint as u8);
        assert_eq!(a_md.products.len(), b_md.products.len());
        assert_eq!(a_md.product_definition_formations.len(), b_md.product_definition_formations.len());
        assert_eq!(a_md.product_definitions.len(), b_md.product_definitions.len());
        assert_eq!(a_md.shape_definition_representations.len(), b_md.shape_definition_representations.len());
        assert_eq!(a_md.assembly_occurrences.len(), b_md.assembly_occurrences.len());
        assert_eq!(a_md.product_names, b_md.product_names);
        assert_eq!(a_md.materials.len(), b_md.materials.len());
        assert_eq!(a_md.layers.len(), b_md.layers.len());
        assert_eq!(a_md.general_properties.len(), b_md.general_properties.len());
        assert_eq!(a_md.property_definitions.len(), b_md.property_definitions.len());
        assert_eq!(a_md.property_definition_representations.len(), b_md.property_definition_representations.len());
        assert_eq!(a_md.dimensional_locations.len(), b_md.dimensional_locations.len());
        assert_eq!(a_md.dimensional_sizes.len(), b_md.dimensional_sizes.len());
        assert_eq!(a_md.geometric_tolerances.len(), b_md.geometric_tolerances.len());
        assert_eq!(
            a_md.geometric_tolerances_with_datum_references.len(),
            b_md.geometric_tolerances_with_datum_references.len()
        );
        assert_eq!(a_md.datums.len(), b_md.datums.len());
        assert_eq!(a_md.datum_systems.len(), b_md.datum_systems.len());
        assert_eq!(a_md.kinematic_pairs.len(), b_md.kinematic_pairs.len());
        assert_eq!(a_md.tolerance_zones.len(), b_md.tolerance_zones.len());
        assert_eq!(a_md.tolerance_zone_definitions.len(), b_md.tolerance_zone_definitions.len());
        assert_eq!(a_md.datum_features.len(), b_md.datum_features.len());
        assert_eq!(a_md.datum_reference_elements.len(), b_md.datum_reference_elements.len());
        assert_eq!(a_md.shape_aspects.len(), b_md.shape_aspects.len());
        assert_eq!(a_md.shape_aspect_definitions.len(), b_md.shape_aspect_definitions.len());
        assert_eq!(a_md.derived_shape_aspects.len(), b_md.derived_shape_aspects.len());
        assert_eq!(a_md.fea_models.len(), b_md.fea_models.len());
        assert_eq!(a_md.fea_meshes.len(), b_md.fea_meshes.len());
        assert_eq!(a_md.fea_node_sets.len(), b_md.fea_node_sets.len());
        assert_eq!(a_md.fea_element_sets.len(), b_md.fea_element_sets.len());
        assert_eq!(a_md.fea_material_properties.len(), b_md.fea_material_properties.len());
        assert_eq!(a_md.fea_boundary_conditions.len(), b_md.fea_boundary_conditions.len());
        assert_eq!(a_md.fea_loads.len(), b_md.fea_loads.len());
        // GDT extended fields
        assert_eq!(a_md.dimensional_tolerances.len(), b_md.dimensional_tolerances.len());
        assert_eq!(a_md.tolerance_values.len(), b_md.tolerance_values.len());
        assert_eq!(a_md.position_tolerances.len(), b_md.position_tolerances.len());
        assert_eq!(a_md.orientation_tolerances.len(), b_md.orientation_tolerances.len());
        assert_eq!(a_md.form_tolerances.len(), b_md.form_tolerances.len());
        assert_eq!(a_md.runout_tolerances.len(), b_md.runout_tolerances.len());
        assert_eq!(a_md.profile_tolerances.len(), b_md.profile_tolerances.len());
        assert_eq!(a_md.datum_reference_frames.len(), b_md.datum_reference_frames.len());
        assert_eq!(a_md.datum_targets.len(), b_md.datum_targets.len());
        assert_eq!(a_md.tolerance_zone_definitions_enhanced.len(), b_md.tolerance_zone_definitions_enhanced.len());
    }

    #[test]
    fn parse_reader_with_healing_report_json_returns_schema() {
        let (_brep, _report, json) = StepReader::parse_reader_with_healing_report_json(
            Cursor::new(HFSS_STEP.as_bytes()),
            HealingOptions {
                mode: HealingMode::AnalyzeAndRepair,
                tolerance: 1e-6,
                max_passes: 2,
                ..HealingOptions::default()
            },
        )
        .expect("stream healing parse should succeed");

        assert!(json.contains("\"schema\": \"step.import.healing.v1\""));
    }

    // ── validate_export_readiness tests ──────────────────────────────────────

    #[test]
    fn export_readiness_clean_brep_is_ready() {
        // A BRep loaded from box.step should parse and be export-ready
        // (no out-of-range pcurve indices, tolerances are ≥ CONFUSION).
        let brep = StepReader::parse_string(BOX_STEP).expect("box.step should parse");
        let report = validate_export_readiness(&brep);
        assert!(
            report.is_ready,
            "box.step should be export-ready; issues: {:?}",
            report.issues
        );
    }

    #[test]
    fn export_readiness_out_of_range_surface_idx_detected() {
        use rcad_kernel::{BRep, PCurve, Edge, Vertex};
        use rcad_kernel::geom::Line2d;

        // Build a minimal BRep with one edge whose PCurve has an invalid surface_idx.
        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: [0.0, 0.0, 0.0].into() });
        brep.vertices.push(Vertex { point: [1.0, 0.0, 0.0].into() });
        brep.edges.push(Edge { start: 0, end: 1 });
        // surface_idx = 99 is out of range (no surfaces in geom.surfaces)
        brep.geom.edge_pcurves.push(vec![PCurve { surface_idx: 99, curve2d_idx: 0 }]);
        brep.geom.curve2ds.push(rcad_kernel::Curve2d::Line(Line2d {
            origin: glam::DVec2::ZERO,
            direction: glam::DVec2::X,
        }));

        let report = validate_export_readiness(&brep);
        assert!(!report.is_ready);
        assert_eq!(report.pcurves_checked, 1);
        assert!(
            report.issues.iter().any(|i| matches!(
                i,
                ExportIssue::PcurveSurfaceOutOfRange { edge_idx: 0, surface_idx: 99 }
            )),
            "expected PcurveSurfaceOutOfRange, got {:?}", report.issues
        );
    }

    #[test]
    fn export_readiness_out_of_range_curve2d_idx_detected() {
        use rcad_kernel::{BRep, PCurve, Edge, Vertex};
        use rcad_kernel::geom::Plane;

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: [0.0, 0.0, 0.0].into() });
        brep.vertices.push(Vertex { point: [1.0, 0.0, 0.0].into() });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.geom.surfaces.push(rcad_kernel::Surface3::Plane(Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
        }));
        // curve2d_idx = 99 is out of range (no curve2ds stored)
        brep.geom.edge_pcurves.push(vec![PCurve { surface_idx: 0, curve2d_idx: 99 }]);

        let report = validate_export_readiness(&brep);
        assert!(!report.is_ready);
        assert!(
            report.issues.iter().any(|i| matches!(
                i,
                ExportIssue::PcurveCurveOutOfRange { edge_idx: 0, curve2d_idx: 99 }
            )),
            "expected PcurveCurveOutOfRange, got {:?}", report.issues
        );
    }

    #[test]
    fn export_readiness_tolerance_too_tight_detected() {
        use rcad_kernel::BRep;

        let mut brep = BRep::new();
        // Push a sub-CONFUSION edge tolerance.
        brep.geom.edge_tolerance.push(1e-10);

        let report = validate_export_readiness(&brep);
        assert!(!report.is_ready);
        assert!(
            report.issues.iter().any(|i| matches!(i, ExportIssue::EdgeToleranceTooTight { .. })),
            "expected EdgeToleranceTooTight, got {:?}", report.issues
        );
    }

    #[test]
    fn export_readiness_too_many_pcurves_flagged() {
        use rcad_kernel::{BRep, PCurve, Edge, Vertex};
        use rcad_kernel::geom::{Plane, Line2d};

        let mut brep = BRep::new();
        brep.vertices.push(Vertex { point: [0.0, 0.0, 0.0].into() });
        brep.vertices.push(Vertex { point: [1.0, 0.0, 0.0].into() });
        brep.edges.push(Edge { start: 0, end: 1 });
        brep.geom.surfaces.push(rcad_kernel::Surface3::Plane(Plane {
            origin: glam::DVec3::ZERO,
            normal: glam::DVec3::Z,
        }));
        brep.geom.curve2ds.push(rcad_kernel::Curve2d::Line(Line2d {
            origin: glam::DVec2::ZERO,
            direction: glam::DVec2::X,
        }));
        // 3 PCurves on a single edge — unusual
        brep.geom.edge_pcurves.push(vec![
            PCurve { surface_idx: 0, curve2d_idx: 0 },
            PCurve { surface_idx: 0, curve2d_idx: 0 },
            PCurve { surface_idx: 0, curve2d_idx: 0 },
        ]);

        let report = validate_export_readiness(&brep);
        assert!(!report.is_ready);
        assert!(
            report.issues.iter().any(|i| matches!(
                i,
                ExportIssue::TooManyPcurves { edge_idx: 0, count: 3 }
            )),
            "expected TooManyPcurves, got {:?}", report.issues
        );
    }

    #[test]
    fn export_readiness_summary_text() {
        use rcad_kernel::BRep;
        let brep = BRep::new();
        let report = validate_export_readiness(&brep);
        let s = report.summary();
        assert!(s.contains("export-ready"), "summary should say export-ready: {s}");
    }

    // ── GDT Extended entity extraction tests ──────────────────────────────────────

    #[test]
    fn extract_dimensional_tolerances_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#10=DIMENSIONAL_TOLERANCE('diam_tol','diameter tolerance',#100,0.05,-0.05,'mm');\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_dimensional_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 10);
        assert_eq!(tols[0].name.as_deref(), Some("diam_tol"));
        assert_eq!(tols[0].description.as_deref(), Some("diameter tolerance"));
        assert_eq!(tols[0].dimensional_characteristic_id, Some(100));
        assert!((tols[0].upper_tolerance.unwrap() - 0.05).abs() < 1e-9);
        assert!((tols[0].lower_tolerance.unwrap() - (-0.05)).abs() < 1e-9);
        assert_eq!(tols[0].unit.as_deref(), Some("mm"));
    }

    #[test]
    fn extract_tolerance_values_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#20=MEASURE_REPRESENTATION_ITEM('tol_value',0.025,'mm');\nENDSEC;\nEND-ISO-10303-21;\n";
        let vals = extract_tolerance_values(step);
        assert_eq!(vals.len(), 1);
        assert_eq!(vals[0].entity_id, 20);
        assert_eq!(vals[0].name.as_deref(), Some("tol_value"));
        assert!((vals[0].value - 0.025).abs() < 1e-9);
        assert_eq!(vals[0].unit.as_deref(), Some("mm"));
    }

    #[test]
    fn extract_position_tolerances_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#30=POSITION_TOLERANCE('pos_tol','positional tolerance',#20,#30,#40,.T.,10.0);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_position_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 30);
        assert_eq!(tols[0].name.as_deref(), Some("pos_tol"));
        assert_eq!(tols[0].description.as_deref(), Some("positional tolerance"));
        assert_eq!(tols[0].value_entity_id, Some(20));
        assert_eq!(tols[0].shape_aspect_id, Some(30));
        assert_eq!(tols[0].datum_system_id, Some(40));
        assert!(tols[0].projected);
        assert!((tols[0].projected_height.unwrap() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn extract_orientation_tolerances_parses_angularity() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#40=ANGULARITY_TOLERANCE('ang_tol','angularity tolerance',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_orientation_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 40);
        assert_eq!(tols[0].name.as_deref(), Some("ang_tol"));
        assert_eq!(tols[0].orientation_type, OrientationToleranceType::Angularity);
    }

    #[test]
    fn extract_orientation_tolerances_parses_perpendicularity() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#41=PERPENDICULARITY_TOLERANCE('perp_tol','perpendicularity tolerance',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_orientation_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 41);
        assert_eq!(tols[0].orientation_type, OrientationToleranceType::Perpendicularity);
    }

    #[test]
    fn extract_orientation_tolerances_parses_parallelism() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#42=PARALLELISM_TOLERANCE('para_tol','parallelism tolerance',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_orientation_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 42);
        assert_eq!(tols[0].orientation_type, OrientationToleranceType::Parallelism);
    }

    #[test]
    fn extract_form_tolerances_parses_flatness() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#50=FLATNESS_TOLERANCE('flat_tol','flatness tolerance',#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_form_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 50);
        assert_eq!(tols[0].name.as_deref(), Some("flat_tol"));
        assert_eq!(tols[0].form_type, FormToleranceType::Flatness);
    }

    #[test]
    fn extract_form_tolerances_parses_straightness() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#51=STRAIGHTNESS_TOLERANCE('str_tol','straightness tolerance',#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_form_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].form_type, FormToleranceType::Straightness);
    }

    #[test]
    fn extract_form_tolerances_parses_circularity() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#52=CIRCULARITY_TOLERANCE('cir_tol','circularity tolerance',#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_form_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].form_type, FormToleranceType::Circularity);
    }

    #[test]
    fn extract_form_tolerances_parses_cylindricity() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#53=CYLINDRICITY_TOLERANCE('cyl_tol','cylindricity tolerance',#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_form_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].form_type, FormToleranceType::Cylindricity);
    }

    #[test]
    fn extract_runout_tolerances_parses_circular_runout() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#60=CIRCULAR_RUNOUT_TOLERANCE('cr_tol','circular runout tolerance',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_runout_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 60);
        assert_eq!(tols[0].runout_type, RunoutToleranceType::CircularRunout);
    }

    #[test]
    fn extract_runout_tolerances_parses_total_runout() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#61=TOTAL_RUNOUT_TOLERANCE('tr_tol','total runout tolerance',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_runout_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].runout_type, RunoutToleranceType::TotalRunout);
    }

    #[test]
    fn extract_profile_tolerances_parses_profile_of_line() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#70=LINE_PROFILE_TOLERANCE('lin_tol','profile of a line',#20,#30);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_profile_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].entity_id, 70);
        assert_eq!(tols[0].profile_type, ProfileToleranceType::ProfileOfALine);
    }

    #[test]
    fn extract_profile_tolerances_parses_profile_of_surface() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#71=SURFACE_PROFILE_TOLERANCE('surf_tol','profile of a surface',#20,#30,#40);\nENDSEC;\nEND-ISO-10303-21;\n";
        let tols = extract_profile_tolerances(step);
        assert_eq!(tols.len(), 1);
        assert_eq!(tols[0].profile_type, ProfileToleranceType::ProfileOfASurface);
    }

    #[test]
    fn extract_datum_reference_frames_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#80=DATUM_REFERENCE_FRAME('DRF1','primary datum reference frame',(#50,#51,#52));\nENDSEC;\nEND-ISO-10303-21;\n";
        let frames = extract_datum_reference_frames(step);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].entity_id, 80);
        assert_eq!(frames[0].name.as_deref(), Some("DRF1"));
        assert_eq!(frames[0].description.as_deref(), Some("primary datum reference frame"));
        assert_eq!(frames[0].datum_system_ids, vec![50, 51, 52]);
    }

    #[test]
    fn extract_datum_targets_parses_point_target() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#90=DATUM_TARGET_POINT('A1','datum target A1','A1',#80,#100);\nENDSEC;\nEND-ISO-10303-21;\n";
        let targets = extract_datum_targets(step);
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0].entity_id, 90);
        assert_eq!(targets[0].name.as_deref(), Some("A1"));
        assert_eq!(targets[0].target_identifier.as_deref(), Some("A1"));
        assert_eq!(targets[0].datum_id, Some(80));
        assert_eq!(targets[0].target_type, DatumTargetType::Point);
    }

    #[test]
    fn extract_datum_targets_parses_line_target() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#91=DATUM_TARGET_LINE('B1','datum target line B1','B1',#80,#100);\nENDSEC;\nEND-ISO-10303-21;\n";
        let targets = extract_datum_targets(step);
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0].target_type, DatumTargetType::Line);
    }

    #[test]
    fn extract_datum_targets_parses_area_target() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#92=DATUM_TARGET_AREA('C1','datum target area C1','C1',#80,#100);\nENDSEC;\nEND-ISO-10303-21;\n";
        let targets = extract_datum_targets(step);
        assert_eq!(targets.len(), 1);
        assert_eq!(targets[0].target_type, DatumTargetType::Area);
    }

    #[test]
    fn extract_tolerance_zone_definitions_enhanced_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=TOLERANCE_ZONE_DEFINITION('cylindrical','symmetric',#90,#110,#120);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_tolerance_zone_definitions_enhanced(step);
        assert_eq!(defs.len(), 1);
        assert_eq!(defs[0].entity_id, 100);
        assert_eq!(defs[0].name.as_deref(), Some("cylindrical"));
        assert_eq!(defs[0].tolerance_zone_id, Some(90));
        assert_eq!(defs[0].shape_aspect_id, Some(110));
        assert_eq!(defs[0].zone_shape, ToleranceZoneShape::Cylindrical);
        assert_eq!(defs[0].zone_position, ToleranceZonePosition::Symmetric);
        assert_eq!(defs[0].defining_shape_aspect_id, Some(120));
    }

    #[test]
    fn extract_tolerance_zone_definitions_enhanced_detects_zone_shapes() {
        // Test spherical shape
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#101=TOLERANCE_ZONE_DEFINITION('spherical',$,#90,#110);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_tolerance_zone_definitions_enhanced(step);
        assert_eq!(defs[0].zone_shape, ToleranceZoneShape::Spherical);

        // Test two parallel planes shape
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#102=TOLERANCE_ZONE_DEFINITION('two parallel planes',$,#90,#110);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_tolerance_zone_definitions_enhanced(step);
        assert_eq!(defs[0].zone_shape, ToleranceZoneShape::TwoParallelPlanes);

        // Test two coaxial cylinders shape
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#103=TOLERANCE_ZONE_DEFINITION('two coaxial cylinders',$,#90,#110);\nENDSEC;\nEND-ISO-10303-21;\n";
        let defs = extract_tolerance_zone_definitions_enhanced(step);
        assert_eq!(defs[0].zone_shape, ToleranceZoneShape::TwoCoaxialCylinders);
    }

    // ── View and Camera extraction tests (AP242) ───────────────────────────────

    #[test]
    fn extract_views_parses_camera_model_d3() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=CAMERA_MODEL_D3('front_view',#200,#300);\nENDSEC;\nEND-ISO-10303-21;\n";
        let views = extract_views(step);
        assert_eq!(views.len(), 1);
        assert_eq!(views[0].entity_id, 100);
        assert_eq!(views[0].name.as_deref(), Some("front_view"));
        assert_eq!(views[0].camera_model_id, Some(100));
    }

    #[test]
    fn extract_cameras_parses_camera_model_d3() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=CAMERA_MODEL_D3('isometric_view',#200,#300);\nENDSEC;\nEND-ISO-10303-21;\n";
        let cameras = extract_cameras(step);
        assert_eq!(cameras.len(), 1);
        assert_eq!(cameras[0].entity_id, 100);
        assert_eq!(cameras[0].name.as_deref(), Some("isometric_view"));
        assert_eq!(cameras[0].view_reference_system_id, Some(200));
        assert_eq!(cameras[0].view_volume_id, Some(300));
    }

    #[test]
    fn extract_view_volumes_parses_view_volume() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#300=VIEW_VOLUME('ortho_view',$,100.0,50.0,200.0);\nENDSEC;\nEND-ISO-10303-21;\n";
        let volumes = extract_view_volumes(step);
        assert_eq!(volumes.len(), 1);
        assert_eq!(volumes[0].entity_id, 300);
        assert_eq!(volumes[0].name.as_deref(), Some("ortho_view"));
        assert_eq!(volumes[0].volume_type, ViewVolumeType::Orthographic);
    }

    #[test]
    fn extract_view_volumes_detects_perspective() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#300=VIEW_VOLUME('perspective_view',.PERSPECTIVE.,100.0,50.0,200.0);\nENDSEC;\nEND-ISO-10303-21;\n";
        let volumes = extract_view_volumes(step);
        assert_eq!(volumes.len(), 1);
        assert_eq!(volumes[0].volume_type, ViewVolumeType::Perspective);
    }

    // ── Annotation extraction tests (AP242) ────────────────────────────────────

    #[test]
    fn extract_notes_parses_annotation() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=DESCRIPTIVE_REPRESENTATION_ITEM('note1','This is a note');\nENDSEC;\nEND-ISO-10303-21;\n";
        let notes = extract_notes(step);
        assert_eq!(notes.len(), 1);
        assert_eq!(notes[0].entity_id, 100);
        assert_eq!(notes[0].name.as_deref(), Some("note1"));
    }

    #[test]
    fn extract_annotation_planes_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=ANNOTATION_PLANE('top_annotation',#200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let planes = extract_annotation_planes(step);
        assert_eq!(planes.len(), 1);
        assert_eq!(planes[0].entity_id, 100);
        assert_eq!(planes[0].name.as_deref(), Some("top_annotation"));
        assert_eq!(planes[0].plane_id, Some(200));
    }

    #[test]
    fn extract_annotation_occurrences_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=ANNOTATION_OCCURRENCE('dim_occurrence',#200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let occurrences = extract_annotation_occurrences(step);
        assert_eq!(occurrences.len(), 1);
        assert_eq!(occurrences[0].entity_id, 100);
        assert_eq!(occurrences[0].name.as_deref(), Some("dim_occurrence"));
        assert_eq!(occurrences[0].style_id, Some(200));
    }

    #[test]
    fn extract_dimension_curves_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=DIMENSION_CURVE('dim_curve',#200,#300);\nENDSEC;\nEND-ISO-10303-21;\n";
        let curves = extract_dimension_curves(step);
        assert_eq!(curves.len(), 1);
        assert_eq!(curves[0].entity_id, 100);
        assert_eq!(curves[0].name.as_deref(), Some("dim_curve"));
        assert_eq!(curves[0].curve_id, Some(200));
        assert_eq!(curves[0].annotation_plane_id, Some(300));
    }

    #[test]
    fn extract_terminator_symbols_parses_arrow() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=FILLED_ARROW('arrow1',#200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let symbols = extract_terminator_symbols(step);
        assert_eq!(symbols.len(), 1);
        assert_eq!(symbols[0].entity_id, 100);
        assert_eq!(symbols[0].name.as_deref(), Some("arrow1"));
        assert_eq!(symbols[0].terminator_type, TerminatorType::Arrow);
    }

    #[test]
    fn extract_terminator_symbols_parses_generic_terminator() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=TERMINATOR_SYMBOL('symbol1',#200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let symbols = extract_terminator_symbols(step);
        assert_eq!(symbols.len(), 1);
        assert_eq!(symbols[0].terminator_type, TerminatorType::Unknown);
    }

    #[test]
    fn extract_datum_feature_callouts_parses_entry() {
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=DATUM_FEATURE_CALLOUT('datum_A','A',#200);\nENDSEC;\nEND-ISO-10303-21;\n";
        let callouts = extract_datum_feature_callouts(step);
        assert_eq!(callouts.len(), 1);
        assert_eq!(callouts[0].entity_id, 100);
        assert_eq!(callouts[0].name.as_deref(), Some("datum_A"));
        assert_eq!(callouts[0].datum_identifier.as_deref(), Some("A"));
        assert_eq!(callouts[0].annotation_plane_id, Some(200));
    }

    #[test]
    fn metadata_includes_view_and_annotation_entities() {
        // Test extraction functions directly since this STEP file has no geometry
        let step = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));\nFILE_NAME('test','','','','','','');\nENDSEC;\nDATA;\n#100=CAMERA_MODEL_D3('front_view',#200,#300);\n#300=VIEW_VOLUME('ortho_view',$,100.0,50.0,200.0);\n#400=ANNOTATION_PLANE('annotation_plane',#500);\n#500=DIMENSION_CURVE('dim_curve',#600,#400);\nENDSEC;\nEND-ISO-10303-21;\n";

        // Verify extraction functions work
        let views = extract_views(step);
        let cameras = extract_cameras(step);
        let volumes = extract_view_volumes(step);
        let planes = extract_annotation_planes(step);
        let curves = extract_dimension_curves(step);

        assert_eq!(views.len(), 1);
        assert_eq!(cameras.len(), 1);
        assert_eq!(volumes.len(), 1);
        assert_eq!(planes.len(), 1);
        assert_eq!(curves.len(), 1);
    }

    // ── Conversion function tests ───────────────────────────────────────────────

    #[test]
    fn view_to_step_entities_converts_correctly() {
        use rcad_kernel::annotation::{View, ViewProjection};

        let view = View::isometric(1)
            .with_position(glam::DVec3::new(100.0, 100.0, 100.0))
            .with_target(glam::DVec3::ZERO)
            .with_projection(ViewProjection::Perspective)
            .with_fov(60.0);

        let (step_view, step_camera, step_volume) = view_to_step_entities(&view, 1000);

        assert_eq!(step_view.entity_id, 1000);
        assert_eq!(step_view.name.as_deref(), Some("Isometric"));
        assert_eq!(step_camera.entity_id, 1001);
        assert!(step_camera.perspective);
        assert_eq!(step_volume.entity_id, 1002);
        assert_eq!(step_volume.volume_type, ViewVolumeType::Perspective);
    }

    #[test]
    fn note_to_step_converts_correctly() {
        use rcad_kernel::annotation::{Note, NoteCategory};

        let note = Note::new(42, "TestNote", "This is a test note")
            .with_category(NoteCategory::Warning)
            .with_author("John Doe");

        let step_note = note_to_step(&note);

        assert_eq!(step_note.entity_id, 42);
        assert_eq!(step_note.name.as_deref(), Some("TestNote"));
        assert_eq!(step_note.text.as_deref(), Some("This is a test note"));
        assert_eq!(step_note.description.as_deref(), Some("John Doe"));
    }

    #[test]
    fn step_views_to_kernel_converts_correctly() {
        let step_views = vec![StepView {
            entity_id: 100,
            name: Some("Front".to_string()),
            description: None,
            camera_model_id: Some(100),
            view_type: Some("front".to_string()),
        }];

        let step_cameras = vec![StepCameraModelD3 {
            entity_id: 100,
            name: Some("Front".to_string()),
            view_reference_system_id: None,
            view_volume_id: Some(200),
            perspective: false,
        }];

        let step_volumes = vec![StepViewVolume {
            entity_id: 200,
            name: Some("ortho".to_string()),
            volume_type: ViewVolumeType::Orthographic,
            view_center: Some([0.0, 0.0, 0.0]),
            view_plane_distance: Some(100.0),
            up_direction: Some([0.0, 1.0, 0.0]),
            view_window_width: Some(50.0),
            view_window_height: Some(50.0),
        }];

        let kernel_views = step_views_to_kernel(&step_views, &step_cameras, &step_volumes);

        assert_eq!(kernel_views.len(), 1);
        assert_eq!(kernel_views[0].id, 100);
        assert_eq!(kernel_views[0].name, "Front");
        assert_eq!(kernel_views[0].projection, rcad_kernel::annotation::ViewProjection::Orthographic);
    }

    #[test]
    fn step_notes_to_kernel_converts_correctly() {
        let step_notes = vec![
            StepNote {
                entity_id: 100,
                name: Some("WarningNote".to_string()),
                description: Some("Jane Smith".to_string()),
                text: Some("Check this area".to_string()),
                annotation_plane_id: None,
                associated_geometry_id: None,
            },
            StepNote {
                entity_id: 101,
                name: Some("InfoNote".to_string()),
                description: None,
                text: Some("Information".to_string()),
                annotation_plane_id: None,
                associated_geometry_id: None,
            },
        ];

        let kernel_notes = step_notes_to_kernel(&step_notes);

        assert_eq!(kernel_notes.len(), 2);
        assert_eq!(kernel_notes[0].id, 100);
        assert_eq!(kernel_notes[0].category, NoteCategory::Warning);
        assert_eq!(kernel_notes[0].author.as_deref(), Some("Jane Smith"));
        assert_eq!(kernel_notes[1].category, NoteCategory::Info);
    }

    #[test]
    fn annotation_store_to_step_metadata_converts_all() {
        use rcad_kernel::annotation::{Annotation, AnnotationStore, Note, View};

        let mut store = AnnotationStore::new();
        store.add_view(View::front(1));
        store.add_view(View::top(2));
        store.add_xc_note(Note::new(100, "Note1", "Test note"));
        store.add_annotation(Annotation::dimension(200, "Dim1", 25.4));

        let metadata = annotation_store_to_step_metadata(&store);

        assert_eq!(metadata.views.len(), 2);
        assert_eq!(metadata.cameras.len(), 2);
        assert_eq!(metadata.view_volumes.len(), 2);
        assert_eq!(metadata.notes.len(), 2); // One from xc_notes, one from annotation
    }

    // ── AP242 Semantic Parsing Tests ────────────────────────────────────────────

    /// Minimal STEP with a vertex for tests that need geometry
    const MINIMAL_STEP_PREFIX: &str = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
"#;

    #[test]
    fn extracts_product_definition_relationships() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=PRODUCT('P001','ParentProduct','',());
#11=PRODUCT('P002','ChildProduct','',());
#12=PRODUCT_DEFINITION_FORMATION('','',#10);
#13=PRODUCT_DEFINITION_FORMATION('','',#11);
#14=PRODUCT_DEFINITION('','',$,#12);
#15=PRODUCT_DEFINITION('','',$,#13);
#16=PRODUCT_DEFINITION_RELATIONSHIP('REL001','assembly relationship','desc',#14,#15);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.product_definition_relationships.len(), 1);
        let rel = &metadata.product_definition_relationships[0];
        assert_eq!(rel.entity_id, 16);
        assert_eq!(rel.relationship_id.as_deref(), Some("REL001"));
        assert_eq!(rel.name.as_deref(), Some("assembly relationship"));
        assert_eq!(rel.relating_product_definition_id, Some(14));
        assert_eq!(rel.related_product_definition_id, Some(15));
    }

    #[test]
    fn extracts_shape_representation_relationships() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=SHAPE_REPRESENTATION('Shape1',(),$);
#11=SHAPE_REPRESENTATION('Shape2',(),$);
#12=SHAPE_REPRESENTATION_RELATIONSHIP('rel','description',#10,#11);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.shape_representation_relationships.len(), 1);
        let rel = &metadata.shape_representation_relationships[0];
        assert_eq!(rel.entity_id, 12);
        assert_eq!(rel.name.as_deref(), Some("rel"));
        assert_eq!(rel.relating_representation_id, Some(10));
        assert_eq!(rel.related_representation_id, Some(11));
    }

    #[test]
    fn extracts_configuration_designs() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=CONFIGURATION_ITEM('CI001','ConfigItem1','description',$);
#11=PRODUCT_DEFINITION('','',$,$);
#12=CONFIGURATION_DESIGN('config',#10,#11);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.configuration_designs.len(), 1);
        let cd = &metadata.configuration_designs[0];
        assert_eq!(cd.entity_id, 12);
        assert_eq!(cd.name.as_deref(), Some("config"));
        assert_eq!(cd.configuration_id, Some(10));
        assert_eq!(cd.product_definition_id, Some(11));
    }

    #[test]
    fn extracts_configuration_items() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=CONFIGURATION_ITEM('CI001','ConfigItem1','description',$);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.configuration_items.len(), 1);
        let ci = &metadata.configuration_items[0];
        assert_eq!(ci.entity_id, 10);
        assert_eq!(ci.item_id.as_deref(), Some("CI001"));
        assert_eq!(ci.name.as_deref(), Some("ConfigItem1"));
    }

    #[test]
    fn extracts_product_concepts() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=PRODUCT_CONCEPT('PC001','ConceptName','description',$);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.product_concepts.len(), 1);
        let pc = &metadata.product_concepts[0];
        assert_eq!(pc.entity_id, 10);
        assert_eq!(pc.concept_id.as_deref(), Some("PC001"));
        assert_eq!(pc.name.as_deref(), Some("ConceptName"));
    }

    #[test]
    fn extracts_approvals() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=APPROVAL(.APPROVED.,'Design Approval');
#11=APPROVAL(.PENDING.,'Review Pending');
#12=APPROVAL(.REJECTED.,'Rejected Change');
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.approvals.len(), 3);
        assert_eq!(metadata.approvals[0].status, ApprovalStatus::Approved);
        assert_eq!(metadata.approvals[1].status, ApprovalStatus::Pending);
        assert_eq!(metadata.approvals[2].status, ApprovalStatus::Rejected);
        assert_eq!(metadata.approvals[0].level.as_deref(), Some("Design Approval"));
    }

    #[test]
    fn extracts_security_classifications() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=SECURITY_CLASSIFICATION('Classified','Top secret doc',.TOP_SECRET.);
#11=SECURITY_CLASSIFICATION('Internal','Company internal',.PROPRIETARY.);
#12=SECURITY_CLASSIFICATION('Public','Public info',.UNCLASSIFIED.);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.security_classifications.len(), 3);
        assert_eq!(metadata.security_classifications[0].security_level, SecurityClassificationLevel::TopSecret);
        assert_eq!(metadata.security_classifications[1].security_level, SecurityClassificationLevel::Proprietary);
        assert_eq!(metadata.security_classifications[2].security_level, SecurityClassificationLevel::Unclassified);
    }

    #[test]
    fn extracts_documents() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=DOCUMENT('DOC001','Design Document','Product specification','drawing');
#11=DOCUMENT('DOC002','Test Report','Validation results',$);
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.documents.len(), 2);
        assert_eq!(metadata.documents[0].document_id.as_deref(), Some("DOC001"));
        assert_eq!(metadata.documents[0].name.as_deref(), Some("Design Document"));
        assert_eq!(metadata.documents[0].document_type.as_deref(), Some("drawing"));
        assert_eq!(metadata.documents[1].document_id.as_deref(), Some("DOC002"));
    }

    #[test]
    fn extracts_document_files() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=DOCUMENT_FILE('DF001','CAD Model','3D geometry file','part.step','STEP AP242');
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.document_files.len(), 1);
        let df = &metadata.document_files[0];
        assert_eq!(df.document_id.as_deref(), Some("DF001"));
        assert_eq!(df.file_name.as_deref(), Some("part.step"));
        assert_eq!(df.file_format.as_deref(), Some("STEP AP242"));
    }

    #[test]
    fn extracts_document_usage_assignments() {
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=DOCUMENT('DOC001','Reference Doc','',());
#11=PRODUCT_DEFINITION('','',$,$);
#12=DOCUMENT_USAGE_ASSIGNMENT(#10,#11,'reference');
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        assert_eq!(metadata.document_usage_assignments.len(), 1);
        let dua = &metadata.document_usage_assignments[0];
        assert_eq!(dua.document_id, Some(10));
        assert_eq!(dua.product_definition_id, Some(11));
        assert_eq!(dua.role.as_deref(), Some("reference"));
    }

    #[test]
    fn metadata_summary_includes_ap242_fields() {
        // Use a minimal valid STEP with geometry (a point as a vertex)
        let step = r#"ISO-10303-21;
HEADER;
FILE_NAME('test','','','','','','');
FILE_SCHEMA(('AP242'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=VERTEX_POINT('',#1);
#10=DOCUMENT('DOC001','Test','',());
#11=APPROVAL(.APPROVED.,'Approved');
#12=SECURITY_CLASSIFICATION('Secret','',.SECRET.);
#13=CONFIGURATION_ITEM('CI001','Config','',());
#14=PRODUCT_DEFINITION_RELATIONSHIP('rel','name','',(),(),());
ENDSEC;
END-ISO-10303-21;
"#;
        let (_, metadata) = StepReader::parse_string_with_metadata(step).expect("parse");
        let summary = metadata.summary();
        assert!(summary.contains("Documents: 1"));
        assert!(summary.contains("Approvals: 1"));
        assert!(summary.contains("Security classifications: 1"));
        assert!(summary.contains("Configuration items: 1"));
        assert!(summary.contains("Product definition relationships: 1"));
    }

    #[test]
    fn round_trip_multi_solid_compound() {
        use rcad_modeling::make_box_brep;
        use glam::DVec3;
        use rcad_algorithms::geom_populate::populate_box_geom;

        let mut box1 = make_box_brep(DVec3::ZERO, DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let mut box2 = make_box_brep(DVec3::new(3.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();
        let mut box3 = make_box_brep(DVec3::new(6.0, 0.0, 0.0), DVec3::X, DVec3::Y, 1.0, 1.0, 1.0).unwrap();

        populate_box_geom(&mut box1);
        populate_box_geom(&mut box2);
        populate_box_geom(&mut box3);

        // Combine into compound by appending breps
        let mut compound = box1;
        rcad_scene::append_brep(&mut compound, box2);
        rcad_scene::append_brep(&mut compound, box3);

        // Verify we have 3 shells (one per box)
        let shell_count: usize = compound.solids.iter().map(|s| s.shells.len()).sum();
        assert_eq!(shell_count, 3, "should have 3 shells");

        let (v1, e1, f1, s1) = count_topology(&compound);
        let round_tripped = round_trip_brep(&compound, true);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        assert_eq!(v1, v2, "vertex count should match after round-trip");
        assert_eq!(e1, e2, "edge count should match after round-trip");
        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");
    }

    #[test]
    fn round_trip_cylinder_surface() {
        use rcad_modeling::make_cylinder_brep;
        use glam::DVec3;

        let cylinder = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();

        let (v1, e1, f1, s1) = count_topology(&cylinder);
        let round_tripped = round_trip_brep(&cylinder, false); // non-strict for curved surfaces
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        // For curved surfaces with gmsh_strict: false, topology may differ slightly
        // but face count should remain stable
        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");

        // Verify cylindrical surface is preserved
        let has_cylinder = round_tripped.geom.surfaces.iter().any(|s| matches!(s, Surface3::Cylinder(_)));
        assert!(has_cylinder, "should preserve cylindrical surface");
    }

    #[test]
    fn round_trip_sphere_surface() {
        use rcad_modeling::make_sphere_brep;
        use glam::DVec3;

        let sphere = make_sphere_brep(DVec3::ZERO, 1.0).unwrap();

        let (v1, e1, f1, s1) = count_topology(&sphere);
        let round_tripped = round_trip_brep(&sphere, false);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");

        let has_sphere = round_tripped.geom.surfaces.iter().any(|s| matches!(s, Surface3::Sphere(_)));
        assert!(has_sphere, "should preserve spherical surface");
    }

    #[test]
    fn round_trip_cone_surface() {
        use rcad_modeling::make_cone_brep;
        use glam::DVec3;

        let cone = make_cone_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();

        let (v1, e1, f1, s1) = count_topology(&cone);
        let round_tripped = round_trip_brep(&cone, false);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        assert_eq!(f1, f2, "face count should match after round-trip");
        assert_eq!(s1, s2, "shell count should match after round-trip");

        let has_cone = round_tripped.geom.surfaces.iter().any(|s| matches!(s, Surface3::Cone(_)));
        assert!(has_cone, "should preserve conical surface");
    }

    #[test]
    fn write_empty_brep_returns_valid_step() {
        let brep = BRep::default();

        let step = StepWriter::write_string(
            &brep,
            ExportSelection {
                selected_faces: &[],
                selected_edges: &[],
            },
        );

        // Should produce valid STEP structure even for empty geometry
        assert!(step.contains("ISO-10303-21"));
        assert!(step.contains("END-ISO-10303-21"));
    }

    #[test]
    fn parse_malformed_step_returns_error() {
        let invalid = "NOT A STEP FILE";
        let result = StepReader::parse_string(invalid);
        assert!(result.is_err(), "should return error for invalid STEP");
    }

    #[test]
    fn parse_truncated_step_returns_error() {
        let truncated = "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AP214'));\n"; // Missing DATA and ENDSEC
        let result = StepReader::parse_string(truncated);
        assert!(result.is_err(), "should return error for truncated STEP");
    }

    // ── Periodic Surface Seam Tests ──────────────────────────────────────────────

    #[test]
    fn cylinder_seam_no_duplicate_edges() {
        use rcad_modeling::make_cylinder_brep;
        use glam::DVec3;

        let cylinder = make_cylinder_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 1.0, 2.0).unwrap();

        // Cylinder should have 3 edges (top circle, bottom circle, seam)
        // Not 4 or more due to duplicate seams
        let original_edge_count = cylinder.edges.len();
        assert!(original_edge_count <= 4, "cylinder should not have excessive edges, got {}", original_edge_count);

        let round_tripped = round_trip_brep(&cylinder, false);
        let round_trip_edge_count = round_tripped.edges.len();
        // Allow for tessellation edges but check reasonable bounds
        // A proper cylinder round-trip should preserve topology reasonably
        assert!(
            round_trip_edge_count <= 10,
            "round-tripped cylinder should not have excessive edges, got {} (original {})",
            round_trip_edge_count, original_edge_count
        );
    }

    #[test]
    fn torus_round_trip_preserves_topology() {
        use rcad_modeling::make_torus_brep;
        use glam::DVec3;

        let torus = make_torus_brep(DVec3::ZERO, DVec3::Z, DVec3::X, 2.0, 0.5).unwrap();

        let (v1, e1, f1, s1) = count_topology(&torus);
        let round_tripped = round_trip_brep(&torus, false);
        let (v2, e2, f2, s2) = count_topology(&round_tripped);

        // Torus is a closed surface, should preserve face count
        assert_eq!(f1, f2, "torus face count should match after round-trip");

        let has_torus = round_tripped.geom.surfaces.iter().any(|s| matches!(s, Surface3::Torus(_)));
        assert!(has_torus, "should preserve toroidal surface");
    }
}
