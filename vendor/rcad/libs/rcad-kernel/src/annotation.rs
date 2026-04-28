//! Annotation object model for CAD annotations.
//!
//! Analogous to OCCT `XCAFNoteObjects` - provides data structures for
//! text annotations, leaders, symbols, and other PMI (Product Manufacturing Information).
//!
//! Also provides view definitions analogous to OCCT `XCAFView` for camera states
//! and drawing views.

use glam::DVec3;
use serde::{Deserialize, Serialize};

// ═══════════════════════════════════════════════════════════════════════════════
// NOTE OBJECTS (XCAFNoteObjects equivalent)
// ═══════════════════════════════════════════════════════════════════════════════

/// Category of a note for classification purposes.
///
/// Analogous to XCAFDoc_Note categories in OCCT.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum NoteCategory {
    /// Informational note.
    #[default]
    Info,
    /// Warning note indicating potential issues.
    Warning,
    /// Comment or general remark.
    Comment,
    /// Requirement specification.
    Requirement,
    /// Approval stamp.
    Approval,
    /// Revision marker.
    Revision,
}

/// A textual note that can be attached to geometry.
///
/// Analogous to OCCT `XCAFDoc_Note` - represents user annotations
/// that can be attached to shapes, faces, edges, or vertices.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Note {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name/label.
    pub name: String,
    /// Text content of the note.
    pub text: String,
    /// Category classification.
    pub category: NoteCategory,
    /// Author who created the note.
    pub author: Option<String>,
    /// Timestamp (ISO 8601 format or Unix timestamp).
    pub timestamp: Option<String>,
    /// 3D position in model space for display.
    pub position: DVec3,
    /// Visibility flag.
    pub visibility: bool,
    /// Links to geometry this note is attached to.
    pub links: Vec<NoteLink>,
}

impl Note {
    /// Create a new note with the given ID, name, and text.
    pub fn new(id: u64, name: impl Into<String>, text: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            text: text.into(),
            category: NoteCategory::default(),
            author: None,
            timestamp: None,
            position: DVec3::ZERO,
            visibility: true,
            links: Vec::new(),
        }
    }

    /// Set the note category.
    pub fn with_category(mut self, category: NoteCategory) -> Self {
        self.category = category;
        self
    }

    /// Set the author.
    pub fn with_author(mut self, author: impl Into<String>) -> Self {
        self.author = Some(author.into());
        self
    }

    /// Set the timestamp.
    pub fn with_timestamp(mut self, timestamp: impl Into<String>) -> Self {
        self.timestamp = Some(timestamp.into());
        self
    }

    /// Set the 3D position.
    pub fn with_position(mut self, position: DVec3) -> Self {
        self.position = position;
        self
    }

    /// Add a link to geometry.
    pub fn add_link(&mut self, link: NoteLink) {
        self.links.push(link);
    }
}

/// Link from a note to specific geometry.
///
/// Defines which geometric entity a note is attached to.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NoteLink {
    /// Type of geometry being linked.
    pub target: NoteTarget,
    /// Optional description of the relationship.
    pub description: Option<String>,
}

impl NoteLink {
    /// Create a link to a shape by index.
    pub fn to_shape(shape_index: usize) -> Self {
        Self {
            target: NoteTarget::Shape { shape_index },
            description: None,
        }
    }

    /// Create a link to a face by shape and face indices.
    pub fn to_face(shape_index: usize, face_index: usize) -> Self {
        Self {
            target: NoteTarget::Face {
                shape_index,
                face_index,
            },
            description: None,
        }
    }

    /// Create a link to an edge by shape and edge indices.
    pub fn to_edge(shape_index: usize, edge_index: usize) -> Self {
        Self {
            target: NoteTarget::Edge {
                shape_index,
                edge_index,
            },
            description: None,
        }
    }

    /// Create a link to a vertex by shape and vertex indices.
    pub fn to_vertex(shape_index: usize, vertex_index: usize) -> Self {
        Self {
            target: NoteTarget::Vertex {
                shape_index,
                vertex_index,
            },
            description: None,
        }
    }

    /// Create a link to a 3D point in model space.
    pub fn to_point(point: DVec3) -> Self {
        Self {
            target: NoteTarget::Point { point },
            description: None,
        }
    }

    /// Add a description to the link.
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }
}

/// Target geometry for a note link.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NoteTarget {
    /// Link to a shape (assembly component or solid).
    Shape {
        /// Index in the assembly/shapes list.
        shape_index: usize,
    },
    /// Link to a face within a shape.
    Face {
        /// Index in the assembly/shapes list.
        shape_index: usize,
        /// Face index within the shape.
        face_index: usize,
    },
    /// Link to an edge within a shape.
    Edge {
        /// Index in the assembly/shapes list.
        shape_index: usize,
        /// Edge index within the shape.
        edge_index: usize,
    },
    /// Link to a vertex within a shape.
    Vertex {
        /// Index in the assembly/shapes list.
        shape_index: usize,
        /// Vertex index within the shape.
        vertex_index: usize,
    },
    /// Link to a 3D point in model space.
    Point {
        /// The 3D point coordinates.
        point: DVec3,
    },
    /// Link to an annotation plane.
    AnnotationPlane {
        /// Index of the annotation plane.
        plane_index: usize,
    },
}

// ═══════════════════════════════════════════════════════════════════════════════
// VIEW DEFINITIONS (XCAFView equivalent)
// ═══════════════════════════════════════════════════════════════════════════════

/// Projection type for a camera view.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ViewProjection {
    /// Orthographic (parallel) projection.
    #[default]
    Orthographic,
    /// Perspective projection.
    Perspective,
}

/// Clipping planes for a view.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewClipping {
    /// Near clipping plane distance.
    pub near: f64,
    /// Far clipping plane distance.
    pub far: f64,
    /// Optional front clipping plane enabled.
    pub front_enabled: bool,
    /// Optional back clipping plane enabled.
    pub back_enabled: bool,
}

impl Default for ViewClipping {
    fn default() -> Self {
        Self {
            near: 0.1,
            far: 10000.0,
            front_enabled: true,
            back_enabled: true,
        }
    }
}

impl ViewClipping {
    /// Create new clipping planes.
    pub fn new(near: f64, far: f64) -> Self {
        Self {
            near,
            far,
            front_enabled: true,
            back_enabled: true,
        }
    }
}

/// A camera view definition.
///
/// Analogous to OCCT `XCAFView` - defines a named camera position
/// for drawings, renderings, or saved viewpoints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct View {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name (e.g., "Front", "Top", "Isometric").
    pub name: String,
    /// Camera position in world coordinates.
    pub camera_position: DVec3,
    /// Target point (look-at point) in world coordinates.
    pub target: DVec3,
    /// Up vector for camera orientation.
    pub up_vector: DVec3,
    /// Projection type (orthographic or perspective).
    pub projection: ViewProjection,
    /// Field of view for perspective projection (in degrees).
    pub fov: f64,
    /// View width for orthographic projection.
    pub view_width: f64,
    /// View height for orthographic projection.
    pub view_height: f64,
    /// Clipping planes.
    pub clipping: ViewClipping,
    /// Whether this view is user-created (vs. standard view).
    pub custom: bool,
    /// Visibility of geometry in this view (layer-like behavior).
    pub visible_shapes: Vec<usize>,
    /// Hidden shapes in this view.
    pub hidden_shapes: Vec<usize>,
}

impl View {
    /// Create a new view with the given ID and name.
    pub fn new(id: u64, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            camera_position: DVec3::new(0.0, 0.0, 100.0),
            target: DVec3::ZERO,
            up_vector: DVec3::Y,
            projection: ViewProjection::Orthographic,
            fov: 45.0,
            view_width: 100.0,
            view_height: 100.0,
            clipping: ViewClipping::default(),
            custom: false,
            visible_shapes: Vec::new(),
            hidden_shapes: Vec::new(),
        }
    }

    /// Create a standard front view.
    pub fn front(id: u64) -> Self {
        Self {
            id,
            name: "Front".to_string(),
            camera_position: DVec3::new(0.0, 0.0, 100.0),
            target: DVec3::ZERO,
            up_vector: DVec3::Y,
            projection: ViewProjection::Orthographic,
            ..Self::new(id, "Front")
        }
    }

    /// Create a standard top view.
    pub fn top(id: u64) -> Self {
        Self {
            id,
            name: "Top".to_string(),
            camera_position: DVec3::new(0.0, 100.0, 0.0),
            target: DVec3::ZERO,
            up_vector: DVec3::NEG_Z,
            projection: ViewProjection::Orthographic,
            ..Self::new(id, "Top")
        }
    }

    /// Create a standard right view.
    pub fn right(id: u64) -> Self {
        Self {
            id,
            name: "Right".to_string(),
            camera_position: DVec3::new(100.0, 0.0, 0.0),
            target: DVec3::ZERO,
            up_vector: DVec3::Z,
            projection: ViewProjection::Orthographic,
            ..Self::new(id, "Right")
        }
    }

    /// Create a standard isometric view.
    pub fn isometric(id: u64) -> Self {
        Self {
            id,
            name: "Isometric".to_string(),
            camera_position: DVec3::new(100.0, 100.0, 100.0),
            target: DVec3::ZERO,
            up_vector: DVec3::Z,
            projection: ViewProjection::Orthographic,
            ..Self::new(id, "Isometric")
        }
    }

    /// Set the camera position.
    pub fn with_position(mut self, position: DVec3) -> Self {
        self.camera_position = position;
        self
    }

    /// Set the target point.
    pub fn with_target(mut self, target: DVec3) -> Self {
        self.target = target;
        self
    }

    /// Set the up vector.
    pub fn with_up(mut self, up: DVec3) -> Self {
        self.up_vector = up;
        self
    }

    /// Set the projection type.
    pub fn with_projection(mut self, projection: ViewProjection) -> Self {
        self.projection = projection;
        self
    }

    /// Set the field of view (for perspective).
    pub fn with_fov(mut self, fov: f64) -> Self {
        self.fov = fov;
        self
    }

    /// Set the view dimensions (for orthographic).
    pub fn with_view_size(mut self, width: f64, height: f64) -> Self {
        self.view_width = width;
        self.view_height = height;
        self
    }

    /// Get the view direction (camera to target).
    pub fn view_direction(&self) -> DVec3 {
        (self.target - self.camera_position).normalize_or_zero()
    }

    /// Get the distance from camera to target.
    pub fn distance(&self) -> f64 {
        (self.target - self.camera_position).length()
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// ANNOTATION TYPES (legacy PMI types)
// ═══════════════════════════════════════════════════════════════════════════════

/// Type of annotation note.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum NoteType {
    /// Plain text note.
    Text,
    /// Symbol (welding, surface finish, etc.).
    Symbol,
    /// Dimension annotation.
    Dimension,
    /// Surface texture specification.
    SurfaceTexture,
    /// Welding symbol.
    WeldSymbol,
    /// Balloon/bubble callout.
    Balloon,
    /// Leader line with arrow.
    Leader,
}

/// Arrow head type for leader lines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ArrowType {
    /// Filled triangle arrow.
    FilledArrow,
    /// Open triangle arrow.
    OpenArrow,
    /// Dot marker.
    Dot,
    /// Origin marker (circle).
    Origin,
    /// No arrow head.
    NoArrow,
}

/// Welding type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum WeldType {
    /// Fillet weld.
    Fillet,
    /// Groove weld (V, U, J, etc.).
    Groove,
    /// Plug weld.
    Plug,
    /// Slot weld.
    Slot,
    /// Spot weld.
    Spot,
    /// Projection weld.
    Projection,
    /// Seam weld.
    Seam,
}

/// General annotation note.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationNote {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name/label.
    pub name: String,
    /// Text content of the note.
    pub content: String,
    /// Type classification.
    pub note_type: NoteType,
    /// 3D position in model space.
    pub position: DVec3,
    /// Indices of geometry this note is attached to.
    pub attached_geometry: Vec<usize>,
    /// Visibility flag.
    pub visibility: bool,
}

impl AnnotationNote {
    /// Create a new annotation note.
    pub fn new(id: u64, name: impl Into<String>, note_type: NoteType, position: DVec3) -> Self {
        Self {
            id,
            name: name.into(),
            content: String::new(),
            note_type,
            position,
            attached_geometry: Vec::new(),
            visibility: true,
        }
    }
}

/// Text annotation with formatting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TextAnnotation {
    /// Unique identifier.
    pub id: u64,
    /// Text content.
    pub text: String,
    /// Font family name.
    pub font: String,
    /// Text height in model units.
    pub height: f64,
    /// Position of text origin.
    pub position: DVec3,
    /// Direction vector for text baseline.
    pub direction: DVec3,
    /// Up vector for text orientation.
    pub up_vector: DVec3,
}

impl TextAnnotation {
    /// Create a new text annotation.
    pub fn new(id: u64, text: impl Into<String>, position: DVec3) -> Self {
        Self {
            id,
            text: text.into(),
            font: "Arial".to_string(),
            height: 5.0,
            position,
            direction: DVec3::X,
            up_vector: DVec3::Z,
        }
    }
}

/// Leader line with arrow.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LeaderLine {
    /// Unique identifier.
    pub id: u64,
    /// Polyline points defining the leader path.
    pub points: Vec<DVec3>,
    /// Arrow head type.
    pub arrow_type: ArrowType,
    /// Arrow head size.
    pub arrow_size: f64,
    /// Index of geometry the leader points to.
    pub attached_geometry: usize,
}

impl LeaderLine {
    /// Create a new leader line.
    pub fn new(id: u64, attached_geometry: usize) -> Self {
        Self {
            id,
            points: Vec::new(),
            arrow_type: ArrowType::FilledArrow,
            arrow_size: 3.0,
            attached_geometry,
        }
    }

    /// Add a point to the leader path.
    pub fn add_point(&mut self, point: DVec3) {
        self.points.push(point);
    }

    /// Get the start point (first in path).
    pub fn start(&self) -> Option<&DVec3> {
        self.points.first()
    }

    /// Get the end point (last in path, where arrow is).
    pub fn end(&self) -> Option<&DVec3> {
        self.points.last()
    }
}

/// Surface texture symbol (roughness specification).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceTextureSymbol {
    /// Unique identifier.
    pub id: u64,
    /// Ra (arithmetic average) roughness value.
    pub roughness_ra: Option<f64>,
    /// Rz (average maximum height) roughness value.
    pub roughness_rz: Option<f64>,
    /// Machining allowance.
    pub machining_allowance: Option<f64>,
    /// Lay direction (e.g., "=", "⊥", "X", "M", "C", "R").
    pub lay_direction: Option<String>,
    /// Index of face this symbol applies to.
    pub attached_face: usize,
}

impl SurfaceTextureSymbol {
    /// Create a new surface texture symbol.
    pub fn new(id: u64, attached_face: usize) -> Self {
        Self {
            id,
            roughness_ra: None,
            roughness_rz: None,
            machining_allowance: None,
            lay_direction: None,
            attached_face,
        }
    }

    /// Create with Ra roughness value.
    pub fn with_ra(id: u64, attached_face: usize, ra: f64) -> Self {
        Self {
            roughness_ra: Some(ra),
            ..Self::new(id, attached_face)
        }
    }
}

/// Welding symbol with specifications.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WeldSymbol {
    /// Unique identifier.
    pub id: u64,
    /// Type of weld.
    pub weld_type: WeldType,
    /// Weld size (leg length for fillet, depth for groove).
    pub size: f64,
    /// Weld length (if partial).
    pub length: Option<f64>,
    /// Pitch (spacing) for intermittent welds.
    pub pitch: Option<f64>,
    /// Whether weld is on arrow side.
    pub arrow_side: bool,
    /// Whether weld is on other side.
    pub other_side: bool,
}

impl WeldSymbol {
    /// Create a new weld symbol.
    pub fn new(id: u64, weld_type: WeldType, size: f64) -> Self {
        Self {
            id,
            weld_type,
            size,
            length: None,
            pitch: None,
            arrow_side: true,
            other_side: false,
        }
    }
}

/// Balloon annotation (numbered callout).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BalloonAnnotation {
    /// Unique identifier.
    pub id: u64,
    /// Balloon number/label.
    pub number: u32,
    /// Position of balloon center.
    pub position: DVec3,
    /// Index of geometry the balloon references.
    pub attached_geometry: usize,
    /// Optional leader line index.
    pub leader_id: Option<u64>,
}

impl BalloonAnnotation {
    /// Create a new balloon annotation.
    pub fn new(id: u64, number: u32, position: DVec3, attached_geometry: usize) -> Self {
        Self {
            id,
            number,
            position,
            attached_geometry,
            leader_id: None,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// UNIFIED ANNOTATION ENTITY
// ═══════════════════════════════════════════════════════════════════════════════

/// Target for an annotation, linking it to geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationTarget {
    /// The target geometry reference.
    pub target: NoteTarget,
    /// Optional point on the target for precise attachment.
    pub point_on_target: Option<DVec3>,
    /// Optional parameter value (e.g., for curves/surfaces).
    pub parameter: Option<f64>,
}

impl AnnotationTarget {
    /// Create a target pointing to a shape.
    pub fn shape(shape_index: usize) -> Self {
        Self {
            target: NoteTarget::Shape { shape_index },
            point_on_target: None,
            parameter: None,
        }
    }

    /// Create a target pointing to a face.
    pub fn face(shape_index: usize, face_index: usize) -> Self {
        Self {
            target: NoteTarget::Face {
                shape_index,
                face_index,
            },
            point_on_target: None,
            parameter: None,
        }
    }

    /// Create a target pointing to a specific point on a face.
    pub fn face_point(shape_index: usize, face_index: usize, point: DVec3) -> Self {
        Self {
            target: NoteTarget::Face {
                shape_index,
                face_index,
            },
            point_on_target: Some(point),
            parameter: None,
        }
    }

    /// Create a target pointing to an edge.
    pub fn edge(shape_index: usize, edge_index: usize) -> Self {
        Self {
            target: NoteTarget::Edge {
                shape_index,
                edge_index,
            },
            point_on_target: None,
            parameter: None,
        }
    }

    /// Create a target pointing to a specific point on an edge.
    pub fn edge_point(shape_index: usize, edge_index: usize, point: DVec3, param: f64) -> Self {
        Self {
            target: NoteTarget::Edge {
                shape_index,
                edge_index,
            },
            point_on_target: Some(point),
            parameter: Some(param),
        }
    }
}

/// Visual style for an annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationStyle {
    /// Color as RGB (0-1 range).
    pub color: [f64; 3],
    /// Line width for lines and curves.
    pub line_width: f64,
    /// Text height for text annotations.
    pub text_height: f64,
    /// Font family name.
    pub font: String,
    /// Whether the annotation is displayed with a leader.
    pub show_leader: bool,
    /// Arrow size for leader lines.
    pub arrow_size: f64,
}

impl Default for AnnotationStyle {
    fn default() -> Self {
        Self {
            color: [0.0, 0.0, 0.0], // Black
            line_width: 0.35,
            text_height: 3.5,
            font: "Arial".to_string(),
            show_leader: true,
            arrow_size: 3.0,
        }
    }
}

impl AnnotationStyle {
    /// Create a new annotation style.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the color.
    pub fn with_color(mut self, r: f64, g: f64, b: f64) -> Self {
        self.color = [r, g, b];
        self
    }

    /// Set the line width.
    pub fn with_line_width(mut self, width: f64) -> Self {
        self.line_width = width;
        self
    }

    /// Set the text height.
    pub fn with_text_height(mut self, height: f64) -> Self {
        self.text_height = height;
        self
    }

    /// Set the font.
    pub fn with_font(mut self, font: impl Into<String>) -> Self {
        self.font = font.into();
        self
    }
}

/// Annotation type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AnnotationKind {
    /// Dimensional annotation (linear, angular, radial, etc.).
    Dimension,
    /// Geometric tolerance (GD&T feature control frame).
    GeometricTolerance,
    /// Surface finish/texture specification.
    SurfaceFinish,
    /// Welding symbol.
    Weld,
    /// Datum reference.
    Datum,
    /// Text note or comment.
    Note,
    /// Balloon/bubble callout.
    Balloon,
    /// Section or detail marker.
    SectionMarker,
    /// Custom/symbolic annotation.
    Custom,
}

/// A unified annotation entity combining notes, dimensions, tolerances.
///
/// This is the primary annotation type for PMI (Product Manufacturing Information)
/// that can represent various kinds of annotations attached to geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Annotation {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name/label.
    pub name: String,
    /// Type of annotation.
    pub kind: AnnotationKind,
    /// Text content (may be multi-line).
    pub text: String,
    /// Targets this annotation is attached to.
    pub targets: Vec<AnnotationTarget>,
    /// 3D position of the annotation in model space.
    pub position: DVec3,
    /// Direction for the annotation (e.g., text baseline direction).
    pub direction: DVec3,
    /// Up vector for orientation.
    pub up_vector: DVec3,
    /// Visual style.
    pub style: AnnotationStyle,
    /// Associated leader line ID (if any).
    pub leader_id: Option<u64>,
    /// Associated dimension curve ID (if any).
    pub dimension_curve_id: Option<u64>,
    /// Visibility flag.
    pub visibility: bool,
    /// View index for view-specific visibility (None = all views).
    pub view_index: Option<usize>,
    /// Additional semantic data (JSON-serializable).
    pub semantic_data: Option<serde_json::Value>,
}

impl Annotation {
    /// Create a new annotation with the given ID, name, and kind.
    pub fn new(id: u64, name: impl Into<String>, kind: AnnotationKind) -> Self {
        Self {
            id,
            name: name.into(),
            kind,
            text: String::new(),
            targets: Vec::new(),
            position: DVec3::ZERO,
            direction: DVec3::X,
            up_vector: DVec3::Z,
            style: AnnotationStyle::default(),
            leader_id: None,
            dimension_curve_id: None,
            visibility: true,
            view_index: None,
            semantic_data: None,
        }
    }

    /// Create a dimensional annotation.
    pub fn dimension(id: u64, name: impl Into<String>, value: f64) -> Self {
        Self {
            text: format!("{:.3}", value),
            ..Self::new(id, name, AnnotationKind::Dimension)
        }
    }

    /// Create a note annotation.
    pub fn note(id: u64, name: impl Into<String>, text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            ..Self::new(id, name, AnnotationKind::Note)
        }
    }

    /// Create a datum annotation.
    pub fn datum(id: u64, name: impl Into<String>, identifier: impl Into<String>) -> Self {
        Self {
            text: identifier.into(),
            ..Self::new(id, name, AnnotationKind::Datum)
        }
    }

    /// Add a target to this annotation.
    pub fn add_target(&mut self, target: AnnotationTarget) {
        self.targets.push(target);
    }

    /// Set the position.
    pub fn with_position(mut self, position: DVec3) -> Self {
        self.position = position;
        self
    }

    /// Set the direction.
    pub fn with_direction(mut self, direction: DVec3) -> Self {
        self.direction = direction;
        self
    }

    /// Set the visibility.
    pub fn with_visibility(mut self, visibility: bool) -> Self {
        self.visibility = visibility;
        self
    }

    /// Set the view index for view-specific visibility.
    pub fn with_view(mut self, view_index: usize) -> Self {
        self.view_index = Some(view_index);
        self
    }

    /// Check if this annotation is visible in the given view.
    pub fn is_visible_in_view(&self, view_index: Option<usize>) -> bool {
        if !self.visibility {
            return false;
        }
        match (self.view_index, view_index) {
            (None, _) => true, // No view restriction
            (Some(av), Some(vv)) => av == vv, // Match specific view
            _ => false,
        }
    }
}

/// Storage container for all annotation data.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AnnotationStore {
    // ── XCAFNoteObjects-style notes ─────────────────────────────────────────────
    /// XCAF-style notes (comment, warning, info, requirement).
    pub xc_notes: Vec<Note>,
    /// XCAF-style views (camera positions).
    pub views: Vec<View>,
    /// Unified annotations (dimensions, tolerances, notes, etc.).
    pub annotations: Vec<Annotation>,

    // ── Legacy PMI annotation types ──────────────────────────────────────────────
    /// All annotation notes (legacy).
    pub notes: Vec<AnnotationNote>,
    /// All text annotations.
    pub text_annotations: Vec<TextAnnotation>,
    /// All leader lines.
    pub leader_lines: Vec<LeaderLine>,
    /// All surface texture symbols.
    pub surface_textures: Vec<SurfaceTextureSymbol>,
    /// All weld symbols.
    pub weld_symbols: Vec<WeldSymbol>,
    /// All balloon annotations.
    pub balloons: Vec<BalloonAnnotation>,
}

impl AnnotationStore {
    /// Create an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an annotation note.
    pub fn add_note(&mut self, note: AnnotationNote) -> usize {
        let idx = self.notes.len();
        self.notes.push(note);
        idx
    }

    /// Add a text annotation.
    pub fn add_text(&mut self, text: TextAnnotation) -> usize {
        let idx = self.text_annotations.len();
        self.text_annotations.push(text);
        idx
    }

    /// Add a leader line.
    pub fn add_leader(&mut self, leader: LeaderLine) -> usize {
        let idx = self.leader_lines.len();
        self.leader_lines.push(leader);
        idx
    }

    /// Add a surface texture symbol.
    pub fn add_surface_texture(&mut self, symbol: SurfaceTextureSymbol) -> usize {
        let idx = self.surface_textures.len();
        self.surface_textures.push(symbol);
        idx
    }

    /// Add a weld symbol.
    pub fn add_weld(&mut self, weld: WeldSymbol) -> usize {
        let idx = self.weld_symbols.len();
        self.weld_symbols.push(weld);
        idx
    }

    /// Add a balloon annotation.
    pub fn add_balloon(&mut self, balloon: BalloonAnnotation) -> usize {
        let idx = self.balloons.len();
        self.balloons.push(balloon);
        idx
    }

    // ── XCAFNoteObjects-style methods ─────────────────────────────────────────

    /// Add an XCAF-style note.
    pub fn add_xc_note(&mut self, note: Note) -> usize {
        let idx = self.xc_notes.len();
        self.xc_notes.push(note);
        idx
    }

    /// Add a view definition.
    pub fn add_view(&mut self, view: View) -> usize {
        let idx = self.views.len();
        self.views.push(view);
        idx
    }

    /// Add a unified annotation.
    pub fn add_annotation(&mut self, annotation: Annotation) -> usize {
        let idx = self.annotations.len();
        self.annotations.push(annotation);
        idx
    }

    /// Get a view by ID.
    pub fn get_view(&self, id: u64) -> Option<&View> {
        self.views.iter().find(|v| v.id == id)
    }

    /// Get a view by name.
    pub fn get_view_by_name(&self, name: &str) -> Option<&View> {
        self.views.iter().find(|v| v.name == name)
    }

    /// Get all notes for a specific shape.
    pub fn notes_for_shape(&self, shape_index: usize) -> Vec<&Note> {
        self.xc_notes
            .iter()
            .filter(|n| {
                n.links.iter().any(|link| match link.target {
                    NoteTarget::Shape { shape_index: si } => si == shape_index,
                    NoteTarget::Face { shape_index: si, .. } => si == shape_index,
                    NoteTarget::Edge { shape_index: si, .. } => si == shape_index,
                    NoteTarget::Vertex { shape_index: si, .. } => si == shape_index,
                    _ => false,
                })
            })
            .collect()
    }

    /// Get all annotations visible in a specific view.
    pub fn annotations_for_view(&self, view_index: Option<usize>) -> Vec<&Annotation> {
        self.annotations
            .iter()
            .filter(|a| a.is_visible_in_view(view_index))
            .collect()
    }

    /// Get all annotations for a specific target.
    pub fn annotations_for_target(&self, shape_index: usize, face_index: Option<usize>) -> Vec<&Annotation> {
        self.annotations
            .iter()
            .filter(|a| {
                a.targets.iter().any(|t| match (&t.target, face_index) {
                    (NoteTarget::Shape { shape_index: si }, _) => *si == shape_index,
                    (NoteTarget::Face { shape_index: si, face_index: fi }, Some(fi2)) => {
                        *si == shape_index && *fi == fi2
                    }
                    (NoteTarget::Face { shape_index: si, .. }, None) => *si == shape_index,
                    (NoteTarget::Edge { shape_index: si, .. }, _) => *si == shape_index,
                    (NoteTarget::Vertex { shape_index: si, .. }, _) => *si == shape_index,
                    _ => false,
                })
            })
            .collect()
    }

    /// Create standard views (front, top, right, isometric).
    pub fn add_standard_views(&mut self, start_id: u64) {
        self.views.push(View::front(start_id));
        self.views.push(View::top(start_id + 1));
        self.views.push(View::right(start_id + 2));
        self.views.push(View::isometric(start_id + 3));
    }

    /// Get all notes attached to specific geometry.
    pub fn notes_for_geometry(&self, geom_idx: usize) -> Vec<&AnnotationNote> {
        self.notes
            .iter()
            .filter(|n| n.attached_geometry.contains(&geom_idx))
            .collect()
    }

    /// Total count of all annotations.
    pub fn total_count(&self) -> usize {
        self.xc_notes.len()
            + self.views.len()
            + self.annotations.len()
            + self.notes.len()
            + self.text_annotations.len()
            + self.leader_lines.len()
            + self.surface_textures.len()
            + self.weld_symbols.len()
            + self.balloons.len()
    }

    /// Count of XCAF-style notes.
    pub fn xc_notes_count(&self) -> usize {
        self.xc_notes.len()
    }

    /// Count of views.
    pub fn views_count(&self) -> usize {
        self.views.len()
    }

    /// Count of unified annotations.
    pub fn annotations_count(&self) -> usize {
        self.annotations.len()
    }

    /// Clear all annotations.
    pub fn clear(&mut self) {
        self.xc_notes.clear();
        self.views.clear();
        self.annotations.clear();
        self.notes.clear();
        self.text_annotations.clear();
        self.leader_lines.clear();
        self.surface_textures.clear();
        self.weld_symbols.clear();
        self.balloons.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_note_creation() {
        let note = AnnotationNote::new(1, "Note1", NoteType::Text, DVec3::ZERO);
        assert_eq!(note.id, 1);
        assert_eq!(note.name, "Note1");
        assert_eq!(note.note_type, NoteType::Text);
        assert!(note.visibility);
    }

    #[test]
    fn text_annotation_creation() {
        let text = TextAnnotation::new(1, "Hello World", DVec3::new(10.0, 0.0, 0.0));
        assert_eq!(text.text, "Hello World");
        assert_eq!(text.font, "Arial");
        assert_eq!(text.position, DVec3::new(10.0, 0.0, 0.0));
    }

    #[test]
    fn leader_line_points() {
        let mut leader = LeaderLine::new(1, 0);
        leader.add_point(DVec3::ZERO);
        leader.add_point(DVec3::new(10.0, 10.0, 0.0));
        assert_eq!(leader.points.len(), 2);
        assert_eq!(*leader.start().unwrap(), DVec3::ZERO);
        assert_eq!(*leader.end().unwrap(), DVec3::new(10.0, 10.0, 0.0));
    }

    #[test]
    fn surface_texture_symbol() {
        let symbol = SurfaceTextureSymbol::with_ra(1, 0, 3.2);
        assert_eq!(symbol.roughness_ra, Some(3.2));
        assert!(symbol.roughness_rz.is_none());
    }

    #[test]
    fn weld_symbol_creation() {
        let weld = WeldSymbol::new(1, WeldType::Fillet, 6.0);
        assert_eq!(weld.weld_type, WeldType::Fillet);
        assert_eq!(weld.size, 6.0);
        assert!(weld.arrow_side);
        assert!(!weld.other_side);
    }

    #[test]
    fn balloon_annotation() {
        let balloon = BalloonAnnotation::new(1, 42, DVec3::new(5.0, 5.0, 0.0), 0);
        assert_eq!(balloon.number, 42);
        assert_eq!(balloon.position, DVec3::new(5.0, 5.0, 0.0));
    }

    #[test]
    fn annotation_store_operations() {
        let mut store = AnnotationStore::new();

        store.add_note(AnnotationNote::new(1, "Note", NoteType::Text, DVec3::ZERO));
        store.add_text(TextAnnotation::new(1, "Text", DVec3::ZERO));
        store.add_leader(LeaderLine::new(1, 0));
        store.add_surface_texture(SurfaceTextureSymbol::new(1, 0));
        store.add_weld(WeldSymbol::new(1, WeldType::Fillet, 5.0));
        store.add_balloon(BalloonAnnotation::new(1, 1, DVec3::ZERO, 0));

        assert_eq!(store.total_count(), 6);
    }

    // ── Tests for XCAFNoteObjects-style notes ───────────────────────────────────

    #[test]
    fn note_creation() {
        let note = Note::new(1, "TestNote", "This is a test note")
            .with_category(NoteCategory::Warning)
            .with_author("John Doe")
            .with_timestamp("2024-01-15T10:30:00Z")
            .with_position(DVec3::new(10.0, 20.0, 30.0));

        assert_eq!(note.id, 1);
        assert_eq!(note.name, "TestNote");
        assert_eq!(note.text, "This is a test note");
        assert_eq!(note.category, NoteCategory::Warning);
        assert_eq!(note.author.as_deref(), Some("John Doe"));
        assert!(note.visibility);
    }

    #[test]
    fn note_link_creation() {
        let link_shape = NoteLink::to_shape(0);
        assert!(matches!(link_shape.target, NoteTarget::Shape { shape_index: 0 }));

        let link_face = NoteLink::to_face(0, 5);
        assert!(matches!(link_face.target, NoteTarget::Face { shape_index: 0, face_index: 5 }));

        let link_edge = NoteLink::to_edge(1, 3);
        assert!(matches!(link_edge.target, NoteTarget::Edge { shape_index: 1, edge_index: 3 }));

        let link_point = NoteLink::to_point(DVec3::new(1.0, 2.0, 3.0));
        assert!(matches!(link_point.target, NoteTarget::Point { .. }));
    }

    #[test]
    fn note_with_links() {
        let mut note = Note::new(1, "LinkedNote", "Note with geometry links");
        note.add_link(NoteLink::to_shape(0));
        note.add_link(NoteLink::to_face(0, 1));

        assert_eq!(note.links.len(), 2);
    }

    // ── Tests for View definitions ───────────────────────────────────────────────

    #[test]
    fn view_creation() {
        let view = View::new(1, "CustomView")
            .with_position(DVec3::new(50.0, 50.0, 50.0))
            .with_target(DVec3::ZERO)
            .with_up(DVec3::Z)
            .with_projection(ViewProjection::Perspective)
            .with_fov(60.0);

        assert_eq!(view.id, 1);
        assert_eq!(view.name, "CustomView");
        assert_eq!(view.camera_position, DVec3::new(50.0, 50.0, 50.0));
        assert_eq!(view.projection, ViewProjection::Perspective);
        assert_eq!(view.fov, 60.0);
    }

    #[test]
    fn standard_views() {
        let front = View::front(1);
        assert_eq!(front.name, "Front");
        assert_eq!(front.projection, ViewProjection::Orthographic);

        let top = View::top(2);
        assert_eq!(top.name, "Top");

        let right = View::right(3);
        assert_eq!(right.name, "Right");

        let iso = View::isometric(4);
        assert_eq!(iso.name, "Isometric");
    }

    #[test]
    fn view_direction_and_distance() {
        let view = View::new(1, "Test")
            .with_position(DVec3::new(0.0, 0.0, 100.0))
            .with_target(DVec3::ZERO);

        let dir = view.view_direction();
        assert!((dir - DVec3::NEG_Z).length() < 1e-10);

        let dist = view.distance();
        assert!((dist - 100.0).abs() < 1e-10);
    }

    #[test]
    fn view_clipping() {
        let clipping = ViewClipping::new(1.0, 1000.0);
        assert!((clipping.near - 1.0).abs() < 1e-10);
        assert!((clipping.far - 1000.0).abs() < 1e-10);
        assert!(clipping.front_enabled);
        assert!(clipping.back_enabled);
    }

    // ── Tests for unified Annotation ─────────────────────────────────────────────

    #[test]
    fn annotation_creation() {
        let annotation = Annotation::note(1, "TestNote", "This is a note")
            .with_position(DVec3::new(10.0, 20.0, 0.0))
            .with_visibility(true);

        assert_eq!(annotation.id, 1);
        assert_eq!(annotation.name, "TestNote");
        assert_eq!(annotation.text, "This is a note");
        assert_eq!(annotation.kind, AnnotationKind::Note);
        assert!(annotation.visibility);
    }

    #[test]
    fn dimension_annotation() {
        let dim = Annotation::dimension(1, "Length", 25.4);
        assert_eq!(dim.kind, AnnotationKind::Dimension);
        assert!(dim.text.contains("25.4"));
    }

    #[test]
    fn datum_annotation() {
        let datum = Annotation::datum(1, "DatumA", "A");
        assert_eq!(datum.kind, AnnotationKind::Datum);
        assert_eq!(datum.text, "A");
    }

    #[test]
    fn annotation_with_target() {
        let mut annotation = Annotation::note(1, "Note", "Test");
        annotation.add_target(AnnotationTarget::face(0, 1));
        annotation.add_target(AnnotationTarget::edge(0, 2));

        assert_eq!(annotation.targets.len(), 2);
    }

    #[test]
    fn annotation_view_visibility() {
        let annotation = Annotation::note(1, "Note", "Test").with_view(2);

        // Should be visible in view 2
        assert!(annotation.is_visible_in_view(Some(2)));

        // Should not be visible in other views
        assert!(!annotation.is_visible_in_view(Some(1)));

        // If no view restriction, visible in all views
        let unrestricted = Annotation::note(2, "Unrestricted", "Test");
        assert!(unrestricted.is_visible_in_view(Some(1)));
        assert!(unrestricted.is_visible_in_view(None));
    }

    #[test]
    fn annotation_style() {
        let style = AnnotationStyle::new()
            .with_color(1.0, 0.0, 0.0)
            .with_line_width(0.5)
            .with_text_height(5.0)
            .with_font("Helvetica");

        assert_eq!(style.color, [1.0, 0.0, 0.0]);
        assert!((style.line_width - 0.5).abs() < 1e-10);
        assert!((style.text_height - 5.0).abs() < 1e-10);
        assert_eq!(style.font, "Helvetica");
    }

    // ── Tests for AnnotationStore with new types ─────────────────────────────────

    #[test]
    fn store_xc_notes() {
        let mut store = AnnotationStore::new();

        store.add_xc_note(Note::new(1, "Note1", "First note"));
        store.add_xc_note(Note::new(2, "Note2", "Second note"));

        assert_eq!(store.xc_notes_count(), 2);
    }

    #[test]
    fn store_views() {
        let mut store = AnnotationStore::new();

        store.add_standard_views(1);
        assert_eq!(store.views_count(), 4);

        // Find view by name
        let front = store.get_view_by_name("Front");
        assert!(front.is_some());
        assert_eq!(front.unwrap().id, 1);

        // Find view by ID
        let top = store.get_view(2);
        assert!(top.is_some());
        assert_eq!(top.unwrap().name, "Top");
    }

    #[test]
    fn store_annotations() {
        let mut store = AnnotationStore::new();

        store.add_annotation(Annotation::dimension(1, "Dim1", 10.0));
        store.add_annotation(Annotation::note(2, "Note1", "Test note"));

        assert_eq!(store.annotations_count(), 2);
    }

    #[test]
    fn store_notes_for_shape() {
        let mut store = AnnotationStore::new();

        let mut note1 = Note::new(1, "ShapeNote", "Note on shape 0");
        note1.add_link(NoteLink::to_shape(0));
        store.add_xc_note(note1);

        let mut note2 = Note::new(2, "FaceNote", "Note on face 1 of shape 0");
        note2.add_link(NoteLink::to_face(0, 1));
        store.add_xc_note(note2);

        let mut note3 = Note::new(3, "OtherNote", "Note on shape 1");
        note3.add_link(NoteLink::to_shape(1));
        store.add_xc_note(note3);

        // Should find notes for shape 0
        let notes_for_0 = store.notes_for_shape(0);
        assert_eq!(notes_for_0.len(), 2);

        // Should find notes for shape 1
        let notes_for_1 = store.notes_for_shape(1);
        assert_eq!(notes_for_1.len(), 1);
    }

    #[test]
    fn store_annotations_for_view() {
        let mut store = AnnotationStore::new();

        let ann1 = Annotation::note(1, "Note1", "Unrestricted").with_visibility(true);
        let ann2 = Annotation::note(2, "Note2", "View2Only").with_view(2);
        let ann3 = Annotation::note(3, "Note3", "View1Only").with_view(1).with_visibility(false);

        store.add_annotation(ann1);
        store.add_annotation(ann2);
        store.add_annotation(ann3);

        // View 1 should see unrestricted + view 1 (but ann3 is hidden)
        let for_view1 = store.annotations_for_view(Some(1));
        assert_eq!(for_view1.len(), 1); // Only ann1

        // View 2 should see unrestricted + view 2
        let for_view2 = store.annotations_for_view(Some(2));
        assert_eq!(for_view2.len(), 2); // ann1 and ann2

        // No view should see only unrestricted
        let for_none = store.annotations_for_view(None);
        assert_eq!(for_none.len(), 1); // Only ann1
    }

    #[test]
    fn store_total_count_includes_new_types() {
        let mut store = AnnotationStore::new();

        store.add_xc_note(Note::new(1, "Note", "Test"));
        store.add_view(View::front(1));
        store.add_annotation(Annotation::note(1, "Ann", "Test"));
        store.add_note(AnnotationNote::new(1, "Legacy", NoteType::Text, DVec3::ZERO));

        assert_eq!(store.xc_notes_count(), 1);
        assert_eq!(store.views_count(), 1);
        assert_eq!(store.annotations_count(), 1);
        assert_eq!(store.notes.len(), 1);
        assert_eq!(store.total_count(), 4);
    }

    #[test]
    fn store_clear() {
        let mut store = AnnotationStore::new();

        store.add_xc_note(Note::new(1, "Note", "Test"));
        store.add_view(View::front(1));
        store.add_annotation(Annotation::note(1, "Ann", "Test"));

        store.clear();

        assert_eq!(store.xc_notes_count(), 0);
        assert_eq!(store.views_count(), 0);
        assert_eq!(store.annotations_count(), 0);
    }
}
