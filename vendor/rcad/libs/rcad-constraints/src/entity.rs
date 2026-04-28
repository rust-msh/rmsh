//! Geometric entity types and their parameter layout.
//!
//! All entity parameters are stored in a flat `Vec<f64>` owned by [`Sketch`].
//! Each entity records its `param_start` offset into that vector.
//!
//! | Kind   | Params (relative offsets)                          |
//! |--------|----------------------------------------------------|
//! | Point  | [0]=x, [1]=y                                       |
//! | Line   | [0]=x1, [1]=y1, [2]=x2, [3]=y2                    |
//! | Circle | [0]=cx, [1]=cy, [2]=r                              |
//! | Arc    | [0]=cx, [1]=cy, [2]=r, [3]=start_angle, [4]=end_angle |

/// Index into the entity list of a [`Sketch`].
pub type EntityId = usize;

/// The geometric kind of an entity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntityKind {
    Point,
    Line,
    Circle,
    Arc,
}

impl EntityKind {
    /// Number of scalar parameters this entity type uses.
    pub fn param_count(self) -> usize {
        match self {
            EntityKind::Point => 2,
            EntityKind::Line => 4,
            EntityKind::Circle => 3,
            EntityKind::Arc => 5,
        }
    }
}

/// Metadata for a single entity stored in a [`Sketch`].
#[derive(Debug, Clone)]
pub struct Entity {
    pub kind: EntityKind,
    /// Index of the first parameter in the sketch's flat `params` vector.
    pub param_start: usize,
}

impl Entity {
    pub fn new(kind: EntityKind, param_start: usize) -> Self {
        Self { kind, param_start }
    }

    /// Absolute parameter index for a relative offset within this entity.
    #[inline]
    pub fn param(&self, offset: usize) -> usize {
        self.param_start + offset
    }
}

/// A reference to a specific 2D point within an entity.
///
/// Used by constraints that operate on points (e.g. Coincident, Fixed).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointRef {
    /// A `Point` entity — the entity itself is the point.
    Point(EntityId),
    /// The start endpoint of a `Line` entity (params [0,1]).
    LineStart(EntityId),
    /// The end endpoint of a `Line` entity (params [2,3]).
    LineEnd(EntityId),
    /// The center of a `Circle` or `Arc` entity (params [0,1]).
    Center(EntityId),
}

impl PointRef {
    /// Return the (x_param_idx, y_param_idx) absolute indices into the sketch
    /// parameter vector.
    pub fn param_indices(&self, entities: &[Entity]) -> (usize, usize) {
        match *self {
            PointRef::Point(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
            PointRef::LineStart(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
            PointRef::LineEnd(id) => {
                let e = &entities[id];
                (e.param(2), e.param(3))
            }
            PointRef::Center(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1))
            }
        }
    }
}
