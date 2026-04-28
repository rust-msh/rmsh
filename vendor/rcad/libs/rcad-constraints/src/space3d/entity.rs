//! 3D geometric entity types and their parameter layout.
//!
//! All entity parameters are stored in a flat `Vec<f64>` owned by [`SpaceSketch`].
//! Each entity records its `param_start` offset into that vector.
//!
//! | Kind       | Params (relative offsets)                              |
//! |------------|--------------------------------------------------------|
//! | SpacePoint | [0]=x, [1]=y, [2]=z                                    |
//! | SpaceLine  | [0]=x1, [1]=y1, [2]=z1, [3]=x2, [4]=y2, [5]=z2       |
//! | Plane      | [0]=nx, [1]=ny, [2]=nz, [3]=d (distance from origin)  |
//! | Sphere     | [0]=cx, [1]=cy, [2]=cz, [3]=r                          |
//! | Cylinder   | [0..3]=axis point, [3..6]=axis direction, [6]=radius  |

/// Index into the entity list of a [`SpaceSketch`].
pub type SpaceEntityId = usize;

/// The geometric kind of a 3D entity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpaceEntityKind {
    SpacePoint,
    SpaceLine,
    Plane,
    Sphere,
    Cylinder,
}

impl SpaceEntityKind {
    /// Number of scalar parameters this entity type uses.
    pub fn param_count(self) -> usize {
        match self {
            SpaceEntityKind::SpacePoint => 3,
            SpaceEntityKind::SpaceLine => 6,
            SpaceEntityKind::Plane => 4,
            SpaceEntityKind::Sphere => 4,
            SpaceEntityKind::Cylinder => 7,
        }
    }
}

/// Metadata for a single 3D entity stored in a [`SpaceSketch`].
#[derive(Debug, Clone)]
pub struct SpaceEntity {
    pub kind: SpaceEntityKind,
    /// Index of the first parameter in the sketch's flat `params` vector.
    pub param_start: usize,
}

impl SpaceEntity {
    pub fn new(kind: SpaceEntityKind, param_start: usize) -> Self {
        Self { kind, param_start }
    }

    /// Absolute parameter index for a relative offset within this entity.
    #[inline]
    pub fn param(&self, offset: usize) -> usize {
        self.param_start + offset
    }
}

/// A reference to a specific 3D point within an entity.
///
/// Used by constraints that operate on points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpacePointRef {
    /// A `SpacePoint` entity — the entity itself is the point.
    Point(SpaceEntityId),
    /// The start endpoint of a `SpaceLine` entity (params [0,1,2]).
    LineStart(SpaceEntityId),
    /// The end endpoint of a `SpaceLine` entity (params [3,4,5]).
    LineEnd(SpaceEntityId),
    /// The center of a `Sphere` entity (params [0,1,2]).
    SphereCenter(SpaceEntityId),
}

impl SpacePointRef {
    /// Return the (x_param_idx, y_param_idx, z_param_idx) absolute indices
    /// into the sketch parameter vector.
    pub fn param_indices(&self, entities: &[SpaceEntity]) -> (usize, usize, usize) {
        match *self {
            SpacePointRef::Point(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::LineStart(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
            SpacePointRef::LineEnd(id) => {
                let e = &entities[id];
                (e.param(3), e.param(4), e.param(5))
            }
            SpacePointRef::SphereCenter(id) => {
                let e = &entities[id];
                (e.param(0), e.param(1), e.param(2))
            }
        }
    }
}

impl From<SpaceEntityId> for SpacePointRef {
    fn from(id: SpaceEntityId) -> Self {
        SpacePointRef::Point(id)
    }
}
