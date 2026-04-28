//! Dimension and tolerance object model.
//!
//! Analogous to OCCT `XCAFDimTolObjects` - provides data structures for
//! dimensional tolerances, geometric tolerances (GDT), and datum systems.

use serde::{Deserialize, Serialize};

/// Dimension type classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DimensionType {
    /// Linear dimension (distance between two points).
    Linear,
    /// Angular dimension (angle between two lines/planes).
    Angular,
    /// Radial dimension (radius of arc/circle).
    Radial,
    /// Diameter dimension (diameter of circle/cylinder).
    Diameter,
    /// Coordinate dimension (X, Y, Z position).
    Coordinate,
    /// Chamfer dimension (45° or angle+distance).
    Chamfer,
    /// Fillet/round radius.
    Fillet,
}

/// Geometric tolerance type (GD&T symbols).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GeometricToleranceType {
    // Location tolerances
    /// Position tolerance (true position).
    Position,
    /// Concentricity tolerance.
    Concentricity,
    /// Symmetry tolerance.
    Symmetry,
    // Orientation tolerances
    /// Perpendicularity tolerance (90°).
    Perpendicularity,
    /// Angularity tolerance (specified angle).
    Angularity,
    /// Parallelism tolerance (0°).
    Parallelism,
    // Form tolerances
    /// Flatness tolerance.
    Flatness,
    /// Circularity (roundness) tolerance.
    Circularity,
    /// Cylindricity tolerance.
    Cylindricity,
    /// Straightness tolerance.
    Straightness,
    // Profile tolerances
    /// Profile of a line tolerance.
    ProfileOfLine,
    /// Profile of a surface tolerance.
    ProfileOfSurface,
    // Runout tolerances
    /// Circular runout tolerance.
    CircularRunout,
    /// Total runout tolerance.
    TotalRunout,
}

/// Tolerance modifier (material condition symbols).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ToleranceModifier {
    /// Maximum material condition (MMC).
    MaximumMaterial,
    /// Least material condition (LMC).
    MinimumMaterial,
    /// Free state condition.
    FreeState,
    /// Projected tolerance zone.
    Projected,
    /// Tangent plane modifier.
    Tangent,
}

/// Dimensional tolerance specification.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DimensionalTolerance {
    /// Unique identifier for this tolerance.
    pub id: u64,
    /// Human-readable name/label.
    pub name: String,
    /// Type of dimension.
    pub tolerance_type: DimensionType,
    /// Nominal (target) value.
    pub nominal_value: f64,
    /// Upper deviation from nominal (positive = larger).
    pub upper_deviation: f64,
    /// Lower deviation from nominal (negative = smaller).
    pub lower_deviation: f64,
    /// Unit of measurement (e.g., "mm", "in", "deg").
    pub unit: String,
    /// Indices of faces this dimension applies to.
    pub attached_faces: Vec<usize>,
    /// Indices of datum references (if applicable).
    pub datum_references: Vec<usize>,
}

impl DimensionalTolerance {
    /// Create a new dimensional tolerance.
    pub fn new(id: u64, name: impl Into<String>, tolerance_type: DimensionType, nominal_value: f64) -> Self {
        Self {
            id,
            name: name.into(),
            tolerance_type,
            nominal_value,
            upper_deviation: 0.0,
            lower_deviation: 0.0,
            unit: "mm".to_string(),
            attached_faces: Vec::new(),
            datum_references: Vec::new(),
        }
    }

    /// Get the upper limit (nominal + upper_deviation).
    pub fn upper_limit(&self) -> f64 {
        self.nominal_value + self.upper_deviation
    }

    /// Get the lower limit (nominal + lower_deviation).
    pub fn lower_limit(&self) -> f64 {
        self.nominal_value + self.lower_deviation
    }

    /// Get the tolerance band width.
    pub fn tolerance_band(&self) -> f64 {
        self.upper_deviation - self.lower_deviation
    }
}

/// Geometric tolerance object (GD&T frame).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeometricToleranceObject {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name/label.
    pub name: String,
    /// Type of geometric tolerance.
    pub tolerance_type: GeometricToleranceType,
    /// Tolerance zone width/diameter.
    pub tolerance_value: f64,
    /// Optional material condition modifier.
    pub modifier: Option<ToleranceModifier>,
    /// Optional datum system this tolerance references.
    pub datum_system_id: Option<usize>,
    /// Indices of geometry (faces/edges) this tolerance applies to.
    pub attached_geometry: Vec<usize>,
}

impl GeometricToleranceObject {
    /// Create a new geometric tolerance.
    pub fn new(
        id: u64,
        name: impl Into<String>,
        tolerance_type: GeometricToleranceType,
        tolerance_value: f64,
    ) -> Self {
        Self {
            id,
            name: name.into(),
            tolerance_type,
            tolerance_value,
            modifier: None,
            datum_system_id: None,
            attached_geometry: Vec::new(),
        }
    }
}

/// Datum reference in a datum system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatumReference {
    /// Index of the datum being referenced.
    pub datum_id: usize,
    /// Optional modifier for this reference.
    pub modifier: Option<ToleranceModifier>,
    /// Precedence order (1 = primary, 2 = secondary, 3 = tertiary).
    pub precedence: u32,
}

/// Datum system (collection of datums with precedence).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatumSystem {
    /// Unique identifier.
    pub id: u64,
    /// Human-readable name (e.g., "A", "A-B", "A-B-C").
    pub name: String,
    /// Ordered list of datum references.
    pub datums: Vec<DatumReference>,
}

impl DatumSystem {
    /// Create a new datum system.
    pub fn new(id: u64, name: impl Into<String>) -> Self {
        Self {
            id,
            name: name.into(),
            datums: Vec::new(),
        }
    }

    /// Add a datum reference.
    pub fn add_datum(&mut self, datum_id: usize, modifier: Option<ToleranceModifier>) {
        let precedence = self.datums.len() as u32 + 1;
        self.datums.push(DatumReference {
            datum_id,
            modifier,
            precedence,
        });
    }

    /// Get the primary datum (first in precedence).
    pub fn primary_datum(&self) -> Option<&DatumReference> {
        self.datums.first()
    }
}

/// Storage container for all dimension and tolerance data.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct DimTolStore {
    /// All dimensional tolerances.
    pub dimensional_tolerances: Vec<DimensionalTolerance>,
    /// All geometric tolerances.
    pub geometric_tolerances: Vec<GeometricToleranceObject>,
    /// All datum systems.
    pub datum_systems: Vec<DatumSystem>,
}

impl DimTolStore {
    /// Create an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a dimensional tolerance.
    pub fn add_dimensional_tolerance(&mut self, tol: DimensionalTolerance) -> usize {
        let idx = self.dimensional_tolerances.len();
        self.dimensional_tolerances.push(tol);
        idx
    }

    /// Add a geometric tolerance.
    pub fn add_geometric_tolerance(&mut self, tol: GeometricToleranceObject) -> usize {
        let idx = self.geometric_tolerances.len();
        self.geometric_tolerances.push(tol);
        idx
    }

    /// Add a datum system.
    pub fn add_datum_system(&mut self, system: DatumSystem) -> usize {
        let idx = self.datum_systems.len();
        self.datum_systems.push(system);
        idx
    }

    /// Get all geometric tolerances attached to a specific face.
    pub fn tolerances_for_face(&self, face_idx: usize) -> Vec<&GeometricToleranceObject> {
        self.geometric_tolerances
            .iter()
            .filter(|t| t.attached_geometry.contains(&face_idx))
            .collect()
    }

    /// Get all dimensional tolerances attached to a specific face.
    pub fn dimensional_tolerances_for_face(&self, face_idx: usize) -> Vec<&DimensionalTolerance> {
        self.dimensional_tolerances
            .iter()
            .filter(|t| t.attached_faces.contains(&face_idx))
            .collect()
    }

    /// Find datum system by ID.
    pub fn get_datum_system(&self, id: u64) -> Option<&DatumSystem> {
        self.datum_systems.iter().find(|d| d.id == id)
    }

    /// Count of all tolerances.
    pub fn total_tolerance_count(&self) -> usize {
        self.dimensional_tolerances.len() + self.geometric_tolerances.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimensional_tolerance_creation() {
        let tol = DimensionalTolerance::new(1, "Length", DimensionType::Linear, 100.0);
        assert_eq!(tol.id, 1);
        assert_eq!(tol.name, "Length");
        assert_eq!(tol.tolerance_type, DimensionType::Linear);
        assert_eq!(tol.nominal_value, 100.0);
        assert_eq!(tol.tolerance_band(), 0.0);
    }

    #[test]
    fn dimensional_tolerance_limits() {
        let mut tol = DimensionalTolerance::new(1, "Length", DimensionType::Linear, 100.0);
        tol.upper_deviation = 0.2;
        tol.lower_deviation = -0.1;
        assert_eq!(tol.upper_limit(), 100.2);
        assert_eq!(tol.lower_limit(), 99.9);
        assert!((tol.tolerance_band() - 0.3).abs() < 1e-10);
    }

    #[test]
    fn geometric_tolerance_creation() {
        let tol = GeometricToleranceObject::new(1, "Flatness", GeometricToleranceType::Flatness, 0.05);
        assert_eq!(tol.id, 1);
        assert_eq!(tol.tolerance_type, GeometricToleranceType::Flatness);
        assert_eq!(tol.tolerance_value, 0.05);
        assert!(tol.modifier.is_none());
    }

    #[test]
    fn datum_system_creation() {
        let mut system = DatumSystem::new(1, "A-B-C");
        system.add_datum(0, None);
        system.add_datum(1, Some(ToleranceModifier::MaximumMaterial));
        system.add_datum(2, None);

        assert_eq!(system.datums.len(), 3);
        assert_eq!(system.primary_datum().unwrap().datum_id, 0);
        assert_eq!(system.primary_datum().unwrap().precedence, 1);
    }

    #[test]
    fn dim_tol_store_operations() {
        let mut store = DimTolStore::new();

        let dim_tol = DimensionalTolerance::new(1, "Diameter", DimensionType::Diameter, 50.0);
        store.add_dimensional_tolerance(dim_tol);

        let geo_tol = GeometricToleranceObject::new(1, "Position", GeometricToleranceType::Position, 0.1);
        store.add_geometric_tolerance(geo_tol);

        let mut datum_system = DatumSystem::new(1, "A");
        datum_system.add_datum(0, None);
        store.add_datum_system(datum_system);

        assert_eq!(store.total_tolerance_count(), 2);
        assert_eq!(store.datum_systems.len(), 1);
    }
}
