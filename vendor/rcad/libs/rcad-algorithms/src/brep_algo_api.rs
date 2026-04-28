//! BRepAlgoAPI-style high-level boolean algorithms.
//!
//! This module provides OCCT BRepAlgoAPI-like API for boolean operations:
//! - **BRepAlgoAPI_Common**: Intersection (common volume)
//! - **BRepAlgoAPI_Fuse**: Union (fuse shapes together)
//! - **BRepAlgoAPI_Cut**: Difference (cut one shape from another)
//! - **BRepAlgoAPI_Section**: Section (intersection curves)
//!
//! # Example
//!
//! ```
//! use rcad_algorithms::brep_algo_api::{BRepAlgoAPI_Fuse, BooleanApiOptions};
//! use rcad_kernel::{BRep, PrimitiveSolid};
//!
//! let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
//! let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 3.0 });
//!
//! let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
//! fuse.set_options(BooleanApiOptions::default().with_fuzzy_value(1e-6));
//!
//! if fuse.build() {
//!     let result = fuse.shape();
//!     println!("Fuse result has {} faces", result.solids[0].shells[0].faces.len());
//! }
//! ```

// Allow OCCT-style naming conventions (e.g., BRepAlgoAPI_Common)
#![allow(non_camel_case_types)]

use rcad_kernel::BRep;
use std::collections::HashMap;

use crate::bopds::ds::DS;
use crate::builder::{BooleanBuilder, BooleanError, BooleanOpType};
use crate::history::{BooleanHistory, EdgeOrigin, FaceOrigin, VertexOrigin, HistoryStatistics};
use crate::pave_filler::PaveFiller;
use crate::section::section;
use crate::tolerance::TOLERANCE_ABS;
use crate::geom_populate;
use crate::bvh;
use crate::{HealingOptions, SimplifyOptions};

/// Options for BRepAlgoAPI boolean operations.
///
/// Provides builder-style configuration for boolean operation parameters.
/// This is a simplified version of the full `BooleanOptions` focused on
/// the common use cases for the BRepAlgoAPI-style interface.
#[derive(Debug, Clone)]
pub struct BooleanApiOptions {
    /// Fuzzy tolerance for near-miss interference detection.
    /// Vertices/edges within this distance are considered coincident.
    pub fuzzy_value: f64,
    /// Use parallel execution where possible.
    pub parallel: bool,
    /// Track history of shape modifications.
    pub history: bool,
    /// Use BVH acceleration for intersection detection.
    pub use_bvh: bool,
    /// Run healing after boolean operation.
    pub run_healing: bool,
    /// Run simplification after boolean operation.
    pub run_simplify: bool,
    /// Simplification options.
    pub simplify_options: SimplifyOptions,
    /// Healing options.
    pub healing_options: HealingOptions,
}

impl Default for BooleanApiOptions {
    fn default() -> Self {
        Self {
            fuzzy_value: TOLERANCE_ABS,
            parallel: false,
            history: false,
            use_bvh: true,
            run_healing: false,
            run_simplify: false,
            simplify_options: SimplifyOptions::default(),
            healing_options: HealingOptions::default(),
        }
    }
}

impl BooleanApiOptions {
    /// Create default options.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set fuzzy tolerance for near-miss interference detection.
    ///
    /// Larger values allow more tolerance for near-coincident geometry.
    pub fn with_fuzzy_value(mut self, tol: f64) -> Self {
        self.fuzzy_value = tol.max(TOLERANCE_ABS);
        self
    }

    /// Enable or disable parallel execution.
    ///
    /// When enabled, uses Rayon for parallel face processing.
    pub fn with_parallel(mut self, enabled: bool) -> Self {
        self.parallel = enabled;
        self
    }

    /// Enable or disable history tracking.
    ///
    /// When enabled, tracks modifications, generations, and deletions.
    pub fn with_history(mut self, enabled: bool) -> Self {
        self.history = enabled;
        self
    }

    /// Enable or disable BVH acceleration.
    pub fn with_bvh(mut self, enabled: bool) -> Self {
        self.use_bvh = enabled;
        self
    }

    /// Enable or disable post-operation healing.
    pub fn with_healing(mut self, enabled: bool) -> Self {
        self.run_healing = enabled;
        self
    }

    /// Enable or disable post-operation simplification.
    pub fn with_simplify(mut self, enabled: bool) -> Self {
        self.run_simplify = enabled;
        self
    }

    /// Set simplification options.
    pub fn with_simplify_options(mut self, options: SimplifyOptions) -> Self {
        self.simplify_options = options;
        self
    }

    /// Set healing options.
    pub fn with_healing_options(mut self, options: HealingOptions) -> Self {
        self.healing_options = options;
        self
    }
}

/// History tracking for boolean operations.
///
/// Provides OCCT BRepAlgoAPI_BuilderShape-like history queries:
/// - Modified shapes (source shape split/modified into result shapes)
/// - Generated shapes (new shapes created during operation)
/// - Deleted shapes (shapes removed during operation)
#[derive(Debug, Clone, Default)]
pub struct BRepHistory {
    /// The underlying boolean history.
    inner: Option<BooleanHistory>,
    /// Map from source face index to result face indices (for shape A).
    modified_a: HashMap<usize, Vec<usize>>,
    /// Map from source face index to result face indices (for shape B).
    modified_b: HashMap<usize, Vec<usize>>,
    /// Generated face indices.
    generated_faces: Vec<usize>,
    /// Generated edge indices.
    generated_edges: Vec<usize>,
    /// Generated vertex indices.
    generated_vertices: Vec<usize>,
    /// Deleted face indices from shape A.
    deleted_a: Vec<usize>,
    /// Deleted face indices from shape B.
    deleted_b: Vec<usize>,
    /// Whether history tracking is enabled.
    is_generated: bool,
}

impl BRepHistory {
    /// Create a new empty history.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create history from a BooleanHistory.
    pub fn from_boolean_history(history: BooleanHistory) -> Self {
        let mut modified_a: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut modified_b: HashMap<usize, Vec<usize>> = HashMap::new();
        let mut generated_faces = Vec::new();
        let mut deleted_a = Vec::new();
        let mut deleted_b = Vec::new();

        // Build modification maps
        for (result_idx, origin) in history.face_origins.iter().enumerate() {
            match origin {
                FaceOrigin::FromA(src_idx) => {
                    modified_a.entry(*src_idx).or_default().push(result_idx);
                }
                FaceOrigin::FromB(src_idx) => {
                    modified_b.entry(*src_idx).or_default().push(result_idx);
                }
                FaceOrigin::Generated => {
                    generated_faces.push(result_idx);
                }
            }
        }

        // Extract generated edges
        let generated_edges: Vec<usize> = history.edge_origins
            .iter()
            .enumerate()
            .filter_map(|(idx, origin)| {
                matches!(origin, EdgeOrigin::Generated).then_some(idx)
            })
            .collect();

        // Extract generated vertices
        let generated_vertices: Vec<usize> = history.vertex_origins
            .iter()
            .enumerate()
            .filter_map(|(idx, origin)| {
                matches!(origin, VertexOrigin::Intersection).then_some(idx)
            })
            .collect();

        // Extract deleted faces
        deleted_a = history.deleted_from_a.clone();
        deleted_b = history.deleted_from_b.clone();

        Self {
            inner: Some(history),
            modified_a,
            modified_b,
            generated_faces,
            generated_edges,
            generated_vertices,
            deleted_a,
            deleted_b,
            is_generated: true,
        }
    }

    /// Check if history tracking was enabled.
    pub fn is_generated(&self) -> bool {
        self.is_generated
    }

    /// Returns true if any shapes were modified.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasModified()`.
    pub fn has_modified(&self) -> bool {
        !self.modified_a.is_empty() || !self.modified_b.is_empty()
    }

    /// Returns true if any shapes were generated.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasGenerated()`.
    pub fn has_generated(&self) -> bool {
        !self.generated_faces.is_empty()
            || !self.generated_edges.is_empty()
            || !self.generated_vertices.is_empty()
    }

    /// Returns true if any shapes were deleted.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::HasDeleted()`.
    pub fn has_deleted(&self) -> bool {
        !self.deleted_a.is_empty() || !self.deleted_b.is_empty()
    }

    /// Get the result faces that came from a source face in shape A.
    /// Returns an empty vector if the face was deleted or not found.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Modified()`.
    pub fn modified_from_a(&self, source_face_idx: usize) -> &[usize] {
        self.modified_a.get(&source_face_idx).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Get the result faces that came from a source face in shape B.
    pub fn modified_from_b(&self, source_face_idx: usize) -> &[usize] {
        self.modified_b.get(&source_face_idx).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Get all generated faces.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::Generated()`.
    pub fn generated_faces(&self) -> &[usize] {
        &self.generated_faces
    }

    /// Get all generated edges.
    pub fn generated_edges(&self) -> &[usize] {
        &self.generated_edges
    }

    /// Get all generated vertices.
    pub fn generated_vertices(&self) -> &[usize] {
        &self.generated_vertices
    }

    /// Check if a face from shape A was deleted.
    /// Analogous to OCCT `BRepAlgoAPI_BuilderShape::IsDeleted()`.
    pub fn is_deleted_from_a(&self, source_face_idx: usize) -> bool {
        self.deleted_a.contains(&source_face_idx)
    }

    /// Check if a face from shape B was deleted.
    pub fn is_deleted_from_b(&self, source_face_idx: usize) -> bool {
        self.deleted_b.contains(&source_face_idx)
    }

    /// Get the underlying BooleanHistory if available.
    pub fn inner(&self) -> Option<&BooleanHistory> {
        self.inner.as_ref()
    }

    /// Get statistics about the history.
    pub fn statistics(&self) -> HistoryStatistics {
        HistoryStatistics {
            modified_faces: self.modified_a.values().map(|v| v.len()).sum::<usize>()
                + self.modified_b.values().map(|v| v.len()).sum::<usize>(),
            modified_edges: 0,
            modified_vertices: 0,
            generated_faces: self.generated_faces.len(),
            generated_edges: self.generated_edges.len(),
            generated_vertices: self.generated_vertices.len(),
            deleted_faces: self.deleted_a.len() + self.deleted_b.len(),
            deleted_edges: 0,
            deleted_vertices: 0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepAlgoAPI_Common (Intersection)
// ─────────────────────────────────────────────────────────────────────────────

/// BRepAlgoAPI_Common - Compute the common (intersection) of two shapes.
///
/// Computes the intersection volume of two shapes.
/// Analogous to OCCT `BRepAlgoAPI_Common`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo_api::BRepAlgoAPI_Common;
/// use rcad_kernel::{BRep, PrimitiveSolid};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 3.0 });
///
/// let mut common = BRepAlgoAPI_Common::new(&box1, &box2);
/// if common.build() {
///     let result = common.shape();
/// }
/// ```
pub struct BRepAlgoAPI_Common<'a> {
    shape1: &'a BRep,
    shape2: &'a BRep,
    options: BooleanApiOptions,
    result: Option<BRep>,
    history: BRepHistory,
    error: Option<BooleanError>,
}

impl<'a> BRepAlgoAPI_Common<'a> {
    /// Create a new Common operation.
    pub fn new(shape1: &'a BRep, shape2: &'a BRep) -> Self {
        Self {
            shape1,
            shape2,
            options: BooleanApiOptions::default(),
            result: None,
            history: BRepHistory::new(),
            error: None,
        }
    }

    /// Set operation options.
    pub fn set_options(&mut self, options: BooleanApiOptions) {
        self.options = options;
    }

    /// Get the current options.
    pub fn options(&self) -> &BooleanApiOptions {
        &self.options
    }

    /// Build the result.
    ///
    /// Returns `true` if the operation succeeded.
    pub fn build(&mut self) -> bool {
        match self.build_internal() {
            Ok(_) => true,
            Err(e) => {
                self.error = Some(e);
                false
            }
        }
    }

    fn build_internal(&mut self) -> Result<(), BooleanError> {
        // Check for empty inputs
        if self.shape1.solids.is_empty() || self.shape2.solids.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        // Ensure geometry is populated for primitives
        let a = self.ensure_geometry(self.shape1);
        let b = self.ensure_geometry(self.shape2);

        // Build DS
        let mut ds = if self.options.fuzzy_value > TOLERANCE_ABS {
            DS::new_with_fuzzy(&a, &b, self.options.fuzzy_value)
        } else {
            DS::new(&a, &b)
        };

        // Run PaveFiller
        let bvh_a = if self.options.use_bvh { Some(bvh::Bvh::build(&a)) } else { None };
        let bvh_b = if self.options.use_bvh { Some(bvh::Bvh::build(&b)) } else { None };

        let mut filler = match (&bvh_a, &bvh_b) {
            (Some(ba), Some(bb)) => PaveFiller::with_bvh(&mut ds, ba, bb),
            _ => PaveFiller::new(&mut ds),
        };
        filler.perform();

        // Build result
        let builder = BooleanBuilder::new(&ds, BooleanOpType::Intersection);
        let (brep, bool_history) = if self.options.parallel {
            builder.build_with_history_par()?
        } else {
            builder.build_with_history()?
        };

        // Store history if requested
        if self.options.history {
            self.history = BRepHistory::from_boolean_history(bool_history);
        }

        // Apply post-processing
        let mut result = brep;
        if self.options.run_healing {
            let (healed, _) = crate::healing::analyze_and_heal(&result, self.options.healing_options.clone());
            result = healed;
        }
        if self.options.run_simplify {
            let (simplified, _) = crate::simplify_brep_post_ops(&result, self.options.simplify_options.clone());
            result = simplified;
        }

        self.result = Some(result);
        Ok(())
    }

    fn ensure_geometry(&self, brep: &'a BRep) -> BRep {
        // Check if geometry needs to be populated
        if brep.geom.surfaces.is_empty() && !brep.solids.is_empty() {
            let mut result = brep.clone();
            geom_populate::populate_box_geom(&mut result);
            result
        } else {
            brep.clone()
        }
    }

    /// Get the result shape.
    ///
    /// Panics if `build()` has not been called or failed.
    pub fn shape(&self) -> &BRep {
        self.result.as_ref().expect("build() must be called before shape()")
    }

    /// Get the result shape, consuming the builder.
    pub fn into_shape(self) -> Option<BRep> {
        self.result
    }

    /// Get the history of the operation.
    pub fn history(&self) -> &BRepHistory {
        &self.history
    }

    /// Get the error if the operation failed.
    pub fn error(&self) -> Option<&BooleanError> {
        self.error.as_ref()
    }

    /// Check if the operation has been built.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepAlgoAPI_Fuse (Union)
// ─────────────────────────────────────────────────────────────────────────────

/// BRepAlgoAPI_Fuse - Compute the union (fuse) of two shapes.
///
/// Computes the union volume of two shapes.
/// Analogous to OCCT `BRepAlgoAPI_Fuse`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo_api::BRepAlgoAPI_Fuse;
/// use rcad_kernel::{BRep, PrimitiveSolid};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 1.0, height: 1.0, depth: 3.0 });
///
/// let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
/// if fuse.build() {
///     let result = fuse.shape();
/// }
/// ```
pub struct BRepAlgoAPI_Fuse<'a> {
    shape1: &'a BRep,
    shape2: &'a BRep,
    options: BooleanApiOptions,
    result: Option<BRep>,
    history: BRepHistory,
    error: Option<BooleanError>,
}

impl<'a> BRepAlgoAPI_Fuse<'a> {
    /// Create a new Fuse operation.
    pub fn new(shape1: &'a BRep, shape2: &'a BRep) -> Self {
        Self {
            shape1,
            shape2,
            options: BooleanApiOptions::default(),
            result: None,
            history: BRepHistory::new(),
            error: None,
        }
    }

    /// Set operation options.
    pub fn set_options(&mut self, options: BooleanApiOptions) {
        self.options = options;
    }

    /// Get the current options.
    pub fn options(&self) -> &BooleanApiOptions {
        &self.options
    }

    /// Build the result.
    ///
    /// Returns `true` if the operation succeeded.
    pub fn build(&mut self) -> bool {
        match self.build_internal() {
            Ok(_) => true,
            Err(e) => {
                self.error = Some(e);
                false
            }
        }
    }

    fn build_internal(&mut self) -> Result<(), BooleanError> {
        // Check for empty inputs
        if self.shape1.solids.is_empty() || self.shape2.solids.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        let a = self.ensure_geometry(self.shape1);
        let b = self.ensure_geometry(self.shape2);

        let mut ds = if self.options.fuzzy_value > TOLERANCE_ABS {
            DS::new_with_fuzzy(&a, &b, self.options.fuzzy_value)
        } else {
            DS::new(&a, &b)
        };

        let bvh_a = if self.options.use_bvh { Some(bvh::Bvh::build(&a)) } else { None };
        let bvh_b = if self.options.use_bvh { Some(bvh::Bvh::build(&b)) } else { None };

        let mut filler = match (&bvh_a, &bvh_b) {
            (Some(ba), Some(bb)) => PaveFiller::with_bvh(&mut ds, ba, bb),
            _ => PaveFiller::new(&mut ds),
        };
        filler.perform();

        let builder = BooleanBuilder::new(&ds, BooleanOpType::Union);
        let (brep, bool_history) = if self.options.parallel {
            builder.build_with_history_par()?
        } else {
            builder.build_with_history()?
        };

        if self.options.history {
            self.history = BRepHistory::from_boolean_history(bool_history);
        }

        let mut result = brep;
        if self.options.run_healing {
            let (healed, _) = crate::healing::analyze_and_heal(&result, self.options.healing_options.clone());
            result = healed;
        }
        if self.options.run_simplify {
            let (simplified, _) = crate::simplify_brep_post_ops(&result, self.options.simplify_options.clone());
            result = simplified;
        }

        self.result = Some(result);
        Ok(())
    }

    fn ensure_geometry(&self, brep: &'a BRep) -> BRep {
        if brep.geom.surfaces.is_empty() && !brep.solids.is_empty() {
            let mut result = brep.clone();
            geom_populate::populate_box_geom(&mut result);
            result
        } else {
            brep.clone()
        }
    }

    /// Get the result shape.
    pub fn shape(&self) -> &BRep {
        self.result.as_ref().expect("build() must be called before shape()")
    }

    /// Get the result shape, consuming the builder.
    pub fn into_shape(self) -> Option<BRep> {
        self.result
    }

    /// Get the history of the operation.
    pub fn history(&self) -> &BRepHistory {
        &self.history
    }

    /// Get the error if the operation failed.
    pub fn error(&self) -> Option<&BooleanError> {
        self.error.as_ref()
    }

    /// Check if the operation has been built.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepAlgoAPI_Cut (Difference)
// ─────────────────────────────────────────────────────────────────────────────

/// BRepAlgoAPI_Cut - Compute the difference (cut) of two shapes.
///
/// Computes shape1 minus shape2.
/// Analogous to OCCT `BRepAlgoAPI_Cut`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo_api::BRepAlgoAPI_Cut;
/// use rcad_kernel::{BRep, PrimitiveSolid};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let cylinder = BRep::from_primitive(PrimitiveSolid::Cylinder { radius: 0.5, height: 3.0 });
///
/// let mut cut = BRepAlgoAPI_Cut::new(&box1, &cylinder);
/// if cut.build() {
///     let result = cut.shape();
/// }
/// ```
pub struct BRepAlgoAPI_Cut<'a> {
    shape1: &'a BRep,
    shape2: &'a BRep,
    options: BooleanApiOptions,
    result: Option<BRep>,
    history: BRepHistory,
    error: Option<BooleanError>,
}

impl<'a> BRepAlgoAPI_Cut<'a> {
    /// Create a new Cut operation.
    /// `shape1` is the shape to cut from, `shape2` is the cutting tool.
    pub fn new(shape1: &'a BRep, shape2: &'a BRep) -> Self {
        Self {
            shape1,
            shape2,
            options: BooleanApiOptions::default(),
            result: None,
            history: BRepHistory::new(),
            error: None,
        }
    }

    /// Set operation options.
    pub fn set_options(&mut self, options: BooleanApiOptions) {
        self.options = options;
    }

    /// Get the current options.
    pub fn options(&self) -> &BooleanApiOptions {
        &self.options
    }

    /// Build the result.
    ///
    /// Returns `true` if the operation succeeded.
    pub fn build(&mut self) -> bool {
        match self.build_internal() {
            Ok(_) => true,
            Err(e) => {
                self.error = Some(e);
                false
            }
        }
    }

    fn build_internal(&mut self) -> Result<(), BooleanError> {
        // Check for empty inputs
        if self.shape1.solids.is_empty() || self.shape2.solids.is_empty() {
            return Err(BooleanError::EmptyInput);
        }

        let a = self.ensure_geometry(self.shape1);
        let b = self.ensure_geometry(self.shape2);

        let mut ds = if self.options.fuzzy_value > TOLERANCE_ABS {
            DS::new_with_fuzzy(&a, &b, self.options.fuzzy_value)
        } else {
            DS::new(&a, &b)
        };

        let bvh_a = if self.options.use_bvh { Some(bvh::Bvh::build(&a)) } else { None };
        let bvh_b = if self.options.use_bvh { Some(bvh::Bvh::build(&b)) } else { None };

        let mut filler = match (&bvh_a, &bvh_b) {
            (Some(ba), Some(bb)) => PaveFiller::with_bvh(&mut ds, ba, bb),
            _ => PaveFiller::new(&mut ds),
        };
        filler.perform();

        let builder = BooleanBuilder::new(&ds, BooleanOpType::Difference);
        let (brep, bool_history) = if self.options.parallel {
            builder.build_with_history_par()?
        } else {
            builder.build_with_history()?
        };

        if self.options.history {
            self.history = BRepHistory::from_boolean_history(bool_history);
        }

        let mut result = brep;
        if self.options.run_healing {
            let (healed, _) = crate::healing::analyze_and_heal(&result, self.options.healing_options.clone());
            result = healed;
        }
        if self.options.run_simplify {
            let (simplified, _) = crate::simplify_brep_post_ops(&result, self.options.simplify_options.clone());
            result = simplified;
        }

        self.result = Some(result);
        Ok(())
    }

    fn ensure_geometry(&self, brep: &'a BRep) -> BRep {
        if brep.geom.surfaces.is_empty() && !brep.solids.is_empty() {
            let mut result = brep.clone();
            geom_populate::populate_box_geom(&mut result);
            result
        } else {
            brep.clone()
        }
    }

    /// Get the result shape.
    pub fn shape(&self) -> &BRep {
        self.result.as_ref().expect("build() must be called before shape()")
    }

    /// Get the result shape, consuming the builder.
    pub fn into_shape(self) -> Option<BRep> {
        self.result
    }

    /// Get the history of the operation.
    pub fn history(&self) -> &BRepHistory {
        &self.history
    }

    /// Get the error if the operation failed.
    pub fn error(&self) -> Option<&BooleanError> {
        self.error.as_ref()
    }

    /// Check if the operation has been built.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BRepAlgoAPI_Section
// ─────────────────────────────────────────────────────────────────────────────

/// BRepAlgoAPI_Section - Compute the section of two shapes.
///
/// Computes the intersection curves/wires between two shapes.
/// Analogous to OCCT `BRepAlgoAPI_Section`.
///
/// # Example
///
/// ```
/// use rcad_algorithms::brep_algo_api::BRepAlgoAPI_Section;
/// use rcad_kernel::{BRep, PrimitiveSolid};
///
/// let box1 = BRep::from_primitive(PrimitiveSolid::Box { width: 2.0, height: 2.0, depth: 2.0 });
/// let box2 = BRep::from_primitive(PrimitiveSolid::Box { width: 3.0, height: 1.0, depth: 1.0 });
///
/// let mut section = BRepAlgoAPI_Section::new(&box1, &box2);
/// if section.build() {
///     let result = section.shape();
/// }
/// ```
pub struct BRepAlgoAPI_Section<'a> {
    shape1: &'a BRep,
    shape2: &'a BRep,
    options: BooleanApiOptions,
    result: Option<BRep>,
    error: Option<BooleanError>,
}

impl<'a> BRepAlgoAPI_Section<'a> {
    /// Create a new Section operation.
    pub fn new(shape1: &'a BRep, shape2: &'a BRep) -> Self {
        Self {
            shape1,
            shape2,
            options: BooleanApiOptions::default(),
            result: None,
            error: None,
        }
    }

    /// Set operation options.
    pub fn set_options(&mut self, options: BooleanApiOptions) {
        self.options = options;
    }

    /// Get the current options.
    pub fn options(&self) -> &BooleanApiOptions {
        &self.options
    }

    /// Build the result.
    ///
    /// Returns `true` if the operation succeeded.
    pub fn build(&mut self) -> bool {
        match self.build_internal() {
            Ok(_) => true,
            Err(e) => {
                self.error = Some(e);
                false
            }
        }
    }

    fn build_internal(&mut self) -> Result<(), BooleanError> {
        // For section, we compute intersection curves between faces
        // Use the section module with a surface derived from shape2

        let a = self.ensure_geometry(self.shape1);
        let b = self.ensure_geometry(self.shape2);

        // Try to get a plane from shape2 for simple section
        if let Some(plane) = self.extract_plane(&b) {
            let result = section(&a, &plane);
            self.result = Some(result);
            return Ok(());
        }

        // For general case, compute face-face intersections
        // This is a simplified implementation that uses the boolean DS
        let mut ds = if self.options.fuzzy_value > TOLERANCE_ABS {
            DS::new_with_fuzzy(&a, &b, self.options.fuzzy_value)
        } else {
            DS::new(&a, &b)
        };

        let bvh_a = if self.options.use_bvh { Some(bvh::Bvh::build(&a)) } else { None };
        let bvh_b = if self.options.use_bvh { Some(bvh::Bvh::build(&b)) } else { None };

        let mut filler = match (&bvh_a, &bvh_b) {
            (Some(ba), Some(bb)) => PaveFiller::with_bvh(&mut ds, ba, bb),
            _ => PaveFiller::new(&mut ds),
        };
        filler.perform();

        // Extract intersection curves from DS
        let result = self.build_section_from_ds(&ds);
        self.result = Some(result);
        Ok(())
    }

    fn extract_plane(&self, brep: &BRep) -> Option<rcad_kernel::geom::Plane> {
        // Check if the BRep is a simple box or plane-based shape
        if brep.solids.len() != 1 {
            return None;
        }

        // Try to find a dominant plane
        for solid in &brep.solids {
            for shell in &solid.shells {
                for (face_idx, _face) in shell.faces.iter().enumerate() {
                    if let Some(surf_idx) = brep.geom.face_surface.get(face_idx).and_then(|o| *o) {
                        if let Some(rcad_kernel::geom::Surface3::Plane(plane)) =
                            brep.geom.surfaces.get(surf_idx)
                        {
                            return Some(*plane);
                        }
                    }
                }
            }
        }

        None
    }

    fn build_section_from_ds(&self, ds: &DS) -> BRep {
        use rcad_kernel::geom::Curve3;
        use rcad_kernel::topology::{Edge, Vertex, Wire, WireEdge};

        let mut result = BRep::new();
        let mut wires: Vec<Wire> = Vec::new();

        // Extract intersection curves
        for ic in &ds.intersection_curves {
            if ic.polyline.len() >= 2 {
                // Build wire from polyline
                let mut wire_edges = Vec::new();

                for i in 0..ic.polyline.len() - 1 {
                    let a = ic.polyline[i];
                    let b = ic.polyline[i + 1];

                    let vi_a = result.vertices.len();
                    result.vertices.push(Vertex { point: a });
                    let vi_b = result.vertices.len();
                    result.vertices.push(Vertex { point: b });

                    let edge_idx = result.edges.len();
                    result.edges.push(Edge {
                        start: vi_a,
                        end: vi_b,
                    });

                    // Store curve geometry
                    let len = (b - a).length();
                    let dir = if len > 1e-10 {
                        (b - a) / len
                    } else {
                        glam::DVec3::X
                    };
                    let curve_idx = result.geom.curves.len();
                    result.geom.curves.push(Curve3::Line(rcad_kernel::geom::Line3 {
                        origin: a,
                        direction: dir,
                    }));

                    while result.geom.edge_curve.len() <= edge_idx {
                        result.geom.edge_curve.push(None);
                    }
                    while result.geom.edge_curve_range.len() <= edge_idx {
                        result.geom.edge_curve_range.push(None);
                    }
                    result.geom.edge_curve[edge_idx] = Some(curve_idx);
                    result.geom.edge_curve_range[edge_idx] = Some([0.0, len]);

                    wire_edges.push(WireEdge::fwd(edge_idx));
                }

                if !wire_edges.is_empty() {
                    wires.push(Wire { edges: wire_edges });
                }
            }
        }

        // Pack wires into result
        if !wires.is_empty() {
            use rcad_kernel::topology::{Face, Shell, Solid};
            let faces: Vec<_> = wires
                .into_iter()
                .map(|w| Face {
                    outer_wire: w,
                    inner_wires: vec![],
                    normal: glam::DVec3::Z,
                    triangles: vec![],
                    mesh_dirty: true,
                })
                .collect();
            result.solids.push(Solid {
                shells: vec![Shell { faces }],
            });
        }

        result
    }

    fn ensure_geometry(&self, brep: &'a BRep) -> BRep {
        if brep.geom.surfaces.is_empty() && !brep.solids.is_empty() {
            let mut result = brep.clone();
            geom_populate::populate_box_geom(&mut result);
            result
        } else {
            brep.clone()
        }
    }

    /// Get the result shape.
    pub fn shape(&self) -> &BRep {
        self.result.as_ref().expect("build() must be called before shape()")
    }

    /// Get the result shape, consuming the builder.
    pub fn into_shape(self) -> Option<BRep> {
        self.result
    }

    /// Get the error if the operation failed.
    pub fn error(&self) -> Option<&BooleanError> {
        self.error.as_ref()
    }

    /// Check if the operation has been built.
    pub fn is_done(&self) -> bool {
        self.result.is_some()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rcad_kernel::PrimitiveSolid;
    use glam::DVec3;

    fn unit_box() -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    fn shifted_box(dx: f64, dy: f64, dz: f64) -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        });
        // Transform the box by shifting
        for vertex in &mut brep.vertices {
            vertex.point.x += dx;
            vertex.point.y += dy;
            vertex.point.z += dz;
        }
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    // ── BooleanApiOptions Tests ─────────────────────────────────────────────────

    #[test]
    fn options_default() {
        let opts = BooleanApiOptions::default();
        assert!(!opts.parallel);
        assert!(!opts.history);
        assert!(opts.use_bvh);
    }

    #[test]
    fn options_builder() {
        let opts = BooleanApiOptions::new()
            .with_fuzzy_value(1e-5)
            .with_parallel(true)
            .with_history(true);

        assert!((opts.fuzzy_value - 1e-5).abs() < 1e-10);
        assert!(opts.parallel);
        assert!(opts.history);
    }

    // ── BRepAlgoAPI_Fuse Tests ───────────────────────────────────────────────────

    #[test]
    fn fuse_two_boxes() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.0, 0.0);

        let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
        assert!(fuse.build());
        assert!(fuse.is_done());

        let result = fuse.shape();
        assert!(!result.solids.is_empty());
        assert!(!result.solids[0].shells.is_empty());
    }

    #[test]
    fn fuse_with_history() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.0, 0.0);

        let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
        fuse.set_options(BooleanApiOptions::default().with_history(true));
        assert!(fuse.build());

        let history = fuse.history();
        assert!(history.is_generated());
    }

    #[test]
    fn fuse_disjoint_boxes() {
        let box1 = unit_box();
        let box2 = shifted_box(5.0, 5.0, 5.0); // Far apart

        let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
        assert!(fuse.build());

        let result = fuse.shape();
        // Result should have both boxes
        assert!(!result.solids.is_empty());
    }

    // ── BRepAlgoAPI_Common Tests ─────────────────────────────────────────────────

    #[test]
    fn common_overlapping_boxes() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.5, 0.5);

        let mut common = BRepAlgoAPI_Common::new(&box1, &box2);
        assert!(common.build());
        assert!(common.is_done());

        let result = common.shape();
        assert!(!result.solids.is_empty());
    }

    #[test]
    fn common_disjoint_boxes() {
        let box1 = unit_box();
        let box2 = shifted_box(5.0, 5.0, 5.0); // Far apart, no intersection

        let mut common = BRepAlgoAPI_Common::new(&box1, &box2);
        // May fail or produce empty result
        let _ = common.build();
    }

    // ── BRepAlgoAPI_Cut Tests ────────────────────────────────────────────────────

    #[test]
    fn cut_box_from_box() {
        let box1 = unit_box();
        let box2 = shifted_box(0.25, 0.25, 0.25);

        let mut cut = BRepAlgoAPI_Cut::new(&box1, &box2);
        assert!(cut.build());
        assert!(cut.is_done());

        let result = cut.shape();
        assert!(!result.solids.is_empty());
    }

    #[test]
    fn cut_with_options() {
        let box1 = unit_box();
        let box2 = shifted_box(0.25, 0.25, 0.25);

        let mut cut = BRepAlgoAPI_Cut::new(&box1, &box2);
        cut.set_options(
            BooleanApiOptions::default()
                .with_history(true)
                .with_bvh(true),
        );
        assert!(cut.build());

        let history = cut.history();
        assert!(history.is_generated());
    }

    // ── BRepAlgoAPI_Section Tests ────────────────────────────────────────────────

    #[test]
    fn section_two_boxes() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.0, 0.0);

        let mut section = BRepAlgoAPI_Section::new(&box1, &box2);
        assert!(section.build());
        assert!(section.is_done());

        let result = section.shape();
        // Section produces wires/edges
        assert!(!result.edges.is_empty() || result.solids.is_empty());
    }

    // ── BRepHistory Tests ────────────────────────────────────────────────────────

    #[test]
    fn history_empty() {
        let history = BRepHistory::new();
        assert!(!history.is_generated());
        assert!(!history.has_modified());
        assert!(!history.has_generated());
        assert!(!history.has_deleted());
    }

    #[test]
    fn history_statistics() {
        let stats = HistoryStatistics::default();
        assert_eq!(stats.total_modified(), 0);
        assert_eq!(stats.total_generated(), 0);
        assert_eq!(stats.total_deleted(), 0);
    }

    // ── Error Handling Tests ─────────────────────────────────────────────────────

    #[test]
    fn error_on_empty_input() {
        let empty = BRep::new();
        let box1 = unit_box();

        let mut fuse = BRepAlgoAPI_Fuse::new(&empty, &box1);
        assert!(!fuse.build());
        assert!(fuse.error().is_some());
    }

    #[test]
    fn into_shape_consumes_builder() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.0, 0.0);

        let mut fuse = BRepAlgoAPI_Fuse::new(&box1, &box2);
        assert!(fuse.build());

        let result = fuse.into_shape();
        assert!(result.is_some());
    }

    // ── Convenience Functions Tests ──────────────────────────────────────────────

    #[test]
    fn common_convenience() {
        let box1 = unit_box();
        let box2 = shifted_box(0.5, 0.5, 0.5);

        let mut common = BRepAlgoAPI_Common::new(&box1, &box2);
        common.set_options(BooleanApiOptions::default().with_parallel(true));
        assert!(common.build());

        let opts = common.options();
        assert!(opts.parallel);
    }
}
