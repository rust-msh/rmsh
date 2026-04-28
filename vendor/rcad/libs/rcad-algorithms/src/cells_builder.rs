//! Reusable split-cell expression builder.
//!
//! This provides a lightweight CellsBuilder API analogous to OCCT
//! `BOPAlgo_CellsBuilder`: callers register reusable cell solids and then
//! evaluate boolean expressions over those cells.

use crate::{BooleanError, BooleanOpType, boolean_op};
use rcad_kernel::BRep;

/// Boolean expression over registered cells.
#[derive(Debug, Clone)]
pub enum CellExpr {
    /// Reference a registered cell by index.
    Cell(usize),
    /// Union of two expressions.
    Union(Box<CellExpr>, Box<CellExpr>),
    /// Intersection of two expressions.
    Intersection(Box<CellExpr>, Box<CellExpr>),
    /// Difference of two expressions: left - right.
    Difference(Box<CellExpr>, Box<CellExpr>),
    /// XOR: symmetric difference (A xor B = (A - B) ∪ (B - A))
    Xor(Box<CellExpr>, Box<CellExpr>),
}

/// Error type for [`CellsBuilder`].
#[derive(Debug)]
pub enum CellsBuilderError {
    /// Referenced cell index does not exist.
    InvalidCellIndex { index: usize, count: usize },
    /// Underlying boolean operation failed.
    Boolean(BooleanError),
}

impl std::fmt::Display for CellsBuilderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCellIndex { index, count } => {
                write!(f, "invalid cell index {index}; available cells: 0..{count}")
            }
            Self::Boolean(e) => write!(f, "boolean operation failed: {e}"),
        }
    }
}

impl std::error::Error for CellsBuilderError {}

impl From<BooleanError> for CellsBuilderError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

/// Reusable cell container and expression evaluator.
#[derive(Debug, Clone, Default)]
pub struct CellsBuilder {
    cells: Vec<BRep>,
}

impl CellsBuilder {
    /// Create an empty builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a builder from precomputed cells.
    pub fn from_cells(cells: Vec<BRep>) -> Self {
        Self { cells }
    }

    /// Add one cell and return its index.
    pub fn add_cell(&mut self, cell: BRep) -> usize {
        self.cells.push(cell);
        self.cells.len() - 1
    }

    /// Number of registered cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Evaluate a boolean expression over registered cells.
    pub fn evaluate(&self, expr: &CellExpr) -> Result<BRep, CellsBuilderError> {
        self.eval_rec(expr)
    }

    fn eval_rec(&self, expr: &CellExpr) -> Result<BRep, CellsBuilderError> {
        match expr {
            CellExpr::Cell(i) => self
                .cells
                .get(*i)
                .cloned()
                .ok_or(CellsBuilderError::InvalidCellIndex {
                    index: *i,
                    count: self.cells.len(),
                }),
            CellExpr::Union(a, b) => self.eval_bin(BooleanOpType::Union, a, b),
            CellExpr::Intersection(a, b) => self.eval_bin(BooleanOpType::Intersection, a, b),
            CellExpr::Difference(a, b) => self.eval_bin(BooleanOpType::Difference, a, b),
            CellExpr::Xor(a, b) => {
                // XOR: (A - B) ∪ (B - A)
                let a_min_b = self.eval_bin(BooleanOpType::Difference, a, b)?;
                let b_min_a = self.eval_bin(BooleanOpType::Difference, b, a)?;
                Ok(boolean_op(BooleanOpType::Union, &a_min_b, &b_min_a)?)
            }
        }
    }

    fn eval_bin(
        &self,
        op: BooleanOpType,
        a: &CellExpr,
        b: &CellExpr,
    ) -> Result<BRep, CellsBuilderError> {
        let left = self.eval_rec(a)?;
        let right = self.eval_rec(b)?;
        Ok(boolean_op(op, &left, &right)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::DVec3;
    use rcad_kernel::PrimitiveSolid;

    fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: w,
            height: h,
            depth: d,
        });
        for v in &mut brep.vertices {
            v.point += DVec3::new(x, y, z);
        }
        brep
    }

    fn face_count_of(brep: &BRep) -> usize {
        brep.solids
            .iter()
            .flat_map(|s| &s.shells)
            .flat_map(|sh| &sh.faces)
            .count()
    }

    #[test]
    fn cells_builder_union_expression_succeeds() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let builder = CellsBuilder::from_cells(vec![a, b]);
        let expr = CellExpr::Union(Box::new(CellExpr::Cell(0)), Box::new(CellExpr::Cell(1)));

        let out = builder.evaluate(&expr).expect("union expression should succeed");
        assert!(face_count_of(&out) > 0);
    }

    #[test]
    fn cells_builder_invalid_index_returns_error() {
        let builder = CellsBuilder::new();
        let expr = CellExpr::Cell(10);
        let err = builder.evaluate(&expr).expect_err("invalid index should fail");
        assert!(matches!(
            err,
            CellsBuilderError::InvalidCellIndex { index: 10, count: 0 }
        ));
    }

    #[test]
    fn cells_builder_xor_expression_succeeds() {
        let a = box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
        let b = box_at(0.5, 0.0, 0.0, 1.0, 1.0, 1.0);

        let builder = CellsBuilder::from_cells(vec![a, b]);

        // XOR: (A ∪ B) - (A ∩ B)
        let expr = CellExpr::Xor(Box::new(CellExpr::Cell(0)), Box::new(CellExpr::Cell(1)));

        let out = builder.evaluate(&expr).expect("XOR expression should succeed");

        // XOR of two overlapping boxes should have more faces than either input
        assert!(face_count_of(&out) > 6);
    }
}
