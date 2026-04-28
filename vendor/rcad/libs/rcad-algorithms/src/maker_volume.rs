//! Build solids from reusable split cells.
//!
//! This is a lightweight baseline analogous to OCCT `BOPAlgo_MakerVolume`:
//! callers register pre-split cell solids and then assemble a final solid from
//! a region mask, an explicit cell index list, or a [`CellExpr`] boolean
//! expression.

use std::collections::HashSet;

use rcad_kernel::BRep;

use crate::cells_builder::{CellExpr, CellsBuilder, CellsBuilderError};
use crate::{BooleanError, GeneralFuseHistory, general_fuse, general_fuse_with_history};

/// Error type for MakerVolume-style solid assembly.
#[derive(Debug)]
pub enum MakerVolumeError {
    /// No cells were provided.
    EmptyInput,
    /// The region selection resolved to no cells.
    EmptySelection,
    /// A region mask did not match the number of registered cells.
    RegionMaskLengthMismatch { expected: usize, got: usize },
    /// A requested cell index was out of bounds.
    InvalidCellIndex { index: usize, count: usize },
    /// Boolean assembly failed.
    Boolean(BooleanError),
    /// Expression-based assembly failed.
    Expression(CellsBuilderError),
}

impl std::fmt::Display for MakerVolumeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "empty input"),
            Self::EmptySelection => write!(f, "empty cell selection"),
            Self::RegionMaskLengthMismatch { expected, got } => {
                write!(f, "region mask length mismatch: expected {expected}, got {got}")
            }
            Self::InvalidCellIndex { index, count } => {
                write!(f, "invalid cell index {index}; available cells: 0..{count}")
            }
            Self::Boolean(source) => write!(f, "boolean assembly failed: {source}"),
            Self::Expression(source) => write!(f, "expression assembly failed: {source}"),
        }
    }
}

impl std::error::Error for MakerVolumeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Boolean(source) => Some(source),
            Self::Expression(source) => Some(source),
            _ => None,
        }
    }
}

impl From<BooleanError> for MakerVolumeError {
    fn from(value: BooleanError) -> Self {
        Self::Boolean(value)
    }
}

impl From<CellsBuilderError> for MakerVolumeError {
    fn from(value: CellsBuilderError) -> Self {
        Self::Expression(value)
    }
}

/// Report describing which cells were selected for a MakerVolume assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MakerVolumeSelection {
    /// Total number of input cells registered on the builder.
    pub input_cell_count: usize,
    /// Unique selected cell indices, in stable encounter order.
    pub selected_cell_indices: Vec<usize>,
}

/// Reusable solid assembler over precomputed split cells.
#[derive(Debug, Clone, Default)]
pub struct MakerVolume {
    cells: Vec<BRep>,
}

impl MakerVolume {
    /// Create an empty MakerVolume builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a MakerVolume builder from precomputed cells.
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

    /// Build a solid from all registered cells.
    pub fn build_all(&self) -> Result<BRep, MakerVolumeError> {
        let indices: Vec<usize> = (0..self.cells.len()).collect();
        self.build_from_indices(&indices)
    }

    /// Build a solid from a boolean region mask.
    pub fn build_from_region_mask(&self, region_mask: &[bool]) -> Result<BRep, MakerVolumeError> {
        let selection = self.selection_from_region_mask(region_mask)?;
        self.build_from_indices(&selection.selected_cell_indices)
    }

    /// Build a solid and per-step history from a region mask.
    pub fn build_from_region_mask_with_history(
        &self,
        region_mask: &[bool],
    ) -> Result<(BRep, GeneralFuseHistory), MakerVolumeError> {
        let selection = self.selection_from_region_mask(region_mask)?;
        self.build_from_indices_with_history(&selection.selected_cell_indices)
    }

    /// Build a solid from an explicit cell index list.
    pub fn build_from_indices(&self, indices: &[usize]) -> Result<BRep, MakerVolumeError> {
        let parts = self.selected_cells(indices)?;
        Ok(general_fuse(&parts)?)
    }

    /// Build a solid and per-step history from an explicit cell index list.
    pub fn build_from_indices_with_history(
        &self,
        indices: &[usize],
    ) -> Result<(BRep, GeneralFuseHistory), MakerVolumeError> {
        let parts = self.selected_cells(indices)?;
        Ok(general_fuse_with_history(&parts)?)
    }

    /// Build a solid from a reusable cell expression.
    pub fn build_from_expr(&self, expr: &CellExpr) -> Result<BRep, MakerVolumeError> {
        if self.cells.is_empty() {
            return Err(MakerVolumeError::EmptyInput);
        }
        let builder = CellsBuilder::from_cells(self.cells.clone());
        Ok(builder.evaluate(expr)?)
    }

    /// Convert a region mask into a validated selection report.
    pub fn selection_from_region_mask(
        &self,
        region_mask: &[bool],
    ) -> Result<MakerVolumeSelection, MakerVolumeError> {
        if self.cells.is_empty() {
            return Err(MakerVolumeError::EmptyInput);
        }
        if region_mask.len() != self.cells.len() {
            return Err(MakerVolumeError::RegionMaskLengthMismatch {
                expected: self.cells.len(),
                got: region_mask.len(),
            });
        }

        let selected_cell_indices: Vec<usize> = region_mask
            .iter()
            .enumerate()
            .filter_map(|(index, enabled)| enabled.then_some(index))
            .collect();

        if selected_cell_indices.is_empty() {
            return Err(MakerVolumeError::EmptySelection);
        }

        Ok(MakerVolumeSelection {
            input_cell_count: self.cells.len(),
            selected_cell_indices,
        })
    }

    fn selected_cells(&self, indices: &[usize]) -> Result<Vec<BRep>, MakerVolumeError> {
        if self.cells.is_empty() {
            return Err(MakerVolumeError::EmptyInput);
        }

        let unique_indices = unique_cell_indices(indices);
        if unique_indices.is_empty() {
            return Err(MakerVolumeError::EmptySelection);
        }

        let mut parts = Vec::with_capacity(unique_indices.len());
        for index in unique_indices {
            let cell = self
                .cells
                .get(index)
                .cloned()
                .ok_or(MakerVolumeError::InvalidCellIndex {
                    index,
                    count: self.cells.len(),
                })?;
            parts.push(cell);
        }
        Ok(parts)
    }
}

/// Convenience helper: assemble a solid from a region mask.
pub fn make_solid_from_region(
    cells: &[BRep],
    region_mask: &[bool],
) -> Result<BRep, MakerVolumeError> {
    MakerVolume::from_cells(cells.to_vec()).build_from_region_mask(region_mask)
}

/// Convenience helper: assemble a solid from a region mask and report history.
pub fn make_solid_from_region_with_history(
    cells: &[BRep],
    region_mask: &[bool],
) -> Result<(BRep, GeneralFuseHistory), MakerVolumeError> {
    MakerVolume::from_cells(cells.to_vec()).build_from_region_mask_with_history(region_mask)
}

/// Convenience helper: assemble a solid from explicit cell indices.
pub fn make_solid_from_cell_indices(
    cells: &[BRep],
    indices: &[usize],
) -> Result<BRep, MakerVolumeError> {
    MakerVolume::from_cells(cells.to_vec()).build_from_indices(indices)
}

fn unique_cell_indices(indices: &[usize]) -> Vec<usize> {
    let mut seen = HashSet::with_capacity(indices.len());
    let mut out = Vec::with_capacity(indices.len());
    for &index in indices {
        if seen.insert(index) {
            out.push(index);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use glam::DVec3;
    use rcad_kernel::{BRep, PrimitiveSolid, properties::volume};

    use super::*;
    use crate::geom_populate;

    fn box_at(x: f64, y: f64, z: f64, w: f64, h: f64, d: f64) -> BRep {
        let mut brep = BRep::from_primitive(PrimitiveSolid::Box {
            width: w,
            height: h,
            depth: d,
        });
        for v in &mut brep.vertices {
            v.point += DVec3::new(x, y, z);
        }
        geom_populate::populate_box_geom(&mut brep);
        brep
    }

    #[test]
    fn maker_volume_build_all_fuses_selected_cells() {
        let volume_builder = MakerVolume::from_cells(vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ]);

        let solid = volume_builder.build_all().expect("build_all should succeed");
        assert!((volume(&solid) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn maker_volume_region_mask_selects_subset() {
        let cells = vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ];

        let solid = make_solid_from_region(&cells, &[true, false, true])
            .expect("region mask should build a two-cell solid");
        assert!((volume(&solid) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn maker_volume_region_mask_reports_selection() {
        let volume_builder = MakerVolume::from_cells(vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ]);

        let selection = volume_builder
            .selection_from_region_mask(&[false, true, true])
            .expect("selection should succeed");
        assert_eq!(selection.input_cell_count, 3);
        assert_eq!(selection.selected_cell_indices, vec![1, 2]);
    }

    #[test]
    fn maker_volume_rejects_invalid_region_mask_length() {
        let volume_builder = MakerVolume::from_cells(vec![box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)]);
        let err = volume_builder
            .build_from_region_mask(&[true, false])
            .expect_err("mask length mismatch should fail");
        assert!(matches!(
            err,
            MakerVolumeError::RegionMaskLengthMismatch { expected: 1, got: 2 }
        ));
    }

    #[test]
    fn maker_volume_rejects_invalid_cell_index() {
        let volume_builder = MakerVolume::from_cells(vec![box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)]);
        let err = volume_builder
            .build_from_indices(&[3])
            .expect_err("invalid cell index should fail");
        assert!(matches!(
            err,
            MakerVolumeError::InvalidCellIndex { index: 3, count: 1 }
        ));
    }

    #[test]
    fn maker_volume_rejects_empty_selection() {
        let volume_builder = MakerVolume::from_cells(vec![box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)]);
        let err = volume_builder
            .build_from_indices(&[])
            .expect_err("empty selection should fail");
        assert!(matches!(err, MakerVolumeError::EmptySelection));
    }

    #[test]
    fn maker_volume_build_from_expr_uses_cells_builder() {
        let volume_builder = MakerVolume::from_cells(vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ]);
        let expr = CellExpr::Union(Box::new(CellExpr::Cell(0)), Box::new(CellExpr::Cell(1)));

        let solid = volume_builder
            .build_from_expr(&expr)
            .expect("expression-based build should succeed");
        assert!((volume(&solid) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn maker_volume_history_tracks_fold_steps() {
        let volume_builder = MakerVolume::from_cells(vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(4.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ]);

        let (_solid, history) = volume_builder
            .build_from_indices_with_history(&[0, 1, 2])
            .expect("history-enabled build should succeed");
        assert_eq!(history.steps.len(), 2);
    }

    #[test]
    fn maker_volume_deduplicates_repeated_indices() {
        let volume_builder = MakerVolume::from_cells(vec![
            box_at(0.0, 0.0, 0.0, 1.0, 1.0, 1.0),
            box_at(2.0, 0.0, 0.0, 1.0, 1.0, 1.0),
        ]);

        let solid = volume_builder
            .build_from_indices(&[0, 1, 1, 0])
            .expect("duplicate selection should still succeed");
        assert!((volume(&solid) - 2.0).abs() < 1e-9);
    }
}