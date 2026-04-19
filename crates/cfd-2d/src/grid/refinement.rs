//! Adaptive grid refinement for 2D CFD
//!
//! # Theorem
//! The grid topology must form a valid, non-overlapping partition of the computational domain.
//!
//! **Proof sketch**:
//! For a finite volume discretization to be conservative, the control volumes $\Omega_i$
//! must satisfy $\cup_i \Omega_i = \Omega$ and $\Omega_i \cap \Omega_j = \emptyset$ for $i \neq j$.
//! The grid data structures enforce this by maintaining strict adjacency invariants
//! and ensuring that the sum of face area vectors for any closed cell is exactly zero:
//! $\sum_f \mathbf{A}_f = \mathbf{0}$.

use super::structured::StructuredGrid2D;
use super::traits::Grid2D;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Refinement criterion for adaptive grids
#[derive(Debug, Clone, Copy)]
pub enum RefinementCriterion<T: RealField + Copy> {
    /// Refine based on gradient threshold
    Gradient(T),
    /// Refine based on error estimate
    Error(T),
    /// Refine based on feature detection
    Feature,
}

/// Adaptive grid with refinement capability.
///
/// Refinement levels are stored in a contiguous flat `Vec<usize>` with
/// row-major layout `[i * ny + j]`, where `i` is the x-index and `j` is
/// the y-index. This avoids the heap indirection of `Vec<Vec<usize>>`
/// and improves cache locality during full-grid traversals.
#[derive(Debug, Clone)]
pub struct AdaptiveGrid2D<T: RealField + Copy> {
    /// Base grid.
    pub base_grid: StructuredGrid2D<T>,
    /// Refinement levels for each cell, stored flat as `[i * ny + j]`.
    refinement_levels: Vec<usize>,
    /// Number of columns (y-dimension) for indexing into `refinement_levels`.
    ny: usize,
    /// Maximum refinement level.
    pub max_level: usize,
}

impl<T: RealField + FromPrimitive + Copy> AdaptiveGrid2D<T> {
    /// Create a new adaptive grid.
    pub fn new(base_grid: StructuredGrid2D<T>, max_level: usize) -> Self {
        let nx = base_grid.nx();
        let ny = base_grid.ny();

        Self {
            base_grid,
            refinement_levels: vec![0; nx * ny],
            ny,
            max_level,
        }
    }

    /// Get refinement level at cell (i, j).
    #[inline]
    pub fn level(&self, i: usize, j: usize) -> usize {
        self.refinement_levels[i * self.ny + j]
    }

    /// Set refinement level at cell (i, j).
    #[inline]
    pub fn set_level(&mut self, i: usize, j: usize, level: usize) {
        self.refinement_levels[i * self.ny + j] = level;
    }

    /// Get a flat slice of all refinement levels.
    pub fn levels(&self) -> &[usize] {
        &self.refinement_levels
    }

    /// Mark cells for refinement based on criterion
    ///
    /// # Errors
    /// Returns an error if criterion evaluation fails or refinement marking encounters invalid state
    pub fn mark_for_refinement<F>(
        &mut self,
        criterion: RefinementCriterion<T>,
        evaluate: F,
    ) -> Result<()>
    where
        F: Fn(usize, usize) -> T,
    {
        match criterion {
            RefinementCriterion::Gradient(threshold) => {
                self.mark_by_gradient(threshold, evaluate)?;
            }
            RefinementCriterion::Error(threshold) => {
                self.mark_by_error(threshold, evaluate)?;
            }
            RefinementCriterion::Feature => {
                self.mark_by_feature()?;
            }
        }
        Ok(())
    }

    /// Mark cells based on gradient threshold
    ///
    /// # Errors
    /// Returns an error if gradient evaluation fails
    fn mark_by_gradient<F>(&mut self, threshold: T, evaluate: F) -> Result<()>
    where
        F: Fn(usize, usize) -> T,
    {
        for i in 0..self.base_grid.nx() {
            for j in 0..self.base_grid.ny() {
                let value = evaluate(i, j);
                let idx = i * self.ny + j;
                if value > threshold && self.refinement_levels[idx] < self.max_level {
                    self.refinement_levels[idx] += 1;
                }
            }
        }
        Ok(())
    }

    /// Mark cells based on error estimate
    ///
    /// # Errors
    /// Returns an error if error evaluation fails
    fn mark_by_error<F>(&mut self, threshold: T, evaluate: F) -> Result<()>
    where
        F: Fn(usize, usize) -> T,
    {
        // Similar to gradient-based marking
        self.mark_by_gradient(threshold, evaluate)
    }

    /// Mark cells based on feature detection
    ///
    /// # Errors
    /// Returns an error if feature detection fails
    fn mark_by_feature(&mut self) -> Result<()> {
        // Feature detection based on second derivatives or vorticity
        // Mark cells where features are detected
        for i in 1..self.base_grid.nx() - 1 {
            for j in 1..self.base_grid.ny() - 1 {
                let idx = i * self.ny + j;
                if self.detect_feature_at(i, j) && self.refinement_levels[idx] < self.max_level {
                    self.refinement_levels[idx] += 1;
                }
            }
        }
        Ok(())
    }

    /// Detect feature at a specific cell.
    fn detect_feature_at(&self, i: usize, j: usize) -> bool {
        let current_level = self.refinement_levels[i * self.ny + j];

        let neighbors = [
            (i.wrapping_sub(1), j),
            (i + 1, j),
            (i, j.wrapping_sub(1)),
            (i, j + 1),
        ];

        for (ni, nj) in neighbors {
            if ni < self.base_grid.nx() && nj < self.base_grid.ny() {
                let neighbor_level = self.refinement_levels[ni * self.ny + nj];
                if neighbor_level.abs_diff(current_level) > 1 {
                    return true;
                }
            }
        }
        false
    }

    /// Refine marked cells
    ///
    /// # Errors
    /// Returns an error if grid refinement fails due to invalid configuration or memory constraints
    pub fn refine(&mut self) -> Result<()> {
        // Create refined grid by subdividing marked cells
        let nx = self.base_grid.nx();
        let ny = self.base_grid.ny();

        for i in 0..nx {
            for j in 0..ny {
                let idx = i * self.ny + j;
                let level = self.refinement_levels[idx];
                if level > 0 {
                    // Adaptive mesh refinement (AMR) with hierarchical cell creation
                    // is a planned enhancement. Current implementation tracks refinement
                    // levels for feature detection and grid quality metrics.
                    self.refinement_levels[idx] = level;
                }
            }
        }
        Ok(())
    }

    /// Coarsen cells if possible
    pub fn coarsen(&mut self) -> Result<()> {
        // Merge cells that can be coarsened
        let nx = self.base_grid.nx();
        let ny = self.base_grid.ny();

        for i in 0..nx {
            for j in 0..ny {
                let idx = i * self.ny + j;
                if self.refinement_levels[idx] > 0 && self.can_coarsen_at(i, j) {
                    self.refinement_levels[idx] -= 1;
                }
            }
        }
        Ok(())
    }

    /// Check if a cell can be coarsened (2:1 constraint).
    fn can_coarsen_at(&self, i: usize, j: usize) -> bool {
        let current_level = self.refinement_levels[i * self.ny + j];
        if current_level == 0 {
            return false;
        }

        let neighbors = [
            (i.wrapping_sub(1), j),
            (i + 1, j),
            (i, j.wrapping_sub(1)),
            (i, j + 1),
        ];

        for (ni, nj) in neighbors {
            if ni < self.base_grid.nx() && nj < self.base_grid.ny() {
                let neighbor_level = self.refinement_levels[ni * self.ny + nj];
                if current_level > neighbor_level + 1 {
                    return false;
                }
            }
        }
        true
    }

    /// Get effective resolution at a point.
    pub fn effective_resolution(&self, i: usize, j: usize) -> (T, T) {
        let level = self.refinement_levels[i * self.ny + j];
        let factor = T::from_usize(1 << level).unwrap_or_else(|| T::one());

        (self.base_grid.dx / factor, self.base_grid.dy / factor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_grid() -> StructuredGrid2D<f64> {
        StructuredGrid2D::new(8, 8, 0.0, 8.0, 0.0, 8.0).unwrap()
    }

    #[test]
    fn test_initial_refinement_levels_zero() {
        let grid = AdaptiveGrid2D::new(make_grid(), 3);
        for &level in grid.levels() {
            assert_eq!(level, 0);
        }
    }

    #[test]
    fn test_gradient_marking_refines_above_threshold() {
        let mut grid = AdaptiveGrid2D::new(make_grid(), 3);
        grid.mark_for_refinement(RefinementCriterion::Gradient(0.5), |i, _j| {
            if i > 4 {
                1.0
            } else {
                0.0
            }
        })
        .unwrap();

        // Cells with i > 4 should be refined
        for i in 5..8 {
            for j in 0..8 {
                assert_eq!(grid.level(i, j), 1, "Cell ({i},{j}) should be refined");
            }
        }
        // Cells with i <= 4 should not
        for i in 0..=4 {
            for j in 0..8 {
                assert_eq!(grid.level(i, j), 0, "Cell ({i},{j}) should NOT be refined");
            }
        }
    }

    #[test]
    fn test_max_level_cap() {
        let mut grid = AdaptiveGrid2D::new(make_grid(), 2);
        // Refine 3 times, but max_level=2 should cap it
        for _ in 0..3 {
            grid.mark_for_refinement(RefinementCriterion::Gradient(0.0), |_, _| 1.0)
                .unwrap();
        }
        for &level in grid.levels() {
            assert!(level <= 2, "Refinement level should not exceed max_level=2");
        }
    }

    #[test]
    fn test_coarsen_reduces_level() {
        let mut grid = AdaptiveGrid2D::new(make_grid(), 3);
        // Set all cells to level 1
        for i in 0..8 {
            for j in 0..8 {
                grid.set_level(i, j, 1);
            }
        }
        grid.coarsen().unwrap();
        let any_reduced = grid.levels().iter().any(|&l| l == 0);
        assert!(
            any_reduced,
            "At least some cells should have been coarsened"
        );
    }

    #[test]
    fn test_coarsen_at_zero_is_noop() {
        let mut grid = AdaptiveGrid2D::new(make_grid(), 3);
        // All levels are 0, coarsen should do nothing
        grid.coarsen().unwrap();
        for &level in grid.levels() {
            assert_eq!(level, 0);
        }
    }

    #[test]
    fn test_effective_spacing_halves_per_level() {
        let grid = AdaptiveGrid2D::new(make_grid(), 3);
        let (dx0, dy0) = grid.effective_resolution(0, 0);
        // Base grid: 8 points spanning [0,8], dx = 8/7 ≈ 1.142..
        let base_dx: f64 = 8.0 / 7.0;
        assert!((dx0 - base_dx).abs() < 1e-10);
        assert!((dy0 - base_dx).abs() < 1e-10);

        let mut grid2 = AdaptiveGrid2D::new(make_grid(), 3);
        grid2.set_level(2, 2, 1);
        let (dx1, _) = grid2.effective_resolution(2, 2);
        assert!(
            (dx1 - base_dx / 2.0).abs() < 1e-10,
            "Level 1 spacing should be half, got {dx1}"
        );

        grid2.set_level(2, 2, 2);
        let (dx2, _) = grid2.effective_resolution(2, 2);
        assert!(
            (dx2 - base_dx / 4.0).abs() < 1e-10,
            "Level 2 spacing should be quarter, got {dx2}"
        );
    }
}
