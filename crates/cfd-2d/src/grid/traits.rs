//! Core traits for 2D grids

use super::boundary::BoundaryType;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};

/// Trait for 2D computational grids
pub trait Grid2D<T: RealField + Copy> {
    /// Get the number of cells in x direction
    fn nx(&self) -> usize;

    /// Get the number of cells in y direction
    fn ny(&self) -> usize;

    /// Get total number of cells
    fn num_cells(&self) -> usize {
        self.nx() * self.ny()
    }

    /// Get cell center coordinates
    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>>;

    /// Get cell volume/area
    fn cell_area(&self, i: usize, j: usize) -> Result<T>;

    /// Get neighboring cell indices (prefer `neighbor_iter` for better performance)
    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)>;

    /// Get neighboring cell indices as iterator (more efficient)
    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)>;

    /// Check if cell is on boundary
    fn is_boundary(&self, i: usize, j: usize) -> bool;

    /// Get boundary type for cell
    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType>;
}
