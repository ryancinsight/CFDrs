//! Grid structure for finite difference operations

use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;
/// 2D grid for finite difference operations
#[derive(Debug, Clone)]
pub struct Grid2D<T: RealField + Copy> {
    /// Grid values
    pub data: DMatrix<T>,
    /// Grid spacing in x-direction
    pub dx: T,
    /// Grid spacing in y-direction
    pub dy: T,
    /// Number of ghost cells
    pub ghost_cells: usize,
}
impl<T: RealField + Copy + FromPrimitive + Copy> Grid2D<T> {
    /// Create a new 2D grid
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, ghost_cells: usize) -> Self {
        let total_nx = nx + 2 * ghost_cells;
        let total_ny = ny + 2 * ghost_cells;
        Self {
            data: DMatrix::zeros(total_nx, total_ny),
            dx,
            dy,
            ghost_cells,
        }
    }
    /// Get interior dimensions (excluding ghost cells)
    pub fn interior_shape(&self) -> (usize, usize) {
        let (total_nx, total_ny) = self.data.shape();
        (
            total_nx - 2 * self.ghost_cells,
            total_ny - 2 * self.ghost_cells,
        )

    }


}
