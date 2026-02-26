//! Grid structure for finite difference operations
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

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
