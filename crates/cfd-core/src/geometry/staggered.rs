//! Staggered grid geometry for finite volume methods.
//!
//! ## Mathematical Invariant
//! For an (nx × ny) staggered grid with spacing (dx, dy):
//! - U-velocity nodes: (nx+1) × ny faces at x = i·dx, y = (j+0.5)·dy
//! - V-velocity nodes: nx × (ny+1) faces at x = (i+0.5)·dx, y = j·dy
//! - Pressure / scalar nodes: nx × ny cells at x = (i+0.5)·dx, y = (j+0.5)·dy
//!
//! This layout satisfies the discrete inf-sup (LBB) condition and prevents
//! checkerboard pressure oscillations without Rhie-Chow correction on the grid itself.
//!
//! ## Reference
//! Harlow, F.H. & Welch, J.E. (1965). Numerical Calculation of Time-Dependent
//! Viscous Incompressible Flow of Fluid with Free Surface. *Physics of Fluids*, 8(12), 2182.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Staggered (MAC) grid descriptor for 2D finite volume methods.
///
/// Stores only the topology and spacing; field data lives in the solver.
/// This is a pure geometry primitive and carries no solver state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StaggeredGrid2D<T: RealField + Copy> {
    /// Number of pressure cells in x
    pub nx: usize,
    /// Number of pressure cells in y
    pub ny: usize,
    /// Uniform cell width [m]
    pub dx: T,
    /// Uniform cell height [m]
    pub dy: T,
    /// Physical domain width (nx * dx) [m]
    pub lx: T,
    /// Physical domain height (ny * dy) [m]
    pub ly: T,
}

impl<T: RealField + Copy + FromPrimitive> StaggeredGrid2D<T> {
    /// Construct from cell counts and total domain size.
    ///
    /// # Panics
    /// Panics if `nx == 0`, `ny == 0`, or if `T` cannot represent the cell spacing.
    #[must_use]
    pub fn new(nx: usize, ny: usize, lx: T, ly: T) -> Self {
        assert!(nx > 0 && ny > 0, "Grid dimensions must be positive");
        let dx = lx / T::from_usize(nx).expect("nx fits in T");
        let dy = ly / T::from_usize(ny).expect("ny fits in T");
        Self { nx, ny, dx, dy, lx, ly }
    }

    // ── Pressure cell centres ──────────────────────────────────────────────

    /// x-coordinate of pressure cell centre (i, j)
    #[inline]
    pub fn x_center(&self, i: usize) -> T {
        let half = T::from_f64(0.5).expect("0.5 fits in T");
        (T::from_usize(i).expect("i fits in T") + half) * self.dx
    }

    /// y-coordinate of pressure cell centre (i, j)
    #[inline]
    pub fn y_center(&self, j: usize) -> T {
        let half = T::from_f64(0.5).expect("0.5 fits in T");
        (T::from_usize(j).expect("j fits in T") + half) * self.dy
    }

    // ── U-face (east/west) x-coordinates ──────────────────────────────────

    /// x-coordinate of the U-velocity face between cells (i-1, j) and (i, j)
    #[inline]
    pub fn x_u_face(&self, i: usize) -> T {
        T::from_usize(i).expect("i fits in T") * self.dx
    }

    // ── V-face (north/south) y-coordinates ────────────────────────────────

    /// y-coordinate of the V-velocity face between cells (i, j-1) and (i, j)
    #[inline]
    pub fn y_v_face(&self, j: usize) -> T {
        T::from_usize(j).expect("j fits in T") * self.dy
    }

    // ── Convenience ───────────────────────────────────────────────────────

    /// Cell-centred index: `j * nx + i`
    #[inline]
    pub fn cell_idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Total number of pressure cells
    #[inline]
    pub fn n_cells(&self) -> usize {
        self.nx * self.ny
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn cell_centres_are_staggered_from_faces() {
        let g = StaggeredGrid2D::<f64>::new(4, 4, 1.0, 1.0);
        // First cell centre at 0.5 * dx = 0.125
        assert_relative_eq!(g.x_center(0), 0.125, epsilon = 1e-12);
        // U-face at i=0 is at x=0.0 (left boundary)
        assert_relative_eq!(g.x_u_face(0), 0.0, epsilon = 1e-12);
        // U-face at i=1 is at x=0.25 (between cells 0 and 1)
        assert_relative_eq!(g.x_u_face(1), 0.25, epsilon = 1e-12);
    }

    #[test]
    fn n_cells_matches_dimensions() {
        let g = StaggeredGrid2D::<f64>::new(8, 6, 2.0, 1.5);
        assert_eq!(g.n_cells(), 48);
    }
}
