//! Staggered grid geometry for finite volume methods.
//!
//! ## Mathematical Invariant
//! For an (nx × ny) staggered grid with spacing (dx, dy):
//! - U-velocity nodes: (nx+1) × ny faces at x = i·dx, y = (j+0.5)·dy
//! - V-velocity nodes: nx × (ny+1) faces at x = (i+0.5)·dx, y = j·dy
//! - Pressure / scalar nodes: nx × ny cells at x = (i+0.5)·dx, y = (j+0.5)·dy
//!
//! Supports non-uniform y-spacing via optional `y_faces` array for resolving
//! extreme contraction ratios (e.g. 45 µm throat in 4 mm channel).
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
///
/// When `y_faces` is `Some`, the grid is non-uniform in y — use [`dy_at`],
/// [`y_center`], and [`y_v_face`] which dispatch correctly for both cases.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StaggeredGrid2D<T: RealField + Copy> {
    /// Number of pressure cells in x
    pub nx: usize,
    /// Number of pressure cells in y
    pub ny: usize,
    /// Uniform cell width [m]
    pub dx: T,
    /// Uniform (or reference) cell height [m].
    /// For non-uniform grids, use [`dy_at`] for the per-cell height.
    pub dy: T,
    /// Physical domain width (nx * dx) [m]
    pub lx: T,
    /// Physical domain height [m]
    pub ly: T,
    /// Non-uniform y-face coordinates (ny+1 entries, from 0 to ly).
    /// When `None`, the grid is uniform in y and `dy` is used everywhere.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y_faces: Option<Vec<T>>,
}

impl<T: RealField + Copy + FromPrimitive> StaggeredGrid2D<T> {
    /// Construct a uniform grid from cell counts and total domain size.
    ///
    /// # Panics
    /// Panics if `nx == 0`, `ny == 0`, or if `T` cannot represent the cell spacing.
    #[must_use]
    pub fn new(nx: usize, ny: usize, lx: T, ly: T) -> Self {
        assert!(nx > 0 && ny > 0, "Grid dimensions must be positive");
        let dx = lx / T::from_usize(nx).expect("nx fits in T");
        let dy = ly / T::from_usize(ny).expect("ny fits in T");
        Self { nx, ny, dx, dy, lx, ly, y_faces: None }
    }

    /// Construct a grid with non-uniform y-spacing from explicit face coordinates.
    ///
    /// `y_face_coords` must have exactly `ny + 1` entries, monotonically increasing
    /// from 0 to `ly`.
    ///
    /// # Panics
    /// Panics if dimensions are invalid or y_face_coords has wrong length.
    #[must_use]
    pub fn new_stretched_y(nx: usize, ny: usize, lx: T, y_face_coords: Vec<T>) -> Self {
        assert!(nx > 0 && ny > 0, "Grid dimensions must be positive");
        assert_eq!(
            y_face_coords.len(),
            ny + 1,
            "y_face_coords must have ny+1 = {} entries, got {}",
            ny + 1,
            y_face_coords.len()
        );
        let dx = lx / T::from_usize(nx).expect("nx fits in T");
        let ly = y_face_coords[ny];
        let dy = ly / T::from_usize(ny).expect("ny fits in T"); // reference (average) dy
        Self {
            nx,
            ny,
            dx,
            dy,
            lx,
            ly,
            y_faces: Some(y_face_coords),
        }
    }

    // ── Per-cell y-spacing ──────────────────────────────────────────────────

    /// Height of pressure cell j.  Returns `y_faces[j+1] - y_faces[j]` for
    /// non-uniform grids, or the uniform `dy` otherwise.
    #[inline]
    pub fn dy_at(&self, j: usize) -> T {
        match &self.y_faces {
            Some(yf) => yf[j + 1] - yf[j],
            None => self.dy,
        }
    }

    /// Distance between the centres of pressure cells j and j+1.
    /// For uniform grids this equals `dy`; for non-uniform grids it is
    /// `(dy[j] + dy[j+1]) / 2`.
    #[inline]
    pub fn dy_face(&self, j: usize) -> T {
        let half = T::from_f64(0.5).expect("0.5 fits in T");
        match &self.y_faces {
            Some(yf) => {
                let dy_j = yf[j + 1] - yf[j];
                let dy_jp1 = yf[j + 2] - yf[j + 1];
                (dy_j + dy_jp1) * half
            }
            None => self.dy,
        }
    }

    // ── Pressure cell centres ──────────────────────────────────────────────

    /// x-coordinate of pressure cell centre (i, j)
    #[inline]
    pub fn x_center(&self, i: usize) -> T {
        let half = T::from_f64(0.5).expect("0.5 fits in T");
        (T::from_usize(i).expect("i fits in T") + half) * self.dx
    }

    /// y-coordinate of pressure cell centre at index j.
    /// For non-uniform grids, returns `(y_faces[j] + y_faces[j+1]) / 2`.
    #[inline]
    pub fn y_center(&self, j: usize) -> T {
        let half = T::from_f64(0.5).expect("0.5 fits in T");
        match &self.y_faces {
            Some(yf) => (yf[j] + yf[j + 1]) * half,
            None => (T::from_usize(j).expect("j fits in T") + half) * self.dy,
        }
    }

    // ── U-face (east/west) x-coordinates ──────────────────────────────────

    /// x-coordinate of the U-velocity face between cells (i-1, j) and (i, j)
    #[inline]
    pub fn x_u_face(&self, i: usize) -> T {
        T::from_usize(i).expect("i fits in T") * self.dx
    }

    // ── V-face (north/south) y-coordinates ────────────────────────────────

    /// y-coordinate of the V-velocity face between cells (i, j-1) and (i, j).
    /// For non-uniform grids, returns `y_faces[j]`.
    #[inline]
    pub fn y_v_face(&self, j: usize) -> T {
        match &self.y_faces {
            Some(yf) => yf[j],
            None => T::from_usize(j).expect("j fits in T") * self.dy,
        }
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

    /// True if this grid uses non-uniform y-spacing.
    #[inline]
    pub fn is_stretched_y(&self) -> bool {
        self.y_faces.is_some()
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

    #[test]
    fn uniform_dy_at_matches_dy() {
        let g = StaggeredGrid2D::<f64>::new(4, 4, 1.0, 1.0);
        for j in 0..4 {
            assert_relative_eq!(g.dy_at(j), g.dy, epsilon = 1e-12);
        }
    }

    #[test]
    fn stretched_grid_preserves_domain_height() {
        let yf = vec![0.0, 0.1, 0.3, 0.6, 1.0];
        let g = StaggeredGrid2D::<f64>::new_stretched_y(4, 4, 1.0, yf);
        assert_relative_eq!(g.ly, 1.0, epsilon = 1e-12);
        assert_relative_eq!(g.dy_at(0), 0.1, epsilon = 1e-12);
        assert_relative_eq!(g.dy_at(1), 0.2, epsilon = 1e-12);
        assert_relative_eq!(g.dy_at(2), 0.3, epsilon = 1e-12);
        assert_relative_eq!(g.dy_at(3), 0.4, epsilon = 1e-12);
    }

    #[test]
    fn stretched_y_center_is_midpoint() {
        let yf = vec![0.0, 0.1, 0.3, 0.6, 1.0];
        let g = StaggeredGrid2D::<f64>::new_stretched_y(4, 4, 1.0, yf);
        assert_relative_eq!(g.y_center(0), 0.05, epsilon = 1e-12);
        assert_relative_eq!(g.y_center(1), 0.20, epsilon = 1e-12);
        assert_relative_eq!(g.y_center(2), 0.45, epsilon = 1e-12);
        assert_relative_eq!(g.y_center(3), 0.80, epsilon = 1e-12);
    }
}
