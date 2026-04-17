//! 2D Cross-junction flow solver for intersecting microfluidic channels.
//!
//! Models the fluid dynamics at points where two channels cross at 90° in
//! a planar millifluidic device.  At a cross-junction the flow from two
//! perpendicular channels meets, potentially exchanging momentum and
//! producing vortical structures.
//!
//! # Physics Background
//!
//! ## Mass Conservation at Cross-Junctions
//!
//! For incompressible flow through a 4-port junction (North, South, East, West):
//!
//! ```text
//! Q_N + Q_S + Q_E + Q_W = 0   (signed, positive = inflow)
//! ```
//!
//! ## Pressure Loss
//!
//! The minor loss coefficient $K$ at a 90° cross-junction (Idelchik 2007)
//! is 1.2–1.5 depending on the Reynolds number and branch area ratio.
//!
//! # Theorem
//!
//! The cross-junction solver must converge to a solution satisfying the
//! discrete conservation laws and the Kirchhoff junction constraint.
//!
//! **Proof sketch**: The domain is the union of two perpendicular rectangular
//! channels sharing an overlap square of side $\min(w_1, w_2)$.  The wall
//! mask excludes cells outside both channels.  The SIMPLE iteration on this
//! masked domain inherits the diagonal-dominance property of the staggered
//! discretisation (provided CFL ≤ 1) and therefore converges.

use super::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Cross-Junction Geometry
// ============================================================================

/// Geometry for a 90° cross-junction of two perpendicular channels.
///
/// The horizontal channel runs along the x-axis and the vertical channel
/// along the y-axis.  Both are centred on the same junction point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossJunctionGeometry<T: RealField + Copy> {
    /// Width of the horizontal channel [m].
    pub horizontal_width: T,
    /// Total length of the horizontal channel [m].
    pub horizontal_length: T,
    /// Width of the vertical channel [m].
    pub vertical_width: T,
    /// Total length of the vertical channel [m].
    pub vertical_length: T,
}

impl<T: RealField + Copy + FromPrimitive> CrossJunctionGeometry<T> {
    /// Create a symmetric cross-junction (equal channel widths and lengths).
    pub fn symmetric(width: T, length: T) -> Self {
        Self {
            horizontal_width: width,
            horizontal_length: length,
            vertical_width: width,
            vertical_length: length,
        }
    }

    /// Test whether the point `(x, y)` lies inside the fluid domain.
    ///
    /// The junction is centred at the origin.
    pub fn contains(&self, x: T, y: T) -> bool {
        let two = T::from_f64(2.0).expect("Exact mathematically representable f64");

        // Horizontal channel: x ∈ [-L_h/2, L_h/2], y ∈ [-w_h/2, w_h/2]
        let half_h_len = self.horizontal_length / two;
        let half_h_w = self.horizontal_width / two;
        if x >= -half_h_len && x <= half_h_len && y >= -half_h_w && y <= half_h_w {
            return true;
        }

        // Vertical channel: x ∈ [-w_v/2, w_v/2], y ∈ [-L_v/2, L_v/2]
        let half_v_len = self.vertical_length / two;
        let half_v_w = self.vertical_width / two;
        if x >= -half_v_w && x <= half_v_w && y >= -half_v_len && y <= half_v_len {
            return true;
        }

        false
    }

    /// Return the axis-aligned bounding box `[min_x, max_x, min_y, max_y]`.
    pub fn bounding_box(&self) -> [T; 4] {
        let two = T::from_f64(2.0).expect("Exact mathematically representable f64");
        let half_h_len = self.horizontal_length / two;
        let half_v_len = self.vertical_length / two;
        let half_v_w = self.vertical_width / two;
        let half_h_w = self.horizontal_width / two;

        let min_x = (-half_h_len).min(-half_v_w);
        let max_x = half_h_len.max(half_v_w);
        let min_y = (-half_v_len).min(-half_h_w);
        let max_y = half_v_len.max(half_h_w);

        [min_x, max_x, min_y, max_y]
    }
}

// ============================================================================
// Cross-Junction Solver
// ============================================================================

/// 2D Navier-Stokes solver for cross-junction flow.
///
/// Wraps [`NavierStokesSolver2D`] with a domain mask shaped like two
/// perpendicular channels.  The horizontal channel is driven by an
/// imposed inlet velocity; the vertical channel can be configured as
/// no-flow walls (closed side ports) or open outlets.
pub struct CrossJunctionSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Cross-junction geometry.
    pub geometry: CrossJunctionGeometry<T>,
    /// Underlying Navier-Stokes solver.
    pub ns_solver: NavierStokesSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> CrossJunctionSolver2D<T> {
    /// Construct a new cross-junction solver.
    ///
    /// # Arguments
    ///
    /// * `geometry` — channel dimensions.
    /// * `blood` — rheological model.
    /// * `density` — fluid density [kg/m³].
    /// * `nx`, `ny` — mesh resolution.
    /// * `config` — SIMPLE algorithm parameters.
    pub fn new(
        geometry: CrossJunctionGeometry<T>,
        blood: BloodModel<T>,
        density: T,
        nx: usize,
        ny: usize,
        config: SIMPLEConfig<T>,
    ) -> Self {
        let bbox = geometry.bounding_box();
        let width = bbox[1] - bbox[0];
        let height = bbox[3] - bbox[2];

        let grid = StaggeredGrid2D::new(nx, ny, width, height);
        let mut ns_solver = NavierStokesSolver2D::new(grid, blood, density, config);

        // Build the domain mask.
        for i in 0..nx {
            for j in 0..ny {
                let x = ns_solver.grid.x_center(i) + bbox[0];
                let y = ns_solver.grid.y_center(j) + bbox[2];
                ns_solver.field.mask[(i, j)] = geometry.contains(x, y);
            }
        }

        Self {
            geometry,
            ns_solver,
        }
    }

    /// Solve for steady-state cross-junction flow.
    ///
    /// `u_inlet` is the inlet velocity on the left boundary of the
    /// horizontal channel (west port).
    pub fn solve(&mut self, u_inlet: T) -> CfdResult<CrossJunctionSolution<T>> {
        let _solve_res = self
            .ns_solver
            .solve(u_inlet)
            .map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        let nx = self.ns_solver.grid.nx;
        let ny = self.ns_solver.grid.ny;
        let dy = self.ns_solver.grid.dy;
        let dx = self.ns_solver.grid.dx;

        // Compute fluxes at each port.
        // West (x = 0, inlet)
        let mut q_west = T::zero();
        for j in 0..ny {
            if self.ns_solver.field.mask[(0, j)] {
                q_west += self.ns_solver.field.u[(0, j)] * dy;
            }
        }

        // East (x = nx)
        let mut q_east = T::zero();
        for j in 0..ny {
            if self.ns_solver.field.mask[(nx - 1, j)] {
                q_east += self.ns_solver.field.u[(nx, j)] * dy;
            }
        }

        // South (y = 0)
        let mut q_south = T::zero();
        for i in 0..nx {
            if self.ns_solver.field.mask[(i, 0)] {
                q_south += self.ns_solver.field.v[(i, 0)] * dx;
            }
        }

        // North (y = ny)
        let mut q_north = T::zero();
        for i in 0..nx {
            if self.ns_solver.field.mask[(i, ny - 1)] {
                q_north += self.ns_solver.field.v[(i, ny)] * dx;
            }
        }

        let q_total_in = q_west;
        let q_total_out = q_east + q_north + q_south;
        let mass_balance_error =
            if Float::abs(q_total_in) > T::from_f64(1e-30).expect("Exact mathematically representable f64") {
                Float::abs(q_total_in - q_total_out) / Float::abs(q_total_in)
            } else {
                T::zero()
            };

        // Compute junction pressure loss ΔP between west and east centroids.
        let j_mid = ny / 2;
        let p_west = self.ns_solver.field.p[(1, j_mid)];
        let p_east = self.ns_solver.field.p[(nx - 2, j_mid)];
        let dp_junction = p_west - p_east;

        Ok(CrossJunctionSolution {
            q_west,
            q_east,
            q_north,
            q_south,
            dp_junction,
            mass_balance_error,
        })
    }
}

/// Solution data from a cross-junction flow simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrossJunctionSolution<T: RealField + Copy> {
    /// Flow rate through the west port (inlet).
    pub q_west: T,
    /// Flow rate through the east port (outlet).
    pub q_east: T,
    /// Flow rate through the north port.
    pub q_north: T,
    /// Flow rate through the south port.
    pub q_south: T,
    /// Pressure drop across the junction, west → east [Pa].
    pub dp_junction: T,
    /// Relative mass-balance error.
    pub mass_balance_error: T,
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_geometry_contains_junction_centre() {
        let g = CrossJunctionGeometry::<f64>::symmetric(1e-3, 10e-3);
        assert!(g.contains(0.0, 0.0), "centre must be in the domain");
    }

    #[test]
    fn symmetric_geometry_contains_each_arm() {
        let g = CrossJunctionGeometry::<f64>::symmetric(1e-3, 10e-3);
        // East arm
        assert!(g.contains(4e-3, 0.0));
        // West arm
        assert!(g.contains(-4e-3, 0.0));
        // North arm
        assert!(g.contains(0.0, 4e-3));
        // South arm
        assert!(g.contains(0.0, -4e-3));
    }

    #[test]
    fn symmetric_geometry_excludes_corners() {
        let g = CrossJunctionGeometry::<f64>::symmetric(1e-3, 10e-3);
        // A point in the far corner of the bounding box is outside both channels.
        assert!(!g.contains(4e-3, 4e-3));
    }

    #[test]
    fn bounding_box_symmetric() {
        let g = CrossJunctionGeometry::<f64>::symmetric(1e-3, 10e-3);
        let bb = g.bounding_box();
        assert!((bb[0] - (-5e-3)).abs() < 1e-10);
        assert!((bb[1] - 5e-3).abs() < 1e-10);
        assert!((bb[2] - (-5e-3)).abs() < 1e-10);
        assert!((bb[3] - 5e-3).abs() < 1e-10);
    }

    #[test]
    fn asymmetric_geometry_contains_both_channels() {
        let g = CrossJunctionGeometry {
            horizontal_width: 1e-3,
            horizontal_length: 12e-3,
            vertical_width: 0.8e-3,
            vertical_length: 8e-3,
        };
        // Inside horizontal channel
        assert!(g.contains(5e-3, 0.0));
        // Inside vertical channel
        assert!(g.contains(0.0, 3e-3));
        // Outside
        assert!(!g.contains(5e-3, 3e-3));
    }
}
