//! 2D Bifurcation flow solver for branching microfluidic junctions
//!
//! This module implements a validated solver for bifurcating flows, which are
//! fundamental for modeling vascular networks and branching millifluidic systems.
//!
//! # Physics Background
//!
//! ## Mass Conservation at Junctions
//!
//! For an incompressible fluid, the total flow rate entering a junction must
//! equal the total flow rate exiting:
//!
//! ```text
//! Q_parent = Q_daughter1 + Q_daughter2
//! ```
//!
//! ## Murray's Law
//!
//! Murray's Law describes the optimal branching architecture in vascular systems:
//! ```text
//! r_p³ = r_d1³ + r_d2³
//! ```
//! where r is the radius (or half-width) of the parent and daughter branches.
//!
//! # Validation Strategy
//!
//! 1. **Mass Conservation**: Verify that flux out equals flux in.
//! 2. **Pressure Drop**: Compare against analytical solutions for branching networks.
//! 3. **Murray's Law**: Validate WSS distribution in optimal vs. non-optimal junctions.

use super::ns_fvm_2d::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Bifurcation Geometry Definition
// ============================================================================

/// Bifurcation geometry configuration for 2D branching channels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationGeometry<T: RealField + Copy> {
    /// Parent channel width [m]
    pub parent_width: T,
    /// Parent channel length [m]
    pub parent_length: T,
    /// Daughter 1 width [m]
    pub daughter1_width: T,
    /// Daughter 1 length [m]
    pub daughter1_length: T,
    /// Daughter 1 angle [radians from x-axis]
    pub daughter1_angle: T,
    /// Daughter 2 width [m]
    pub daughter2_width: T,
    /// Daughter 2 length [m]
    pub daughter2_length: T,
    /// Daughter 2 angle [radians from x-axis]
    pub daughter2_angle: T,
}

impl<T: RealField + Copy + FromPrimitive> BifurcationGeometry<T> {
    /// Create a new symmetric bifurcation
    pub fn new_symmetric(
        parent_width: T,
        parent_length: T,
        daughter_width: T,
        daughter_length: T,
        angle: T,
    ) -> Self {
        Self {
            parent_width,
            parent_length,
            daughter1_width: daughter_width,
            daughter1_length: daughter_length,
            daughter1_angle: angle,
            daughter2_width: daughter_width,
            daughter2_length: daughter_length,
            daughter2_angle: -angle,
        }
    }

    /// Check if a point (x, y) is within the fluid domain
    pub fn contains(&self, x: T, y: T) -> bool {
        // Parent branch: Horizontal from x=0 to parent_length, centered at y=0
        let half_pw = self.parent_width / T::from_f64(2.0).unwrap();
        if x >= T::zero() && x <= self.parent_length && y >= -half_pw && y <= half_pw {
            return true;
        }

        // Daughter 1
        if self.in_segment(x, y, self.parent_length, T::zero(), self.daughter1_angle, self.daughter1_length, self.daughter1_width) {
            return true;
        }

        // Daughter 2
        if self.in_segment(x, y, self.parent_length, T::zero(), self.daughter2_angle, self.daughter2_length, self.daughter2_width) {
            return true;
        }

        false
    }

    fn in_segment(
        &self,
        x: T,
        y: T,
        start_x: T,
        start_y: T,
        angle: T,
        length: T,
        width: T,
    ) -> bool {
        let dx = x - start_x;
        let dy = y - start_y;

        let cos_a = angle.cos();
        let sin_a = angle.sin();

        // Local coordinates (aligned with branch)
        let lx = dx * cos_a + dy * sin_a;
        let ly = -dx * sin_a + dy * cos_a;

        let half_w = width / T::from_f64(2.0).unwrap();

        lx >= T::zero() && lx <= length && ly >= -half_w && ly <= half_w
    }

    /// Get bounding box [min_x, max_x, min_y, max_y]
    pub fn bounding_box(&self) -> [T; 4] {
        let d1_end_x = self.parent_length + self.daughter1_length * self.daughter1_angle.cos();
        let d1_end_y = self.daughter1_length * self.daughter1_angle.sin();
        let d2_end_x = self.parent_length + self.daughter2_length * self.daughter2_angle.cos();
        let d2_end_y = self.daughter2_length * self.daughter2_angle.sin();

        let half_pw = self.parent_width / T::from_f64(2.0).unwrap();
        let half_d1w = self.daughter1_width / T::from_f64(2.0).unwrap();
        let half_d2w = self.daughter2_width / T::from_f64(2.0).unwrap();

        let min_x = T::zero();
        let max_x = d1_end_x.max(d2_end_x).max(self.parent_length);
        
        let min_y = d1_end_y.min(d2_end_y).min(-half_pw).min(d1_end_y - half_d1w).min(d2_end_y - half_d2w);
        let max_y = d1_end_y.max(d2_end_y).max(half_pw).max(d1_end_y + half_d1w).max(d2_end_y + half_d2w);

        [min_x, max_x, min_y, max_y]
    }
}

// ============================================================================
// Bifurcation Solver
// ============================================================================

/// Solver for branching flow junctions
pub struct BifurcationSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    pub geometry: BifurcationGeometry<T>,
    pub ns_solver: NavierStokesSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> BifurcationSolver2D<T> {
    /// Create a new bifurcation solver
    pub fn new(
        geometry: BifurcationGeometry<T>,
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

        // Populate mask
        for i in 0..nx {
            for j in 0..ny {
                let x = ns_solver.grid.x_center(i) + bbox[0];
                let y = ns_solver.grid.y_center(j) + bbox[2];
                ns_solver.field.mask[i][j] = geometry.contains(x, y);
            }
        }

        Self {
            geometry,
            ns_solver,
        }
    }

    /// Solve the bifurcation flow
    pub fn solve(&mut self, u_inlet: T) -> CfdResult<BifurcationSolution<T>> {
        let _solve_res = self.ns_solver.solve(u_inlet).map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        // Extract outlet flow rates
        let nx = self.ns_solver.grid.nx;
        let ny = self.ns_solver.grid.ny;
        let dy = self.ns_solver.grid.dy;

        let mut q_d1 = T::zero();
        let mut q_d2 = T::zero();

        // Find mid-y of outlet regions
        let mut mid_y = T::zero();
        let mut y_coords = Vec::new();
        let mut n_outlet_mask = 0;
        for j in 0..ny {
            if self.ns_solver.field.mask[nx-1][j] {
                y_coords.push(self.ns_solver.grid.y_center(j));
                n_outlet_mask += 1;
            }
        }

        if !y_coords.is_empty() {
            let y_min = y_coords.iter().cloned().fold(y_coords[0], Float::min);
            let y_max = y_coords.iter().cloned().fold(y_coords[0], Float::max);
            mid_y = (y_min + y_max) / T::from_f64(2.0).unwrap();
        }

        for j in 0..ny {
            if self.ns_solver.field.mask[nx-1][j] {
                let uj = self.ns_solver.field.u[nx][j];
                let q_local = uj * dy;
                if self.ns_solver.grid.y_center(j) > mid_y {
                    q_d1 = q_d1 + q_local;
                } else {
                    q_d2 = q_d2 + q_local;
                }
            }
        }

        let mut q_in_actual = T::zero();
        for j in 0..ny {
            if self.ns_solver.field.mask[0][j] {
                q_in_actual = q_in_actual + self.ns_solver.field.u[0][j] * dy;
            }
        }
        let q_parent = q_in_actual;

        Ok(BifurcationSolution {
            q_parent,
            q_daughter1: q_d1,
            q_daughter2: q_d2,
            mass_balance_error: Float::abs(q_parent - (q_d1 + q_d2)) / q_parent,
        })
    }
}

/// Results from a bifurcation flow simulation
#[derive(Debug, Serialize, Deserialize)]
pub struct BifurcationSolution<T: RealField + Copy> {
    pub q_parent: T,
    pub q_daughter1: T,
    pub q_daughter2: T,
    pub mass_balance_error: T,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_symmetric_bifurcation_mass_conservation() {
        let geom = BifurcationGeometry::new_symmetric(
            0.001, // 1mm parent
            0.002, // 2mm parent length
            0.0007, // ~0.7mm daughters
            0.002, // 2mm daughters length
            0.5,   // ~28 degrees
        );

        let blood = BloodModel::Newtonian(0.0035);
        let density = 1060.0;
        let nx = 50;
        let ny = 30;

        let mut config = SIMPLEConfig::default();
        config.max_iterations = 5000;
        config.tolerance = 1e-5;
        config.alpha_u = 0.5;
        config.alpha_p = 0.2;
        
        let mut solver = BifurcationSolver2D::new(geom, blood, density, nx, ny, config);
        let result = solver.solve(0.1);

        assert!(result.is_ok());
        let sol = result.unwrap();
        
        println!("Bifurcation Mass Balance Error: {:?}", sol.mass_balance_error);
        // Error should be small if converging
        assert!(sol.mass_balance_error < 0.05);

        // Symmetric case: daughters should have roughly equal flow
        let flow_diff = Float::abs(sol.q_daughter1 - sol.q_daughter2) / sol.q_parent;
        assert!(flow_diff < 0.05);
    }

    #[test]
    fn test_non_newtonian_bifurcation_mass_conservation() {
        let mu_inf = 0.0035;    // Pa·s
        let yield_stress = 0.005; // Pa
        let hematocrit = 0.45;
        let density = 1060.0;   // kg/m³
        let casson = CassonBlood::new(density, yield_stress, mu_inf, hematocrit);
        let blood = BloodModel::Casson(casson);

        let geom = BifurcationGeometry {
            parent_width: 0.001,
            parent_length: 0.002,
            daughter1_width: 0.0007,
            daughter1_length: 0.002,
            daughter1_angle: 0.5,
            daughter2_width: 0.0007,
            daughter2_length: 0.002,
            daughter2_angle: -0.5,
        };

        let mut config = SIMPLEConfig::default();
        config.max_iterations = 2000;
        config.tolerance = 1e-4;
        config.alpha_u = 0.5;
        config.alpha_p = 0.2;

        let mut solver = BifurcationSolver2D::new(geom, blood, density, 50, 30, config);
        let u_inlet = 0.1;
        
        let sol = solver.solve(u_inlet).unwrap();
        
        println!("Non-Newtonian Bifurcation Q_parent: {:?}", sol.q_parent);
        println!("Non-Newtonian Bifurcation Q_d1: {:?}, Q_d2: {:?}", sol.q_daughter1, sol.q_daughter2);
        println!("Non-Newtonian Mass Balance Error: {:?}", sol.mass_balance_error);

        // Mass balance should still be perfect regardless of the viscosity model
        assert!(sol.mass_balance_error < 0.01);
    }
}
