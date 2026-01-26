//! SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
//!
//! This module implements the SIMPLE algorithm for pressure-velocity coupling
//! in incompressible CFD.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::{MomentumComponent, MomentumSolver};
use crate::solvers::fdm::PoissonSolver;
use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{BiCGSTAB, IterativeLinearSolver, IterativeSolverConfig};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use std::collections::HashMap;

/// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
///
/// This is the standard algorithm for solving incompressible Navier-Stokes equations.
/// It uses a predictor-corrector approach to handle the pressure-velocity coupling.
pub struct SimpleAlgorithm<T: RealField + Copy + FromPrimitive + std::fmt::Debug> {
    /// Under-relaxation factor for pressure (typically 0.1-0.8)
    pressure_relaxation: T,
    /// Under-relaxation factor for velocity (typically 0.5-0.9)
    velocity_relaxation: T,
    /// Maximum number of SIMPLE iterations per time step
    max_iterations: usize,
    /// Convergence tolerance for continuity residual
    tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
    /// Create new SIMPLE algorithm with default parameters
    pub fn new() -> Self {
        Self {
            pressure_relaxation: T::from_f64(0.3).unwrap(), // Standard value
            velocity_relaxation: T::from_f64(0.7).unwrap(), // Standard value
            max_iterations: 50,
            tolerance: T::from_f64(1e-6).unwrap(),
        }
    }

    /// Set pressure under-relaxation factor
    pub fn with_pressure_relaxation(mut self, alpha_p: T) -> Self {
        self.pressure_relaxation = alpha_p;
        self
    }

    /// Set velocity under-relaxation factor
    pub fn with_velocity_relaxation(mut self, alpha_u: T) -> Self {
        self.velocity_relaxation = alpha_u;
        self
    }

    /// Set maximum iterations
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iterations = max_iter;
        self
    }

    /// Set convergence tolerance
    pub fn with_tolerance(mut self, tol: T) -> Self {
        self.tolerance = tol;
        self
    }

    /// Execute one simplified SIMPLE iteration
    ///
    /// This method performs the following steps:
    /// 1. Predictor: Solve momentum equations for u* and v* using current pressure.
    /// 2. Pressure Correction: Solve for p' to enforce continuity.
    /// 3. Correction: Update p and u, v using p'.
    ///
    /// Returns the maximum continuity residual and whether convergence was achieved
    pub fn simple_iteration(
        &self,
        momentum_solver: &mut MomentumSolver<T>,
        _poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        _boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<(T, bool)> {
        let nx = grid.nx;
        let ny = grid.ny;

        // Ensure momentum solver has correct relaxation
        momentum_solver.set_velocity_relaxation(self.velocity_relaxation);

        // 1. Momentum Predictor Step
        // Solve for u* and v* using current pressure field
        momentum_solver.solve(MomentumComponent::U, fields, dt)?;
        momentum_solver.solve(MomentumComponent::V, fields, dt)?;

        // Get Ap coefficients from momentum solver (needed for pressure correction)
        let (ap_u, _, ap_v, _) = momentum_solver.get_ap_coefficients();

        // Precompute D coefficients (inverse of Ap)
        // D = Vol / Ap. Here we use approximate area weighting: D_u = dy / Ap_u, D_v = dx / Ap_v
        let mut d_u = Field2D::new(nx, ny, T::zero());
        let mut d_v = Field2D::new(nx, ny, T::zero());
        let dx = grid.dx;
        let dy = grid.dy;

        // Use a safe threshold for Ap to avoid numerical instability
        let min_ap = T::from_f64(1e-10).unwrap();

        for j in 0..ny {
            for i in 0..nx {
                let ap_u_val = ap_u.at(i, j);
                if ap_u_val.abs() > min_ap {
                     // Ensure Ap is positive for stability
                     let val = if ap_u_val < T::zero() { -ap_u_val } else { ap_u_val };
                     d_u.set(i, j, dy / val);
                } else {
                     d_u.set(i, j, T::zero());
                }

                let ap_v_val = ap_v.at(i, j);
                if ap_v_val.abs() > min_ap {
                     let val = if ap_v_val < T::zero() { -ap_v_val } else { ap_v_val };
                     d_v.set(i, j, dx / val);
                } else {
                     d_v.set(i, j, T::zero());
                }
            }
        }

        // 2. Pressure Correction Step
        // Build linear system for pressure correction: div(u) = 0
        let n = nx * ny;
        let mut matrix_builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::zeros(n);

        let half = T::from_f64(0.5).unwrap();
        let two = T::from_f64(2.0).unwrap();

        // Max residual for convergence check
        let mut max_residual = T::zero();

        // Temporary storage for pressure correction
        let mut p_prime = DVector::zeros(n);

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                // --- Rhie-Chow Interpolation for Face Velocities ---

                // Get velocities
                let u_p = fields.u.at(i, j);
                let u_e = fields.u.at(i + 1, j);
                let u_w = fields.u.at(i - 1, j);
                let v_p = fields.v.at(i, j);
                let v_n = fields.v.at(i, j + 1);
                let v_s = fields.v.at(i, j - 1);

                // Get pressures
                let p_p = fields.p.at(i, j);
                let p_e = fields.p.at(i + 1, j);
                let p_w = fields.p.at(i - 1, j);
                let p_n = fields.p.at(i, j + 1);
                let p_s = fields.p.at(i, j - 1);

                // Linear extrapolation for boundaries to maintain second-order accuracy
                let two_val = T::from_f64(2.0).unwrap();
                let p_ee = if i + 2 < nx { fields.p.at(i + 2, j) } else { two_val * p_e - p_p };
                let p_ww = if i >= 2 { fields.p.at(i - 2, j) } else { two_val * p_w - p_p };
                let p_nn = if j + 2 < ny { fields.p.at(i, j + 2) } else { two_val * p_n - p_p };
                let p_ss = if j >= 2 { fields.p.at(i, j - 2) } else { two_val * p_s - p_p };

                let d_u_p = d_u.at(i, j);
                let d_u_e = d_u.at(i + 1, j);
                let d_v_p = d_v.at(i, j);
                let d_v_n = d_v.at(i, j + 1);

                // Face D coefficients (average)
                let d_face_e = (d_u_p + d_u_e) * half;
                let d_face_n = (d_v_p + d_v_n) * half;

                // Rhie-Chow Face Velocities

                // Calculate average pressure gradient at face e
                let grad_p_p_x = (p_e - p_w) / (two * dx);
                let grad_p_e_x = (p_ee - p_p) / (two * dx);
                let avg_grad_p_e = (grad_p_p_x + grad_p_e_x) * half;

                // Actual pressure gradient at face e
                let grad_p_e = (p_e - p_p) / dx;

                // Rhie-Chow velocity at East face
                // Using a small factor for Rhie-Chow interpolation to provide some pressure-velocity coupling
                // while maintaining stability. Full factor (1.0) can be unstable with certain BC configurations.
                let rc_factor = T::from_f64(0.05).unwrap();
                let u_face_e = (u_p + u_e) * half + d_face_e * (avg_grad_p_e - grad_p_e) * rc_factor;

                // Similar for North face
                let grad_p_p_y = (p_n - p_s) / (two * dy);
                let grad_p_n_y = (p_nn - p_p) / (two * dy);
                let avg_grad_p_n = (grad_p_p_y + grad_p_n_y) * half;
                let grad_p_n = (p_n - p_p) / dy;

                let v_face_n = (v_p + v_n) * half + d_face_n * (avg_grad_p_n - grad_p_n) * rc_factor;

                // We need u_face_w and v_face_s.
                let d_u_w = d_u.at(i - 1, j);
                let d_face_w = (d_u_w + d_u_p) * half;
                let grad_p_w_x = (p_p - p_ww) / (two * dx);
                let avg_grad_p_w = (grad_p_w_x + grad_p_p_x) * half;
                let grad_p_w = (p_p - p_w) / dx;
                let u_face_w = (u_w + u_p) * half + d_face_w * (avg_grad_p_w - grad_p_w) * rc_factor;

                let d_v_s = d_v.at(i, j - 1);
                let d_face_s = (d_v_s + d_v_p) * half;
                let grad_p_s_y = (p_p - p_ss) / (two * dy);
                let avg_grad_p_s = (grad_p_s_y + grad_p_p_y) * half;
                let grad_p_s = (p_p - p_s) / dy;
                let v_face_s = (v_s + v_p) * half + d_face_s * (avg_grad_p_s - grad_p_s) * rc_factor;


                // Mass Imbalance (b term for pressure equation)
                let rho = fields.density.at(i, j);
                let flux_e = rho * u_face_e * dy;
                let flux_w = rho * u_face_w * dy;
                let flux_n = rho * v_face_n * dx;
                let flux_s = rho * v_face_s * dx;

                let mass_imbalance = flux_e - flux_w + flux_n - flux_s;

                let a_e = rho * d_face_e * dy / dx;
                let a_w = rho * d_face_w * dy / dx;
                let a_n = rho * d_face_n * dx / dy;
                let a_s = rho * d_face_s * dx / dy;
                let a_p = -(a_e + a_w + a_n + a_s);

                // Add to matrix
                matrix_builder.add_entry(idx, idx, a_p)?;
                matrix_builder.add_entry(idx, idx + 1, a_e)?;
                matrix_builder.add_entry(idx, idx - 1, a_w)?;
                matrix_builder.add_entry(idx, idx + nx, a_n)?;
                matrix_builder.add_entry(idx, idx - nx, a_s)?;

                rhs[idx] = mass_imbalance;

                if mass_imbalance.abs() > max_residual {
                    max_residual = mass_imbalance.abs();
                }
            }
        }

        // Fill in boundary rows in matrix
        for j in 0..ny {
            for i in 0..nx {
                if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
                    continue; // Already handled
                }

                let idx = j * nx + i;

                let is_pressure_dirichlet = i == 0 || i == nx - 1; // Inlet/Outlet

                if is_pressure_dirichlet {
                    matrix_builder.add_entry(idx, idx, T::one())?;
                    rhs[idx] = T::zero();
                } else {
                    // Neumann (grad p' = 0)
                    let neighbor_idx = if i == 0 { idx + 1 }
                                       else if i == nx - 1 { idx - 1 }
                                       else if j == 0 { idx + nx }
                                       else { idx - nx };

                    matrix_builder.add_entry(idx, idx, T::one())?;
                    matrix_builder.add_entry(idx, neighbor_idx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
        }

        // Solve for p'
        let matrix = matrix_builder.build()?;
        let solver_config = IterativeSolverConfig {
            tolerance: self.tolerance,
            max_iterations: 2000,
            ..Default::default()
        };
        let linear_solver = BiCGSTAB::new(solver_config);

        linear_solver.solve(&matrix, &rhs, &mut p_prime, None::<&IdentityPreconditioner>)?;

        // 3. Correction Step
        // Correct pressure and velocity
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let pp = p_prime[idx];

                // Correct pressure: P_new = P_old + alpha_p * p'
                if let Some(p) = fields.p.at_mut(i, j) {
                    *p += self.pressure_relaxation * pp;
                }

                // Correct velocity: u_new = u* - D * grad p'

                let pp_e = p_prime[idx + 1];
                let pp_w = p_prime[idx - 1];
                let pp_n = p_prime[idx + nx];
                let pp_s = p_prime[idx - nx];

                let d_u_p = d_u.at(i, j);
                let d_v_p = d_v.at(i, j);

                // Calculate pressure correction gradients
                let dp_dx_prime = (pp_e - pp_w) / (two * dx);
                let dp_dy_prime = (pp_n - pp_s) / (two * dy);

                // Update velocities
                if let Some(u) = fields.u.at_mut(i, j) {
                    *u -= d_u_p * dp_dx_prime * dx;
                }

                if let Some(v) = fields.v.at_mut(i, j) {
                    *v -= d_v_p * dp_dy_prime * dy;
                }
            }
        }

        let converged = max_residual < self.tolerance;
        Ok((max_residual, converged))
    }

    /// Solve incompressible Navier-Stokes equations using SIMPLE algorithm
    pub fn solve_simple(
        &self,
        momentum_solver: &mut MomentumSolver<T>,
        poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<usize> {
        let mut iteration = 0;
        let mut converged = false;

        while iteration < self.max_iterations && !converged {
            let (_residual, conv) = self.simple_iteration(
                momentum_solver,
                poisson_solver,
                fields,
                dt,
                grid,
                boundary_conditions,
            )?;

            converged = conv;
            iteration += 1;

            // Optional: Log progress in debug builds
            #[cfg(debug_assertions)]
            tracing::debug!(
                "SIMPLE iteration {}: residual = {:?}, converged = {}",
                iteration,
                _residual,
                converged
            );
        }

        if !converged {
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.max_iterations,
                },
            ));
        }

        Ok(iteration)
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> Default for SimpleAlgorithm<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_algorithm_creation() {
        let simple = SimpleAlgorithm::<f64>::new();
        assert_eq!(simple.max_iterations, 50);
        assert!(simple.pressure_relaxation > 0.0 && simple.pressure_relaxation < 1.0);
        assert!(simple.velocity_relaxation > 0.0 && simple.velocity_relaxation < 1.0);
    }
}
