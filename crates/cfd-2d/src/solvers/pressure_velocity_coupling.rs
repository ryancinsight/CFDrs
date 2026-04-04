//! Pressure-Velocity Coupling Algorithms for Incompressible Flow
//!
//! This module implements algorithms for solving the coupled momentum and continuity
//! equations in incompressible CFD. The primary algorithm is SIMPLE (Semi-Implicit
//! Method for Pressure-Linked Equations).
//!
//! # Theorem (SIMPLE Pressure Correction — Patankar 1980)
//!
//! The pressure correction equation $\nabla \cdot (\mathbf{d}\,\nabla p') = \nabla \cdot \mathbf{u}^*$
//! is a Poisson equation with an M-matrix coefficient structure. Under appropriate
//! velocity relaxation ($0 < \alpha_u < 1$) the coupled outer iteration reduces
//! the continuity residual monotonically.
//!
//! **Proof sketch**:
//! The face velocity $\mathbf{u}_f^* = \bar{\mathbf{u}}_f + \mathbf{d}_f(\nabla p^* - \nabla p_f^*)$
//! yields a pressure Poisson equation with positive definite $\mathbf{d}_f = V_f/A_{P,f}$.
//! The resulting matrix is symmetric negative-definite (M-matrix), and its spectral
//! radius under relaxation is bounded by $\max(1 - \alpha_p, \alpha_u) < 1$.


use crate::fields::{SimulationFields, Field2D};
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::MomentumSolver;
use crate::physics::turbulence::TurbulenceModel;
use crate::solvers::continuity::{
    max_forward_continuity_residual, pointwise_forward_continuity_residual,
};
use crate::solvers::fdm::PoissonSolver;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::error::{Error, Result};
use cfd_math::linear_solver::{ConjugateGradient, IterativeSolverConfig, preconditioners::IncompleteLU};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField, Vector2};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
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
    /// Convergence tolerance for pressure residual
    tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
    /// Create new SIMPLE algorithm with default parameters
    pub fn new() -> Self {
        Self {
            pressure_relaxation: T::from_f64(0.3).unwrap_or_else(num_traits::Zero::zero), // Standard value
            velocity_relaxation: T::from_f64(0.7).unwrap_or_else(num_traits::Zero::zero), // Standard value
            max_iterations: 50,
            tolerance: T::from_f64(1e-6).unwrap_or_else(num_traits::Zero::zero),
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



    /// Construct the pressure correction equation RHS
    fn construct_pressure_correction_equation(
        &self,
        u_star: &Field2D<T>,
        v_star: &Field2D<T>,
        rhs: &mut DVector<T>,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        dt: T,
    ) -> Result<()> {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;

        let mut u_at = |ii: usize, jj: usize| u_star.at(ii, jj);
        let mut v_at = |ii: usize, jj: usize| v_star.at(ii, jj);

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Continuity residual: ∇·u* = 0 (incompressibility)
                let continuity_residual = pointwise_forward_continuity_residual(
                    i,
                    j,
                    nx,
                    ny,
                    dx,
                    dy,
                    &mut u_at,
                    &mut v_at,
                );

                // Pressure Poisson equation: ∇²p' = (ρ/dt) * ∇·u*
                let rho = fields.density.at(i, j);
                rhs[idx] = (rho / dt) * continuity_residual;
            }
        }

        Ok(())
    }

    /// Build sparse matrix for pressure Poisson equation: ∇²p' = rhs
    ///
    /// Uses 5-point Laplacian stencil on structured grid:
    /// ∇²p = (p_{i+1,j} + p_{i-1,j} - 2p_{i,j})/dx² + (p_{i,j+1} + p_{i,j-1} - 2p_{i,j})/dy²
    ///
    /// **Theorem (Discrete Laplacian)**: The 5-point stencil provides second-order
    /// accurate discretization of the Laplacian operator on uniform Cartesian grids.
    ///
    /// **Reference**: Numerical Recipes in C, Chapter 19: Partial Differential Equations
    fn build_pressure_poisson_matrix(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<()> {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;

        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let two = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

        // Build matrix for interior points (Dirichlet boundaries assumed p' = 0)
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Check if this is a boundary point (assume Dirichlet p' = 0)
                let is_boundary = i == 0 || i == (nx as usize) - 1 ||
                                j == 0 || j == (ny as usize) - 1;

                if is_boundary {
                    // Dirichlet boundary condition: p' = 0
                    matrix_builder.add_entry(idx, idx, T::one())?;
                } else {
                    // Interior point: 5-point Laplacian stencil
                    // Center coefficient: -2/dx² - 2/dy²
                    let center_coeff = -two / dx2 - two / dy2;
                    matrix_builder.add_entry(idx, idx, center_coeff)?;

                    // Left neighbor (i-1): +1/dx²
                    let left_idx = j * nx + (i - 1);
                    matrix_builder.add_entry(idx, left_idx, T::one() / dx2)?;

                    // Right neighbor (i+1): +1/dx²
                    let right_idx = j * nx + (i + 1);
                    matrix_builder.add_entry(idx, right_idx, T::one() / dx2)?;

                    // Bottom neighbor (j-1): +1/dy²
                    let bottom_idx = (j - 1) * nx + i;
                    matrix_builder.add_entry(idx, bottom_idx, T::one() / dy2)?;

                    // Top neighbor (j+1): +1/dy²
                    let top_idx = (j + 1) * nx + i;
                    matrix_builder.add_entry(idx, top_idx, T::one() / dy2)?;
                }
            }
        }

        Ok(())
    }

    /// Compute continuity residual at a grid point
    fn compute_continuity_residual(
        &self,
        i: usize,
        j: usize,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> T {
        let mut u_at = |ii: usize, jj: usize| fields.u.at(ii, jj);
        let mut v_at = |ii: usize, jj: usize| fields.v.at(ii, jj);

        pointwise_forward_continuity_residual(
            i,
            j,
            grid.nx,
            grid.ny,
            grid.dx,
            grid.dy,
            &mut u_at,
            &mut v_at,
        )
    }

    /// Solve incompressible Navier-Stokes equations using SIMPLE algorithm
    ///
    /// This is the main driver function that iterates SIMPLE until convergence.
    pub fn solve_simple(
        &self,
        momentum_solver: &mut MomentumSolver<T>,
        poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<usize> {
        let nx = grid.nx;
        let ny = grid.ny;

        // Workspace
        let mut u_old = fields.u.clone();
        let mut v_old = fields.v.clone();
        let mut p_old = fields.p.clone();
        let mut pressure_correction = DVector::zeros(nx * ny);

        // Precompute Sparse Matrix
        let mut matrix_builder = SparseMatrixBuilder::new(nx * ny, nx * ny);
        self.build_pressure_poisson_matrix(&mut matrix_builder, grid)?;
        let poisson_matrix = matrix_builder.build()?;
        let cg_config = IterativeSolverConfig {
            tolerance: T::from_f64(1e-8).unwrap_or_else(num_traits::Zero::zero),
            max_iterations: 1000,
            ..Default::default()
        };
        let cg_solver = ConjugateGradient::new(cg_config);
        let ilu_preconditioner = IncompleteLU::new(&poisson_matrix)?;

        let mut iteration = 0;
        let mut converged = false;

        while iteration < self.max_iterations && !converged {
            // Save state
            for i in 0..nx {
                for j in 0..ny {
                    u_old[(i, j)] = fields.u[(i, j)];
                    v_old[(i, j)] = fields.v[(i, j)];
                    p_old[(i, j)] = fields.p[(i, j)];
                }
            }

            // Momentum solve
            let _coeffs_u = momentum_solver.solve_with_coefficients(
                crate::physics::momentum::MomentumComponent::U,
                fields,
                dt,
            )?;
            let _coeffs_v = momentum_solver.solve_with_coefficients(
                crate::physics::momentum::MomentumComponent::V,
                fields,
                dt,
            )?;

            // Pressure correction eqn
            self.construct_pressure_correction_equation(
                &fields.u, &fields.v, &mut pressure_correction,
                fields, grid, dt,
            )?;

            // CG Solve
            let p_prime_vec = cg_solver.solve_preconditioned(
                &poisson_matrix,
                &pressure_correction,
                &ilu_preconditioner,
                None,
            )?;

            let rho = fields.density.at(0, 0); 
            let correction_factor = dt / rho;

            for j in 0..ny {
                for i in 0..nx {
                    let idx = j * nx + i;
                    let p_prime_val = p_prime_vec[idx];
                    fields.p[(i, j)] = p_old[(i, j)] + self.pressure_relaxation * p_prime_val;

                    let dp_dx = if i > 0 && i < nx - 1 {
                        (p_prime_vec[idx + 1] - p_prime_vec[idx - 1]) / (grid.dx + grid.dx)
                    } else {
                        T::zero()
                    };
                    let dp_dy = if j > 0 && j < ny - 1 {
                        (p_prime_vec[idx + nx] - p_prime_vec[idx - nx]) / (grid.dy + grid.dy)
                    } else {
                        T::zero()
                    };

                    let base_u = fields.u[(i, j)];
                    let base_v = fields.v[(i, j)];
                    fields.u[(i, j)] = base_u - correction_factor * dp_dx;
                    fields.v[(i, j)] = base_v - correction_factor * dp_dy;
                }
            }

            for j in 0..ny {
                for i in 0..nx {
                    fields.u[(i, j)] = self.velocity_relaxation * fields.u[(i, j)] +
                                      (T::one() - self.velocity_relaxation) * u_old[(i, j)];
                    fields.v[(i, j)] = self.velocity_relaxation * fields.v[(i, j)] +
                                      (T::one() - self.velocity_relaxation) * v_old[(i, j)];
                }
            }

            let max_residual = max_forward_continuity_residual(
                nx,
                ny,
                grid.dx,
                grid.dy,
                |i, j| fields.u.at(i, j),
                |i, j| fields.v.at(i, j),
            );

            converged = max_residual < self.tolerance;
            iteration += 1;

            #[cfg(debug_assertions)]
            tracing::debug!(
                "SIMPLE iteration {}: residual = {:?}, converged = {}",
                iteration,
                max_residual,
                converged
            );
        }

        if !converged {
            return Err(Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.max_iterations,
                },
            ));
        }

        Ok(iteration)
    }
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::Debug> Default for SimpleAlgorithm<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::{SimulationFields, Field2D};
    use crate::grid::StructuredGrid2D;
    use crate::physics::momentum::MomentumSolver;
    use crate::solvers::fdm::PoissonSolver;
    use nalgebra::DVector;
    use approx::assert_relative_eq;

    #[test]
    fn test_simple_algorithm_creation() {
        let simple = SimpleAlgorithm::<f64>::new();
        assert_eq!(simple.max_iterations, 50);
        assert!(simple.pressure_relaxation > 0.0 && simple.pressure_relaxation < 1.0);
        assert!(simple.velocity_relaxation > 0.0 && simple.velocity_relaxation < 1.0);
    }

    #[test]
    fn test_simple_algorithm_configuration() {
        let simple = SimpleAlgorithm::<f64>::new()
            .with_pressure_relaxation(0.5)
            .with_velocity_relaxation(0.8)
            .with_max_iterations(100)
            .with_tolerance(1e-8);

        assert_eq!(simple.pressure_relaxation, 0.5);
        assert_eq!(simple.velocity_relaxation, 0.8);
        assert_eq!(simple.max_iterations, 100);
        assert_eq!(simple.tolerance, 1e-8);
    }

    #[test]
    fn test_continuity_residual_computation() {
        let simple = SimpleAlgorithm::<f64>::new();

        // Create simple 2x2 grid
        let grid = StructuredGrid2D::<f64>::new(2, 2, 0.0, 1.0, 0.0, 1.0).unwrap();
        let mut fields = SimulationFields::<f64>::new(2, 2);

        // Set up a divergence-free velocity field: u = y, v = -x
        // ∇·u = ∂u/∂x + ∂v/∂y = 0 + 0 = 0
        for i in 0..2 {
            for j in 0..2 {
                fields.u[(i, j)] = j as f64 * grid.dy; // u = y
                fields.v[(i, j)] = -(i as f64) * grid.dx; // v = -x
            }
        }

        // Check continuity residual at interior point (0,0)
        let residual = simple.compute_continuity_residual(0, 0, &fields, &grid);
        assert!(residual.abs() < 1e-10, "Divergence should be zero for this field");

        // Set up a divergent field: u = x, v = y
        // ∇·u = ∂u/∂x + ∂v/∂y = 1 + 1 = 2
        for i in 0..2 {
            for j in 0..2 {
                fields.u[(i, j)] = i as f64 * grid.dx; // u = x
                fields.v[(i, j)] = j as f64 * grid.dy; // v = y
            }
        }

        let residual = simple.compute_continuity_residual(0, 0, &fields, &grid);
        assert!(residual > 1.0, "Divergence should be positive for this field");
    }

    #[test]
    fn test_pressure_correction_equation_construction() {
        let simple = SimpleAlgorithm::<f64>::new();

        // Create 3x3 grid for better testing
        let grid = StructuredGrid2D::<f64>::new(3, 3, 0.0, 1.0, 0.0, 1.0).unwrap();
        let mut fields = SimulationFields::<f64>::new(3, 3);
        let mut rhs = DVector::zeros(9);

        // Set up a simple divergent velocity field
        let mut u_star = Field2D::new(3, 3, 0.0f64);
        let mut v_star = Field2D::new(3, 3, 0.0f64);
        u_star[(1, 1)] = 1.0;
        v_star[(1, 1)] = 1.0;

        let dt = 0.01;

        simple.construct_pressure_correction_equation(
            &u_star, &v_star, &mut rhs, &fields, &grid, dt,
        ).unwrap();

        // Check that RHS has non-zero values where divergence occurs
        // The center point should have significant RHS contribution
        let center_idx = 4; // (1,1) in row-major order
        assert!(rhs[center_idx].abs() > 0.0, "Pressure correction RHS should be non-zero at divergent points");
    }


}
