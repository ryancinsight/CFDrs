//! Pressure-Velocity Coupling Algorithms for Incompressible Flow
//!
//! This module implements algorithms for solving the coupled momentum and continuity
//! equations in incompressible CFD. The primary algorithm is SIMPLE (Semi-Implicit
//! Method for Pressure-Linked Equations).

use crate::fields::{SimulationFields, Field2D};
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::MomentumSolver;
use crate::physics::turbulence::TurbulenceModel;
use crate::pressure_velocity::rhie_chow::RhieChowInterpolation;
use crate::solvers::fdm::PoissonSolver;
use cfd_core::boundary::BoundaryCondition;
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

    /// Execute one SIMPLE iteration
    ///
    /// Returns the maximum residual and whether convergence was achieved
    pub fn simple_iteration(
        &self,
        momentum_solver: &mut MomentumSolver<T>,
        poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<(T, bool)> {
        let nx = grid.nx;
        let ny = grid.ny;

        // Store old velocities for under-relaxation
        let u_old = fields.u.clone();
        let v_old = fields.v.clone();
        let p_old = fields.p.clone();

        // 1. MOMENTUM PREDICTION STEP
        // Solve momentum equations with current pressure guess

        // Solve U-momentum equation
        let mut coeffs_u = momentum_solver.solve_with_coefficients(
            crate::physics::momentum::MomentumComponent::U,
            fields,
            dt,
        )?;

        // Solve V-momentum equation
        let mut coeffs_v = momentum_solver.solve_with_coefficients(
            crate::physics::momentum::MomentumComponent::V,
            fields,
            dt,
        )?;

        // Extract predicted velocities from momentum solver
        // (In practice, these would be extracted from the linear system solution)
        let u_star = fields.u.clone(); // Predicted U-velocity
        let v_star = fields.v.clone(); // Predicted V-velocity

        // Build previous velocity field as Vector2 for transient Rhie–Chow term
        let mut u_old_vec = Field2D::new(grid.nx, grid.ny, Vector2::new(T::zero(), T::zero()));
        for j in 0..grid.ny {
            for i in 0..grid.nx {
                u_old_vec[(i, j)] = Vector2::new(u_old[(i, j)], v_old[(i, j)]);
            }
        }

        // Apply Rhie-Chow interpolation using actual A_p coefficients from momentum solver
        self.apply_rhie_chow_interpolation(
            &u_star,
            &v_star,
            fields,
            grid,
            dt,
            &coeffs_u.ap,
            &coeffs_v.ap,
            &u_old_vec,
        )?;

        // 2. PRESSURE CORRECTION STEP
        // Solve pressure Poisson equation: ∇²p' = ∇·u*

        let mut pressure_correction = DVector::zeros(nx * ny);

        // Construct pressure correction equation
        self.construct_pressure_correction_equation(
            &u_star, &v_star, &mut pressure_correction,
            fields, grid, dt,
        )?;

        // Build sparse matrix for pressure Poisson equation: ∇²p' = rhs
        let mut matrix_builder = SparseMatrixBuilder::new(nx * ny, nx * ny);
        self.build_pressure_poisson_matrix(&mut matrix_builder, grid)?;

        // Solve using conjugate gradient with ILU preconditioning
        let poisson_matrix = matrix_builder.build()?;
        let cg_config = IterativeSolverConfig {
            tolerance: T::from_f64(1e-8).unwrap(),
            max_iterations: 1000,
            ..Default::default()
        };
        let cg_solver = ConjugateGradient::new(cg_config);
        let ilu_preconditioner = IncompleteLU::new(&poisson_matrix)?;

        // Solve ∇²p' = rhs for pressure correction p'
        let p_prime_vec = cg_solver.solve_preconditioned(
            &poisson_matrix,
            &pressure_correction,
            &ilu_preconditioner,
            None,
        )?;

        // Convert solution vector back to grid format
        let mut p_prime = fields.p.clone();
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                p_prime[(i, j)] = p_prime_vec[idx];
            }
        }

        // Apply pressure under-relaxation
        for i in 0..nx {
            for j in 0..ny {
                let idx = j * nx + i;
                fields.p[(i, j)] = p_old[(i, j)] + self.pressure_relaxation * p_prime[(i, j)];
            }
        }

        // 3. VELOCITY CORRECTION STEP
        // Correct velocities using Rhie–Chow corrected prediction as base:
        // u = u_RC* - (dt/ρ)∇p'

        let rho = fields.density.at(0, 0); // Assume constant density
        let correction_factor = dt / rho;

        for i in 0..nx {
            for j in 0..ny {
                // Compute pressure gradient
                let dp_dx = if i > 0 && i < (nx as usize) - 1 {
                    (fields.p[(i + 1, j)] - fields.p[(i - 1, j)]) / (grid.dx + grid.dx)
                } else {
                    T::zero()
                };

                let dp_dy = if j > 0 && j < (ny as usize) - 1 {
                    (fields.p[(i, j + 1)] - fields.p[(i, j - 1)]) / (grid.dy + grid.dy)
                } else {
                    T::zero()
                };

                // Correct velocities on top of Rhie–Chow interpolated prediction
                let base_u = fields.u[(i, j)];
                let base_v = fields.v[(i, j)];
                fields.u[(i, j)] = base_u - correction_factor * dp_dx;
                fields.v[(i, j)] = base_v - correction_factor * dp_dy;
            }
        }

        // Apply velocity under-relaxation
        for i in 0..nx {
            for j in 0..ny {
                fields.u[(i, j)] = self.velocity_relaxation * fields.u[(i, j)] +
                                  (T::one() - self.velocity_relaxation) * u_old[(i, j)];
                fields.v[(i, j)] = self.velocity_relaxation * fields.v[(i, j)] +
                                  (T::one() - self.velocity_relaxation) * v_old[(i, j)];
            }
        }

        // 4. CHECK CONVERGENCE
        // Compute maximum residual of continuity equation
        let mut max_residual = T::zero();
        for i in 0..nx {
            for j in 0..ny {
                let residual = self.compute_continuity_residual(i, j, fields, grid);
                let abs_residual = if residual >= T::zero() { residual } else { -residual };
                if abs_residual > max_residual {
                    max_residual = abs_residual;
                }
            }
        }

        let converged = max_residual < self.tolerance;

        Ok((max_residual, converged))
    }

    /// Apply Rhie-Chow interpolation to prevent pressure oscillations
    ///
    /// This prevents checkerboard pressure oscillations by including pressure
    /// gradients in the velocity interpolation at cell faces.
    ///
    /// **Theorem (Rhie-Chow Consistency)**: For colocated grids, the pressure correction equation
    /// becomes ill-conditioned, leading to checkerboard oscillations in pressure and velocity fields.
    /// The face velocity is computed as: u_f = ū_f + d_f * [(∇p)_cells - (∇p)_face]
    ///
    /// **Reference**: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent
    /// flow past an airfoil with trailing edge separation." AIAA Journal, 21(11), 1525-1532.
    fn apply_rhie_chow_interpolation(
        &self,
        u_star: &Field2D<T>,
        v_star: &Field2D<T>,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        dt: T,
        ap_u: &Field2D<T>,
        ap_v: &Field2D<T>,
        u_old_vec: &Field2D<Vector2<T>>,
    ) -> Result<()> {
        // Create Rhie-Chow interpolator
        let mut rhie_chow = RhieChowInterpolation::new(grid);

        // Update momentum coefficients (A_p from discretized momentum equations)
        // Use actual coefficients provided by momentum solver
        rhie_chow.update_u_coefficients(ap_u);
        rhie_chow.update_v_coefficients(ap_v);

        // Provide previous velocity field for transient correction term
        rhie_chow.update_old_velocity(u_old_vec);

        // Convert predicted velocities to Vector2 field format for Rhie–Chow interpolation
        let mut velocity_field = Field2D::new(grid.nx, grid.ny, Vector2::new(T::zero(), T::zero()));

        // Copy velocities to field format (interior cells only)
        for i in 0..grid.nx {
            for j in 0..grid.ny {
                velocity_field[(i, j)] = Vector2::new(u_star[(i, j)], v_star[(i, j)]);
            }
        }

        // Apply Rhie-Chow interpolation to prevent pressure-velocity oscillations
        // This ensures that the interpolated face velocities satisfy the momentum equations
        // and prevent checkerboard pressure patterns

        // For each interior cell, apply Rhie-Chow correction to velocities
        for i in 1..grid.nx - 1 {
            for j in 1..grid.ny - 1 {
                // Compute Rhie-Chow corrected face velocities and apply correction to cell velocities
                // This prevents the decoupling that causes pressure oscillations

                // East face correction (affects u-velocity at cell i,j)
                if i < grid.nx - 1 {
                    let u_face_corrected = rhie_chow.face_velocity_x(
                        &velocity_field,
                        &fields.p,
                        grid.dx,
                        grid.dy,
                        Some(dt),
                        i,
                        j,
                    );

                    // Apply Rhie-Chow correction to u-velocity (full correction, no blending)
                    fields.u[(i, j)] = u_face_corrected;
                }

                // North face correction (affects v-velocity at cell i,j)
                if j < grid.ny - 1 {
                    let v_face_corrected = rhie_chow.face_velocity_y(
                        &velocity_field,
                        &fields.p,
                        grid.dx,
                        grid.dy,
                        Some(dt),
                        i,
                        j,
                    );

                    // Apply Rhie-Chow correction to v-velocity (full correction, no blending)
                    fields.v[(i, j)] = v_face_corrected;
                }
            }
        }

        Ok(())
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

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute divergence of predicted velocity field
                let du_dx = if i < (nx as usize) - 1 {
                    (u_star.at(i + 1, j) - u_star.at(i, j)) / dx
                } else {
                    T::zero()
                };

                let dv_dy = if j < (ny as usize) - 1 {
                    (v_star.at(i, j + 1) - v_star.at(i, j)) / dy
                } else {
                    T::zero()
                };

                // Continuity residual: ∇·u* = 0 (incompressibility)
                let continuity_residual = du_dx + dv_dy;

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
        let two = T::from_f64(2.0).unwrap();

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
        let dx = grid.dx;
        let dy = grid.dy;
        let nx = grid.nx;
        let ny = grid.ny;

        // Compute velocity divergence
        let du_dx = if i < (nx as usize) - 1 {
            (fields.u[(i + 1, j)] - fields.u[(i, j)]) / dx
        } else {
            T::zero()
        };

        let dv_dy = if j < (ny as usize) - 1 {
            (fields.v[(i, j + 1)] - fields.v[(i, j)]) / dy
        } else {
            T::zero()
        };

        du_dx + dv_dy
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
        let mut iteration = 0;
        let mut converged = false;

        while iteration < self.max_iterations && !converged {
            let (residual, conv) = self.simple_iteration(
                momentum_solver,
                poisson_solver,
                fields,
                dt,
                grid,
                boundary_conditions,
            )?;

            converged = conv;
            iteration += 1;

            // Log progress (in debug builds)
            #[cfg(debug_assertions)]
            tracing::debug!(
                "SIMPLE iteration {}: residual = {:?}, converged = {}",
                iteration,
                residual,
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

    #[test]
    fn test_rhie_chow_interpolation_uses_ap_coefficients() {
        let simple = SimpleAlgorithm::<f64>::new();

        // 3x3 grid
        let grid = StructuredGrid2D::<f64>::new(3, 3, 0.0, 1.0, 0.0, 1.0).unwrap();

        // Two field states to compare effects of different A_p
        let mut fields_high = SimulationFields::<f64>::new(3, 3);
        let mut fields_low = SimulationFields::<f64>::new(3, 3);

        // Impose a quadratic pressure field in x-direction to ensure
        // (∇p)_cells ≠ (∇p)_face and produce a non-zero correction
        for j in 0..3 {
            for i in 0..3 {
                let p_val = (i as f64).powi(2);
                fields_high.p[(i, j)] = p_val;
                fields_low.p[(i, j)] = p_val;
            }
        }

        // Predicted velocities (zeros)
        let u_star = Field2D::new(3, 3, 0.0f64);
        let v_star = Field2D::new(3, 3, 0.0f64);

        // Momentum coefficients A_p
        let ap_high = Field2D::new(3, 3, 10.0f64);
        let ap_low = Field2D::new(3, 3, 0.1f64);

        // Old velocity buffer (zeros)
        let u_old_vec = Field2D::new(3, 3, nalgebra::Vector2::new(0.0f64, 0.0f64));

        let dt = 0.01f64;

        // Apply Rhie–Chow with high A_p (smaller correction expected)
        simple
            .apply_rhie_chow_interpolation(
                &u_star,
                &v_star,
                &mut fields_high,
                &grid,
                dt,
                &ap_high,
                &ap_high,
                &u_old_vec,
            )
            .unwrap();

        // Apply Rhie–Chow with low A_p (larger correction expected)
        simple
            .apply_rhie_chow_interpolation(
                &u_star,
                &v_star,
                &mut fields_low,
                &grid,
                dt,
                &ap_low,
                &ap_low,
                &u_old_vec,
            )
            .unwrap();

        let u_center_high = fields_high.u[(1, 1)];
        let u_center_low = fields_low.u[(1, 1)];

        // Non-zero correction under a pressure gradient
        assert!(u_center_high.abs() > 0.0, "Rhie–Chow must produce non-zero correction for pressure gradient");
        // Verify correction magnitude depends on A_p via d_f = V / A_p
        assert!(
            u_center_low.abs() > u_center_high.abs(),
            "Smaller A_p should increase Rhie–Chow correction magnitude",
        );
    }
}
