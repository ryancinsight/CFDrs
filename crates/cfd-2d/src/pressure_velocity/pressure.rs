//! Pressure correction solver for STANDARD algorithm

use crate::grid::StructuredGrid2D;
use cfd_math::{ConjugateGradient, LinearSolver, SparseMatrixBuilder};
use nalgebra::{DVector, RealField, Vector2};
use num_traits::FromPrimitive;

/// Pressure correction solver
pub struct PressureCorrectionSolver<T: RealField + Copy> {
    /// Grid
    grid: StructuredGrid2D<T>,
    /// Linear solver for pressure Poisson equation
    linear_solver: ConjugateGradient<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> PressureCorrectionSolver<T> {
    /// Create new pressure correction solver
    pub fn new(grid: StructuredGrid2D<T>) -> cfd_core::error::Result<Self> {
        let config = cfd_core::solver::LinearSolverConfig {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-8).unwrap_or_else(|| T::zero()),
            preconditioning: true,
        };

        Ok(Self {
            grid,
            linear_solver: ConjugateGradient::new(config),
        })
    }

    /// Solve pressure correction equation
    ///
    /// The pressure correction equation is derived from the continuity equation:
    /// ∇²p' = ρ/Δt * ∇·u*
    ///
    /// where p' is the pressure correction and u* is the predicted velocity
    pub fn solve_pressure_correction(
        &self,
        u_star: &Vec<Vec<Vector2<T>>>,
        dt: T,
        rho: T,
    ) -> cfd_core::error::Result<Vec<Vec<T>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        // Build discrete Laplacian matrix and divergence source term
        let n = (nx - 2) * (ny - 2); // Interior points only
        let mut builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::zeros(n);

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);
        let coeff = rho / dt;

        // Build system for interior points
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);

                // Laplacian stencil
                let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());
                let ap = -two * (dx2_inv + dy2_inv);
                builder.add_entry(idx, idx, ap)?;

                // Neighbors
                if i > 1 {
                    builder.add_entry(idx, idx - (ny - 2), dx2_inv)?;
                }
                if i < nx - 2 {
                    builder.add_entry(idx, idx + (ny - 2), dx2_inv)?;
                }
                if j > 1 {
                    builder.add_entry(idx, idx - 1, dy2_inv)?;
                }
                if j < ny - 2 {
                    builder.add_entry(idx, idx + 1, dy2_inv)?;
                }

                // Divergence of predicted velocity
                let div_u = (u_star[i + 1][j].x - u_star[i - 1][j].x)
                    / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx)
                    + (u_star[i][j + 1].y - u_star[i][j - 1].y)
                        / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);

                rhs[idx] = coeff * div_u;
            }
        }

        // Solve the linear system
        let matrix = builder.build()?;
        let p_correction_vec = self.linear_solver.solve(&matrix, &rhs, None)?;

        // Convert back to 2D grid
        let mut p_correction = vec![vec![T::zero(); ny]; nx];
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                p_correction[i][j] = p_correction_vec[idx];
            }
        }

        // Apply Neumann boundary conditions (zero gradient)
        for i in 0..nx {
            p_correction[i][0] = p_correction[i][1];
            p_correction[i][ny - 1] = p_correction[i][ny - 2];
        }
        for j in 0..ny {
            p_correction[0][j] = p_correction[1][j];
            p_correction[nx - 1][j] = p_correction[nx - 2][j];
        }

        Ok(p_correction)
    }

    /// Correct velocity field using pressure correction
    pub fn correct_velocity(
        &self,
        u_star: &mut Vec<Vec<Vector2<T>>>,
        p_correction: &Vec<Vec<T>>,
        dt: T,
        rho: T,
        alpha: T,
    ) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        let factor = dt / rho;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                // Velocity correction from pressure gradient
                let dp_dx = (p_correction[i + 1][j] - p_correction[i - 1][j])
                    / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
                let dp_dy = (p_correction[i][j + 1] - p_correction[i][j - 1])
                    / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);

                // Apply velocity correction with relaxation
                u_star[i][j].x -= alpha * factor * dp_dx;
                u_star[i][j].y -= alpha * factor * dp_dy;
            }
        }
    }

    /// Correct pressure field
    pub fn correct_pressure(&self, p: &mut Vec<Vec<T>>, p_correction: &Vec<Vec<T>>, alpha: T) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        for i in 0..nx {
            for j in 0..ny {
                p[i][j] += alpha * p_correction[i][j];
            }
        }
    }
}
