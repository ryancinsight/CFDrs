//! Pressure correction solver for STANDARD algorithm

use super::config::PressureLinearSolver;
use crate::grid::StructuredGrid2D;
use cfd_math::linear_solver::preconditioners::{
    AMGConfig, AlgebraicMultigrid, CoarseningStrategy, IdentityPreconditioner, MultigridCycle,
};
use cfd_math::linear_solver::{BiCGSTAB, ConjugateGradient, IterativeLinearSolver, GMRES};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField, Vector2};
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// Pressure correction solver supporting multiple linear solver backends
pub struct PressureCorrectionSolver<T: RealField + Copy> {
    /// Grid
    grid: StructuredGrid2D<T>,
    /// Linear solver choice
    solver_type: PressureLinearSolver,
    /// Conjugate Gradient solver (for symmetric systems)
    cg_solver: ConjugateGradient<T>,
    /// BiCGSTAB solver (for non-symmetric systems)
    bicgstab_solver: BiCGSTAB<T>,
    /// GMRES solver (industry standard for SIMPLE/PISO)
    gmres_solver: Option<GMRES<T>>,
}

impl<T: RealField + Copy + FromPrimitive + Debug> PressureCorrectionSolver<T> {
    /// Create new pressure correction solver with specified linear solver
    ///
    /// # Arguments
    ///
    /// * `grid` - Computational grid
    /// * `solver_type` - Linear solver selection (CG, BiCGSTAB, or GMRES)
    ///
    /// # Recommended Solver Selection
    ///
    /// - **GMRES**: Industry standard for SIMPLE/PISO (default)
    ///   - Reference: Saad (2003), Patankar (1980)
    ///   - Best for non-symmetric pressure Poisson equations
    ///   - restart_dim=30 optimal for most CFD applications
    /// - **BiCGSTAB**: Good alternative for non-symmetric systems
    /// - **CG**: Only for symmetric systems (not typical in CFD)
    pub fn new(
        grid: StructuredGrid2D<T>,
        solver_type: PressureLinearSolver,
    ) -> cfd_core::error::Result<Self> {
        let config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: crate::constants::solver::DEFAULT_MAX_ITERATIONS,
            tolerance: T::from_f64(crate::constants::solver::DEFAULT_TOLERANCE)
                .expect("Failed to convert pressure solver tolerance"),
            use_preconditioner: true,
            use_parallel_spmv: false,
        };

        let gmres_solver = match solver_type {
            PressureLinearSolver::GMRES { restart_dim } => Some(GMRES::new(config, restart_dim)),
            _ => None,
        };

        Ok(Self {
            grid,
            solver_type,
            cg_solver: ConjugateGradient::new(config),
            bicgstab_solver: BiCGSTAB::new(config),
            gmres_solver,
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

        if n == 0 {
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }
        if n == 1 {
            // Single interior point: fix pressure correction to zero
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }

        let reference_idx = 0usize;
        let system_size = n - 1;

        let map_index = |idx: usize| -> Option<usize> {
            if idx == reference_idx {
                None
            } else if idx < reference_idx {
                Some(idx)
            } else {
                Some(idx - 1)
            }
        };

        let mut builder = SparseMatrixBuilder::new(system_size, system_size);
        let mut rhs = DVector::zeros(system_size);

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);
        let coeff = rho / dt;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);

                if idx == reference_idx {
                    continue;
                }

                let row_idx = map_index(idx).expect("row index must exist");

                // Laplacian stencil - diagonal is negative sum of neighbor coefficients
                let two = T::from_f64(cfd_core::constants::numerical::common::TWO)
                    .unwrap_or_else(|| T::one() + T::one());
                let ap = -two * (dx2_inv + dy2_inv);

                builder.add_entry(row_idx, row_idx, ap)?;

                // Neighbors
                if i > 1 {
                    let neighbor_idx = idx - (ny - 2);
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dx2_inv)?;
                    }
                }
                if i < nx - 2 {
                    let neighbor_idx = idx + (ny - 2);
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dx2_inv)?;
                    }
                }
                if j > 1 {
                    let neighbor_idx = idx - 1;
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dy2_inv)?;
                    }
                }
                if j < ny - 2 {
                    let neighbor_idx = idx + 1;
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dy2_inv)?;
                    }
                }

                // Divergence of predicted velocity
                let div_u = (u_star[i + 1][j].x - u_star[i - 1][j].x)
                    / (T::from_f64(cfd_core::constants::numerical::common::TWO)
                        .unwrap_or_else(|| T::zero())
                        * dx)
                    + (u_star[i][j + 1].y - u_star[i][j - 1].y)
                        / (T::from_f64(cfd_core::constants::numerical::common::TWO)
                            .unwrap_or_else(|| T::zero())
                            * dy);

                rhs[row_idx] = coeff * div_u;
            }
        }

        // Solve the linear system using selected solver
        let matrix = builder.build()?;
        let mut p_correction_vec = DVector::zeros(matrix.nrows());

        // Create AMG preconditioner for this solve
        // AMG is particularly effective for the pressure Poisson equation
        let amg_config = AMGConfig {
            cycle_type: MultigridCycle::V,
            nu1: 2,
            nu2: 2,
            max_levels: 5,
            min_coarse_size: 10,
            coarsening: CoarseningStrategy::RugeStueben,
        };

        let amg_preconditioner = match AlgebraicMultigrid::with_config(&matrix, amg_config) {
            Ok(amg) => Some(amg),
            Err(_) => {
                tracing::debug!("AMG preconditioner construction failed, using identity");
                None
            }
        };

        match self.solver_type {
            PressureLinearSolver::ConjugateGradient => {
                // Use AMG preconditioner if available, otherwise identity
                if let Some(ref amg) = amg_preconditioner {
                    self.cg_solver
                        .solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                } else {
                    self.cg_solver.solve(
                        &matrix,
                        &rhs,
                        &mut p_correction_vec,
                        None::<&IdentityPreconditioner>,
                    )?;
                }
            }
            PressureLinearSolver::BiCGSTAB => {
                // Use AMG preconditioner if available, otherwise identity
                if let Some(ref amg) = amg_preconditioner {
                    self.bicgstab_solver
                        .solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                } else {
                    self.bicgstab_solver.solve(
                        &matrix,
                        &rhs,
                        &mut p_correction_vec,
                        None::<&IdentityPreconditioner>,
                    )?;
                }
            }
            PressureLinearSolver::GMRES { .. } => {
                if let Some(ref solver) = self.gmres_solver {
                    // Use AMG preconditioner if available, otherwise identity
                    if let Some(ref amg) = amg_preconditioner {
                        solver.solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                    } else {
                        solver.solve(
                            &matrix,
                            &rhs,
                            &mut p_correction_vec,
                            None::<&IdentityPreconditioner>,
                        )?;
                    }
                } else {
                    return Err(cfd_core::error::Error::InvalidConfiguration(
                        "GMRES solver not initialized".to_string(),
                    ));
                }
            }
        }

        // Convert back to 2D grid
        let mut p_correction = vec![vec![T::zero(); ny]; nx];
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                let value = if idx == reference_idx {
                    T::zero()
                } else if let Some(col_idx) = map_index(idx) {
                    p_correction_vec[col_idx]
                } else {
                    T::zero()
                };
                p_correction[i][j] = value;
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

    /// Solve pressure correction equation using face velocities (Rhie-Chow)
    ///
    /// The pressure correction equation is derived from the continuity equation:
    /// ∇²p' = ρ/Δt * ∇·u_face
    ///
    /// # Arguments
    /// * `u_face` - Face velocity u-component at east face (i+1/2, j)
    /// * `v_face` - Face velocity v-component at north face (i, j+1/2)
    /// * `dt` - Time step
    /// * `rho` - Density
    pub fn solve_pressure_correction_from_faces(
        &self,
        u_face: &Vec<Vec<T>>,
        v_face: &Vec<Vec<T>>,
        dt: T,
        rho: T,
    ) -> cfd_core::error::Result<Vec<Vec<T>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        // Build discrete Laplacian matrix and divergence source term
        let n = (nx - 2) * (ny - 2); // Interior points only

        if n == 0 {
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }
        if n == 1 {
            // Single interior point: fix pressure correction to zero
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }

        let reference_idx = 0usize;
        let system_size = n - 1;

        let map_index = |idx: usize| -> Option<usize> {
            if idx == reference_idx {
                None
            } else if idx < reference_idx {
                Some(idx)
            } else {
                Some(idx - 1)
            }
        };

        let mut builder = SparseMatrixBuilder::new(system_size, system_size);
        let mut rhs = DVector::zeros(system_size);

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);
        let coeff = rho / dt;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);

                if idx == reference_idx {
                    continue;
                }

                let row_idx = map_index(idx).expect("row index must exist");

                // Laplacian stencil - diagonal is negative sum of neighbor coefficients
                let two = T::from_f64(cfd_core::constants::numerical::common::TWO)
                    .unwrap_or_else(|| T::one() + T::one());
                let ap = -two * (dx2_inv + dy2_inv);

                builder.add_entry(row_idx, row_idx, ap)?;

                // Neighbors
                if i > 1 {
                    let neighbor_idx = idx - (ny - 2);
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dx2_inv)?;
                    }
                }
                if i < nx - 2 {
                    let neighbor_idx = idx + (ny - 2);
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dx2_inv)?;
                    }
                }
                if j > 1 {
                    let neighbor_idx = idx - 1;
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dy2_inv)?;
                    }
                }
                if j < ny - 2 {
                    let neighbor_idx = idx + 1;
                    if let Some(col_idx) = map_index(neighbor_idx) {
                        builder.add_entry(row_idx, col_idx, dy2_inv)?;
                    }
                }

                // Divergence of face velocity (Rhie-Chow)
                // u_face[i][j] is u_{i+1/2, j} (East)
                // u_face[i-1][j] is u_{i-1/2, j} (West)
                // v_face[i][j] is v_{i, j+1/2} (North)
                // v_face[i][j-1] is v_{i, j-1/2} (South)
                let div_u = (u_face[i][j] - u_face[i - 1][j]) / dx
                    + (v_face[i][j] - v_face[i][j - 1]) / dy;

                rhs[row_idx] = coeff * div_u;
            }
        }

        // Solve the linear system using selected solver
        let matrix = builder.build()?;
        let mut p_correction_vec = DVector::zeros(matrix.nrows());

        // Create AMG preconditioner for this solve
        let amg_config = AMGConfig {
            cycle_type: MultigridCycle::V,
            nu1: 2,
            nu2: 2,
            max_levels: 5,
            min_coarse_size: 10,
            coarsening: CoarseningStrategy::RugeStueben,
        };

        let amg_preconditioner = match AlgebraicMultigrid::with_config(&matrix, amg_config) {
            Ok(amg) => Some(amg),
            Err(_) => {
                tracing::debug!("AMG preconditioner construction failed, using identity");
                None
            }
        };

        match self.solver_type {
            PressureLinearSolver::ConjugateGradient => {
                if let Some(ref amg) = amg_preconditioner {
                    self.cg_solver
                        .solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                } else {
                    self.cg_solver.solve(
                        &matrix,
                        &rhs,
                        &mut p_correction_vec,
                        None::<&IdentityPreconditioner>,
                    )?;
                }
            }
            PressureLinearSolver::BiCGSTAB => {
                if let Some(ref amg) = amg_preconditioner {
                    self.bicgstab_solver
                        .solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                } else {
                    self.bicgstab_solver.solve(
                        &matrix,
                        &rhs,
                        &mut p_correction_vec,
                        None::<&IdentityPreconditioner>,
                    )?;
                }
            }
            PressureLinearSolver::GMRES { .. } => {
                if let Some(ref solver) = self.gmres_solver {
                    if let Some(ref amg) = amg_preconditioner {
                        solver.solve(&matrix, &rhs, &mut p_correction_vec, Some(amg))?;
                    } else {
                        solver.solve(
                            &matrix,
                            &rhs,
                            &mut p_correction_vec,
                            None::<&IdentityPreconditioner>,
                        )?;
                    }
                } else {
                    return Err(cfd_core::error::Error::InvalidConfiguration(
                        "GMRES solver not initialized".to_string(),
                    ));
                }
            }
        }

        // Convert back to 2D grid
        let mut p_correction = vec![vec![T::zero(); ny]; nx];
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                let value = if idx == reference_idx {
                    T::zero()
                } else if let Some(col_idx) = map_index(idx) {
                    p_correction_vec[col_idx]
                } else {
                    T::zero()
                };
                p_correction[i][j] = value;
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
                    / (T::from_f64(cfd_core::constants::numerical::common::TWO)
                        .unwrap_or_else(|| T::zero())
                        * dx);
                let dp_dy = (p_correction[i][j + 1] - p_correction[i][j - 1])
                    / (T::from_f64(cfd_core::constants::numerical::common::TWO)
                        .unwrap_or_else(|| T::zero())
                        * dy);

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
