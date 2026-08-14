use super::ops::{from_usize, l2_norm, vector_add, vector_add_assign, vector_sub};
use super::{GmgMatrix, GmgVector, NonlinearOperator};
use crate::error::Result;
use cfd_core::error::Error;
use eunomia::{FloatElement, NumericElement, RealField};

/// Geometric multigrid hierarchy for structured grids
#[derive(Debug, Clone)]
pub struct GeometricMultigrid<T: RealField> {
    /// Grid dimensions for each level (finest to coarsest)
    pub(super) grid_sizes: Vec<(usize, usize)>,
    /// System matrices for each level
    pub(super) matrices: Vec<GmgMatrix<T>>,
    /// Relaxation parameter (ω for weighted Jacobi/SOR)
    pub(super) relaxation_param: T,
    /// Number of pre-smoothing iterations
    pub(super) nu1: usize,
    /// Number of post-smoothing iterations
    pub(super) nu2: usize,
}

impl<T: RealField + FloatElement> GeometricMultigrid<T> {
    /// Create geometric multigrid hierarchy for 2D Poisson equation
    ///
    /// # Arguments
    ///
    /// * `nx`, `ny` - Dimensions of finest grid
    /// * `max_levels` - Maximum number of multigrid levels
    ///
    /// # Returns
    ///
    /// Geometric multigrid solver configured for 2D Poisson equation
    pub fn new(nx: usize, ny: usize, max_levels: usize) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConfiguration(
                "Grid dimensions must be positive".to_string(),
            ));
        }

        let mut grid_sizes = Vec::new();
        let mut matrices = Vec::new();

        // Build hierarchy from finest to coarsest
        let mut current_nx = nx;
        let mut current_ny = ny;
        let mut current_h = <T as NumericElement>::ONE / from_usize::<T>(nx.max(ny));

        for _level in 0..max_levels {
            grid_sizes.push((current_nx, current_ny));

            // Create discrete Laplacian matrix for this level
            let matrix = Self::create_poisson_matrix(current_nx, current_ny, current_h)?;
            matrices.push(matrix);

            // Stop if grid is too small for further coarsening
            if current_nx <= 3 || current_ny <= 3 {
                break;
            }

            // Coarsen grid (simple 2:1 coarsening)
            current_nx = current_nx.div_ceil(2);
            current_ny = current_ny.div_ceil(2);
            current_h *= <T as FloatElement>::from_f64(2.0);
        }

        Ok(Self {
            grid_sizes,
            matrices,
            relaxation_param: <T as FloatElement>::from_f64(0.8), // Weighted Jacobi
            nu1: 2,                                               // Pre-smoothing iterations
            nu2: 2,                                               // Post-smoothing iterations
        })
    }

    /// Create discrete Poisson matrix (-Δu = f) for structured grid
    pub(super) fn create_poisson_matrix(nx: usize, ny: usize, h: T) -> Result<GmgMatrix<T>> {
        let n = nx * ny;
        let mut matrix = GmgMatrix::zeros([n, n]);

        let h_squared = h * h;
        let four = <T as FloatElement>::from_f64(4.0);
        let minus_one = -<T as NumericElement>::ONE;

        for i in 0..nx {
            for j in 0..ny {
                let idx = j * nx + i;

                // Diagonal element
                matrix[[idx, idx]] = four / h_squared;

                // Off-diagonal elements (5-point stencil)
                if i > 0 {
                    let left_idx = j * nx + (i - 1);
                    matrix[[idx, left_idx]] = minus_one / h_squared;
                }
                if i < nx - 1 {
                    let right_idx = j * nx + (i + 1);
                    matrix[[idx, right_idx]] = minus_one / h_squared;
                }
                if j > 0 {
                    let down_idx = (j - 1) * nx + i;
                    matrix[[idx, down_idx]] = minus_one / h_squared;
                }
                if j < ny - 1 {
                    let up_idx = (j + 1) * nx + i;
                    matrix[[idx, up_idx]] = minus_one / h_squared;
                }
            }
        }

        Ok(matrix)
    }

    /// Solve nonlinear system using Full Approximation Scheme (FAS) multigrid
    ///
    /// ## Full Approximation Scheme (FAS) Algorithm
    ///
    /// **Mathematical Foundation**: For nonlinear problems F(u) = f, FAS solves
    /// the full nonlinear problem on all grid levels rather than just the residual equation.
    ///
    /// **Key Difference from Linear MG**: Instead of restricting the residual r = f - Au,
    /// FAS restricts the solution u and computes the coarse-level right-hand side as:
    ///
    /// f² = R(r) + A²(R(u)) where R is the restriction operator
    ///
    /// **Algorithm Steps**:
    /// 1. Relax on fine grid: Apply nonlinear smoother to reduce high-frequency errors
    /// 2. Restrict solution: u² = R(u)
    /// 3. Compute FAS right-hand side: f² = R(f - F(u) + F²(u²))
    /// 4. Recursively solve coarse problem: F²(u²) = f²
    /// 5. Prolongate and correct: u += P(u² - R(u))
    /// 6. Post-relax on fine grid
    ///
    /// **Convergence Theory**: FAS converges for problems where the nonlinear operator
    /// satisfies appropriate smoothing and approximation properties.
    ///
    /// **Literature**: Brandt (1977), Trottenberg et al. (2001) Section 9.3
    pub fn solve_fas<Op: NonlinearOperator<T>>(
        &self,
        operator: &Op,
        rhs: &GmgVector<T>,
        tolerance: T,
        max_iterations: usize,
    ) -> Result<(GmgVector<T>, usize, T)> {
        if rhs.shape()[0] != self.matrices[0].shape()[0] {
            return Err(Error::InvalidConfiguration(
                "RHS vector size mismatch".to_string(),
            ));
        }

        let n = rhs.shape()[0];
        let mut u = GmgVector::zeros([n]);

        // Initial residual computation
        let initial_residual = operator.residual(&u, rhs, 0);
        let mut residual_norm = l2_norm(&initial_residual);

        for iteration in 0..max_iterations {
            // Perform one FAS V-cycle
            self.fas_v_cycle(operator, &mut u, rhs, 0);

            // Check convergence using nonlinear residual
            let new_residual = operator.residual(&u, rhs, 0);
            let new_residual_norm = l2_norm(&new_residual);

            if new_residual_norm < tolerance {
                return Ok((u, iteration + 1, new_residual_norm));
            }

            residual_norm = new_residual_norm;
        }

        Ok((u, max_iterations, residual_norm))
    }

    /// Solve linear system using standard geometric multigrid
    pub fn solve(
        &self,
        rhs: &GmgVector<T>,
        tolerance: T,
        max_iterations: usize,
    ) -> Result<(GmgVector<T>, usize, T)> {
        if rhs.shape()[0] != self.matrices[0].shape()[0] {
            return Err(Error::InvalidConfiguration(
                "RHS vector size mismatch".to_string(),
            ));
        }

        let n = rhs.shape()[0];
        let mut u = GmgVector::zeros([n]);
        let mut residual_norm = self.compute_residual_norm(&self.matrices[0], &u, rhs);

        for iteration in 0..max_iterations {
            // Perform one V-cycle
            self.v_cycle(&mut u, rhs, 0);

            // Check convergence
            let new_residual_norm = self.compute_residual_norm(&self.matrices[0], &u, rhs);

            if new_residual_norm < tolerance {
                return Ok((u, iteration + 1, new_residual_norm));
            }

            residual_norm = new_residual_norm;
        }

        Ok((u, max_iterations, residual_norm))
    }

    /// Perform one V-cycle iteration
    fn v_cycle(&self, u: &mut GmgVector<T>, f: &GmgVector<T>, level: usize) {
        let current_matrix = &self.matrices[level];
        let (nx, ny) = self.grid_sizes[level];

        // Pre-smoothing
        self.jacobi_relaxation(current_matrix, u, f, self.relaxation_param, self.nu1);

        // Compute residual
        let residual = self.compute_residual(current_matrix, u, f);

        // Restrict to coarse grid (if not coarsest level)
        if level < self.matrices.len() - 1 {
            let coarse_residual = self.restrict_residual(
                &residual,
                nx,
                ny,
                self.grid_sizes[level + 1].0,
                self.grid_sizes[level + 1].1,
            );

            // Recursively solve on coarse grid
            let mut coarse_correction = GmgVector::zeros([coarse_residual.shape()[0]]);
            self.v_cycle(&mut coarse_correction, &coarse_residual, level + 1);

            // Prolongate correction back to fine grid
            let fine_correction = self.prolongate_correction(
                &coarse_correction,
                self.grid_sizes[level + 1].0,
                self.grid_sizes[level + 1].1,
                nx,
                ny,
            );

            // Add correction
            vector_add_assign(u, &fine_correction);
        } else {
            // Coarsest level: solve directly (or use iterative method)
            // For now, use additional Jacobi iterations
            self.jacobi_relaxation(current_matrix, u, f, self.relaxation_param, 10);
        }

        // Post-smoothing
        self.jacobi_relaxation(current_matrix, u, f, self.relaxation_param, self.nu2);
    }

    /// Perform one FAS V-cycle for nonlinear problems
    fn fas_v_cycle<Op: NonlinearOperator<T>>(
        &self,
        operator: &Op,
        u: &mut GmgVector<T>,
        f: &GmgVector<T>,
        level: usize,
    ) {
        // Pre-smoothing: Apply nonlinear relaxation
        operator.smooth(u, f, level, self.nu1);

        // If not the coarsest level, restrict and solve on coarse grid
        if level < self.matrices.len() - 1 {
            // Restrict solution to coarse grid
            let coarse_u_restricted = operator.restrict_solution(u, level);

            // Compute fine grid residual: r = f - F(u)
            let fine_residual = operator.residual(u, f, level);

            // Restrict the residual
            let coarse_residual_restricted = operator.restrict_residual(&fine_residual, level);

            // Compute F_coarse(u_coarse_restricted)
            let coarse_operator_applied = operator.apply(&coarse_u_restricted, level + 1);

            // Compute FAS right-hand side: f_coarse = R(r) + F_coarse(R(u))
            let coarse_rhs = vector_add(&coarse_residual_restricted, &coarse_operator_applied);

            // Recursively solve on coarse grid
            let mut coarse_correction = coarse_u_restricted.clone(); // Start with restricted solution
            self.fas_v_cycle(operator, &mut coarse_correction, &coarse_rhs, level + 1);

            // Compute coarse grid correction: e_coarse = u_coarse - R(u)
            let coarse_correction_delta = vector_sub(&coarse_correction, &coarse_u_restricted);

            // Prolongate correction to fine grid
            let fine_correction = operator.prolongate(&coarse_correction_delta, level);

            // Add correction to fine grid solution
            vector_add_assign(u, &fine_correction);
        } else {
            // Coarsest level: solve nonlinear system directly
            let _ = operator.coarsest_solve(u, f, level);
        }

        // Post-smoothing
        operator.smooth(u, f, level, self.nu2);
    }
}
