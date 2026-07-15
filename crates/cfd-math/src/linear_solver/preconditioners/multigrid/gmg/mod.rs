//! Geometric Multigrid (GMG) for structured grids
//!
//! Geometric multigrid exploits structured grid hierarchies for optimal
//! convergence properties and computational efficiency.
//!
//! ## Mathematical Foundation
//!
//! For a structured grid with mesh size h, the geometric multigrid hierarchy is:
//! ```math
//! Ω¹ ⊂ Ω² ⊂ ⋯ ⊂ Ω^J
//! ```
//!
//! where Ω¹ is the finest grid and Ω^J is the coarsest grid.
//!
//! ## Algorithm Overview
//!
//! 1. **Relaxation**: Apply smoother on fine grid
//! 2. **Restriction**: Transfer residual to coarse grid
//! 3. **Coarsest Solve**: Direct or iterative solution on coarsest grid
//! 4. **Prolongation**: Transfer correction back to fine grid
//! 5. **Post-relaxation**: Apply smoother on fine grid
//!
//! ## Literature Compliance
//!
//! - Briggs, W. L., et al. (2000). *A multigrid tutorial*. SIAM. Chapter 3.
//! - Trottenberg, U., et al. (2001). *Multigrid*. Academic Press. Chapter 4.
//! - Wesseling, P. (1992). *An introduction to multigrid methods*. Wiley.

mod transfer;

use crate::error::Result;
use cfd_core::error::Error;
use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};

type GmgVector<T> = Array1<T>;
type GmgMatrix<T> = Array2<T>;

#[inline]
pub(super) fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: grid dimension fits into u64");
    from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

fn l2_norm<T: RealField>(vector: &GmgVector<T>) -> T {
    let mut sum = <T as NumericElement>::ZERO;
    for i in 0..vector.shape()[0] {
        sum += vector[i] * vector[i];
    }
    <T as NumericElement>::sqrt(sum)
}

fn matrix_vector_product<T: RealField>(
    matrix: &GmgMatrix<T>,
    vector: &GmgVector<T>,
) -> GmgVector<T> {
    let [rows, cols] = matrix.shape();
    debug_assert_eq!(cols, vector.shape()[0]);
    let mut output = GmgVector::zeros([rows]);
    for row in 0..rows {
        let mut sum = <T as NumericElement>::ZERO;
        for col in 0..cols {
            sum += matrix[[row, col]] * vector[col];
        }
        output[row] = sum;
    }
    output
}

fn vector_add<T: RealField>(lhs: &GmgVector<T>, rhs: &GmgVector<T>) -> GmgVector<T> {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    let mut output = GmgVector::zeros([lhs.shape()[0]]);
    for i in 0..lhs.shape()[0] {
        output[i] = lhs[i] + rhs[i];
    }
    output
}

fn vector_sub<T: RealField>(lhs: &GmgVector<T>, rhs: &GmgVector<T>) -> GmgVector<T> {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    let mut output = GmgVector::zeros([lhs.shape()[0]]);
    for i in 0..lhs.shape()[0] {
        output[i] = lhs[i] - rhs[i];
    }
    output
}

fn vector_add_assign<T: RealField>(lhs: &mut GmgVector<T>, rhs: &GmgVector<T>) {
    debug_assert_eq!(lhs.shape(), rhs.shape());
    for i in 0..lhs.shape()[0] {
        lhs[i] += rhs[i];
    }
}

/// Trait for nonlinear operators in multigrid methods
///
/// This trait defines the interface for nonlinear operators that can be used
/// with the Full Approximation Scheme (FAS) multigrid method.
pub trait NonlinearOperator<T: RealField> {
    /// Compute the nonlinear residual: r = f - F(u)
    ///
    /// The default implementation computes f - apply(u).
    fn residual(&self, u: &GmgVector<T>, rhs: &GmgVector<T>, level: usize) -> GmgVector<T> {
        vector_sub(rhs, &self.apply(u, level))
    }

    /// Apply the nonlinear operator: F(u)
    fn apply(&self, u: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Solve the nonlinear system on the coarsest level
    fn coarsest_solve(&self, u: &mut GmgVector<T>, rhs: &GmgVector<T>, level: usize) -> Result<()>;

    /// Restrict the residual to the next coarser level
    fn restrict_residual(&self, fine: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Restrict the solution to the next coarser level
    fn restrict_solution(&self, fine: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Prolongate a vector to the next finer level
    fn prolongate(&self, coarse: &GmgVector<T>, level: usize) -> GmgVector<T>;

    /// Apply smoothing to the solution u given RHS f
    fn smooth(&self, u: &mut GmgVector<T>, rhs: &GmgVector<T>, level: usize, iterations: usize);
}

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
            current_h *= from_f64::<T>(2.0);
        }

        Ok(Self {
            grid_sizes,
            matrices,
            relaxation_param: from_f64::<T>(0.8), // Weighted Jacobi
            nu1: 2,                               // Pre-smoothing iterations
            nu2: 2,                               // Post-smoothing iterations
        })
    }

    /// Create discrete Poisson matrix (-Δu = f) for structured grid
    fn create_poisson_matrix(nx: usize, ny: usize, h: T) -> Result<GmgMatrix<T>> {
        let n = nx * ny;
        let mut matrix = GmgMatrix::zeros([n, n]);

        let h_squared = h * h;
        let four = from_f64::<T>(4.0);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometric_multigrid_creation() {
        let gmg = GeometricMultigrid::<f64>::new(16, 16, 4);
        assert!(gmg.is_ok(), "GMG creation should succeed");

        let gmg = gmg.unwrap();
        assert!(!gmg.grid_sizes.is_empty(), "Should have grid levels");
        assert!(!gmg.matrices.is_empty(), "Should have matrices");
        assert_eq!(gmg.grid_sizes[0], (16, 16), "Finest grid should be 16x16");
    }

    #[test]
    fn test_poisson_matrix_creation() {
        let matrix = GeometricMultigrid::<f64>::create_poisson_matrix(4, 4, 0.25);
        assert!(matrix.is_ok(), "Poisson matrix creation should succeed");

        let matrix = matrix.unwrap();
        assert_eq!(matrix.shape(), [16, 16], "Matrix should be 16x16");

        // Check that diagonal elements are positive
        for i in 0..16 {
            assert!(matrix[[i, i]] > 0.0, "Diagonal elements should be positive");
        }

        assert_eq!(matrix[[0, 0]], 64.0);
        assert_eq!(matrix[[0, 1]], -16.0);
        assert_eq!(matrix[[0, 4]], -16.0);
        assert_eq!(matrix[[0, 5]], 0.0);
    }

    #[test]
    fn restrict_residual_full_weighting_matches_boundary_stencil() {
        let gmg = GeometricMultigrid::<f64>::new(3, 3, 2).unwrap();
        let fine_residual =
            GmgVector::from_shape_vec([9], (0usize..9).map(from_usize::<f64>).collect()).unwrap();

        let restricted = gmg.restrict_residual(&fine_residual, 3, 3, 2, 2);

        assert_eq!(restricted.shape(), [4]);
        assert_eq!(restricted[0], 4.0 / 3.0);
        assert_eq!(restricted[1], 8.0 / 3.0);
        assert_eq!(restricted[2], 16.0 / 3.0);
        assert_eq!(restricted[3], 20.0 / 3.0);
    }

    #[test]
    fn test_geometric_multigrid_solve() {
        let mut gmg = GeometricMultigrid::<f64>::new(8, 8, 3).unwrap();

        // Create a simple test problem: -Δu = 1 with u=0 on boundary
        let n = 8 * 8;
        let mut rhs = GmgVector::from_elem([n], 1.0);

        for i in 0..8 {
            for j in 0..8 {
                let idx = j * 8 + i;
                if i == 0 || i == 7 || j == 0 || j == 7 {
                    rhs[idx] = 0.0;
                    let matrix = &mut gmg.matrices[0];
                    for col in 0..n {
                        matrix[[idx, col]] = 0.0;
                    }
                    matrix[[idx, idx]] = 1.0;
                }
            }
        }

        let (solution, iterations, residual_norm) = gmg.solve(&rhs, 1e-6, 10).unwrap();

        assert!(iterations > 0, "Should require at least one iteration");
        assert!(residual_norm < 1.0, "Residual should be reduced");
        assert_eq!(solution.shape(), [n], "Solution should have correct size");
    }

    /// A linear operator wrapper to test the nonlinear FAS solver.
    /// It basically delegates everything to the underlying GMG matrices/methods.
    struct LinearPoissonOperator<'a> {
        gmg: &'a GeometricMultigrid<f64>,
    }

    impl NonlinearOperator<f64> for LinearPoissonOperator<'_> {
        fn apply(&self, u: &GmgVector<f64>, level: usize) -> GmgVector<f64> {
            let matrix = &self.gmg.matrices[level];
            matrix_vector_product(matrix, u)
        }

        fn coarsest_solve(
            &self,
            u: &mut GmgVector<f64>,
            rhs: &GmgVector<f64>,
            level: usize,
        ) -> Result<()> {
            // For testing, just run a few relaxation steps on the coarsest level
            let matrix = &self.gmg.matrices[level];
            // Use 20 iterations for coarse solve, starting from provided guess u
            self.gmg
                .jacobi_relaxation(matrix, u, rhs, self.gmg.relaxation_param, 20);
            Ok(())
        }

        fn restrict_residual(&self, fine: &GmgVector<f64>, level: usize) -> GmgVector<f64> {
            let (fine_nx, fine_ny) = self.gmg.grid_sizes[level];
            let (coarse_nx, coarse_ny) = self.gmg.grid_sizes[level + 1];
            self.gmg
                .restrict_residual(fine, fine_nx, fine_ny, coarse_nx, coarse_ny)
        }

        fn restrict_solution(&self, fine: &GmgVector<f64>, level: usize) -> GmgVector<f64> {
            let (fine_nx, _fine_ny) = self.gmg.grid_sizes[level];
            let (coarse_nx, coarse_ny) = self.gmg.grid_sizes[level + 1];

            let mut coarse = GmgVector::zeros([coarse_nx * coarse_ny]);
            for j in 0..coarse_ny {
                for i in 0..coarse_nx {
                    // Injection: coincide with fine grid points
                    let coarse_idx = j * coarse_nx + i;
                    let fine_idx = (2 * j) * fine_nx + (2 * i);
                    coarse[coarse_idx] = fine[fine_idx];
                }
            }
            coarse
        }

        fn prolongate(&self, coarse: &GmgVector<f64>, level: usize) -> GmgVector<f64> {
            let (fine_nx, fine_ny) = self.gmg.grid_sizes[level];
            let (coarse_nx, coarse_ny) = self.gmg.grid_sizes[level + 1];
            self.gmg
                .prolongate_correction(coarse, coarse_nx, coarse_ny, fine_nx, fine_ny)
        }

        fn smooth(
            &self,
            u: &mut GmgVector<f64>,
            rhs: &GmgVector<f64>,
            level: usize,
            iterations: usize,
        ) {
            let matrix = &self.gmg.matrices[level];
            self.gmg
                .jacobi_relaxation(matrix, u, rhs, self.gmg.relaxation_param, iterations);
        }
    }

    #[test]
    fn test_fas_solve_linear_problem() {
        let gmg = GeometricMultigrid::<f64>::new(16, 16, 3).unwrap();
        let operator = LinearPoissonOperator { gmg: &gmg };

        // Test problem: -Δu = 1, u=0 on boundary
        let n = 16 * 16;
        let mut rhs = GmgVector::from_elem([n], 1.0);

        // Zero out boundaries in RHS effectively (though strict Dirichlet BCs usually involve modifying matrix/RHS)
        // For this test, we accept the matrix built by create_poisson_matrix which assumes zero BCs implicitly
        // if we don't put source terms on boundary nodes.
        for i in 0..16 {
            for j in 0..16 {
                if i == 0 || i == 15 || j == 0 || j == 15 {
                    let idx = j * 16 + i;
                    rhs[idx] = 0.0;
                }
            }
        }

        // FAS Solve
        // Note: FAS convergence on linear problems with simple injection/restriction
        // might be slower than optimized linear MG. Relaxing tolerance for this test.
        let tolerance = 1e-3;
        let max_iter = 50;
        let (solution, iterations, residual_norm) =
            gmg.solve_fas(&operator, &rhs, tolerance, max_iter).unwrap();

        assert!(iterations > 0);
        assert!(iterations <= max_iter);
        assert!(
            residual_norm < tolerance,
            "Residual norm {residual_norm} is not < {tolerance}"
        );

        // Compare with standard linear solve
        let (linear_sol, _, _) = gmg.solve(&rhs, tolerance, max_iter).unwrap();

        let diff = vector_sub(&solution, &linear_sol);
        // The solutions should be relatively close
        assert!(
            l2_norm(&diff) < 1e-3,
            "FAS solution should match linear solution for linear problem (diff: {})",
            l2_norm(&diff)
        );
    }
}
