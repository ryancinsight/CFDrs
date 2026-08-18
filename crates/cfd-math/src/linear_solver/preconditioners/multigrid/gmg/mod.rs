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

mod multigrid;
mod ops;
mod transfer;

pub use multigrid::GeometricMultigrid;
use ops::vector_sub;
#[cfg(test)]
use ops::{from_usize, l2_norm, matrix_vector_product};

use crate::error::Result;
use eunomia::RealField;
use leto::{Array1, Array2};

type GmgVector<T> = Array1<T>;
type GmgMatrix<T> = Array2<T>;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_geometric_multigrid_creation() {
        let gmg = GeometricMultigrid::<f64>::new(16, 16, 4);
        let gmg = gmg.expect("invariant: positive grid dimensions create GMG");
        assert!(!gmg.grid_sizes.is_empty(), "Should have grid levels");
        assert!(!gmg.matrices.is_empty(), "Should have matrices");
        assert_eq!(gmg.grid_sizes[0], (16, 16), "Finest grid should be 16x16");
    }

    #[test]
    fn test_poisson_matrix_creation() {
        let matrix = GeometricMultigrid::<f64>::create_poisson_matrix(4, 4, 0.25);
        let matrix = matrix.expect("invariant: positive grid dimensions create Poisson matrix");
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
        let gmg = GeometricMultigrid::<f64>::new(3, 3, 2)
            .expect("invariant: positive grid dimensions create GMG");
        let fine_residual =
            GmgVector::from_shape_vec([9], (0usize..9).map(from_usize::<f64>).collect())
                .expect("invariant: residual data matches its shape");

        let restricted = gmg.restrict_residual(&fine_residual, 3, 3, 2, 2);

        assert_eq!(restricted.shape(), [4]);
        assert_eq!(restricted[0], 4.0 / 3.0);
        assert_eq!(restricted[1], 8.0 / 3.0);
        assert_eq!(restricted[2], 16.0 / 3.0);
        assert_eq!(restricted[3], 20.0 / 3.0);
    }

    #[test]
    fn test_geometric_multigrid_solve() {
        let mut gmg = GeometricMultigrid::<f64>::new(8, 8, 3)
            .expect("invariant: positive grid dimensions create GMG");

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

        let (solution, iterations, residual_norm) = gmg
            .solve(&rhs, 1e-6, 10)
            .expect("invariant: configured Poisson solve converges");

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
        let gmg = GeometricMultigrid::<f64>::new(16, 16, 3)
            .expect("invariant: positive grid dimensions create GMG");
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
        let (solution, iterations, residual_norm) = gmg
            .solve_fas(&operator, &rhs, tolerance, max_iter)
            .expect("invariant: configured FAS solve converges");

        assert!(iterations > 0);
        assert!(iterations <= max_iter);
        assert!(
            residual_norm < tolerance,
            "Residual norm {residual_norm} is not < {tolerance}"
        );

        // Compare with standard linear solve
        let (linear_sol, _, _) = gmg
            .solve(&rhs, tolerance, max_iter)
            .expect("invariant: configured linear solve converges");

        let diff = vector_sub(&solution, &linear_sol);
        // The solutions should be relatively close
        assert!(
            l2_norm(&diff) < 1e-3,
            "FAS solution should match linear solution for linear problem (diff: {})",
            l2_norm(&diff)
        );
    }
}
