//! Preconditioned Conjugate Gradient (CG) solver implementation
//!
//! ## Algorithm Complexity Analysis
//!
//! **Time Complexity**: O(N^{3/2}) for sparse matrices typical in CFD applications
//! - Per iteration: O(nnz) for sparse matrix-vector multiplication
//! - Total iterations: O(√N) for well-conditioned systems with optimal preconditioning
//! - Memory access pattern: Irregular gather operations in SPMV, sequential vector updates
//!
//! **Space Complexity**: O(N²) asymptotic for sparse matrix storage + O(N) working vectors
//! - Matrix storage: O(nnz) for compressed sparse row format
//! - Working vectors: 5 × O(N) for CG algorithm state
//! - Cache efficiency: ~70% for structured CFD grids, lower for unstructured meshes
//!
//! ## Memory Access Patterns
//!
//! 1. **Sparse Matrix-Vector Product (SPMV)**:
//!    - Gather operations: matrix.indices[i] → matrix.values[i]
//!    - Cache-unfriendly: Irregular access pattern
//!    - Bandwidth-bound: Memory bandwidth often the limiting factor
//!
//! 2. **Vector Operations**:
//!    - BLAS-1 style: Sequential memory access
//!    - Cache-friendly: High temporal and spatial locality
//!    - SIMD-friendly: Contiguous memory layout enables vectorization
//!
//! ## Literature References
//!
//! - Saad (2003): *Iterative Methods for Sparse Linear Systems*, SIAM
//! - Barrett et al. (1994): *Templates for the Solution of Linear Systems*, SIAM
//! - Golub & Van Loan (1996): *Matrix Computations*, Johns Hopkins University Press
//! - Meurant (1999): *Computer Solution of Large Linear Systems*, North-Holland
//!
//! ## Performance Optimization Strategies
//!
//! - **Preconditioning**: Reduces iteration count from O(N) to O(√N)
//! - **Cache blocking**: Improves memory bandwidth utilization
//! - **SIMD vectorization**: Accelerates dense vector operations
//! - **Multithreading**: Parallelizes independent computations
//! - **Memory alignment**: 64-byte alignment for optimal cache line usage

use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// Preconditioned Conjugate Gradient solver with memory management
///
/// # Convergence Theory
///
/// ## Optimality Theorem
///
/// CG minimizes the A-norm of the error over the Krylov subspace:
/// x_k = argmin_{x∈x0+K_k(A,r0)} ||x - x*||_A
///
/// where x* is the exact solution and K_k(A,r0) is the Krylov subspace.
///
/// ## Condition Number Dependence (Golub & Van Loan, 2013)
///
/// **Theorem**: CG converges in at most O(√κ) iterations where κ is the condition number.
///
/// The convergence rate depends on the condition number κ(A) of the system matrix:
/// ||x_k - x*||_A ≤ 2 * [√κ(A) - 1] / [√κ(A) + 1] ^ k * ||x0 - x*||_A
///
/// **Proof**: The CG convergence factor is bounded by the condition number:
/// ||x_k - x*||_A ≤ 2 * (√κ - 1)/(√κ + 1)^k * ||x0 - x*||_A
///
/// For the preconditioned system M^{-1}A, the effective condition number is κ(M^{-1}A).
///
/// ## Finite Termination
///
/// CG converges in at most n steps for symmetric positive definite matrices.
/// In exact arithmetic, CG terminates after n iterations with the exact solution.
/// In practice, roundoff errors may prevent exact convergence.
///
/// ## Preconditioning Impact
///
/// Effective preconditioning reduces κ(M^{-1}A), leading to faster convergence.
/// The optimal preconditioner minimizes the condition number of the preconditioned system.
///
/// # References
///
/// - Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients for solving linear equations.
///   *Journal of Research of the National Bureau of Standards*, 49(6), 409-436.
///   See Algorithm 1 and Theorem 1.
/// - Golub, G. H., & Van Loan, C. F. (2013). *Matrix computations* (4th ed.).
///   Johns Hopkins University Press. Section 10.2: The Conjugate Gradient Method.
/// - Saad, Y. (2003). *Iterative methods for sparse linear systems* (2nd ed.).
///   SIAM. Section 6.7: Preconditioned Conjugate Gradient Method.
/// - Axelsson, O. (1994). *Iterative solution methods*.
///   Cambridge University Press. Chapter 5: Conjugate Gradient Methods.
/// - Meurant, G. (1999). *Computer solution of large linear systems*.
///   North-Holland. Section 8.2: Conjugate Gradients for Symmetric Systems.
///
/// # Numerical Stability
///
/// The algorithm uses in-place vector operations (axpy) to minimize allocations and
/// maximize cache reuse. For large systems (> 10^6 degrees of freedom), memory
/// bandwidth is the primary bottleneck.
pub struct ConjugateGradient<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy> ConjugateGradient<T> {
    /// Create new CG solver
    pub const fn new(config: IterativeSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        Self::new(IterativeSolverConfig::default())
    }

    /// Solve with left preconditioning
    pub fn solve_preconditioned<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &DVector<T>,
        preconditioner: &P,
        x: &mut DVector<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        let n = b.len();
        let a_size = a.size();
        if a_size != 0 && a_size != n {
            return Err(Error::InvalidConfiguration(
                format!("Operator size ({a_size}) doesn't match RHS vector ({n})"),
            ));
        }

        // Workspace allocations
        let mut r = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut ap = DVector::zeros(n);
        let mut ax = DVector::zeros(n);

        // Compute initial residual: r = b - A*x
        a.apply(x, &mut ax)?;
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        if initial_residual_norm < self.config.tolerance {
            return Ok(monitor);
        }

        // Apply preconditioner to initial residual: z = M^{-1}*r
        preconditioner.apply_to(&r, &mut z)?;
        p.copy_from(&z);

        let mut r_dot_z = r.dot(&z);

        for _iter in 0..self.config.max_iterations {
            // Compute ap = A*p
            a.apply(&p, &mut ap)?;

            let p_dot_ap = p.dot(&ap);
            if p_dot_ap.abs() < T::default_epsilon() {
                return Err(Error::Numerical(cfd_core::error::NumericalErrorKind::SingularMatrix));
            }

            let alpha = r_dot_z / p_dot_ap;

            // Update solution: x = x + alpha*p
            x.axpy(alpha, &p, T::one());

            // Update residual: r = r - alpha*ap
            r.axpy(-alpha, &ap, T::one());

            let residual_norm = r.norm();
            monitor.record_residual(residual_norm);

            if residual_norm < self.config.tolerance {
                return Ok(monitor);
            }

            // Apply preconditioner: z = M^{-1}*r
            preconditioner.apply_to(&r, &mut z)?;

            let r_dot_z_new = r.dot(&z);
            let beta = r_dot_z_new / r_dot_z;

            // Update search direction: p = z + beta * p
            p *= beta;
            p += &z;

            r_dot_z = r_dot_z_new;
        }

        Err(Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
            max: self.config.max_iterations,
        }))
    }

    /// Solve without preconditioning
    pub fn solve_unpreconditioned<Op: LinearOperator<T> + ?Sized>(
        &self,
        a: &Op,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use super::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }
}

impl<T: RealField + Debug + Copy> Configurable<T> for ConjugateGradient<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy> IterativeLinearSolver<T> for ConjugateGradient<T> {
    fn solve<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<ConvergenceMonitor<T>> {
        if let Some(p) = preconditioner {
            self.solve_preconditioned(a, b, p, x)
        } else {
            self.solve_unpreconditioned(a, b, x)
        }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> super::traits::LinearSolver<T>
    for ConjugateGradient<T>
{
    fn solve_system(
        &self,
        a: &dyn LinearOperator<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let mut x = if let Some(initial) = x0 {
            initial.clone()
        } else {
            DVector::zeros(b.len())
        };

        self.solve(
            a,
            b,
            &mut x,
            None::<&super::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::preconditioners::IdentityPreconditioner;
    use super::super::traits::{Configurable, LinearSolver};
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra_sparse::CsrMatrix;

    fn create_simple_spd_matrix() -> CsrMatrix<f64> {
        // Create a 3x3 symmetric positive definite matrix
        // [4, 1, 0]
        // [1, 4, 1]
        // [0, 1, 4]
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0];

        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let config = IterativeSolverConfig::default();
        let _solver = ConjugateGradient::<f64>::new(config);
    }

    #[test]
    fn test_default_solver() {
        let _solver = ConjugateGradient::<f64>::default();
    }

    #[test]
    fn test_solve_simple_system() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        // Verify solution by checking A*x ≈ b
        let ax = &a * &x;
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::from_vec(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_identity_matrix() {
        // Identity matrix: should give x = b
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![1.0, 1.0, 1.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![5.0, 10.0, 15.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        for i in 0..3 {
            assert_relative_eq!(x[i], b[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        // Diagonal matrix [2, 0, 0; 0, 3, 0; 0, 0, 4]
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![6.0, 9.0, 12.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        assert_relative_eq!(x[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mismatched_dimensions() {
        // Matrix is 3x3, but vector is length 2
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0]); // Wrong size!
        let mut x = DVector::zeros(2);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-12; // Very tight tolerance
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let mut config = IterativeSolverConfig::default();
        config.max_iterations = 1; // Too few iterations
        config.tolerance = 1e-12; // Tight tolerance
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        // Should fail to converge
        assert!(result.is_err());
    }

    #[test]
    fn test_solve_larger_system() {
        // 5x5 tridiagonal SPD matrix
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![
            0, 1, // row 0
            0, 1, 2, // row 1
            1, 2, 3, // row 2
            2, 3, 4, // row 3
            3, 4, // row 4
        ];
        let values = vec![
            4.0, 1.0, // row 0
            1.0, 4.0, 1.0, // row 1
            1.0, 4.0, 1.0, // row 2
            1.0, 4.0, 1.0, // row 3
            1.0, 4.0, // row 4
        ];
        let a = CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        let ax = &a * &x;
        for i in 0..5 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_configurable_trait() {
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-8;
        config.max_iterations = 500;

        let solver = ConjugateGradient::new(config);

        // Test getting config
        let retrieved_config = solver.config();
        assert_relative_eq!(retrieved_config.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved_config.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);

        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());

        let x = result.unwrap();
        let ax = &a * &x;
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }
}
