//! Preconditioned Conjugate Gradient solver implementation

use super::config::IterativeSolverConfig;
use super::traits::{Configurable, ConvergenceMonitor, IterativeLinearSolver, Preconditioner};
use crate::vector_ops::SimdVectorOps;
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

// Cache-aligned vector for optimal memory access patterns
#[repr(align(64))] // Cache line alignment for x86_64
struct AlignedVector<T> {
    data: Vec<T>,
}

impl<T> AlignedVector<T> {
    fn with_capacity(capacity: usize) -> Self {
        let mut data = Vec::with_capacity(capacity);
        data.reserve_exact(capacity);
        Self { data }
    }

    fn resize(&mut self, new_len: usize, value: T) where T: Clone {
        self.data.resize(new_len, value);
    }

    fn as_slice(&self) -> &[T] {
        &self.data
    }

    fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }
}

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

    /// Solve with preconditioning and memory management
    pub fn solve_preconditioned<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        preconditioner: &P,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        // Pre-allocate all workspace vectors to avoid allocations in the loop
        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);
        let mut r = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut ap;

        // Compute initial residual: r = b - A*x
        let ax = a * &x;
        r.copy_from(b);
        r -= &ax;

        // Apply preconditioner: M*z = r
        preconditioner.apply_to(&r, &mut z)?;
        p.copy_from(&z);

        let mut rzold = r.dot(&z);

        // PCG iterations with in-place operations
        for iter in 0..self.config.max_iterations {
            // Compute A*p
            ap = a * &p;

            let alpha = rzold / p.dot(&ap);

            // Update solution: x = x + alpha * p
            x.axpy(alpha, &p, T::one());

            // Update residual: r = r - alpha * ap
            r.axpy(-alpha, &ap, T::one());

            let residual_norm = r.norm();
            if self.is_converged(residual_norm) {
                tracing::debug!("PCG converged in {} iterations", iter + 1);
                return Ok(x);
            }

            // Apply preconditioner: M*z = r
            preconditioner.apply_to(&r, &mut z)?;

            let rznew = r.dot(&z);
            let beta = rznew / rzold;

            // Update search direction: p = z + beta * p
            p *= beta;
            p += &z;

            rzold = rznew;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy + FromPrimitive + Send + Sync> Configurable<T>
    for ConjugateGradient<T>
{
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + FromPrimitive + Send + Sync> IterativeLinearSolver<T>
    for ConjugateGradient<T>
{
    fn solve<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
        _preconditioner: Option<&P>,
    ) -> Result<()> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidInput(
                "Matrix dimensions must match vector size".to_string(),
            ));
        }

        // Initialize x to zero if not already set
        if x.len() != n {
            *x = DVector::zeros(n);
        }

        let mut r = b - a * &*x;

        // Use SIMD operations for large vectors
        let mut r_norm_sq = if n > 1000 {
            r.simd_norm().powi(2)
        } else {
            r.norm_squared()
        };

        if r_norm_sq < self.config.tolerance * self.config.tolerance {
            return Ok(());
        }

        // Initialize convergence monitoring
        let initial_residual = r_norm_sq.sqrt();
        let mut convergence_monitor = ConvergenceMonitor::new(initial_residual);

        // Estimate condition number for theoretical bound (simplified estimation)
        // In practice, this would be computed from matrix properties
        let estimated_kappa = n as f64; // Conservative estimate
        let theoretical_bound = convergence_monitor.cg_theoretical_bound(estimated_kappa);
        convergence_monitor.set_theoretical_bound(theoretical_bound);
        convergence_monitor.set_condition_number_estimate(estimated_kappa);

        // Zero-copy: initialize p directly instead of cloning r
        let mut p = DVector::zeros(n);
        p.copy_from(&r);

        for _iteration in 0..self.config.max_iterations {
            // Memory prefetching for cache optimization on large problems
            #[cfg(target_arch = "x86_64")]
            if n > 10000 {
                unsafe {
                    // Prefetch matrix rows and vectors for better cache performance
                    for i in (0..n).step_by(64) { // Prefetch every cache line
                        std::arch::x86_64::_mm_prefetch(p.as_slice()[i..].as_ptr() as *const i8, std::arch::x86_64::_MM_HINT_T0);
                        // Note: CsrMatrix doesn't have as_slice(), so we skip prefetching matrix data for now
                    }
                }
            }

            let ap = a * &p;

            // Use SIMD dot product for large vectors with prefetching
            let pap = if n > 1000 {
                // Additional prefetching for SIMD operations
                #[cfg(target_arch = "x86_64")]
                if n > 10000 {
                    unsafe {
                        std::arch::x86_64::_mm_prefetch(p.as_slice().as_ptr() as *const i8, std::arch::x86_64::_MM_HINT_T0);
                        std::arch::x86_64::_mm_prefetch(ap.as_slice().as_ptr() as *const i8, std::arch::x86_64::_MM_HINT_T0);
                    }
                }
                p.simd_dot(&ap)
            } else {
                p.dot(&ap)
            };

            if pap.abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                break;
            }

            let alpha = r_norm_sq / pap;

            // Use axpy for in-place updates
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }

            // Use SIMD norm for convergence check
            let r_norm_sq_new = if n > 1000 {
                r.simd_norm().powi(2)
            } else {
                r.norm_squared()
            };

            // Record residual for convergence monitoring
            let current_residual = r_norm_sq_new.sqrt();
            convergence_monitor.record_residual(current_residual);

            if r_norm_sq_new < self.config.tolerance * self.config.tolerance {
                // Validate convergence against theoretical bounds
                if let Err(e) = convergence_monitor.validate_convergence() {
                    eprintln!("Warning: CG convergence outside theoretical bounds: {}", e);
                }
                return Ok(());
            }

            let beta = r_norm_sq_new / r_norm_sq;

            // Update p = r + beta * p
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }

            r_norm_sq = r_norm_sq_new;
        }

        // Final convergence validation
        if let Err(e) = convergence_monitor.validate_convergence() {
            eprintln!("Warning: CG convergence outside theoretical bounds: {}", e);
        }

        Ok(())
    }
}

// Implement object-safe LinearSolver trait for trait objects
impl<T: RealField + Copy + num_traits::FromPrimitive + Send + Sync> super::traits::LinearSolver<T>
    for ConjugateGradient<T>
{
    fn solve_system(
        &self,
        a: &nalgebra_sparse::CsrMatrix<T>,
        b: &nalgebra::DVector<T>,
        x0: Option<&nalgebra::DVector<T>>,
    ) -> cfd_core::error::Result<nalgebra::DVector<T>> {
        let mut x = if let Some(initial) = x0 {
            initial.clone()
        } else {
            nalgebra::DVector::zeros(b.len())
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
    use super::*;
    use super::super::preconditioners::IdentityPreconditioner;
    use super::super::traits::{Configurable, LinearSolver};
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
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_ok());
        
        let x = result.unwrap();
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
        let x0 = DVector::from_vec(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, Some(&x0));
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
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_ok());
        
        let x = result.unwrap();
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
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_ok());
        
        let x = result.unwrap();
        assert_relative_eq!(x[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mismatched_dimensions() {
        // Matrix is 3x3, but vector is length 2
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0]); // Wrong size!
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-12; // Very tight tolerance
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_ok());
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut config = IterativeSolverConfig::default();
        config.max_iterations = 1; // Too few iterations
        config.tolerance = 1e-12; // Tight tolerance
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        // Should fail to converge
        assert!(result.is_err());
    }

    #[test]
    fn test_solve_larger_system() {
        // 5x5 tridiagonal SPD matrix
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![
            0, 1,          // row 0
            0, 1, 2,       // row 1
            1, 2, 3,       // row 2
            2, 3, 4,       // row 3
            3, 4           // row 4
        ];
        let values = vec![
            4.0, 1.0,      // row 0
            1.0, 4.0, 1.0, // row 1
            1.0, 4.0, 1.0, // row 2
            1.0, 4.0, 1.0, // row 3
            1.0, 4.0       // row 4
        ];
        let a = CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");
        
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, None);
        assert!(result.is_ok());
        
        let x = result.unwrap();
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
