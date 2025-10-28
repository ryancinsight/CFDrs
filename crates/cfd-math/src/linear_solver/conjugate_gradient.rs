//! Preconditioned Conjugate Gradient solver implementation

use super::config::IterativeSolverConfig;
use super::traits::{Configurable, IterativeLinearSolver, Preconditioner};
use crate::vector_ops::SimdVectorOps;
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// Preconditioned Conjugate Gradient solver with memory management
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

        // Zero-copy: initialize p directly instead of cloning r
        let mut p = DVector::zeros(n);
        p.copy_from(&r);

        for _ in 0..self.config.max_iterations {
            let ap = a * &p;

            // Use SIMD dot product for large vectors
            let pap = if n > 1000 {
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

            if r_norm_sq_new < self.config.tolerance * self.config.tolerance {
                return Ok(());
            }

            let beta = r_norm_sq_new / r_norm_sq;

            // Update p = r + beta * p
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }

            r_norm_sq = r_norm_sq_new;
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
