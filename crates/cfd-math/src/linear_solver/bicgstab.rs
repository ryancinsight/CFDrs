//! BiCGSTAB solver implementation

use super::config::IterativeSolverConfig;
use super::preconditioners::IdentityPreconditioner;
use super::traits::{Configurable, IterativeLinearSolver, Preconditioner};
use crate::sparse::{spmv, spmv_parallel};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// BiCGSTAB solver with efficient memory management
///
/// # Convergence Theory
///
/// ## Acceleration over BiCG
///
/// BiCGSTAB (BiConjugate Gradient Stabilized) addresses instability issues in the
/// standard BiCG method by introducing additional stabilization terms. The algorithm
/// combines BiCG with GMRES-like stabilization to prevent divergence.
///
/// Unlike BiCG, which may exhibit erratic convergence or breakdown, BiCGSTAB
/// provides smoother convergence behavior and better numerical stability.
///
/// ## Convergence Properties
///
/// BiCGSTAB converges for non-symmetric matrices where the field of values
/// lies in the right half-plane. The convergence rate depends on:
/// - Condition number of the system matrix
/// - Effectiveness of preconditioning
/// - Initial residual smoothness
///
/// ## Breakdown Prevention
///
/// BiCGSTAB handles two types of breakdown:
/// 1. **Primary breakdown**: ρ_new = 0 (division by zero)
/// 2. **Secondary breakdown**: α = 0 or ω = 0 (stagnation)
///
/// Robust implementations include checks for near-breakdown conditions.
///
/// # References
///
/// - Van der Vorst, H. A. (1992). Bi-CGSTAB: A fast and smoothly converging variant
///   of Bi-CG for the solution of nonsymmetric linear systems. *SIAM Journal on
///   Scientific and Statistical Computing*, 13(2), 631-644.
///   See Algorithm 1 and convergence analysis in Section 3.
/// - Sleijpen, G. L. G., & Fokkema, D. R. (1993). BiCGstab(l) and other hybrid
///   Bi-CG methods. *Numerical Algorithms*, 7(1), 75-109.
///   See Section 2.1 for BiCGSTAB algorithm derivation.
/// - Saad, Y. (2003). *Iterative methods for sparse linear systems* (2nd ed.).
///   SIAM. Section 7.3: BiCG and Variants, including BiCGSTAB.
/// - Sonneveld, P. (1989). CGS, a fast Lanczos-type solver for nonsymmetric linear systems.
///   *SIAM Journal on Scientific and Statistical Computing*, 10(1), 36-52.
///   Original CGS method, precursor to BiCGSTAB.
/// - Fletcher, R. (1976). Conjugate gradient methods for indefinite systems.
///   In *Numerical Analysis* (pp. 73-89). Springer.
///   BiCG foundation and convergence theory.
pub struct BiCGSTAB<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy> BiCGSTAB<T> {
    /// Perform matrix-vector multiplication, choosing serial or parallel based on config
    fn spmv(&self, a: &CsrMatrix<T>, x: &DVector<T>, y: &mut DVector<T>)
    where
        T: Send + Sync,
    {
        if self.config.use_parallel_spmv {
            spmv_parallel(a, x, y);
        } else {
            spmv(a, x, y);
        }
    }
    /// Create new BiCGSTAB solver
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

    /// Solve with preconditioning using efficient memory management
    ///
    /// # Arguments
    /// * `a` - The coefficient matrix (must be square)
    /// * `b` - The right-hand side vector
    /// * `preconditioner` - The preconditioner to use
    /// * `x` - On entry: initial guess; On exit: solution vector
    ///
    /// # Returns
    /// * `Ok(())` if converged successfully
    /// * `Err(...)` if failed to converge or numerical breakdown occurred
    pub fn solve_preconditioned<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        preconditioner: &P,
        x: &mut DVector<T>, // Changed: mutable reference instead of Option
    ) -> Result<()> {
        // Changed: returns Result<()> instead of Result<DVector<T>>
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        if x.len() != n {
            return Err(Error::InvalidConfiguration(
                "Solution vector dimension doesn't match system size".to_string(),
            ));
        }

        // Pre-allocate ALL workspace vectors outside the loop
        let mut r = DVector::zeros(n);
        let mut r0_hat = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n); // Pre-allocated for matrix-vector product
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n); // Pre-allocated for matrix-vector product
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);
        let mut ax = DVector::zeros(n); // Pre-allocated for initial residual

        // Compute initial residual: r = b - A*x
        self.spmv(a, x, &mut ax); // Allocation-free multiplication
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();

        // Check if already converged
        if self.is_converged(initial_residual_norm) {
            tracing::debug!("BiCGSTAB converged at initial guess");
            return Ok(());
        }

        // Use a robust breakdown tolerance based on machine epsilon
        let epsilon = T::default_epsilon();
        let breakdown_tolerance = epsilon * epsilon; // Simple, robust choice

        r0_hat.copy_from(&r); // Shadow residual

        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for iter in 0..self.config.max_iterations {
            let rho_new = r0_hat.dot(&r);

            // Primary breakdown condition: rho approaching zero
            if rho_new.abs() < breakdown_tolerance {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            let beta = (rho_new / rho) * (alpha / omega);

            // p = r + beta * (p - omega * v) - using in-place operations
            p.axpy(-omega, &v, T::one());
            p *= beta;
            p += &r;

            // Solve M*z = p and compute v = A*z
            preconditioner.apply_to(&p, &mut z)?;
            self.spmv(a, &z, &mut v); // Allocation-free multiplication

            alpha = rho_new / r0_hat.dot(&v);

            // s = r - alpha * v
            s.copy_from(&r);
            s.axpy(-alpha, &v, T::one());

            // REMOVED: Redundant convergence check on s
            // The standard BiCGSTAB algorithm only checks at the end of iteration

            // Solve M*z2 = s and compute t = A*z2
            preconditioner.apply_to(&s, &mut z2)?;
            self.spmv(a, &z2, &mut t); // Allocation-free multiplication

            let t_dot_t = t.dot(&t);
            if t_dot_t.abs() < breakdown_tolerance {
                // If t is nearly zero, we can't compute omega
                // But we can still update x with what we have
                x.axpy(alpha, &z, T::one());
                // Check if this partial update achieved convergence
                self.spmv(a, x, &mut ax);
                r.copy_from(b);
                r -= &ax;
                if self.is_converged(r.norm()) {
                    tracing::debug!(
                        "BiCGSTAB converged in {} iterations (t breakdown)",
                        iter + 1
                    );
                    return Ok(());
                }
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            omega = s.dot(&t) / t_dot_t;

            // Update solution: x = x + alpha*z + omega*z2
            x.axpy(alpha, &z, T::one());
            x.axpy(omega, &z2, T::one());

            // Update residual efficiently using swap
            // Instead of r = s - omega * t, we swap and update
            std::mem::swap(&mut r, &mut s); // O(1) swap instead of O(n) copy
            r.axpy(-omega, &t, T::one()); // Now r contains the new residual

            // Single convergence check per iteration (standard BiCGSTAB)
            if self.is_converged(r.norm()) {
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(());
            }

            // REMOVED: Check on omega - not a standard breakdown condition
            // The primary breakdown is caught by the rho check at loop start

            rho = rho_new;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Solve without preconditioning
    pub fn solve_unpreconditioned(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<()> {
        let identity = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &identity, x)
    }

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy> Configurable<T> for BiCGSTAB<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy> IterativeLinearSolver<T> for BiCGSTAB<T> {
    fn solve<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<()> {
        if let Some(p) = preconditioner {
            self.solve_preconditioned(a, b, p, x)
        } else {
            self.solve_unpreconditioned(a, b, x)
        }
    }
}

// Implement object-safe LinearSolver trait for trait objects
impl<T: RealField + Copy + num_traits::FromPrimitive> super::traits::LinearSolver<T>
    for BiCGSTAB<T>
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

    fn create_nonsymmetric_matrix() -> CsrMatrix<f64> {
        // Create a 3x3 nonsymmetric matrix (for testing BiCGSTAB)
        // [5, 1, 0]
        // [2, 4, 1]
        // [0, 2, 3]
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![5.0, 1.0, 2.0, 4.0, 1.0, 2.0, 3.0];
        
        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    fn create_diagonal_matrix() -> CsrMatrix<f64> {
        // Diagonal matrix [2, 0, 0; 0, 3, 0; 0, 0, 4]
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];
        
        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let config = IterativeSolverConfig::default();
        let _solver = BiCGSTAB::<f64>::new(config);
    }

    #[test]
    fn test_default_solver() {
        let _solver = BiCGSTAB::<f64>::default();
    }

    #[test]
    fn test_solve_nonsymmetric_system() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
        
        // Verify solution by checking A*x ≈ b
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![6.0, 9.0, 12.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
        
        // Solution should be [3.0, 3.0, 3.0]
        assert_relative_eq!(x[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::from_vec(vec![0.5, 0.5, 0.5]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
        
        // Verify solution
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_mismatched_dimensions_matrix() {
        let a = create_nonsymmetric_matrix(); // 3x3
        let b = DVector::from_vec(vec![1.0, 2.0]); // Wrong size!
        let mut x = DVector::zeros(2);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_mismatched_dimensions_solution() {
        let a = create_nonsymmetric_matrix(); // 3x3
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(2); // Wrong size!
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-12; // Very tight tolerance
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let mut config = IterativeSolverConfig::default();
        config.max_iterations = 1; // Too few iterations
        config.tolerance = 1e-12; // Tight tolerance
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        // Should fail to converge
        assert!(result.is_err());
    }

    #[test]
    fn test_already_converged() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![2.0, 3.0, 4.0]);
        let mut x = DVector::from_vec(vec![1.0, 1.0, 1.0]); // Exact solution
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
        
        // Solution should remain unchanged
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_configurable_trait() {
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-8;
        config.max_iterations = 500;
        
        let solver = BiCGSTAB::new(config);
        
        // Test getting config
        let retrieved_config = solver.config();
        assert_relative_eq!(retrieved_config.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved_config.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        
        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());
        
        let x = result.unwrap();
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_larger_nonsymmetric_system() {
        // 5x5 nonsymmetric tridiagonal matrix
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![
            0, 1,          // row 0
            0, 1, 2,       // row 1
            1, 2, 3,       // row 2
            2, 3, 4,       // row 3
            3, 4           // row 4
        ];
        let values = vec![
            5.0, 1.0,      // row 0
            2.0, 4.0, 1.0, // row 1
            2.0, 4.0, 1.0, // row 2
            2.0, 4.0, 1.0, // row 3
            2.0, 3.0       // row 4
        ];
        let a = CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");
        
        let b = DVector::from_vec(vec![6.0, 11.0, 11.0, 11.0, 8.0]);
        let mut x = DVector::zeros(5);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;
        
        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
        
        let mut ax = DVector::zeros(5);
        spmv(&a, &x, &mut ax);
        for i in 0..5 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }
}
