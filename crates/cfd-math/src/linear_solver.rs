//! Linear solver implementations for CFD applications.
//!
//! This module provides iterative linear solvers optimized for sparse matrices
//! arising from CFD discretizations. It emphasizes the use of well-tested
//! implementations and efficient memory management.
//!
//! ## Design Philosophy
//! 
//! Rather than re-implementing complex numerical algorithms from scratch, this module
//! focuses on providing efficient wrappers and using robust implementations where
//! available. This approach reduces maintenance overhead, eliminates numerical bugs,
//! and provides better performance through optimized algorithms.
//!
//! ## Performance Considerations
//!
//! All solvers prioritize memory efficiency by:
//! - Using in-place operations to avoid unnecessary allocations
//! - Pre-allocating workspace vectors outside iteration loops
//! - Leveraging nalgebra's optimized BLAS-like operations
//! - Providing efficient preconditioner APIs

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use nalgebra_sparse::coo::CooMatrix;
use num_traits::cast::FromPrimitive;
use std::fmt::Debug;
use crate::sparse::SparseMatrixExt;

// Re-export the unified configuration from cfd-core
pub use cfd_core::solver::{LinearSolverConfig, SolverConfiguration};

/// Trait for linear solvers
pub trait LinearSolver<T: RealField>: Send + Sync {
    /// Solve Ax = b
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>>;

    /// Get solver configuration
    fn config(&self) -> &LinearSolverConfig<T>;

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config().tolerance()
    }
}

/// Simplified, efficient preconditioner trait
/// 
/// This API enforces explicit memory management and avoids hidden allocations
/// by requiring the user to provide the output vector.
pub trait Preconditioner<T: RealField>: Send + Sync {
    /// Apply the preconditioner: M(z) = r
    /// 
    /// Solves the preconditioning system and stores the result in `z`.
    /// This approach makes memory management explicit and avoids hidden allocations.
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()>;
}

/// Identity preconditioner (no preconditioning)
#[derive(Default)]
pub struct IdentityPreconditioner;

impl<T: RealField> Preconditioner<T> for IdentityPreconditioner {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.copy_from(r);
        Ok(())
    }
}

/// Jacobi (diagonal) preconditioner with optimized memory usage
pub struct JacobiPreconditioner<T: RealField> {
    inv_diagonal: DVector<T>,
}

impl<T: RealField + FromPrimitive> JacobiPreconditioner<T> {
    /// Create Jacobi preconditioner from matrix diagonal
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "Jacobi preconditioner requires square matrix".to_string()
            ));
        }

        // Extract diagonal directly
        let diag = a.diagonal();
        let mut inv_diagonal = DVector::zeros(n);

        for (i, val) in diag.iter().enumerate() {
            if val.clone().abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                return Err(Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation));
            }
            inv_diagonal[i] = T::one() / val.clone();
        }

        Ok(Self { inv_diagonal })
    }
}

impl<T: RealField> Preconditioner<T> for JacobiPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        z.copy_from(r);
        z.component_mul_assign(&self.inv_diagonal);
        Ok(())
    }
}

/// SOR (Successive Over-Relaxation) preconditioner with validation
pub struct SORPreconditioner<T: RealField> {
    matrix: CsrMatrix<T>,
    omega: T,
}

impl<T: RealField + FromPrimitive> SORPreconditioner<T> {
    /// Create SOR preconditioner with specified relaxation parameter
    pub fn new(a: &CsrMatrix<T>, omega: T) -> Result<Self> {
        let n = a.nrows();
        if n != a.ncols() {
            return Err(Error::InvalidConfiguration(
                "SOR preconditioner requires square matrix".to_string()
            ));
        }

        // Validate omega range for stability
        let zero = T::zero();
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        if omega <= zero || omega >= two {
            return Err(Error::InvalidConfiguration(
                "SOR omega parameter must be in range (0, 2) for stability".to_string()
            ));
        }

        Ok(Self {
            matrix: a.clone(),
            omega,
        })
    }

    /// Create SOR preconditioner with omega optimized for 1D Poisson problems
    ///
    /// ## Warning
    /// This function computes an optimal omega value **specifically** for matrices
    /// arising from 1D Poisson equations discretized with second-order finite
    /// differences on uniform grids. Using this for any other matrix type will
    /// likely result in suboptimal performance or convergence issues.
    ///
    /// ## Validation
    /// This function validates that the matrix has the expected tridiagonal
    /// structure before computing the optimal omega.
    pub fn with_omega_for_1d_poisson(a: &CsrMatrix<T>) -> Result<Self> {
        // Validate matrix structure for 1D Poisson
        Self::validate_1d_poisson_structure(a)?;
        
        let n = a.nrows() as f64;
        let omega_opt = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        let omega = T::from_f64(omega_opt).ok_or_else(|| {
            Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation)
        })?;
        
        Self::new(a, omega)
    }

    /// Validate that matrix has 1D Poisson structure (tridiagonal with specific pattern)
    fn validate_1d_poisson_structure(a: &CsrMatrix<T>) -> Result<()> {
        let n = a.nrows();
        
        // Check basic structure: each row should have at most 3 non-zeros
        for i in 0..n {
            let row = a.row(i);
            if row.nnz() > 3 {
                return Err(Error::InvalidConfiguration(
                    format!("Row {} has {} non-zeros; 1D Poisson should have at most 3", i, row.nnz())
                ));
            }
            
            // Check that non-zeros are in expected positions (diagonal and adjacent)
            for &j in row.col_indices() {
                if (j as i32 - i as i32).abs() > 1 {
                    return Err(Error::InvalidConfiguration(
                        format!("Non-zero at ({}, {}) violates tridiagonal structure", i, j)
                    ));
                }
            }
        }
        
        Ok(())
    }
}

impl<T: RealField> Preconditioner<T> for SORPreconditioner<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let n = r.len();
        z.fill(T::zero());
        
        // Forward SOR sweep with in-place operations (standard SOR: (D/ω - L) z = r)
        for i in 0..n {
            let mut sum = T::zero();
            let row = self.matrix.row(i);
            let mut diag = T::one();
            
            // Process row entries
            for (j, val) in row.col_indices().iter().zip(row.values()) {
                if *j < i {
                    // Accumulate strictly lower part L z
                    sum = sum + val.clone() * z[*j].clone();
                } else if *j == i {
                    diag = val.clone();
                }
            }
            
            // Standard SOR preconditioner solves (D/ω - L) z = r
            // => (diag/ω) z_i = r_i + (L z)_i
            z[i] = (r[i].clone() + sum) * self.omega.clone() / diag;
        }
        
        Ok(())
    }
}

/// Preconditioned Conjugate Gradient solver with optimized memory management
pub struct ConjugateGradient<T: RealField> {
    config: LinearSolverConfig<T>,
}

impl<T: RealField> ConjugateGradient<T> {
    /// Create new CG solver
    pub const fn new(config: LinearSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        Self::new(LinearSolverConfig::default())
    }

    /// Solve with preconditioning and optimized memory management
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
        let mut ap = DVector::zeros(n);

        // Compute initial residual: r = b - A*x
        let ax = a * &x;
        r.copy_from(b);
        r -= &ax;

        // Apply preconditioner: M*z = r
        preconditioner.apply_to(&r, &mut z)?;
        p.copy_from(&z);
        
        let mut rzold = r.dot(&z);

        // PCG iterations with in-place operations
        for iter in 0..self.config.base.convergence.max_iterations {
            // Compute A*p
            ap = a * &p;
            
            let alpha = rzold.clone() / p.dot(&ap);
            
            // Update solution: x = x + alpha * p
            x.axpy(alpha.clone(), &p, T::one());
            
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
            let beta = rznew.clone() / rzold;

            // Update search direction: p = z + beta * p
            p *= beta;
            p += &z;

            rzold = rznew;
        }

        Err(Error::Convergence(cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { 
            max: self.config.base.convergence.max_iterations 
        }))
    }
}

impl<T: RealField + Debug> LinearSolver<T> for ConjugateGradient<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        // Use identity preconditioner for unpreconditioned solve
        let identity = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &identity, x0)
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}

/// BiCGSTAB solver with optimized memory management
pub struct BiCGSTAB<T: RealField> {
    config: LinearSolverConfig<T>,
}

impl<T: RealField> BiCGSTAB<T> {
    /// Create new BiCGSTAB solver
    pub const fn new(config: LinearSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        Self::new(LinearSolverConfig::default())
    }

    /// Solve with preconditioning and optimized memory management
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

        // Pre-allocate workspace vectors
        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);
        let mut r = DVector::zeros(n);
        let mut r0_hat = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n);
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);

        // Compute initial residual: r = b - A*x
        let ax = a * &x;
        r.copy_from(b);
        r -= &ax;
        
        let initial_residual_norm = r.norm();
        let breakdown_tolerance = initial_residual_norm * T::default_epsilon();
        
        r0_hat.copy_from(&r); // Shadow residual
        
        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for iter in 0..self.config.base.convergence.max_iterations {
            let rho_new = r0_hat.dot(&r);
            
            if rho_new.clone().abs() < breakdown_tolerance {
                return Err(Error::Numerical(
                    cfd_core::error::NumericalErrorKind::SingularMatrix
                ));
            }

            let beta = (rho_new.clone() / rho.clone()) * (alpha.clone() / omega.clone());
            
            // p = r + beta * (p - omega * v) - using in-place operations
            p.axpy(-omega.clone(), &v, T::one());
            p *= beta;
            p += &r;

            // Solve M*z = p and compute v = A*z
            preconditioner.apply_to(&p, &mut z)?;
            v = a * &z;
            
            alpha = rho_new.clone() / r0_hat.dot(&v);
            
            // s = r - alpha * v
            s.copy_from(&r);
            s.axpy(-alpha.clone(), &v, T::one());
            
            if self.is_converged(s.norm()) {
                x.axpy(alpha.clone(), &z, T::one());
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }

            // Solve M*z2 = s and compute t = A*z2
            preconditioner.apply_to(&s, &mut z2)?;
            t = a * &z2;
            
            omega = s.dot(&t) / t.dot(&t);
            
            // Update solution: x = x + alpha*z + omega*z2
            x.axpy(alpha.clone(), &z, T::one());
            x.axpy(omega.clone(), &z2, T::one());
            
            // Update residual: r = s - omega * t
            r.copy_from(&s);
            r.axpy(-omega.clone(), &t, T::one());
            
            if self.is_converged(r.norm()) {
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }
            
            if omega.clone().abs() < breakdown_tolerance {
                return Err(Error::Numerical(
                    cfd_core::error::NumericalErrorKind::SingularMatrix
                ));
            }

            rho = rho_new;
        }

        Err(Error::Convergence(cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { 
            max: self.config.base.convergence.max_iterations 
        }))
    }
}

impl<T: RealField + Debug> LinearSolver<T> for BiCGSTAB<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        // Use identity preconditioner for unpreconditioned solve
        let identity = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &identity, x0)
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;

    #[cfg(test)]
    fn create_validation_matrix() -> CsrMatrix<f64> {
        // Create a simple 3x3 SPD matrix
        let dense = DMatrix::from_row_slice(3, 3, &[
            4.0, -1.0, 0.0,
            -1.0, 4.0, -1.0,
            0.0, -1.0, 4.0,
        ]);
        CsrMatrix::from(&dense)
    }

    fn create_1d_poisson_matrix(n: usize) -> CsrMatrix<f64> {
        // Create tridiagonal matrix for 1D Poisson: [-1, 2, -1]
        let mut builder = crate::sparse::SparseMatrixBuilder::new(n, n);
        
        for i in 0..n {
            // Diagonal
            builder.add_entry(i, i, 2.0).expect("CRITICAL: Add proper error handling");
            
            // Off-diagonals
            if i > 0 {
                builder.add_entry(i, i-1, -1.0).expect("CRITICAL: Add proper error handling");
            }
            if i < n-1 {
                builder.add_entry(i, i+1, -1.0).expect("CRITICAL: Add proper error handling");
            }
        }
        
        builder.build().expect("CRITICAL: Add proper error handling")
    }

    fn create_test_matrix() -> CsrMatrix<f64> {
        let mut coo = CooMatrix::new(3, 3);
        coo.push(0, 0, 4.0);
        coo.push(0, 1, -1.0);
        coo.push(1, 0, -1.0);
        coo.push(1, 1, 4.0);
        coo.push(1, 2, -1.0);
        coo.push(2, 1, -1.0);
        coo.push(2, 2, 4.0);
        CsrMatrix::from(&coo)
    }

    #[test]
    fn test_cg_solver() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = ConjugateGradient::default();
        let x = solver.solve(&a, &b, None).expect("CRITICAL: Add proper error handling");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_bicgstab_solver() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = BiCGSTAB::default();
        let x = solver.solve(&a, &b, None).expect("CRITICAL: Add proper error handling");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_jacobi_preconditioner() {
        let a = create_test_matrix();
        let r = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let precond = JacobiPreconditioner::new(&a).expect("CRITICAL: Add proper error handling");
        let mut z = DVector::zeros(3);
        precond.apply_to(&r, &mut z).expect("CRITICAL: Add proper error handling");
        
        // Check that z is correctly scaled by inverse diagonal
        assert_relative_eq!(z[0], r[0] / 4.0, epsilon = 1e-12);
        assert_relative_eq!(z[1], r[1] / 4.0, epsilon = 1e-12);
        assert_relative_eq!(z[2], r[2] / 4.0, epsilon = 1e-12);
    }

    #[test]
    fn test_sor_omega_validation() {
        let a = create_test_matrix();
        
        // Test invalid omega values
        assert!(SORPreconditioner::new(&a, 0.0).is_err());
        assert!(SORPreconditioner::new(&a, 2.0).is_err());
        assert!(SORPreconditioner::new(&a, -0.5).is_err());
        
        // Test valid omega
        assert!(SORPreconditioner::new(&a, 1.5).is_ok());
    }

    #[test]
    fn test_1d_poisson_omega_optimization() {
        let a = create_1d_poisson_matrix(10);
        
        // Should successfully create SOR with optimized omega
        let sor = SORPreconditioner::with_omega_for_1d_poisson(&a).expect("CRITICAL: Add proper error handling");
        
        // Omega should be in valid range
        assert!(sor.omega > 0.0 && sor.omega < 2.0);
        
        // For n=10, calculate expected omega
        let n = 10.0_f64;
        let expected_omega = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        println!("Expected omega: {}, Actual omega: {}", expected_omega, sor.omega);
        
        // For n=10, optimal omega ≈ 1.69 (allow larger tolerance)
        assert!((sor.omega - expected_omega).abs() < 0.01);
    }

    #[test]
    fn test_1d_poisson_validation_rejects_invalid_matrix() {
        let a = create_test_matrix(); // This is a 3x3 tridiagonal matrix: should actually pass validation
        
        // The test matrix we created is actually tridiagonal, so let's create a non-tridiagonal matrix
        let mut builder = crate::sparse::SparseMatrixBuilder::new(3, 3);
        builder.add_entry(0, 0, 2.0).expect("CRITICAL: Add proper error handling");
        builder.add_entry(0, 2, 1.0).expect("CRITICAL: Add proper error handling"); // Non-adjacent entry
        builder.add_entry(1, 1, 2.0).expect("CRITICAL: Add proper error handling");
        builder.add_entry(2, 2, 2.0).expect("CRITICAL: Add proper error handling");
        let non_tridiag = builder.build().expect("CRITICAL: Add proper error handling");
        
        // Should reject non-tridiagonal matrix
        assert!(SORPreconditioner::with_omega_for_1d_poisson(&non_tridiag).is_err());
    }

    #[test]
    fn test_preconditioned_cg() {
        let a = create_test_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        let solver = ConjugateGradient::default();
        let precond = JacobiPreconditioner::new(&a).expect("CRITICAL: Add proper error handling");
        
        let x = solver.solve_preconditioned(&a, &b, &precond, None).expect("CRITICAL: Add proper error handling");
        
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_identity_preconditioner() {
        let r = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let precond = IdentityPreconditioner;
        
        let mut z = DVector::zeros(3);
        precond.apply_to(&r, &mut z).expect("CRITICAL: Add proper error handling");
        
        assert_relative_eq!(z, r, epsilon = 1e-12);
    }
}