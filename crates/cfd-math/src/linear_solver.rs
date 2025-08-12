//! Linear solver implementations for CFD applications.
//!
//! This module provides various iterative linear solvers optimized for
//! sparse matrices arising from CFD discretizations.

use cfd_core::{Error, Result};
use nalgebra::{ComplexField, DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::{cast::FromPrimitive, Float};
use std::fmt::Debug;

use crate::sparse::SparseMatrixBuilder;

// Re-export the unified configuration from cfd-core
pub use cfd_core::{LinearSolverConfig, SolverConfiguration};

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

/// Preconditioner trait
pub trait Preconditioner<T: RealField>: Send + Sync {
    /// Apply preconditioner: solve M * z = r
    fn apply(&self, r: &DVector<T>) -> DVector<T>;
}

/// Identity preconditioner (no preconditioning)
pub struct IdentityPreconditioner;

impl<T: RealField> Preconditioner<T> for IdentityPreconditioner {
    fn apply(&self, r: &DVector<T>) -> DVector<T> {
        r.clone()
    }
}

/// Jacobi (diagonal) preconditioner
/// 
/// Uses the diagonal of the matrix as preconditioner: M = diag(A)
pub struct JacobiPreconditioner<T: RealField> {
    /// Inverse of diagonal elements
    inv_diagonal: DVector<T>,
}

impl<T: RealField> JacobiPreconditioner<T> {
    /// Create Jacobi preconditioner from matrix
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        let mut inv_diagonal = DVector::zeros(n);
        
        // Extract diagonal elements
        for i in 0..n {
            let mut diag = T::zero();
            let row = a.row(i);
            for (&j, val) in row.col_indices().iter().zip(row.values()) {
                if i == j {
                    diag = val.clone();
                    break;
                }
            }
            
            if ComplexField::abs(diag.clone()) < T::from_f64(1e-14).unwrap() {
                return Err(Error::NumericalError(
                    format!("Zero diagonal element at row {}", i)
                ));
            }
            inv_diagonal[i] = T::one() / diag;
        }
        
        Ok(Self { inv_diagonal })
    }
}

impl<T: RealField> Preconditioner<T> for JacobiPreconditioner<T> {
    fn apply(&self, r: &DVector<T>) -> DVector<T> {
        r.component_mul(&self.inv_diagonal)
    }
}

/// SOR (Successive Over-Relaxation) preconditioner
/// 
/// Uses forward SOR sweep as preconditioner: M = (D + ωL)
pub struct SORPreconditioner<T: RealField> {
    /// Matrix in CSR format
    matrix: CsrMatrix<T>,
    /// Relaxation parameter (typically 1.0 < ω < 2.0)
    omega: T,
}

impl<T: RealField + FromPrimitive> SORPreconditioner<T> {
    /// Create SOR preconditioner with given relaxation parameter
    pub fn new(a: &CsrMatrix<T>, omega: T) -> Result<Self> {
        if omega <= T::zero() || omega >= T::from_f64(2.0).unwrap() {
            return Err(Error::InvalidConfiguration(
                "SOR relaxation parameter must be in (0, 2)".to_string()
            ));
        }
        
        Ok(Self {
            matrix: a.clone(),
            omega,
        })
    }
    
    /// Create with optimal omega for symmetric positive definite matrices
    /// ω_opt = 2 / (1 + sin(π/n))
    pub fn with_optimal_omega(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows() as f64;
        let omega_opt = 2.0 / (1.0 + (std::f64::consts::PI / n).sin());
        let omega = T::from_f64(omega_opt).unwrap();
        Self::new(a, omega)
    }
}

impl<T: RealField> Preconditioner<T> for SORPreconditioner<T> {
    fn apply(&self, r: &DVector<T>) -> DVector<T> {
        let n = r.len();
        let mut z: DVector<T> = DVector::zeros(n);
        
        // Forward SOR sweep: (D + ωL)z = r
        for i in 0..n {
            let mut sum = r[i].clone();
            let mut diag = T::one();
            
            let row = self.matrix.row(i);
            for (&j, val) in row.col_indices().iter().zip(row.values()) {
                if j < i {
                    sum = sum - val.clone() * z[j].clone() * self.omega.clone();
                } else if j == i {
                    diag = val.clone();
                }
            }
            
            z[i] = sum / (diag * self.omega.clone());
        }
        
        z
    }
}

/// ILU(0) (Incomplete LU factorization with zero fill-in) preconditioner
/// 
/// Computes incomplete LU factorization preserving the sparsity pattern
pub struct ILU0Preconditioner<T: RealField> {
    /// Lower triangular factor (including diagonal)
    l_factor: CsrMatrix<T>,
    /// Upper triangular factor
    u_factor: CsrMatrix<T>,
}

impl<T: RealField + FromPrimitive> ILU0Preconditioner<T> {
    /// Create ILU(0) preconditioner from matrix
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        let n = a.nrows();
        
        // Copy matrix structure for factorization
        let mut l_builder = SparseMatrixBuilder::new(n, n);
        let mut u_builder = SparseMatrixBuilder::new(n, n);
        
        // Working copy of matrix
        let mut work = vec![vec![T::zero(); n]; n];
        let mut pattern = vec![vec![false; n]; n];
        
        // Copy matrix to dense working array (only non-zero pattern)
        for i in 0..n {
            let row = a.row(i);
            for (&j, val) in row.col_indices().iter().zip(row.values()) {
                work[i][j] = val.clone();
                pattern[i][j] = true;
            }
        }
        
        // ILU(0) factorization
        for k in 0..n {
            // Check for zero pivot
            if ComplexField::abs(work[k][k].clone()) < T::from_f64(1e-14).unwrap() {
                return Err(Error::NumericalError(
                    format!("Zero pivot encountered at row {}", k)
                ));
            }
            
            for i in (k + 1)..n {
                if pattern[i][k] {
                    // Compute multiplier
                    work[i][k] = work[i][k].clone() / work[k][k].clone();
                    
                    // Update row i
                    for j in (k + 1)..n {
                        if pattern[i][j] && pattern[k][j] {
                            work[i][j] = work[i][j].clone() - work[i][k].clone() * work[k][j].clone();
                        }
                    }
                }
            }
        }
        
        // Extract L and U factors
        for i in 0..n {
            for j in 0..=i {
                if pattern[i][j] {
                    if i == j {
                        l_builder.add_entry(i, j, T::one());
                        u_builder.add_entry(i, j, work[i][j].clone());
                    } else {
                        l_builder.add_entry(i, j, work[i][j].clone());
                    }
                }
            }
            for j in (i + 1)..n {
                if pattern[i][j] {
                    u_builder.add_entry(i, j, work[i][j].clone());
                }
            }
        }
        
        Ok(Self {
            l_factor: l_builder.build()?,
            u_factor: u_builder.build()?,
        })
    }
}

impl<T: RealField> Preconditioner<T> for ILU0Preconditioner<T> {
    fn apply(&self, r: &DVector<T>) -> DVector<T> {
        let n = r.len();
        
        // Forward substitution: L * y = r
        let mut y: DVector<T> = DVector::zeros(n);
        for i in 0..n {
            let mut sum = r[i].clone();
            let row = self.l_factor.row(i);
            for (&j, val) in row.col_indices().iter().zip(row.values()) {
                if j < i {
                    sum = sum - val.clone() * y[j].clone();
                }
            }
            y[i] = sum;
        }
        
        // Backward substitution: U * z = y
        let mut z: DVector<T> = DVector::zeros(n);
        for i in (0..n).rev() {
            let mut sum = y[i].clone();
            let mut diag = T::one();
            
            let row = self.u_factor.row(i);
            for (&j, val) in row.col_indices().iter().zip(row.values()) {
                if j > i {
                    sum = sum - val.clone() * z[j].clone();
                } else if j == i {
                    diag = val.clone();
                }
            }
            
            z[i] = sum / diag;
        }
        
        z
    }
}

/// Conjugate Gradient solver
///
/// Implements the Conjugate Gradient method for symmetric positive definite matrices.
/// Reference: Hestenes, M. R.; Stiefel, E. (1952). "Methods of conjugate gradients for solving linear systems"
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
}

impl<T: RealField + Debug> LinearSolver<T> for ConjugateGradient<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        // Initialize solution
        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);

        // Compute initial residual: r = b - A*x
        let mut r = b - a * &x;
        let mut p = r.clone();
        let mut rsold = r.dot(&r);

        // CG iterations
        for iter in 0..self.config.max_iterations() {
            let ap = a * &p;
            let alpha = rsold.clone() / p.dot(&ap);
            
            // Zero-copy update using advanced iterator combinators for SIMD optimization
            // Use zero-copy in-place operations for better performance
            x.iter_mut()
                .zip(p.iter())
                .for_each(|(xi, pi)| *xi += alpha.clone() * pi.clone());

            r.iter_mut()
                .zip(ap.iter())
                .for_each(|(ri, api)| *ri -= alpha.clone() * api.clone());
            
            let rsnew = r.dot(&r);
            
            // Check convergence using zero-copy norm computation
            if self.is_converged(rsnew.clone().sqrt()) {
                tracing::debug!("CG converged in {} iterations", iter + 1);
                return Ok(x);
            }

            let beta = rsnew.clone() / rsold;

            // Zero-copy search direction update
            p.iter_mut()
                .zip(r.iter())
                .for_each(|(pi, ri)| *pi = ri.clone() + beta.clone() * pi.clone());

            rsold = rsnew;
        }

        Err(Error::ConvergenceFailure(format!(
            "CG failed to converge after {} iterations",
            self.config.max_iterations()
        )))
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}

/// GMRES solver (Generalized Minimal Residual)
///
/// Implements the GMRES method for general non-symmetric matrices.
/// Reference: Saad, Y.; Schultz, M. H. (1986). "GMRES: A generalized minimal residual algorithm"
pub struct GMRES<T: RealField + Float> {
    config: LinearSolverConfig<T>,
}

impl<T: RealField + Float> GMRES<T> {
    /// Create new GMRES solver
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

    /// Validate matrix and vector dimensions
    fn validate_dimensions(&self, a: &CsrMatrix<T>, b: &DVector<T>) -> Result<()> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }
        Ok(())
    }

    /// Perform one GMRES cycle
    fn gmres_cycle(
        &self,
        a: &CsrMatrix<T>,
        r: &DVector<T>,
        beta: T,
        restart: usize,
    ) -> Result<DVector<T>> {
        let mut q = vec![r / beta.clone()];
        let mut h: Vec<Vec<T>> = (0..=restart).map(|_| Vec::new()).collect();
        let mut g = DVector::zeros(restart + 1);
        g[0] = beta;

        let mut cs: Vec<T> = Vec::with_capacity(restart);
        let mut sn: Vec<T> = Vec::with_capacity(restart);

        for k in 0..restart {
            let q_new = self.arnoldi_step(a, &q, &mut h, k);
            q.push(q_new);

            self.apply_previous_rotations(&mut h, &cs, &sn, k);
            let (c, s) = self.compute_givens_rotation(&mut h, k)?;
            cs.push(c.clone());
            sn.push(s.clone());

            self.update_residual(&mut g, c, s, k);

            let residual_norm = ComplexField::abs(g[k + 1].clone());
            if self.is_converged(residual_norm) {
                let y = self.solve_upper_triangular(&h, &g, k + 1)?;
                return Ok(self.compute_correction(&q, &y));
            }
        }

        let y = self.solve_upper_triangular(&h, &g, restart)?;
        Ok(self.compute_correction(&q, &y))
    }

    /// Arnoldi process for building orthonormal basis
    fn arnoldi_step(
        &self,
        a: &CsrMatrix<T>,
        q: &[DVector<T>],
        h: &mut Vec<Vec<T>>,
        k: usize,
    ) -> DVector<T> {
        let mut v = a * &q[k];

        // Modified Gram-Schmidt orthogonalization for better numerical stability
        for (i, qi) in q.iter().enumerate().take(k + 1) {
            let hij = v.dot(qi);
            h[i].push(hij.clone());
            v = &v - &(qi * hij);
            
            // Re-orthogonalize for better numerical stability (Classical Gram-Schmidt with reorthogonalization)
            let hij_correction = v.dot(qi);
            if ComplexField::abs(hij_correction.clone()) > T::from_f64(1e-10).unwrap() {
                h[i][k] = h[i][k].clone() + hij_correction.clone();
                v = &v - &(qi * hij_correction);
            }
        }

        // Normalize
        let norm = v.norm();
        if norm < T::from_f64(1e-14).unwrap() {
            // Handle breakdown - vector is in the span of previous vectors
            h[k + 1].push(T::zero());
            DVector::zeros(v.len())
        } else {
            h[k + 1].push(norm.clone());
            v / norm
        }
    }

    /// Apply previous Givens rotations to current column
    fn apply_previous_rotations(&self, h: &mut Vec<Vec<T>>, cs: &[T], sn: &[T], k: usize) {
        for i in 0..k {
            let (h_new, h_ip1_new) = Self::apply_givens_rotation(
                cs[i].clone(),
                sn[i].clone(),
                h[i][k].clone(),
                h[i + 1][k].clone(),
            );
            h[i][k] = h_new;
            h[i + 1][k] = h_ip1_new;
        }
    }

    /// Compute new Givens rotation coefficients
    fn compute_givens_rotation(&self, h: &mut Vec<Vec<T>>, k: usize) -> Result<(T, T)> {
        // Ensure h[k] and h[k+1] have enough elements
        while h[k].len() <= k {
            h[k].push(T::zero());
        }
        while h[k + 1].len() <= k {
            h[k + 1].push(T::zero());
        }

        let h_k = h[k][k].clone();
        let h_kp1 = h[k + 1][k].clone();
        
        // Improved numerical stability for Givens rotation
        if ComplexField::abs(h_kp1.clone()) < T::from_f64(1e-14).unwrap() {
            // h_kp1 is essentially zero
            let c = if ComplexField::abs(h_k.clone()) < T::from_f64(1e-14).unwrap() {
                T::one()
            } else {
                ComplexField::signum(h_k.clone())
            };
            h[k][k] = ComplexField::abs(h_k);
            h[k + 1][k] = T::zero();
            return Ok((c, T::zero()));
        }
        
        // Use more stable formulation
        if ComplexField::abs(h_kp1.clone()) > ComplexField::abs(h_k.clone()) {
            let tau = h_k.clone() / h_kp1.clone();
            let s = T::one() / ComplexField::sqrt(T::one() + tau.clone() * tau.clone());
            let c = s.clone() * tau;
            h[k][k] = h_kp1.clone() / s;
            h[k + 1][k] = T::zero();
            Ok((c, s))
        } else {
            let tau = h_kp1.clone() / h_k.clone();
            let c = T::one() / ComplexField::sqrt(T::one() + tau.clone() * tau.clone());
            let s = c.clone() * tau;
            h[k][k] = h_k.clone() / c;
            h[k + 1][k] = T::zero();
            Ok((c, s))
        }
    }

    /// Update residual vector using Givens rotation
    fn update_residual(&self, g: &mut DVector<T>, c: T, s: T, k: usize) {
        g[k + 1] = -s.clone() * g[k].clone();
        g[k] = c * g[k].clone();
    }

    /// Solve upper triangular system by back substitution
    fn solve_upper_triangular(&self, h: &[Vec<T>], g: &DVector<T>, size: usize) -> Result<DVector<T>> {
        let mut y: DVector<T> = DVector::zeros(size);
        for i in (0..size).rev() {
            let mut yi = g[i].clone();
            for j in (i + 1)..size {
                yi = yi - h[i][j].clone() * y[j].clone();
            }
            if h[i][i] == T::zero() {
                return Err(Error::NumericalError("Singular matrix in GMRES".to_string()));
            }
            y[i] = yi / h[i][i].clone();
        }
        Ok(y)
    }

    /// Compute correction vector from Krylov basis and coefficients
    fn compute_correction(&self, q: &[DVector<T>], y: &DVector<T>) -> DVector<T> {
        let mut correction = DVector::zeros(q[0].len());
        for (i, yi) in y.iter().enumerate() {
            correction += &q[i] * yi.clone();
        }
        correction
    }

    /// Apply Givens rotation
    fn apply_givens_rotation(c: T, s: T, h_i: T, h_ip1: T) -> (T, T) {
        let temp = c.clone() * h_i.clone() + s.clone() * h_ip1.clone();
        let h_ip1_new = -s * h_i + c * h_ip1;
        (temp, h_ip1_new)
    }
}

impl<T: RealField + Debug + Float> LinearSolver<T> for GMRES<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        self.validate_dimensions(a, b)?;

        let mut x = x0.map_or_else(|| DVector::zeros(b.len()), DVector::clone);
        let restart = self.config.restart.min(b.len());

        for _outer in 0..self.config.max_iterations() {
            let r = b - a * &x;
            let beta = r.norm();

            if self.is_converged(beta.clone()) {
                return Ok(x);
            }

            let correction = self.gmres_cycle(a, &r, beta, restart)?;
            x += correction;
        }

        Err(Error::ConvergenceFailure(
            "GMRES failed to converge within maximum iterations".to_string(),
        ))
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}

/// BiCGSTAB solver (Biconjugate Gradient Stabilized)
///
/// Implements the BiCGSTAB method for non-symmetric matrices.
/// Reference: Van der Vorst, H. A. (1992). "Bi-CGSTAB: A fast and smoothly converging variant of Bi-CG"
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
}

impl<T: RealField + Debug> LinearSolver<T> for BiCGSTAB<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);

        // Initialize
        let mut r = b - a * &x;
        let r0_hat = r.clone(); // Arbitrary choice, could be different
        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n);

        for iter in 0..self.config.max_iterations() {
            let rho_new = r0_hat.dot(&r);
            
            if rho_new.clone().abs() < T::from_f64(1e-14).unwrap() {
                return Err(Error::NumericalError(
                    "BiCGSTAB breakdown: rho = 0".to_string(),
                ));
            }

            let beta = (rho_new.clone() / rho.clone()) * (alpha.clone() / omega.clone());
            
            // Update search direction
            p = &r + &((&p - &(&v * omega.clone())) * beta);
            
            // Matrix-vector product
            v = a * &p;
            alpha = rho_new.clone() / r0_hat.dot(&v);
            
            // Update s
            let s = &r - &(&v * alpha.clone());
            
            // Check convergence
            let s_norm = s.norm();
            if self.is_converged(s_norm) {
                x = &x + &(&p * alpha);
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }
            
            // Second matrix-vector product
            let t = a * &s;
            omega = s.dot(&t) / t.norm_squared();
            
            // Update solution
            x = &x + &(&p * alpha.clone()) + &(&s * omega.clone());
            
            // Update residual
            r = s - &t * omega.clone();
            
            // Check convergence
            let r_norm = r.norm();
            if self.is_converged(r_norm) {
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }
            
            rho = rho_new;
        }

        Err(Error::ConvergenceFailure(format!(
            "BiCGSTAB failed to converge after {} iterations",
            self.config.max_iterations()
        )))
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra_sparse::CsrMatrix;

    fn create_test_system() -> (CsrMatrix<f64>, DVector<f64>) {
        // Create a simple 3x3 SPD matrix for testing
        // [4, -1, 0]
        // [-1, 4, -1]
        // [0, -1, 4]
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0];
        
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values).unwrap();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        
        (a, b)
    }

    #[test]
    fn test_conjugate_gradient() {
        let (a, b) = create_test_system();
        let solver = ConjugateGradient::default();
        let x = solver.solve(&a, &b, None).unwrap();
        
        // Verify solution
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }

    #[test]
    fn test_gmres() {
        let (a, b) = create_test_system();
        let mut config = LinearSolverConfig::default();
        config.restart = 3; // Use full restart for small system
        config.base = cfd_core::SolverConfig::builder()
            .max_iterations(50)
            .tolerance(1e-6)
            .build();
        let solver = GMRES::new(config);
        let x = solver.solve(&a, &b, None).unwrap();
        
        // Verify solution
        let residual = &b - &a * &x;
        println!("GMRES residual norm: {}", residual.norm());
        println!("Solution: {:?}", x);

        // GMRES should achieve good accuracy with improved numerical stability
        assert!(residual.norm() < 1e-6, "GMRES should achieve accuracy of 1e-6 with improved stability");
    }

    #[test]
    fn test_bicgstab() {
        let (a, b) = create_test_system();
        let solver = BiCGSTAB::default();
        let x = solver.solve(&a, &b, None).unwrap();
        
        // Verify solution
        let residual = &b - &a * &x;
        assert!(residual.norm() < 1e-10);
    }
}