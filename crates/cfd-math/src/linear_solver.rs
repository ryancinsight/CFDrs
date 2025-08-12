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
    /// Apply preconditioner: solve M * z = r, storing result in z
    /// 
    /// # Arguments
    /// * `r` - Right-hand side vector
    /// * `z` - Output vector (pre-allocated)
    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>);
}

/// Identity preconditioner (no preconditioning)
pub struct IdentityPreconditioner;

impl<T: RealField> Preconditioner<T> for IdentityPreconditioner {
    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>) {
        z.copy_from(r);
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
    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>) {
        // z = inv_diagonal .* r (element-wise multiplication)
        z.copy_from(r);
        z.component_mul_assign(&self.inv_diagonal);
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
    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>) {
        let n = r.len();
        z.fill(T::zero());
        
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
    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>) {
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
        z.fill(T::zero());
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
            
            // Use in-place nalgebra operations for efficiency
            x.axpy(alpha.clone(), &p, T::one());  // x = x + alpha * p
            r.axpy(-alpha.clone(), &ap, T::one()); // r = r - alpha * ap
            
            let rsnew = r.dot(&r);
            
            // Check convergence
            if self.is_converged(rsnew.clone().sqrt()) {
                tracing::debug!("CG converged in {} iterations", iter + 1);
                return Ok(x);
            }

            let beta = rsnew.clone() / rsold;

            // Update search direction: p = r + beta * p
            // We need to compute p_new = r + beta * p_old
            p *= beta;  // First scale p by beta
            p += &r;    // Then add r

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

    /// GMRES iteration (one restart cycle)
    fn gmres_cycle(
        &self,
        a: &CsrMatrix<T>,
        r: &DVector<T>,
        beta: T,
        restart: usize,
    ) -> Result<DVector<T>> {
        let n = r.len();
        
        // Initialize Krylov subspace basis
        let mut q = Vec::with_capacity(restart + 1);
        q.push(r / beta.clone());
        
        // Hessenberg matrix stored as flat vector in column-major order
        // H is (restart+1) x restart
        let mut h = vec![T::zero(); (restart + 1) * restart];
        
        // Givens rotation coefficients
        let mut cs = Vec::with_capacity(restart);
        let mut sn = Vec::with_capacity(restart);
        
        // Right-hand side for least squares problem
        let mut g = DVector::zeros(restart + 1);
        g[0] = beta;
        
        for k in 0..restart {
            // Arnoldi step
            let qk_plus_1 = self.arnoldi_step(a, &q, &mut h, k, restart)?;
            q.push(qk_plus_1);
            
            // Apply previous Givens rotations to new column
            self.apply_previous_rotations(&mut h, &cs, &sn, k, restart);
            
            // Compute new Givens rotation
            let (c, s) = self.compute_givens_rotation(&mut h, k, restart)?;
            cs.push(c.clone());
            sn.push(s.clone());
            
            // Update residual
            self.update_residual(&mut g, c, s, k);
            
            // Check convergence
            let residual_norm = ComplexField::abs(g[k + 1].clone());
            if self.is_converged(residual_norm) {
                let y = self.solve_upper_triangular(&h, &g, k + 1, restart)?;
                return Ok(self.compute_correction(&q, &y));
            }
        }

        let y = self.solve_upper_triangular(&h, &g, restart, restart)?;
        Ok(self.compute_correction(&q, &y))
    }

    /// Arnoldi process for building orthonormal basis with in-place operations
    fn arnoldi_step(
        &self,
        a: &CsrMatrix<T>,
        q: &[DVector<T>],
        h: &mut Vec<T>,
        k: usize,
        restart: usize,
    ) -> Result<DVector<T>> {
        let mut v = a * &q[k];

        // Modified Gram-Schmidt orthogonalization with in-place operations
        for i in 0..=k {
            let hij = v.dot(&q[i]);
            
            // Store in column-major order: h[i][k] = h[i + k * (restart + 1)]
            h[i + k * (restart + 1)] = hij.clone();
            
            // In-place subtraction: v = v - hij * q[i]
            v.axpy(-hij.clone(), &q[i], T::one());
            
            // Re-orthogonalize for better numerical stability
            let hij_correction = v.dot(&q[i]);
            if ComplexField::abs(hij_correction.clone()) > T::from_f64(1e-10).unwrap() {
                h[i + k * (restart + 1)] = h[i + k * (restart + 1)].clone() + hij_correction.clone();
                v.axpy(-hij_correction, &q[i], T::one());
            }
        }

        // Normalize
        let norm = v.norm();
        if norm < T::from_f64(1e-14).unwrap() {
            // Handle breakdown - vector is in the span of previous vectors
            h[(k + 1) + k * (restart + 1)] = T::zero();
            Ok(DVector::zeros(v.len()))
        } else {
            h[(k + 1) + k * (restart + 1)] = norm.clone();
            Ok(v / norm)
        }
    }

    /// Apply previous Givens rotations to current column
    fn apply_previous_rotations(&self, h: &mut Vec<T>, cs: &[T], sn: &[T], k: usize, restart: usize) {
        for i in 0..k {
            let idx_i = i + k * (restart + 1);
            let idx_ip1 = (i + 1) + k * (restart + 1);
            
            let (h_new, h_ip1_new) = Self::apply_givens_rotation(
                cs[i].clone(),
                sn[i].clone(),
                h[idx_i].clone(),
                h[idx_ip1].clone(),
            );
            h[idx_i] = h_new;
            h[idx_ip1] = h_ip1_new;
        }
    }

    /// Compute new Givens rotation coefficients
    fn compute_givens_rotation(&self, h: &mut Vec<T>, k: usize, restart: usize) -> Result<(T, T)> {
        let idx_k = k + k * (restart + 1);
        let idx_kp1 = (k + 1) + k * (restart + 1);
        
        let h_k = h[idx_k].clone();
        let h_kp1 = h[idx_kp1].clone();
        
        // Improved numerical stability for Givens rotation
        if ComplexField::abs(h_kp1.clone()) < T::from_f64(1e-14).unwrap() {
            // h_kp1 is essentially zero
            let c = if ComplexField::abs(h_k.clone()) < T::from_f64(1e-14).unwrap() {
                T::one()
            } else {
                ComplexField::signum(h_k.clone())
            };
            h[idx_k] = ComplexField::abs(h_k);
            h[idx_kp1] = T::zero();
            return Ok((c, T::zero()));
        }
        
        // Use more stable formulation
        if ComplexField::abs(h_kp1.clone()) > ComplexField::abs(h_k.clone()) {
            let tau = h_k.clone() / h_kp1.clone();
            let s = T::one() / ComplexField::sqrt(T::one() + tau.clone() * tau.clone());
            let c = s.clone() * tau;
            h[idx_k] = h_kp1.clone() / s;
            h[idx_kp1] = T::zero();
            Ok((c, s))
        } else {
            let tau = h_kp1.clone() / h_k.clone();
            let c = T::one() / ComplexField::sqrt(T::one() + tau.clone() * tau.clone());
            let s = c.clone() * tau;
            h[idx_k] = h_k.clone() / c;
            h[idx_kp1] = T::zero();
            Ok((c, s))
        }
    }

    /// Update residual vector using Givens rotation
    fn update_residual(&self, g: &mut DVector<T>, c: T, s: T, k: usize) {
        g[k + 1] = -s.clone() * g[k].clone();
        g[k] = c * g[k].clone();
    }

    /// Solve upper triangular system by back substitution
    fn solve_upper_triangular(&self, h: &[T], g: &DVector<T>, size: usize, restart: usize) -> Result<DVector<T>> {
        let mut y: DVector<T> = DVector::zeros(size);
        for i in (0..size).rev() {
            let mut yi = g[i].clone();
            for j in (i + 1)..size {
                let h_ij = h[i + j * (restart + 1)];
                yi = yi - h_ij.clone() * y[j].clone();
            }
            let h_ii = h[i + i * (restart + 1)];
            if h_ii == T::zero() {
                return Err(Error::NumericalError("Singular matrix in GMRES".to_string()));
            }
            y[i] = yi / h_ii.clone();
        }
        Ok(y)
    }

    /// Compute correction vector from Krylov basis and coefficients
    fn compute_correction(&self, q: &[DVector<T>], y: &DVector<T>) -> DVector<T> {
        let mut correction = DVector::zeros(q[0].len());
        for (i, yi) in y.iter().enumerate() {
            correction.axpy(yi.clone(), &q[i], T::one());
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
        let initial_residual_norm = r.norm();
        
        // Compute breakdown tolerance relative to initial residual
        let breakdown_tolerance = initial_residual_norm * T::from_f64(1e-14).unwrap();
        
        let r0_hat = r.clone(); // Arbitrary choice, could be different
        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n);
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n);

        for iter in 0..self.config.max_iterations() {
            let rho_new = r0_hat.dot(&r);
            
            // Use relative tolerance for breakdown detection
            if rho_new.clone().abs() < breakdown_tolerance {
                return Err(Error::NumericalError(
                    "BiCGSTAB breakdown: rho is close to zero relative to initial residual".to_string(),
                ));
            }

            let beta = (rho_new.clone() / rho.clone()) * (alpha.clone() / omega.clone());
            
            // Update search direction using in-place operations
            // p = r + beta * (p - omega * v)
            p.axpy(-omega.clone(), &v, T::one()); // p = p - omega * v
            p *= beta;                      // p = beta * p
            p += &r;                        // p = p + r
            
            // Matrix-vector product
            v = a * &p;
            alpha = rho_new.clone() / r0_hat.dot(&v);
            
            // Update s using in-place operations
            s.copy_from(&r);
            s.axpy(-alpha.clone(), &v, T::one()); // s = r - alpha * v
            
            // Check convergence
            let s_norm = s.norm();
            if self.is_converged(s_norm) {
                x.axpy(alpha.clone(), &p, T::one()); // x = x + alpha * p
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }
            
            // Second matrix-vector product
            t = a * &s;
            omega = s.dot(&t) / t.dot(&t);
            
            // Update solution using in-place operations
            x.axpy(alpha.clone(), &p, T::one());  // x = x + alpha * p
            x.axpy(omega.clone(), &s, T::one());  // x = x + omega * s
            
            // Update residual using in-place operations
            r.copy_from(&s);
            r.axpy(-omega.clone(), &t, T::one()); // r = s - omega * t
            
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