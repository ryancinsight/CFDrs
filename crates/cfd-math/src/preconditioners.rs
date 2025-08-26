//! Extended preconditioners for iterative linear solvers
//!
//! This module provides state-of-the-art preconditioning techniques
//! essential for practical CFD applications.

use nalgebra::{DVector, RealField};
use cfd_core::numeric;
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use crate::linear_solver::{LinearSolverError, Result, Preconditioner};
/// Incomplete Cholesky factorization preconditioner (IC(0))
/// 
/// For symmetric positive definite matrices, computes L such that L*L^T ≈ A
/// where L has the same sparsity pattern as the lower triangular part of A.
pub struct IncompleteCholesky<T: RealField + Copy> {
    /// Lower triangular factor stored in CSR format
    l_factor: CsrMatrix<T>,
    /// Workspace for forward/backward substitution
    workspace: Vec<T>,
}
impl<T: RealField + Copy> IncompleteCholesky<T> {
    /// Construct IC(0) preconditioner from symmetric positive definite matrix
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        // Validate matrix is square
        if a.nrows() != a.ncols() {
            return Err(LinearSolverError::NonSquareMatrix {
                rows: a.nrows(),
                cols: a.ncols(),
            });
        }
        
        let n = a.nrows();
        // Check symmetry (with tolerance)
        let tol = cfd_core::numeric::from_f64(1e-10)?;
        let mut max_asymmetry = T::zero();
        let mut max_i = 0;
        let mut max_j = 0;
        for i in 0..n {
            for j in i+1..n {
                let a_ij = a.get_entry(i, j).unwrap_or(T::zero());
                let a_ji = a.get_entry(j, i).unwrap_or(T::zero());
                let diff = (a_ij - a_ji).abs();
                
                if diff > max_asymmetry {
                    max_asymmetry = diff;
                    max_i = i;
                    max_j = j;
                }
            }
        if max_asymmetry > tol {
            return Err(LinearSolverError::NonSymmetric {
                max_diff: max_asymmetry.to_subset().unwrap_or(0.0),
                i: max_i,
                j: max_j,
        // Build L factor with same sparsity as lower triangular A
        let mut l_vals = Vec::new();
        let mut l_indices = Vec::new();
        let mut l_offsets = vec![0];
        // Temporary storage for rows of L
        let mut l_rows: Vec<Vec<(usize, T)>> = vec![Vec::new(); n];
        // IC(0) factorization algorithm
            let mut l_ii = a.get_entry(i, i).unwrap_or(T::zero());
            
            // Subtract contributions from previous columns
            for &(k, l_ik) in &l_rows[i] {
                if k < i {
                    l_ii -= l_ik * l_ik;
            // Check for positive definiteness
            if l_ii <= T::zero() {
                return Err(LinearSolverError::SingularMatrix);
            l_ii = l_ii.sqrt();
            l_rows[i].push((i, l_ii));
            // Compute entries in column i below diagonal
            let row_start = a.row_offsets()[i];
            let row_end = a.row_offsets()[i + 1];
            for idx in row_start..row_end {
                let j = a.col_indices()[idx];
                if j > i {  // Only process lower triangular part
                    let mut l_ji = a.get_entry(j, i).unwrap_or(T::zero());
                    
                    // Subtract contributions from previous columns
                    for &(k, l_ik) in &l_rows[i] {
                        if k < i {
                            // Find l_jk
                            let l_jk = l_rows[j].iter()
                                .find(|(col, _)| *col == k)
                                .map(|(_, val)| *val)
                                .unwrap_or(T::zero());
                            l_ji -= l_ik * l_jk;
                        }
                    }
                    l_ji = l_ji / l_ii;
                    l_rows[j].push((i, l_ji));
        // Convert to CSR format
            l_rows[i].sort_by_key(|(j, _)| *j);
            for (j, val) in &l_rows[i] {
                l_indices.push(*j);
                l_vals.push(*val);
            l_offsets.push(l_indices.len());
        let l_factor = CsrMatrix::try_from_csr_data(
            n, n, l_offsets, l_indices, l_vals
        ).map_err(|e| LinearSolverError::FactorizationFailed {
            reason: format!("Failed to create L factor: {}", e),
        })?;
        Ok(Self {
            l_factor,
            workspace: vec![T::zero(); n],
        })
    }
    
    /// Forward substitution: solve L*y = b
    fn forward_substitution(&self, b: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        let n = self.l_factor.nrows();
            let mut sum = b[i];
            let row_start = self.l_factor.row_offsets()[i];
            let row_end = self.l_factor.row_offsets()[i + 1];
            let mut diag = T::one();
            for k in row_start..row_end {
                let j = self.l_factor.col_indices()[k];
                let val = self.l_factor.values()[k];
                if j < i {
                    sum -= val * y[j];
                } else if j == i {
                    diag = val;
            if diag.abs() < cfd_core::numeric::from_f64(1e-14)? {
                return Err(LinearSolverError::ZeroDiagonalElement {
                    index: i,
                    value: diag.to_subset().unwrap_or(0.0),
                });
            y[i] = sum / diag;
        Ok(())
    /// Backward substitution: solve L^T*x = y
    fn backward_substitution(&self, y: &DVector<T>, x: &mut DVector<T>) -> Result<()> {
        // Process rows in reverse order for L^T
        for i in (0..n).rev() {
            let mut sum = y[i];
            // For L^T, we need to access column i of L
            // This requires iterating through all rows looking for entries in column i
                let row_start = self.l_factor.row_offsets()[j];
                let row_end = self.l_factor.row_offsets()[j + 1];
                for k in row_start..row_end {
                    if self.l_factor.col_indices()[k] == i {
                        sum -= self.l_factor.values()[k] * x[j];
                        break;
            // Get diagonal element
                if self.l_factor.col_indices()[k] == i {
                    diag = self.l_factor.values()[k];
                    break;
            x[i] = sum / diag;
impl<T: RealField + Copy> Preconditioner<T> for IncompleteCholesky<T> {
    }

    fn apply(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let mut y = DVector::zeros(r.len());
        // Solve L*y = r
        self.forward_substitution(r, &mut y)?;
        // Solve L^T*z = y
        self.backward_substitution(&y, z)?;
    }

    fn setup(&mut self, a: &CsrMatrix<T>) -> Result<()> {
        *self = Self::new(a)?;
/// SSOR (Symmetric Successive Over-Relaxation) preconditioner
///
/// For symmetric matrices, applies forward and backward SOR sweeps
    }

}

pub struct SSORPreconditioner<T: RealField + Copy> {
    matrix: CsrMatrix<T>,
    omega: T,
    diagonal: Vec<T>,
impl<T: RealField + Copy + FromPrimitive> SSORPreconditioner<T> {
    /// Create SSOR preconditioner with specified relaxation parameter
    ///
    /// # Optimal omega values
    /// - For Poisson equation: ω ≈ 2/(1 + sin(π/n))
    /// - For general problems: typically ω ∈ [1.0, 1.5]
    /// - ω = 1.0 gives symmetric Gauss-Seidel
    pub fn new(a: &CsrMatrix<T>, omega: T) -> Result<Self> {
        // Validate omega parameter
        if omega <= T::zero() || omega >= cfd_core::numeric::from_f64(2.0)? {
            return Err(LinearSolverError::InvalidOmega {
                omega: omega.to_subset().unwrap_or(0.0),
        // Extract diagonal
        let mut diagonal = Vec::with_capacity(a.nrows());
        for i in 0..a.nrows() {
            let diag = a.get_entry(i, i).unwrap_or(T::zero());
            diagonal.push(diag);
            matrix: a,
            omega,
            diagonal,
    /// Forward SOR sweep
    }

    fn forward_sweep(&self, b: &DVector<T>, x: &mut DVector<T>) {
        let n = self.matrix.nrows();
            let row_start = self.matrix.row_offsets()[i];
            let row_end = self.matrix.row_offsets()[i + 1];
                let j = self.matrix.col_indices()[k];
                if j != i {
                    sum -= self.matrix.values()[k] * x[j];
            x[i] = (T::one() - self.omega) * x[i] + self.omega * sum / self.diagonal[i];
    /// Backward SOR sweep
    }

    fn backward_sweep(&self, b: &DVector<T>, x: &mut DVector<T>) {
impl<T: RealField + Copy + FromPrimitive> Preconditioner<T> for SSORPreconditioner<T> {
        // Initialize z with zeros
        z.fill(T::zero());
        // Forward sweep
        self.forward_sweep(r, z);
        // Backward sweep
        self.backward_sweep(r, z);
        *self = Self::new(a, self.omega)?;
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;
    fn create_spd_matrix() -> CsrMatrix<f64> {
        // Create a simple SPD matrix: tridiagonal [-1, 2, -1]
        let n = 5;
        let mut triplets = Vec::new();
            triplets.push((i, i, 2.0));
            if i > 0 {
                triplets.push((i, i-1, -1.0));
            if i < n-1 {
                triplets.push((i, i+1, -1.0));
        CsrMatrix::try_from_triplets(n, n, triplets).expect("CRITICAL: Add proper error handling")
    #[test]
    }

    fn test_incomplete_cholesky() {
        let a = create_spd_matrix();
        let ic = IncompleteCholesky::new(&a);
        assert!(ic.is_ok());
    }

    fn test_ssor_preconditioner() {
        let ssor = SSORPreconditioner::new(&a, 1.2);
        assert!(ssor.is_ok());
        // Test invalid omega
        let invalid = SSORPreconditioner::new(&a, 2.5);
        assert!(matches!(invalid, Err(LinearSolverError::InvalidOmega { .. })));


}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
}
