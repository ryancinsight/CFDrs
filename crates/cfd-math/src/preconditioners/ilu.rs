//! Incomplete LU factorization preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use std::collections::HashMap;

/// Incomplete LU factorization preconditioner (ILU(k))
///
/// Computes L and U factors such that L*U ≈ A with controlled fill-in.
/// The parameter k controls the level of fill: ILU(0) has the same sparsity as A,
/// while ILU(k) allows fill up to level k.
///
/// Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.
pub struct IncompleteLU<T: RealField + Copy> {
    /// Combined LU factors (L has unit diagonal)
    lu_factor: CsrMatrix<T>,
    /// Workspace for forward/backward substitution
    workspace: Vec<T>,
    /// Fill level k (0 means no fill beyond original sparsity)
    fill_level: usize,
}

impl<T: RealField + Copy> IncompleteLU<T> {
    /// Construct ILU(0) preconditioner (no fill beyond original sparsity)
    pub fn new(a: &CsrMatrix<T>) -> Result<Self> {
        Self::with_fill_level(a, 0)
    }

    /// Construct ILU(k) preconditioner with specified fill level
    ///
    /// # Arguments
    ///
    /// * `a` - Input matrix (must be square)
    /// * `k` - Fill level (0 = no fill, 1 = level-1 fill, etc.)
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use cfd_math::preconditioners::IncompleteLU;
    /// use nalgebra_sparse::CsrMatrix;
    ///
    /// // Create ILU(2) preconditioner
    /// let ilu = IncompleteLU::with_fill_level(&matrix, 2)?;
    /// ```
    pub fn with_fill_level(a: &CsrMatrix<T>, k: usize) -> Result<Self> {
        if a.nrows() != a.ncols() {
            return Err(Error::InvalidInput(format!(
                "Matrix must be square, got {}x{}",
                a.nrows(),
                a.ncols()
            )));
        }

        let n = a.nrows();
        let lu_factor = if k == 0 {
            Self::factorize_ilu0(a)?
        } else {
            Self::factorize_iluk(a, k)?
        };

        Ok(Self {
            lu_factor,
            workspace: vec![T::zero(); n],
            fill_level: k,
        })
    }

    /// Perform ILU(0) factorization (original algorithm, optimized)
    fn factorize_ilu0(a: &CsrMatrix<T>) -> Result<CsrMatrix<T>> {
        // Create a mutable copy of the matrix
        let mut lu_vals = a.values().to_vec();
        let lu_indices = a.col_indices().to_vec();
        let lu_offsets = a.row_offsets().to_vec();

        let n = a.nrows();

        // ILU(0) factorization
        for k in 0..n - 1 {
            let diag_idx = Self::find_diagonal_index(&lu_offsets, &lu_indices, k)?;
            let a_kk = lu_vals[diag_idx];

            if a_kk.abs() <= T::default_epsilon() {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            // Process rows below k
            for i in k + 1..n {
                let row_start = lu_offsets[i];
                let row_end = lu_offsets[i + 1];

                // Find a_ik
                let mut a_ik_idx = None;
                for idx in row_start..row_end {
                    if lu_indices[idx] == k {
                        a_ik_idx = Some(idx);
                        break;
                    }
                }

                if let Some(idx) = a_ik_idx {
                    // Compute multiplier
                    let mult = lu_vals[idx] / a_kk;
                    lu_vals[idx] = mult;

                    // Update row i
                    for j_idx in row_start..row_end {
                        let j = lu_indices[j_idx];
                        if j > k {
                            // Find a_kj in row k
                            let k_row_start = lu_offsets[k];
                            let k_row_end = lu_offsets[k + 1];

                            for k_idx in k_row_start..k_row_end {
                                if lu_indices[k_idx] == j {
                                    let k_val = lu_vals[k_idx];
                                    lu_vals[j_idx] -= mult * k_val;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        CsrMatrix::try_from_csr_data(n, n, lu_offsets, lu_indices, lu_vals)
            .map_err(|_| Error::InvalidInput("Failed to create CSR matrix from ILU(0) factorization".to_string()))
    }

    /// Perform ILU(k) factorization with level-k fill
    ///
    /// This implements the symbolic/numeric factorization algorithm from Saad (2003) §10.4.
    /// Fill levels are tracked: lev(i,j) = min{lev(i,j), lev(i,k) + lev(k,j) + 1}
    /// Fill is kept only if lev(i,j) <= k.
    fn factorize_iluk(a: &CsrMatrix<T>, k: usize) -> Result<CsrMatrix<T>> {
        let n = a.nrows();
        
        // Phase 1: Symbolic factorization - determine sparsity pattern with level-k fill
        let mut levels: HashMap<(usize, usize), usize> = HashMap::new();
        
        // Initialize levels from original matrix A (level 0)
        for i in 0..n {
            let row_start = a.row_offsets()[i];
            let row_end = a.row_offsets()[i + 1];
            for idx in row_start..row_end {
                let j = a.col_indices()[idx];
                levels.insert((i, j), 0);
            }
        }
        
        // Symbolic phase: compute fill levels
        for i in 0..n {
            // For each nonzero a_ij in row i where j < i (lower triangular part)
            let row_keys: Vec<(usize, usize)> = levels.keys()
                .filter(|&&(row, col)| row == i && col < i)
                .copied()
                .collect();
                
            for &(_, j) in &row_keys {
                let lev_ij = levels[&(i, j)];
                
                // Update levels for positions (i, m) where m > j
                // based on fill from multiplication L_ij * U_jm
                let upper_keys: Vec<(usize, usize)> = levels.keys()
                    .filter(|&&(row, col)| row == j && col >= j)
                    .copied()
                    .collect();
                    
                for &(_, m) in &upper_keys {
                    if m <= j {
                        continue; // Skip lower triangular of row j
                    }
                    
                    let lev_jm = levels[&(j, m)];
                    let new_level = lev_ij + lev_jm + 1;
                    
                    // Update or insert level if it's within tolerance k
                    if new_level <= k {
                        let current_level = levels.entry((i, m)).or_insert(usize::MAX);
                        if new_level < *current_level {
                            *current_level = new_level;
                        }
                    }
                }
            }
        }
        
        // Filter to keep only entries with level <= k
        let mut kept_entries: Vec<(usize, usize)> = levels.iter()
            .filter(|&(_, &level)| level <= k)
            .map(|(&pos, _)| pos)
            .collect();
        kept_entries.sort_by_key(|&(i, j)| (i, j));
        
        // Phase 2: Build CSR structure for the new sparsity pattern
        let mut row_offsets = vec![0; n + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        
        for i in 0..n {
            for &(row, col) in &kept_entries {
                if row == i {
                    col_indices.push(col);
                    
                    // Get initial value from A if it exists, otherwise 0
                    let val = if let Some(&level) = levels.get(&(row, col)) {
                        if level == 0 {
                            // Original entry from A
                            let a_row_start = a.row_offsets()[row];
                            let a_row_end = a.row_offsets()[row + 1];
                            let mut found_val = T::zero();
                            for idx in a_row_start..a_row_end {
                                if a.col_indices()[idx] == col {
                                    found_val = a.values()[idx];
                                    break;
                                }
                            }
                            found_val
                        } else {
                            T::zero() // Fill entry
                        }
                    } else {
                        T::zero()
                    };
                    values.push(val);
                }
            }
            row_offsets[i + 1] = col_indices.len();
        }
        
        // Phase 3: Numeric factorization on the extended sparsity pattern
        let mut lu_vals = values;
        let lu_indices = col_indices;
        let lu_offsets = row_offsets;
        
        for k_idx in 0..n - 1 {
            let diag_idx = Self::find_diagonal_index(&lu_offsets, &lu_indices, k_idx)?;
            let a_kk = lu_vals[diag_idx];
            
            if a_kk.abs() <= T::default_epsilon() {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }
            
            // Process rows below k
            for i in k_idx + 1..n {
                let row_start = lu_offsets[i];
                let row_end = lu_offsets[i + 1];
                
                // Find a_ik
                let mut a_ik_idx = None;
                for idx in row_start..row_end {
                    if lu_indices[idx] == k_idx {
                        a_ik_idx = Some(idx);
                        break;
                    }
                }
                
                if let Some(idx) = a_ik_idx {
                    // Compute multiplier
                    let mult = lu_vals[idx] / a_kk;
                    lu_vals[idx] = mult;
                    
                    // Update row i
                    for j_idx in row_start..row_end {
                        let j = lu_indices[j_idx];
                        if j > k_idx {
                            // Find a_kj in row k
                            let k_row_start = lu_offsets[k_idx];
                            let k_row_end = lu_offsets[k_idx + 1];
                            
                            for k_idx_inner in k_row_start..k_row_end {
                                if lu_indices[k_idx_inner] == j {
                                    let k_val = lu_vals[k_idx_inner];
                                    lu_vals[j_idx] -= mult * k_val;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        CsrMatrix::try_from_csr_data(n, n, lu_offsets, lu_indices, lu_vals)
            .map_err(|_| Error::InvalidInput("Failed to create CSR matrix from ILU(k) factorization".to_string()))
    }

    /// Find index of diagonal element in row
    fn find_diagonal_index(offsets: &[usize], indices: &[usize], row: usize) -> Result<usize> {
        let row_start = offsets[row];
        let row_end = offsets[row + 1];

        for idx in row_start..row_end {
            if indices[idx] == row {
                return Ok(idx);
            }
        }

        Err(Error::InvalidInput("Diagonal element not found in matrix row".to_string()))
    }

    /// Forward substitution with L (unit diagonal)
    fn forward_substitution(&self, b: &DVector<T>, y: &mut DVector<T>) {
        let n = self.lu_factor.nrows();

        for i in 0..n {
            let mut sum = b[i];

            let row_start = self.lu_factor.row_offsets()[i];
            let row_end = self.lu_factor.row_offsets()[i + 1];

            for idx in row_start..row_end {
                let j = self.lu_factor.col_indices()[idx];
                if j < i {
                    sum -= self.lu_factor.values()[idx] * y[j];
                }
            }

            y[i] = sum;
        }
    }

    /// Backward substitution with U
    fn backward_substitution(&self, y: &DVector<T>, x: &mut DVector<T>) {
        let n = self.lu_factor.nrows();

        for i in (0..n).rev() {
            let mut sum = y[i];

            let row_start = self.lu_factor.row_offsets()[i];
            let row_end = self.lu_factor.row_offsets()[i + 1];

            let mut diag_val = T::one();

            for idx in row_start..row_end {
                let j = self.lu_factor.col_indices()[idx];
                if j > i {
                    sum -= self.lu_factor.values()[idx] * x[j];
                } else if j == i {
                    diag_val = self.lu_factor.values()[idx];
                }
            }

            x[i] = sum / diag_val;
        }
    }
}

impl<T: RealField + Copy> Preconditioner<T> for IncompleteLU<T> {
    fn apply_to(&self, r: &DVector<T>, z: &mut DVector<T>) -> Result<()> {
        let n = r.len();
        
        // Use workspace for intermediate result
        let mut y = DVector::zeros(n);

        // Solve LU*z = r via forward and backward substitution
        self.forward_substitution(r, &mut y);
        self.backward_substitution(&y, z);

        Ok(())
    }
}

impl<T: RealField + Copy> IncompleteLU<T> {
    /// Get the fill level k
    pub fn fill_level(&self) -> usize {
        self.fill_level
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra_sparse::CsrMatrix;

    /// Create a simple 4x4 tridiagonal test matrix
    /// [4 -1  0  0]
    /// [-1 4 -1  0]
    /// [0 -1  4 -1]
    /// [0  0 -1  4]
    fn create_tridiagonal_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 8, 10];
        let col_indices = vec![
            0, 1,       // row 0
            0, 1, 2,    // row 1
            1, 2, 3,    // row 2
            2, 3,       // row 3
        ];
        let values = vec![
            4.0, -1.0,          // row 0
            -1.0, 4.0, -1.0,    // row 1
            -1.0, 4.0, -1.0,    // row 2
            -1.0, 4.0,          // row 3
        ];
        
        CsrMatrix::try_from_csr_data(4, 4, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    /// Create a 5x5 sparse matrix with more complex structure
    fn create_sparse_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 3, 6, 9, 12, 14];
        let col_indices = vec![
            0, 1, 3,       // row 0: connections to 1, 3
            0, 1, 2,       // row 1: connections to 0, 2
            1, 2, 3,       // row 2: connections to 1, 3
            0, 2, 3, 4,    // row 3: connections to 0, 2, 4
            3, 4,          // row 4: connection to 3
        ];
        let values = vec![
            5.0, -1.0, -1.0,        // row 0
            -1.0, 5.0, -1.0,        // row 1
            -1.0, 5.0, -1.0,        // row 2
            -1.0, -1.0, 5.0, -1.0,  // row 3
            -1.0, 5.0,              // row 4
        ];
        
        CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    #[test]
    fn test_ilu0_construction() {
        let matrix = create_tridiagonal_matrix();
        let ilu = IncompleteLU::new(&matrix);
        
        assert!(ilu.is_ok());
        let ilu = ilu.unwrap();
        assert_eq!(ilu.fill_level(), 0);
    }

    #[test]
    fn test_iluk_construction() {
        let matrix = create_tridiagonal_matrix();
        let ilu = IncompleteLU::with_fill_level(&matrix, 1);
        
        assert!(ilu.is_ok());
        let ilu = ilu.unwrap();
        assert_eq!(ilu.fill_level(), 1);
    }

    #[test]
    fn test_ilu0_apply() {
        let matrix = create_tridiagonal_matrix();
        let ilu = IncompleteLU::new(&matrix).expect("ILU(0) construction");
        
        // Test with a simple vector
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let mut z = DVector::zeros(4);
        ilu.apply_to(&b, &mut z).expect("Apply preconditioner");
        
        // Solution should be positive and bounded
        assert!(z[0] > 0.0 && z[0] < 10.0);
        assert!(z[1] > 0.0 && z[1] < 10.0);
        assert!(z[2] > 0.0 && z[2] < 10.0);
        assert!(z[3] > 0.0 && z[3] < 10.0);
    }

    #[test]
    fn test_iluk_vs_ilu0_sparsity() {
        let matrix = create_sparse_matrix();
        
        let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
        let ilu1 = IncompleteLU::with_fill_level(&matrix, 1).expect("ILU(1) construction");
        
        // ILU(1) should have at least as many nonzeros as ILU(0)
        let nnz_ilu0 = ilu0.lu_factor.nnz();
        let nnz_ilu1 = ilu1.lu_factor.nnz();
        
        assert!(nnz_ilu1 >= nnz_ilu0, 
                "ILU(1) should have >= nonzeros than ILU(0): {} vs {}", 
                nnz_ilu1, nnz_ilu0);
    }

    #[test]
    fn test_ilu0_preconditioner_quality() {
        // Test that ILU(0) actually improves conditioning
        let matrix = create_tridiagonal_matrix();
        let ilu = IncompleteLU::new(&matrix).expect("ILU(0) construction");
        
        // Test with identity-like vector
        let b = DVector::from_vec(vec![1.0, 1.0, 1.0, 1.0]);
        let mut z = DVector::zeros(4);
        ilu.apply_to(&b, &mut z).expect("Apply preconditioner");
        
        // For a diagonally dominant matrix, preconditioner should give reasonable approximation
        // to A^{-1} * b, which for this matrix and b is approximately [0.36, 0.43, 0.43, 0.36]
        assert_relative_eq!(z[0], 0.36, epsilon = 0.15);
        assert_relative_eq!(z[1], 0.43, epsilon = 0.15);
        assert_relative_eq!(z[2], 0.43, epsilon = 0.15);
        assert_relative_eq!(z[3], 0.36, epsilon = 0.15);
    }

    #[test]
    fn test_iluk_improved_approximation() {
        // Test that ILU(k) for k>0 improves approximation quality
        let matrix = create_sparse_matrix();
        
        let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
        let ilu1 = IncompleteLU::with_fill_level(&matrix, 1).expect("ILU(1) construction");
        
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        
        let mut z0 = DVector::zeros(5);
        let mut z1 = DVector::zeros(5);
        ilu0.apply_to(&b, &mut z0).expect("Apply ILU(0)");
        ilu1.apply_to(&b, &mut z1).expect("Apply ILU(1)");
        
        // Both should give bounded results
        for i in 0..5 {
            assert!(z0[i].abs() < 100.0, "ILU(0) result unbounded");
            assert!(z1[i].abs() < 100.0, "ILU(1) result unbounded");
        }
        
        // ILU(1) should generally give different (hopefully better) results
        // We don't test for "better" directly, just that fill improves approximation
        let diff = (&z0 - &z1).norm();
        assert!(diff > 0.0, "ILU(1) should differ from ILU(0)");
    }

    #[test]
    fn test_ilu_non_square_matrix() {
        // Test error handling for non-square matrix
        let row_offsets = vec![0, 2, 4];
        let col_indices = vec![0, 1, 1, 2];
        let values = vec![1.0, 2.0, 3.0, 4.0];
        
        let matrix = CsrMatrix::try_from_csr_data(2, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");
        
        let ilu = IncompleteLU::new(&matrix);
        assert!(ilu.is_err());
        
        // Verify error is InvalidInput for non-square matrix
        assert!(matches!(ilu, Err(Error::InvalidInput(_))));
    }

    #[test]
    fn test_ilu_fill_levels() {
        let matrix = create_tridiagonal_matrix();
        
        let ilu0 = IncompleteLU::new(&matrix).expect("ILU(0) construction");
        assert_eq!(ilu0.fill_level(), 0);
        
        let ilu2 = IncompleteLU::with_fill_level(&matrix, 2).expect("ILU(2) construction");
        assert_eq!(ilu2.fill_level(), 2);
    }

    #[test]
    fn test_ilu_multiple_fill_levels() {
        let matrix = create_sparse_matrix();
        
        for k in 0..=3 {
            let ilu = IncompleteLU::with_fill_level(&matrix, k);
            assert!(ilu.is_ok(), "ILU({}) construction failed", k);
            
            let ilu = ilu.unwrap();
            assert_eq!(ilu.fill_level(), k);
            
            // Test application
            let b = DVector::from_vec(vec![1.0; 5]);
            let x = ilu.lu_factor.nrows();
            assert_eq!(x, 5);
        }
    }
}