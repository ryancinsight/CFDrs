//! ILU(k) factorization with level-k fill
//!
//! This implements the symbolic/numeric factorization algorithm from Saad (2003) §10.4.
//! Fill levels are tracked: lev(i,j) = min{lev(i,j), lev(i,k) + lev(k,j) + 1}
//! Fill is kept only if lev(i,j) <= k.
//!
//! Reference: Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). SIAM, §10.4.

use cfd_core::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use nalgebra_sparse::CsrMatrix;
use std::collections::HashMap;

use super::utils;

/// Perform ILU(k) factorization with level-k fill
///
/// # Arguments
///
/// * `a` - Input sparse matrix in CSR format
/// * `k` - Fill level (controls sparsity: 0 = no fill, higher k = more fill)
///
/// # Returns
///
/// LU factors as a single CSR matrix (L below diagonal with unit diagonal, U on and above diagonal)
pub fn factorize<T: RealField + Copy>(a: &CsrMatrix<T>, k: usize) -> Result<CsrMatrix<T>> {
    let n = a.nrows();

    // Phase 1: Symbolic factorization - determine sparsity pattern with level-k fill
    let levels = symbolic_phase(a, k);

    // Phase 2: Build CSR structure for the new sparsity pattern
    let (lu_offsets, lu_indices, mut lu_vals) = build_csr_structure(a, &levels, k, n);

    // Phase 3: Numeric factorization on the extended sparsity pattern
    numeric_phase(&mut lu_vals, &lu_offsets, &lu_indices, n)?;

    CsrMatrix::try_from_csr_data(n, n, lu_offsets, lu_indices, lu_vals).map_err(|_| {
        Error::InvalidInput("Failed to create CSR matrix from ILU(k) factorization".to_string())
    })
}

/// Phase 1: Symbolic factorization - determine sparsity pattern
fn symbolic_phase<T: RealField + Copy>(
    a: &CsrMatrix<T>,
    k: usize,
) -> HashMap<(usize, usize), usize> {
    let n = a.nrows();
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
        let row_keys: Vec<(usize, usize)> = levels
            .keys()
            .filter(|&&(row, col)| row == i && col < i)
            .copied()
            .collect();

        for &(_, j) in &row_keys {
            let lev_ij = levels[&(i, j)];

            // Update levels for positions (i, m) where m > j
            // based on fill from multiplication L_ij * U_jm
            let upper_keys: Vec<(usize, usize)> = levels
                .keys()
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

    levels
}

/// Phase 2: Build CSR structure for new sparsity pattern
#[allow(clippy::type_complexity)]
fn build_csr_structure<T: RealField + Copy>(
    a: &CsrMatrix<T>,
    levels: &HashMap<(usize, usize), usize>,
    k: usize,
    n: usize,
) -> (Vec<usize>, Vec<usize>, Vec<T>) {
    // Filter to keep only entries with level <= k
    let mut kept_entries: Vec<(usize, usize)> = levels
        .iter()
        .filter(|&(_, &level)| level <= k)
        .map(|(&pos, _)| pos)
        .collect();
    kept_entries.sort_by_key(|&(i, j)| (i, j));

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

    (row_offsets, col_indices, values)
}

/// Phase 3: Numeric factorization on extended sparsity pattern
fn numeric_phase<T: RealField + Copy>(
    lu_vals: &mut [T],
    lu_offsets: &[usize],
    lu_indices: &[usize],
    n: usize,
) -> Result<()> {
    for k_idx in 0..n - 1 {
        let diag_idx = utils::find_diagonal_index(lu_offsets, lu_indices, k_idx)?;
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

    Ok(())
}
