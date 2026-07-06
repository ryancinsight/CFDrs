//! Incomplete Cholesky factorization preconditioner

use crate::linear_solver::Preconditioner;
use cfd_core::error::{Error, NumericalErrorKind, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, Scalar as LetoScalar};

// Tolerance for symmetry checking
const SYMMETRY_TOLERANCE: f64 = 1e-10;

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
fn vector_len<T>(vector: &Array1<T>) -> usize {
    vector.shape()[0]
}

fn validate_vector_len<T>(name: &str, vector: &Array1<T>, expected: usize) -> Result<()> {
    let actual = vector_len(vector);
    if actual != expected {
        return Err(Error::InvalidConfiguration(format!(
            "{name} length mismatch: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

fn csr_value<T: LetoScalar>(matrix: &CsrMatrix<T>, row: usize, col: usize) -> T {
    let matrix_row = matrix.row(row);
    for (&candidate_col, &value) in matrix_row.col_indices().iter().zip(matrix_row.values()) {
        if candidate_col == col {
            return value;
        }
    }
    <T as NumericElement>::ZERO
}

/// Incomplete Cholesky factorization preconditioner (IC(0))
///
/// For symmetric positive definite matrices, computes L such that L*L^T ≈ A
/// where L has the same sparsity pattern as the lower triangular part of A.
pub struct IncompleteCholesky<T: RealField + Copy + LetoScalar> {
    /// Lower triangular factor stored in CSR format
    l_factor: CsrMatrix<T>,
}

impl<T: RealField + Copy + LetoScalar> IncompleteCholesky<T> {
    /// Construct IC(0) preconditioner from symmetric positive definite matrix
    pub fn new(a: &CsrMatrix<T>) -> Result<Self>
    where
        T: FloatElement,
    {
        // Validate matrix is square
        if a.nrows() != a.ncols() {
            return Err(Error::InvalidInput(format!(
                "Matrix must be square, got {}x{}",
                a.nrows(),
                a.ncols()
            )));
        }

        // Check symmetry
        Self::check_symmetry(a)?;

        // Perform IC(0) factorization
        let l_factor = Self::factorize(a)?;

        Ok(Self { l_factor })
    }

    /// Check matrix symmetry
    fn check_symmetry(a: &CsrMatrix<T>) -> Result<()>
    where
        T: FloatElement,
    {
        let tol = from_f64(SYMMETRY_TOLERANCE);
        let n = a.nrows();
        let mut max_asymmetry = <T as NumericElement>::ZERO;
        let mut max_i = 0;
        let mut max_j = 0;

        for i in 0..n {
            for j in i + 1..n {
                let a_ij = csr_value(a, i, j);
                let a_ji = csr_value(a, j, i);
                let diff = NumericElement::abs(a_ij - a_ji);

                if diff > max_asymmetry {
                    max_asymmetry = diff;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        if max_asymmetry > tol {
            return Err(Error::InvalidInput(format!(
                "Matrix is not symmetric: max asymmetry {:.2e} at ({}, {})",
                NumericElement::to_f64(max_asymmetry),
                max_i,
                max_j
            )));
        }

        Ok(())
    }

    /// Perform IC(0) factorization
    fn factorize(a: &CsrMatrix<T>) -> Result<CsrMatrix<T>>
    where
        T: FloatElement,
    {
        let n = a.nrows();
        let mut l_vals = Vec::new();
        let mut l_indices = Vec::new();
        let mut l_offsets = vec![0];

        // Row buffer for the incomplete lower triangular factor.
        let mut l_rows: Vec<Vec<(usize, T)>> = vec![Vec::new(); n];

        // IC(0) factorization algorithm
        for i in 0..n {
            let mut l_ii = csr_value(a, i, i);

            // Subtract contributions from previous columns
            for &(k, l_ik) in &l_rows[i] {
                if k < i {
                    l_ii -= l_ik * l_ik;
                }
            }

            // Check for positive definiteness
            if l_ii <= <T as NumericElement>::ZERO {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            l_ii = NumericElement::sqrt(l_ii);
            l_rows[i].push((i, l_ii));

            // Compute entries in column i below diagonal
            let row = a.row(i);
            for &j in row.col_indices() {
                if j > i {
                    // Only process lower triangular part
                    let mut l_ji = csr_value(a, j, i);

                    // Subtract contributions from previous columns
                    for &(k, l_ik) in &l_rows[i] {
                        if k < i {
                            // Find l_jk in row j
                            if let Some(&(_, l_jk)) = l_rows[j].iter().find(|(col, _)| *col == k) {
                                l_ji -= l_ik * l_jk;
                            }
                        }
                    }

                    l_ji = l_ji / l_ii;
                    l_rows[j].push((i, l_ji));
                }
            }
        }

        // Convert to CSR format
        for row in &l_rows {
            let mut sorted_row = row.clone();
            sorted_row.sort_by_key(|(col, _)| *col);

            for &(col, val) in &sorted_row {
                l_indices.push(col);
                l_vals.push(val);
            }
            l_offsets.push(l_indices.len());
        }

        CsrMatrix::from_parts(l_vals, l_indices, l_offsets, n, n).map_err(|error| {
            Error::InvalidInput(format!(
                "Failed to create Leto CSR matrix from IC(0) factorization: {error}"
            ))
        })
    }

    /// Forward substitution: solve L*y = b
    fn forward_substitution(&self, b: &Array1<T>, y: &mut Array1<T>) {
        let n = self.l_factor.nrows();

        for i in 0..n {
            let mut sum = b[i];

            let row = self.l_factor.row(i);
            for (&j, &value) in row.col_indices().iter().zip(row.values()) {
                if j < i {
                    sum -= value * y[j];
                } else if j == i {
                    sum = sum / value;
                    break;
                }
            }

            y[i] = sum;
        }
    }

    /// Backward substitution: solve L^T*x = y
    fn backward_substitution(&self, y: &Array1<T>, x: &mut Array1<T>) {
        let n = self.l_factor.nrows();

        for i in (0..n).rev() {
            let mut sum = y[i];

            // Process contributions from later rows
            for j in i + 1..n {
                let row = self.l_factor.row(j);
                for (&col, &value) in row.col_indices().iter().zip(row.values()) {
                    if col == i {
                        sum -= value * x[j];
                        break;
                    }
                }
            }

            // Divide by diagonal element
            let row = self.l_factor.row(i);
            for (&col, &value) in row.col_indices().iter().zip(row.values()) {
                if col == i {
                    x[i] = sum / value;
                    break;
                }
            }
        }
    }
}

impl<T: RealField + Copy + FloatElement + LetoScalar> Preconditioner<T> for IncompleteCholesky<T> {
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> Result<()> {
        let n = self.l_factor.nrows();
        validate_vector_len("IncompleteCholesky residual", r, n)?;
        validate_vector_len("IncompleteCholesky output", z, n)?;

        // Allocate intermediate result for two-stage solve
        // NOTE: Cannot use pre-allocated workspace with immutable &self API
        let mut y = Array1::zeros([n]);
        let mut solution = Array1::zeros([n]);

        // Solve L*L^T*z = r via forward and backward substitution
        self.forward_substitution(r, &mut y);
        self.backward_substitution(&y, &mut solution);

        for idx in 0..vector_len(z) {
            z[idx] = solution[idx];
        }

        Ok(())
    }
}
