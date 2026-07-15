//! Sparse matrix operations and extensions

use crate::linear_solver::LinearOperator;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{
    spgemm, spmv_into as leto_spmv_into, CsrMatrix, RealScalar as LetoRealScalar,
    Scalar as LetoScalar,
};

impl<T: RealField + Copy + Send + Sync + LetoScalar> LinearOperator<T> for CsrMatrix<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        try_spmv(self, x, y)
    }

    fn size(&self) -> usize {
        self.nrows()
    }
}

/// Sparse matrix-vector multiplication (SpMV): y = A * x
///
/// Uses the Atlas-owned Leto CSR provider. Vectors are Leto arrays, so the
/// sparse operation surface does not expose legacy matrix or vector storage.
///
/// # Arguments
/// * `a` - Sparse matrix in CSR format
/// * `x` - Input vector (must have length = a.ncols())
/// * `y` - Output vector (must have length = a.nrows(), will be overwritten)
///
/// # Panics
/// Panics if vector dimensions don't match matrix dimensions
pub fn spmv<T>(a: &CsrMatrix<T>, x: &Array1<T>, y: &mut Array1<T>)
where
    T: Copy + LetoScalar,
{
    try_spmv(a, x, y).expect("invariant: sparse matrix-vector product inputs are valid");
}

fn map_leto_sparse_error(operation: &str, error: leto::LetoError) -> Error {
    Error::InvalidConfiguration(format!("Leto CSR {operation} failed: {error}"))
}

/// Fallible sparse matrix-vector multiplication (SpMV): `y = A * x`.
///
pub fn try_spmv<T>(a: &CsrMatrix<T>, x: &Array1<T>, y: &mut Array1<T>) -> Result<()>
where
    T: Copy + LetoScalar,
{
    try_leto_spmv(a, x, y)
}

/// Fallible Leto CSR sparse matrix-vector multiplication: `y = A * x`.
///
/// This is the Atlas-native sparse operator path. It lets iterative solvers
/// consume `leto_ops::CsrMatrix` directly through [`try_spmv`].
pub fn try_leto_spmv<T>(a: &CsrMatrix<T>, x: &Array1<T>, y: &mut Array1<T>) -> Result<()>
where
    T: Copy + LetoScalar,
{
    if x.shape()[0] != a.ncols() {
        return Err(Error::InvalidConfiguration(format!(
            "Input vector dimension mismatch for SpMV: matrix is {}x{}, vector length is {}",
            a.nrows(),
            a.ncols(),
            x.shape()[0]
        )));
    }
    if y.shape()[0] != a.nrows() {
        return Err(Error::InvalidConfiguration(format!(
            "Output vector dimension mismatch for SpMV: matrix is {}x{}, vector length is {}",
            a.nrows(),
            a.ncols(),
            y.shape()[0]
        )));
    }

    let mut output = vec![<T as NumericElement>::ZERO; a.nrows()];
    leto_spmv_into(a, &x.view(), &mut output)
        .map_err(|e| Error::InvalidConfiguration(format!("Leto SpMV failed: {e}")))?;
    for idx in 0..output.len() {
        y[idx] = output[idx];
    }
    Ok(())
}

/// Fallible sparse matrix-matrix multiplication (SpMM): `C = A * B`.
///
pub fn try_sparse_sparse_mul<T>(a: &CsrMatrix<T>, b: &CsrMatrix<T>) -> Result<CsrMatrix<T>>
where
    T: RealField + Copy + LetoScalar,
{
    if a.ncols() != b.nrows() {
        return Err(Error::InvalidConfiguration(format!(
            "Matrix dimension mismatch for multiplication: lhs is {}x{}, rhs is {}x{}",
            a.nrows(),
            a.ncols(),
            b.nrows(),
            b.ncols()
        )));
    }

    spgemm(a, b).map_err(|e| Error::InvalidConfiguration(format!("Leto SpGEMM failed: {e}")))
}

/// Fallible sparse matrix transpose.
///
pub fn try_sparse_transpose<T>(matrix: &CsrMatrix<T>) -> Result<CsrMatrix<T>>
where
    T: RealField + Copy + LetoScalar,
{
    Ok(matrix.transpose())
}

/// Sparse matrix-matrix multiplication (SpMM): C = A * B
///
/// Multiplies two CSR matrices and returns the result in CSR format.
/// Delegates to [`try_sparse_sparse_mul`] and panics on invalid inputs, matching
/// the pre-existing infallible API contract.
pub fn sparse_sparse_mul<T>(a: &CsrMatrix<T>, b: &CsrMatrix<T>) -> CsrMatrix<T>
where
    T: RealField + Copy + LetoScalar,
{
    try_sparse_sparse_mul(a, b).expect("invariant: sparse matrix product inputs are valid")
}

/// Extension trait for sparse matrix operations
pub trait SparseMatrixExt<T: RealField + Copy> {
    /// Extract diagonal elements
    fn diagonal(&self) -> Array1<T>;

    /// Set diagonal elements
    fn set_diagonal(&mut self, diag: &Array1<T>) -> Result<()>;

    /// Scale matrix by a scalar
    fn scale(&mut self, factor: T);

    /// Add identity matrix scaled by factor
    fn add_identity(&mut self, factor: T) -> Result<()>;

    /// Compute Frobenius norm
    fn frobenius_norm(&self) -> T
    where
        T: NumericElement;

    /// Compute condition number estimate
    fn condition_estimate(&self) -> Result<T>
    where
        T: FloatElement + LetoRealScalar;

    /// Check if matrix is diagonally dominant
    fn is_diagonally_dominant(&self) -> bool
    where
        T: NumericElement;

    /// Apply row scaling
    fn scale_rows(&mut self, scaling: &Array1<T>) -> Result<()>;

    /// Apply column scaling  
    fn scale_columns(&mut self, scaling: &Array1<T>) -> Result<()>;
}

impl<T: RealField + Copy + LetoScalar> SparseMatrixExt<T> for CsrMatrix<T> {
    fn diagonal(&self) -> Array1<T> {
        Array1::from_shape_vec([self.nrows().min(self.ncols())], self.diagonal())
            .expect("invariant: CSR diagonal length matches matrix diagonal shape")
    }

    fn set_diagonal(&mut self, diag: &Array1<T>) -> Result<()> {
        if diag.shape()[0] != self.nrows().min(self.ncols()) {
            return Err(Error::InvalidConfiguration(
                "Diagonal size mismatch".to_string(),
            ));
        }

        // CSR format limitation: Diagonal modification not supported due to
        // immutable sparse structure. Use set_diagonal during matrix construction.
        Err(Error::InvalidConfiguration(
            "Direct diagonal modification not supported for CSR format".to_string(),
        ))
    }

    fn scale(&mut self, factor: T) {
        self.scale_values(factor);
    }

    fn add_identity(&mut self, _factor: T) -> Result<()> {
        if self.nrows() != self.ncols() {
            return Err(Error::InvalidConfiguration(
                "Matrix must be square to add identity".to_string(),
            ));
        }

        // This operation requires rebuilding the matrix structure
        // For CSR format, this is complex and would require conversion
        Err(Error::InvalidConfiguration(
            "Adding identity not directly supported for CSR format".to_string(),
        ))
    }

    fn frobenius_norm(&self) -> T
    where
        T: NumericElement,
    {
        self.frobenius_norm()
    }

    fn condition_estimate(&self) -> Result<T>
    where
        T: FloatElement + LetoRealScalar,
    {
        self.condition_estimate()
            .map_err(|error| map_leto_sparse_error("condition estimate", error))
    }

    fn is_diagonally_dominant(&self) -> bool
    where
        T: NumericElement,
    {
        self.is_strictly_diagonally_dominant()
    }

    fn scale_rows(&mut self, scaling: &Array1<T>) -> Result<()> {
        let scaling_values: Vec<T> = (0..scaling.shape()[0]).map(|idx| scaling[idx]).collect();
        self.scale_rows(&scaling_values)
            .map_err(|error| map_leto_sparse_error("row scaling", error))?;
        Ok(())
    }

    fn scale_columns(&mut self, scaling: &Array1<T>) -> Result<()> {
        let scaling_values: Vec<T> = (0..scaling.shape()[0]).map(|idx| scaling[idx]).collect();
        self.scale_columns(&scaling_values)
            .map_err(|error| map_leto_sparse_error("column scaling", error))?;
        Ok(())
    }
}
