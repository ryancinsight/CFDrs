//! FEM-local dense-vector helpers for Leto-backed sparse assembly.

use cfd_core::error::Result;
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use leto::Array1;

use cfd_core::CfdScalar;

pub(super) fn build_with_vector_rhs<T>(
    builder: SparseMatrixBuilder<T>,
    mut rhs: Array1<T>,
    _context: &str,
) -> Result<(SparseMatrix<T>, Array1<T>)>
where
    T: CfdScalar,
{
    let matrix = builder.build_with_rhs(&mut rhs)?;
    Ok((matrix, rhs))
}
