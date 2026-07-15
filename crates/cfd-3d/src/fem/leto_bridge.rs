//! FEM-local dense-vector helpers for Leto-backed sparse assembly.

use cfd_core::error::Result;
use cfd_math::linear_solver::{IterativeLinearSolver, Preconditioner};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use leto::Array1;

use crate::scalar::Cfd3dScalar;

pub(super) fn iterative_solve<T, S, P>(
    solver: &S,
    matrix: &SparseMatrix<T>,
    rhs: &Array1<T>,
    solution: &mut Array1<T>,
    preconditioner: Option<&P>,
    _context: &str,
) -> Result<(usize, Option<T>)>
where
    T: Cfd3dScalar,
    S: IterativeLinearSolver<T>,
    P: Preconditioner<T>,
{
    let monitor = solver.solve(matrix, rhs, solution, preconditioner)?;
    Ok((monitor.iteration, monitor.residual_history.last().copied()))
}

pub(super) fn build_with_vector_rhs<T>(
    builder: SparseMatrixBuilder<T>,
    mut rhs: Array1<T>,
    _context: &str,
) -> Result<(SparseMatrix<T>, Array1<T>)>
where
    T: Cfd3dScalar,
{
    let matrix = builder.build_with_rhs(&mut rhs)?;
    Ok((matrix, rhs))
}
