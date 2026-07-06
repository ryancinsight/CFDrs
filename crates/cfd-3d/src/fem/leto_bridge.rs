//! FEM-local conversion boundary for nalgebra work vectors and Leto-backed
//! sparse assembly.

use cfd_core::error::{Error, Result};
use cfd_math::linear_solver::{IterativeLinearSolver, Preconditioner};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use leto::{Array1, Storage};
use nalgebra::DVector;

use crate::scalar::Cfd3dScalar;

pub(super) fn vector_to_array<T: Cfd3dScalar>(
    vector: &DVector<T>,
    context: &str,
) -> Result<Array1<T>> {
    Array1::from_shape_vec([vector.len()], vector.as_slice().to_vec()).map_err(|error| {
        Error::InvalidConfiguration(format!("Invalid {context} Leto vector: {error}"))
    })
}

pub(super) fn array_to_vector<T: Cfd3dScalar>(array: Array1<T>) -> DVector<T> {
    DVector::from_vec(array.storage().as_slice().to_vec())
}

pub(super) fn iterative_solve<T, S, P>(
    solver: &S,
    matrix: &SparseMatrix<T>,
    rhs: &DVector<T>,
    solution: &mut DVector<T>,
    preconditioner: Option<&P>,
    context: &str,
) -> Result<(usize, Option<T>)>
where
    T: Cfd3dScalar,
    S: IterativeLinearSolver<T>,
    P: Preconditioner<T>,
{
    let rhs_array = vector_to_array(rhs, context)?;
    let mut solution_array = vector_to_array(solution, context)?;
    let result = solver.solve(matrix, &rhs_array, &mut solution_array, preconditioner);
    *solution = DVector::from_vec(solution_array.storage().as_slice().to_vec());
    let monitor = result?;
    Ok((monitor.iteration, monitor.residual_history.last().copied()))
}

pub(super) fn build_with_vector_rhs<T>(
    builder: SparseMatrixBuilder<T>,
    rhs: DVector<T>,
    context: &str,
) -> Result<(SparseMatrix<T>, DVector<T>)>
where
    T: Cfd3dScalar,
{
    let mut rhs_array = vector_to_array(&rhs, context)?;
    let matrix = builder.build_with_rhs(&mut rhs_array)?;
    Ok((matrix, array_to_vector(rhs_array)))
}
