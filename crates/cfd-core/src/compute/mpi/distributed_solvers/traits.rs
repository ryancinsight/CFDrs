use super::vector::DistributedVector;
use crate::compute::mpi::error::MpiResult;
use nalgebra::{DVector, RealField};

/// Distributed linear operator trait for matrix-free operations
pub trait DistributedLinearOperator<T: RealField> {
    /// Apply the operator to a distributed vector
    fn apply(&self, x: &DistributedVector<T>, y: &mut DistributedVector<T>) -> MpiResult<()>;

    /// Get local dimension owned by this process
    fn local_dimension(&self) -> usize;

    /// Get global dimension across all processes
    fn global_dimension(&self) -> usize;

    /// Extract the diagonal of the local operator
    fn extract_diagonal(&self) -> DVector<T>;

    /// Assemble the local matrix (optional, for Schwarz/ILU)
    fn assemble_local_matrix(&self) -> Option<nalgebra::DMatrix<T>> {
        None
    }
}

/// Preconditioner trait for distributed solvers
pub trait Preconditioner<T: RealField> {
    /// Apply preconditioner: z = M^-1 * r
    fn apply(&self, r: &DistributedVector<T>, z: &mut DistributedVector<T>) -> MpiResult<()>;
}
