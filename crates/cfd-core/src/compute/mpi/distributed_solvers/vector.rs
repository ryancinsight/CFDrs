use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::decomposition::LocalSubdomain;
use crate::compute::mpi::error::MpiResult;
use crate::compute::mpi::ghost_cells::GhostCellManager;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// Distributed vector with MPI-aware data distribution
#[derive(Debug, Clone)]
pub struct DistributedVector<T: RealField> {
    /// Local data owned by this process (excluding ghost cells)
    pub local_data: DVector<T>,
    /// MPI communicator
    pub(super) communicator: MpiCommunicator,
    /// Domain decomposition information
    pub(super) subdomain: LocalSubdomain,
    /// Ghost cell manager for boundary exchange
    pub(super) ghost_manager: Option<GhostCellManager<T>>,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedVector<T> {
    /// Create a new distributed vector
    pub fn new(
        local_size: usize,
        communicator: &MpiCommunicator,
        subdomain: LocalSubdomain,
        ghost_manager: Option<GhostCellManager<T>>,
    ) -> Self {
        let local_data = DVector::zeros(local_size);
        Self {
            local_data,
            communicator: communicator.clone(),
            subdomain,
            ghost_manager,
        }
    }

    /// Create distributed vector from local data
    pub fn from_local_data(
        data: DVector<T>,
        communicator: &MpiCommunicator,
        subdomain: LocalSubdomain,
        ghost_manager: Option<GhostCellManager<T>>,
    ) -> Self {
        Self {
            local_data: data,
            communicator: communicator.clone(),
            subdomain,
            ghost_manager,
        }
    }

    /// Global dot product across all processes
    pub fn dot(&self, other: &DistributedVector<T>) -> MpiResult<T> {
        let local_dot = self.local_data.dot(&other.local_data);
        let mut global_dot = T::zero();
        self.communicator.all_reduce_sum(&mut global_dot, local_dot);
        Ok(global_dot)
    }

    /// Global L2 norm across all processes
    pub fn norm(&self) -> MpiResult<T> {
        let local_norm_sq = self.local_data.norm_squared();
        let mut global_norm_sq = T::zero();
        self.communicator
            .all_reduce_sum(&mut global_norm_sq, local_norm_sq);
        Ok(global_norm_sq.sqrt())
    }

    /// Scale vector by scalar
    pub fn scale(&mut self, alpha: T) {
        self.local_data.scale_mut(alpha);
    }

    /// Add scaled vector: y = y + alpha * x
    pub fn axpy(&mut self, alpha: T, x: &DistributedVector<T>) {
        self.local_data.axpy(alpha, &x.local_data, T::one());
    }

    /// Copy from another distributed vector
    pub fn copy_from(&mut self, other: &DistributedVector<T>) {
        self.local_data.copy_from(&other.local_data);
    }
}
