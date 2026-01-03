//! Matrix-free linear algebra and compute acceleration for CFD.
//!
//! This module provides the infrastructure for matrix-free computations,
//! including GPU compute contexts and specialized operators.
//!
//! The iterative solvers themselves are unified in the parent `linear_solver` module
//! and work with any implementation of the `LinearOperator` trait.

#[cfg(feature = "gpu")]
pub use cfd_core::compute::gpu::GpuContext;

#[cfg(feature = "mpi")]
pub mod parallel_solvers;

// Re-export core types from the unified locations
pub use crate::linear_solver::traits::LinearOperator;
#[cfg(feature = "gpu")]
pub use crate::linear_solver::traits::GpuLinearOperator;
pub use crate::linear_solver::operators::{
    IdentityOperator, ScaledOperator,
    LaplacianOperator2D, PoissonOperator3D,
    MomentumOperator1D, MomentumOperator2D, EnergyOperator2D
};

#[cfg(feature = "gpu")]
pub use crate::linear_solver::operators::gpu::{BoundaryType, GpuLaplacianOperator2D, DispatchMetrics};

// Conditionally re-export MPI-related types when the 'mpi' feature is enabled
#[cfg(feature = "mpi")]
pub use parallel_solvers::{
    CommunicationOptimization, CommunicationOptimizer, CommunicationOverlap,
    LoadBalancingRecommendations, LoadBalancingStrategy, ParallelLoadBalancer,
    ParallelMatrixFreeBiCGSTAB,
};

#[cfg(test)]
mod tests;
