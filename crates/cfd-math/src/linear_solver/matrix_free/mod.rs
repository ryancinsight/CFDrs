//! Matrix-free linear solvers for CFD applications.
//!
//! This module provides matrix-free implementations of iterative linear solvers
//! that avoid explicit matrix storage and assembly. Matrix-free methods are
//! particularly advantageous for CFD applications because:
//!
//! - **Memory Efficiency**: No need to store large sparse matrices (80-90% reduction)
//! - **Scalability**: Enables solution of very large problems limited only by physics
//! - **Performance**: Better cache locality and natural parallelization
//! - **Flexibility**: Easy integration with complex physics and boundary conditions
//!
//! ## Core Components
//!
//! - [`LinearOperator`]: Trait for matrix-vector products without matrix storage
//! - [`MatrixFreeCG`]: Conjugate Gradient solver using operator abstraction
//! - [`MatrixFreeGMRES`]: GMRES solver using operator abstraction
//!
//! ## Usage Example
//!
//! ```rust,ignore
//! // Define a linear operator (e.g., discretized Laplacian)
//! struct LaplacianOperator {
//!     nx: usize,
//!     ny: usize,
//!     dx: f64,
//!     dy: f64,
//! }
//!
//! impl LinearOperator<f64> for LaplacianOperator {
//!     fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
//!         // Compute y = A*x without storing A
//!         // ... implementation ...
//!         Ok(())
//!     }
//!
//!     fn size(&self) -> usize {
//!         self.nx * self.ny
//!     }
//! }
//!
//! // Solve using matrix-free CG
//! let operator = LaplacianOperator { nx: 100, ny: 100, dx: 0.01, dy: 0.01 };
//! let solver = MatrixFreeCG::new(config);
//! let mut x = vec![0.0; operator.size()];
//! let b = vec![1.0; operator.size()];
//!
//! solver.solve(&operator, &b, &mut x)?;
//! ```

mod operator;
mod cg;
mod gmres;
mod bicgstab;
mod cfd_operators;
mod traits;
mod gpu_compute;
#[cfg(feature = "mpi")]
mod parallel_solvers;

// GPU modules are conditionally compiled only when the 'gpu' feature is enabled
#[cfg(feature = "gpu")]
mod gpu_operators;
#[cfg(feature = "gpu")]
mod gpu_solvers;

pub use operator::{LinearOperator, IdentityOperator, ScaledOperator};
pub use cg::MatrixFreeCG;
pub use gmres::MatrixFreeGMRES;
pub use bicgstab::MatrixFreeBiCGSTAB;
pub use cfd_operators::{LaplacianOperator2D, MomentumOperator1D, PoissonOperator3D, MomentumOperator2D, EnergyOperator2D};
pub use traits::MatrixFreeSolver;

// Conditionally re-export GPU-related types when the 'gpu' feature is enabled
#[cfg(feature = "gpu")]
pub use gpu_solvers::{GpuMatrixFreeGMRES, GpuMatrixFreeBiCGSTAB};
#[cfg(feature = "gpu")]
pub use operator::GpuLinearOperator;
#[cfg(feature = "gpu")]
pub use gpu_operators::{GpuLaplacianOperator2D, GpuPoissonOperator3D, GpuMomentumOperator2D};

// Conditionally re-export MPI-related types when the 'mpi' feature is enabled
#[cfg(feature = "mpi")]
pub use parallel_solvers::{
    ParallelMatrixFreeBiCGSTAB,
    ParallelLoadBalancer,
    CommunicationOptimizer,
    CommunicationOverlap,
    LoadBalancingStrategy,
    LoadBalancingRecommendations,
    CommunicationOptimization,
};

#[cfg(test)]
mod tests;
