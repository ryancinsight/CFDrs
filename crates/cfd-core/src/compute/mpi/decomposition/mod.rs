//! Domain decomposition for MPI parallelization.
//!
//! Manages spatial partitioning of simulation domains across MPI ranks,
//! including 1D/2D/3D Cartesian decomposition, recursive coordinate
//! bisection (RCB), load balancing, and adaptive mesh refinement support.

/// Adaptive mesh refinement support for distributed grids.
mod amr;
/// Domain decomposition manager and core algorithms.
mod domain;
/// Load balancing metrics and dynamic repartitioning.
mod load_balancer;
/// Neighbor information for ghost cell exchange.
mod neighbors;
/// Recursive coordinate bisection and weighted partitioning.
mod partition;
/// Domain decomposition strategy enumeration.
mod strategy;

pub use amr::{AdaptiveMeshRefinement, RefinementCriteria};
pub use domain::DomainDecomposition;
pub use load_balancer::{LoadBalanceMetrics, LoadBalancer};
pub use neighbors::{NeighborDirection, NeighborInfo};
pub use strategy::DecompositionStrategy;
