//! MPI-based distributed memory parallelization for CFD
//!
//! This module provides MPI-based domain decomposition and distributed computing
//! capabilities for large-scale CFD simulations. It enables scaling beyond single-node
//! shared memory limitations through:
//!
//! - Domain decomposition with load balancing
//! - Ghost cell communication patterns
//! - Parallel I/O and data distribution
//! - Scalable linear algebra operations
//!
//! ## System Requirements
//!
//! **MPI Library Installation Required:**
//!
//! ### Linux (Ubuntu/Debian):
//! ```bash
//! sudo apt-get install libopenmpi-dev openmpi-bin
//! ```
//!
//! ### macOS (with Homebrew):
//! ```bash
//! brew install open-mpi
//! ```
//!
//! ### Windows (Microsoft MPI):
//! Download and install MSMPI from Microsoft:
//! https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi
//!
//! ## Architecture
//!
//! The MPI implementation follows a hierarchical design:
//!
//! 1. **Communicator Management**: MPI universe and process groups
//! 2. **Domain Decomposition**: Grid partitioning and load balancing
//! 3. **Data Distribution**: Field data layout and ownership
//! 4. **Communication**: Ghost cell exchanges and reductions
//! 5. **Collective Operations**: Parallel I/O and global reductions
//!
//! ## Usage
//!
//! ### Basic MPI Setup and Communication
//! ```no_run
//! use cfd_core::compute::mpi::*;
//!
//! // Initialize MPI environment
//! let universe = MpiUniverse::new()?;
//! let world = universe.world();
//!
//! // Create domain decomposition
//! let global_extents = GlobalExtents::new_2d(100, 100, (0.0, 1.0, 0.0, 1.0));
//! let decomp = DomainDecomposition::new(global_extents, &world, DecompositionStrategy::Cartesian2D)?;
//!
//! // Create distributed grid with ghost cells
//! let dist_grid = DistributedGrid::new(decomp.local_subdomain().clone())?;
//!
//! // Initialize CFD field data (velocity components and pressure)
//! let mut velocity_u = vec![vec![nalgebra::Vector2::zeros(); dist_grid.ny_total()]; dist_grid.nx_total()];
//! let mut velocity_v = vec![vec![nalgebra::Vector2::zeros(); dist_grid.ny_total()]; dist_grid.nx_total()];
//! let mut pressure = vec![vec![0.0f64; dist_grid.ny_total()]; dist_grid.nx_total()];
//!
//! // Create ghost cell manager for communication
//! let ghost_manager = GhostCellManager::new(
//!     world,
//!     decomp.neighbors().clone(),
//!     1, // ghost layers
//! )?;
//!
//! // Main CFD time stepping loop
//! for time_step in 0..num_steps {
//!     // Update ghost cells before computation
//!     ghost_manager.update_ghost_cells(
//!         &mut velocity_u,
//!         &mut velocity_v,
//!         &mut pressure,
//!         decomp.local_subdomain(),
//!     )?;
//!
//!     // Perform CFD computation on owned cells only
//!     // (ghost cells provide boundary conditions from neighboring processes)
//!     compute_cfd_step(&mut velocity_u, &mut velocity_v, &mut pressure, &dist_grid)?;
//!
//!     // Optional: Use async communication for overlap with next computation
//!     // let requests = ghost_manager.update_ghost_cells_async(...)?;
//!     // compute_next_step(...);
//!     // ghost_manager.complete_async_updates(...)?;
//! }
//! ```
//!
//! ### Distributed Linear Solvers
//! ```no_run
//! use cfd_core::compute::mpi::*;
//!
//! // Create distributed Laplacian operator
//! let operator = DistributedLaplacian2D::new(&decomp, &world, dx, dy)?;
//!
//! // Setup parallel preconditioner
//! let preconditioner = BlockJacobiPreconditioner::new(&operator, &decomp, &world)?;
//!
//! // Create distributed GMRES solver
//! let solver = DistributedGMRES::new(operator, preconditioner, &world, krylov_dim);
//!
//! // Solve distributed system
//! let solution = solver.solve(&rhs, &initial_guess, tolerance, max_iter)?;
//! ```
//!
//! ### Parallel I/O Operations
//! ```no_run
//! use cfd_io::vtk::ParallelVtkWriter;
//! use cfd_io::hdf5_module::ParallelHdf5Writer;
//!
//! // Parallel VTK output
//! let vtk_writer = ParallelVtkWriter::new(&world)?;
//! vtk_writer.write_vtk_file("output.vtk", &points, &cells, &cell_types, &point_data, &cell_data)?;
//!
//! // Parallel HDF5 checkpointing
//! let hdf5_writer = ParallelHdf5Writer::new(&world)?;
//! hdf5_writer.write_checkpoint("checkpoint.h5", time_step, simulation_time, &datasets)?;
//! ```
//!
//! ### Dynamic Load Balancing
//! ```no_run
//! use cfd_core::compute::mpi::*;
//!
//! // Create load balancer
//! let mut load_balancer = LoadBalancer::new(&world, decomp, 1.2, 100)?;
//!
//! // Assess current load balance
//! let local_workload = compute_local_workload();
//! let metrics = load_balancer.assess_load_balance(local_workload)?;
//!
//! // Trigger repartitioning if needed
//! if load_balancer.should_repartition(&metrics) {
//!     let new_workloads = gather_workload_distribution();
//!     let new_decomp = load_balancer.repartition(&new_workloads)?;
//!     // Update local subdomain and redistribute data
//! }
//! ```
//!
//! ### Adaptive Mesh Refinement
//! ```no_run
//! use cfd_core::compute::mpi::*;
//!
//! // Setup refinement criteria
//! let criteria = RefinementCriteria {
//!     error_threshold: 1e-3,
//!     coarsening_threshold: 1e-5,
//!     max_refinement_ratio: 4,
//! };
//!
//! // Create AMR system with load balancing
//! let load_balancer = LoadBalancer::new(&world, decomp, 1.2, 100)?;
//! let mut amr = AdaptiveMeshRefinement::new(5, criteria, Some(load_balancer));
//!
//! // During simulation, adapt mesh based on error estimates
//! let error_estimates = compute_error_estimates();
//! let new_decomp = amr.adapt_mesh(&error_estimates, &current_decomp)?;
//!
//! if new_decomp.local_subdomain() != current_decomp.local_subdomain() {
//!     // Mesh was repartitioned - redistribute data
//!     redistribute_data_to_new_decomposition(&new_decomp)?;
//! }
//! ```
//!
//! ## Communication Patterns
//!
//! ### Blocking Communication (Synchronous)
//! ```no_run
//! # use cfd_core::compute::mpi::*;
//! # let ghost_manager = todo!();
//! # let mut velocity_u = todo!();
//! # let mut velocity_v = todo!();
//! # let mut pressure = todo!();
//! # let subdomain = todo!();
//! // Standard blocking exchange - waits for completion
//! ghost_manager.update_ghost_cells(
//!     &mut velocity_u,
//!     &mut velocity_v,
//!     &mut pressure,
//!     &subdomain,
//! )?;
//! ```
//!
//! ### Non-Blocking Communication (Asynchronous)
//! ```no_run
//! # use cfd_core::compute::mpi::*;
//! # let ghost_manager = todo!();
//! # let mut velocity_u = todo!();
//! # let mut velocity_v = todo!();
//! # let mut pressure = todo!();
//! # let subdomain = todo!();
//! // Start asynchronous communication
//! let requests = ghost_manager.update_ghost_cells_async(
//!     &mut velocity_u,
//!     &mut velocity_v,
//!     &mut pressure,
//!     &subdomain,
//! )?;
//!
//! // Perform computation while communication happens in background
//! perform_computation_on_interior_cells()?;
//!
//! // Wait for communication to complete
//! ghost_manager.complete_async_updates(
//!     &mut velocity_u,
//!     &mut velocity_v,
//!     &mut pressure,
//!     &subdomain,
//!     requests,
//! )?;
//! ```
//!
//! ### Performance Validation & Scaling Tests
//! ```no_run
//! use cfd_core::compute::mpi::performance_validation::*;
//!
//! // Create performance validator
//! let validator = PerformanceValidator::<f64>::new(&world)?;
//!
//! // Define test core counts for scaling analysis
//! let core_counts = vec![1, 2, 4, 8, 16, 32, 64];
//!
//! // Run strong scaling test (fixed problem size)
//! let strong_results = validator.run_strong_scaling_test(
//!     |comm, cores, tol| {
//!         // Run your CFD simulation here
//!         // Return performance metrics
//!         Ok(PerformanceMetrics::new())
//!     },
//!     &core_counts,
//!     1e-6
//! )?;
//!
//! // Run weak scaling test (problem size scales with cores)
//! let weak_results = validator.run_weak_scaling_test(
//!     |comm, cores, tol| {
//!         // Run your CFD simulation with scaled problem size
//!         // Return performance metrics
//!         Ok(PerformanceMetrics::new())
//!     },
//!     &core_counts,
//!     1e-6
//! )?;
//!
//! // Assess production readiness
//! let readiness_report = validator.assess_production_readiness()?;
//!
//! // Print scaling assessment
//! println!("Strong scaling grade: {:?}", strong_results.assessment.grade);
//! println!("Production readiness: {}%", readiness_report.overall_score);
//! ```
//!
//! ## Compilation
//!
//! Enable MPI support with the `mpi` feature:
//!
//! ```bash
//! cargo build --features mpi
//! ```
//!
//! Without MPI installed, the feature will fail to compile with clear error messages.

mod communicator;
mod decomposition;
mod distributed_grid;
mod distributed_solvers;
mod error;
mod ghost_cells;
#[cfg(feature = "mpi")]
mod performance_validation;
#[cfg(test)]
mod tests;

pub use communicator::*;
pub use decomposition::{AdaptiveMeshRefinement, LoadBalanceMetrics, LoadBalancer, RefinementCriteria};
pub use distributed_grid::*;
pub use distributed_solvers::*;
pub use error::*;
pub use ghost_cells::*;
#[cfg(feature = "mpi")]
pub use performance_validation::*;

/// MPI rank identifier (process ID within communicator)
pub type Rank = i32;

/// MPI process count
pub type Size = i32;

/// Global domain extents for distributed grids
#[derive(Debug, Clone, Copy)]
pub struct GlobalExtents {
    /// Global number of cells in x direction
    pub nx_global: usize,
    /// Global number of cells in y direction
    pub ny_global: usize,
    /// Global number of cells in z direction (for 3D)
    pub nz_global: usize,
    /// Global domain bounds
    pub bounds: (f64, f64, f64, f64, f64, f64), // (x_min, x_max, y_min, y_max, z_min, z_max)
}

impl GlobalExtents {
    /// Create 2D global extents
    pub fn new_2d(nx: usize, ny: usize, bounds: (f64, f64, f64, f64)) -> Self {
        Self {
            nx_global: nx,
            ny_global: ny,
            nz_global: 1,
            bounds: (bounds.0, bounds.1, bounds.2, bounds.3, 0.0, 0.0),
        }
    }

    /// Create 3D global extents
    pub fn new_3d(nx: usize, ny: usize, nz: usize, bounds: (f64, f64, f64, f64, f64, f64)) -> Self {
        Self {
            nx_global: nx,
            ny_global: ny,
            nz_global: nz,
            bounds,
        }
    }
}

/// Local subdomain information for each MPI rank
#[derive(Debug, Clone)]
pub struct LocalSubdomain {
    /// MPI rank owning this subdomain
    pub rank: Rank,
    /// Local grid dimensions (owned cells only, excluding ghost cells)
    pub nx_local: usize,
    pub ny_local: usize,
    pub nz_local: usize,
    /// Global indices of the first owned cell in this subdomain
    pub i_start_global: usize,
    pub j_start_global: usize,
    pub k_start_global: usize,
    /// Number of ghost cells in each direction
    pub ghost_layers: usize,
}

impl LocalSubdomain {
    /// Total local grid size including ghost cells
    pub fn total_nx(&self) -> usize {
        self.nx_local + 2 * self.ghost_layers
    }

    /// Total local grid size including ghost cells
    pub fn total_ny(&self) -> usize {
        self.ny_local + 2 * self.ghost_layers
    }

    /// Total local grid size including ghost cells
    pub fn total_nz(&self) -> usize {
        self.nz_local + 2 * self.ghost_layers
    }

    /// Check if local indices are within owned region (excluding ghost cells)
    pub fn is_owned_cell(&self, i: usize, j: usize, k: usize) -> bool {
        i >= self.ghost_layers &&
        i < self.ghost_layers + self.nx_local &&
        j >= self.ghost_layers &&
        j < self.ghost_layers + self.ny_local &&
        k >= self.ghost_layers &&
        k < self.ghost_layers + self.nz_local
    }

    /// Convert local index (including ghosts) to global index
    pub fn local_to_global(&self, i_local: usize, j_local: usize, k_local: usize) -> (usize, usize, usize) {
        let i_global = self.i_start_global + (i_local - self.ghost_layers);
        let j_global = self.j_start_global + (j_local - self.ghost_layers);
        let k_global = self.k_start_global + (k_local - self.ghost_layers);
        (i_global, j_global, k_global)
    }

    /// Check if this subdomain owns a global cell
    pub fn owns_global_cell(&self, i_global: usize, j_global: usize, k_global: usize) -> bool {
        i_global >= self.i_start_global &&
        i_global < self.i_start_global + self.nx_local &&
        j_global >= self.j_start_global &&
        j_global < self.j_start_global + self.ny_local &&
        k_global >= self.k_start_global &&
        k_global < self.k_start_global + self.nz_local
    }
}
