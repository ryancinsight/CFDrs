//! Load balancing metrics and dynamic repartitioning.

use super::domain::DomainDecomposition;
use super::strategy::DecompositionStrategy;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::error::{MpiError, MpiResult};
use crate::compute::mpi::LocalSubdomain;

/// Load balancing metrics and algorithms
pub struct LoadBalancer {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Current domain decomposition
    current_decomp: DomainDecomposition,
    /// Load imbalance threshold for triggering repartitioning
    imbalance_threshold: f64,
    /// Minimum cells per process to avoid excessive overhead
    min_cells_per_process: usize,
}

impl LoadBalancer {
    /// Create new load balancer
    pub fn new(
        communicator: &MpiCommunicator,
        initial_decomp: DomainDecomposition,
        imbalance_threshold: f64,
        min_cells_per_process: usize,
    ) -> Self {
        Self {
            communicator: communicator.clone(),
            current_decomp: initial_decomp,
            imbalance_threshold,
            min_cells_per_process,
        }
    }

    /// Assess current load balance across processes
    pub fn assess_load_balance(&self, local_workload: usize) -> MpiResult<LoadBalanceMetrics> {
        let size = self.communicator.size() as usize;

        // Gather workload from all processes
        let mut all_workloads = vec![0i32; size];
        let local_workload_i32 = local_workload as i32;

        // All-gather workloads
        self.communicator
            .all_gather(&local_workload_i32, &mut all_workloads);

        // Calculate load balance metrics
        let workloads: Vec<usize> = all_workloads.iter().map(|&x| x as usize).collect();
        let total_work: usize = workloads.iter().sum();
        let avg_work = total_work as f64 / size as f64;
        let max_work = *workloads.iter().max().unwrap() as f64;
        let min_work = *workloads.iter().min().unwrap() as f64;

        let imbalance_ratio = max_work / avg_work;
        let efficiency = avg_work / max_work;

        Ok(LoadBalanceMetrics {
            max_load: max_work as usize,
            min_load: min_work as usize,
            avg_load: avg_work,
            imbalance_ratio,
            efficiency,
            needs_rebalancing: imbalance_ratio > self.imbalance_threshold,
        })
    }

    /// Check if repartitioning is needed based on current metrics
    pub fn should_repartition(&self, metrics: &LoadBalanceMetrics) -> bool {
        metrics.needs_rebalancing
    }

    /// Perform dynamic repartitioning
    pub fn repartition(&mut self, new_workloads: &[usize]) -> MpiResult<DomainDecomposition> {
        let rank = self.communicator.rank();
        let size = self.communicator.size();

        // Calculate optimal new distribution
        let total_work: usize = new_workloads.iter().sum();
        let target_work_per_proc = total_work / size as usize;

        // Create new subdomain assignment
        let new_subdomain =
            self.compute_new_subdomain(rank, target_work_per_proc, new_workloads)?;

        // Update neighbor relationships
        let neighbors = DomainDecomposition::compute_neighbors(
            &new_subdomain,
            self.current_decomp.global_extents(),
            size,
            DecompositionStrategy::RecursiveBisection,
        )?;

        // Create new decomposition
        let new_decomp = DomainDecomposition {
            global_extents: self.current_decomp.global_extents().clone(),
            strategy: DecompositionStrategy::RecursiveBisection,
            communicator: self.communicator.clone(),
            local_subdomain: new_subdomain,
            neighbors,
        };

        self.current_decomp = new_decomp.clone();
        Ok(new_decomp)
    }

    /// Compute new subdomain for this process
    fn compute_new_subdomain(
        &self,
        rank: i32,
        _target_work: usize,
        workloads: &[usize],
    ) -> MpiResult<LocalSubdomain> {
        let size = self.communicator.size() as usize;
        if workloads.len() != size {
            return Err(MpiError::DecompositionError(
                "Workloads length must match number of ranks".to_string(),
            ));
        }

        // Gather current volumes from all ranks
        let local_sub = &self.current_decomp.local_subdomain;
        let local_volume =
            (local_sub.nx_local as u64) * (local_sub.ny_local as u64) * (local_sub.nz_local as u64);

        let mut all_volumes = vec![0u64; size];
        self.communicator
            .all_gather(&local_volume, &mut all_volumes);

        // Compute weights proportional to Volume / Load (inverse density)
        // to equalize load in the new partition.
        // weight = volume / (load + epsilon)
        let weights: Vec<f64> = workloads
            .iter()
            .zip(all_volumes.iter())
            .map(|(&load, &vol)| vol as f64 / (load as f64 + 1.0))
            .collect();

        let partitions =
            super::partition::partition_rcb(self.current_decomp.global_extents(), &weights)?;
        let min_cells = self.min_cells_per_process;
        for subdomain in &partitions {
            let cell_count = subdomain.nx_local * subdomain.ny_local * subdomain.nz_local;
            if cell_count < min_cells {
                return Err(MpiError::DecompositionError(
                    "Repartitioning produced subdomain below minimum cell threshold".to_string(),
                ));
            }
        }
        partitions
            .into_iter()
            .find(|subdomain| subdomain.rank == rank)
            .ok_or_else(|| {
                MpiError::DecompositionError("Rank not found in repartitioned domain".to_string())
            })
    }

    /// Get current decomposition
    pub fn current_decomposition(&self) -> &DomainDecomposition {
        &self.current_decomp
    }
}

/// Load balance assessment metrics
#[derive(Debug, Clone)]
pub struct LoadBalanceMetrics {
    /// Maximum load across all processes
    pub max_load: usize,
    /// Minimum load across all processes
    pub min_load: usize,
    /// Average load per process
    pub avg_load: f64,
    /// Load imbalance ratio (max/avg)
    pub imbalance_ratio: f64,
    /// Parallel efficiency (avg/max)
    pub efficiency: f64,
    /// Whether rebalancing is recommended
    pub needs_rebalancing: bool,
}
