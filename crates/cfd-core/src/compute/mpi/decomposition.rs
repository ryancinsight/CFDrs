//! Domain decomposition for MPI parallelization

use super::communicator::MpiCommunicator;
use super::error::{MpiError, MpiResult};
use super::{GlobalExtents, LocalSubdomain, Rank, Size};
use std::collections::HashMap;

/// Domain decomposition strategy
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DecompositionStrategy {
    /// Simple 1D decomposition along x-direction
    Simple1D,
    /// 2D Cartesian decomposition (x and y directions)
    Cartesian2D,
    /// Recursive coordinate bisection (for load balancing)
    RecursiveBisection,
    /// METIS-based graph partitioning (if available)
    Metis,
}

/// Domain decomposition manager
#[derive(Debug)]
pub struct DomainDecomposition {
    /// Global domain information
    global_extents: GlobalExtents,
    /// Decomposition strategy
    strategy: DecompositionStrategy,
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Local subdomain for this rank
    local_subdomain: LocalSubdomain,
    /// Neighbor information for ghost cell exchange
    neighbors: HashMap<Rank, NeighborInfo>,
}

impl DomainDecomposition {
    /// Create a new domain decomposition
    pub fn new(
        global_extents: GlobalExtents,
        communicator: &MpiCommunicator,
        strategy: DecompositionStrategy,
    ) -> MpiResult<Self> {
        let rank = communicator.rank();
        let size = communicator.size();

        if size == 1 {
            // Single process - no decomposition needed
            let local_subdomain = LocalSubdomain {
                rank,
                nx_local: global_extents.nx_global,
                ny_local: global_extents.ny_global,
                nz_local: global_extents.nz_global,
                i_start_global: 0,
                j_start_global: 0,
                k_start_global: 0,
                ghost_layers: 0,
            };

            return Ok(Self {
                global_extents,
                strategy,
                communicator: communicator.clone(),
                local_subdomain,
                neighbors: HashMap::new(),
            });
        }

        // Perform domain decomposition based on strategy
        let local_subdomain = match strategy {
            DecompositionStrategy::Simple1D => Self::decompose_1d(&global_extents, rank, size)?,
            DecompositionStrategy::Cartesian2D => Self::decompose_2d(&global_extents, rank, size)?,
            DecompositionStrategy::RecursiveBisection => {
                Self::decompose_bisection(&global_extents, rank, size)?
            }
            DecompositionStrategy::Metis => Self::decompose_metis(&global_extents, rank, size)?,
        };

        // Determine neighbor relationships
        let neighbors = Self::compute_neighbors(&local_subdomain, &global_extents, size)?;

        Ok(Self {
            global_extents,
            strategy,
            communicator: communicator.clone(),
            local_subdomain,
            neighbors,
        })
    }

    /// Get the local subdomain for this rank
    pub fn local_subdomain(&self) -> &LocalSubdomain {
        &self.local_subdomain
    }

    /// Get neighbor information
    pub fn neighbors(&self) -> &HashMap<Rank, NeighborInfo> {
        &self.neighbors
    }

    /// Get global extents
    pub fn global_extents(&self) -> &GlobalExtents {
        &self.global_extents
    }

    /// Simple 1D decomposition along x-direction
    // TODO: Implement load-balanced decomposition with consideration of cell weights/complexity
    fn decompose_1d(global: &GlobalExtents, rank: Rank, size: Size) -> MpiResult<LocalSubdomain> {
        let nx_per_process = global.nx_global / size as usize;
        let remainder = global.nx_global % size as usize;

        // Distribute remainder cells to first few processes
        let mut i_start = 0;
        let mut nx_local = nx_per_process;

        for r in 0..rank {
            let extra = if r < remainder as i32 { 1 } else { 0 };
            i_start += nx_per_process + extra;
        }

        if rank < remainder as i32 {
            nx_local += 1;
        }

        Ok(LocalSubdomain {
            rank,
            nx_local,
            ny_local: global.ny_global,
            nz_local: global.nz_global,
            i_start_global: i_start,
            j_start_global: 0,
            k_start_global: 0,
            ghost_layers: 1, // Single layer of ghost cells
        })
    }

    /// 2D Cartesian decomposition
    fn decompose_2d(global: &GlobalExtents, rank: Rank, size: Size) -> MpiResult<LocalSubdomain> {
        // Find optimal 2D processor grid
        let (px, py) = Self::find_optimal_2d_grid(size as usize);

        if px * py != size as usize {
            return Err(MpiError::DecompositionError(format!(
                "Cannot create {} x {} grid for {} processes",
                px, py, size
            )));
        }

        // Convert rank to 2D coordinates
        let rx = rank as usize % px;
        let ry = rank as usize / px;

        // Compute local dimensions
        let nx_local = Self::distribute_dimension(global.nx_global, px, rx);
        let ny_local = Self::distribute_dimension(global.ny_global, py, ry);

        // Compute global starting indices
        let i_start = Self::compute_global_start(global.nx_global, px, rx);
        let j_start = Self::compute_global_start(global.ny_global, py, ry);

        Ok(LocalSubdomain {
            rank,
            nx_local,
            ny_local,
            nz_local: global.nz_global,
            i_start_global: i_start,
            j_start_global: j_start,
            k_start_global: 0,
            ghost_layers: 1,
        })
    }

    /// TODO: Implement recursive bisection decomposition (Zoltan/ParMETIS integration).
    fn decompose_bisection(
        _global: &GlobalExtents,
        _rank: Rank,
        _size: Size,
    ) -> MpiResult<LocalSubdomain> {
        // TODO: Implement graph partitioning based decomposition; current path is unavailable.
        Err(MpiError::NotAvailable(
            "Recursive bisection not implemented".to_string(),
        ))
    }

    /// TODO: Implement METIS-based decomposition.
    fn decompose_metis(
        _global: &GlobalExtents,
        _rank: Rank,
        _size: Size,
    ) -> MpiResult<LocalSubdomain> {
        // TODO: Integrate METIS and implement mesh/graph partitioning based decomposition.
        Err(MpiError::NotAvailable(
            "METIS decomposition not implemented".to_string(),
        ))
    }

    /// Find optimal 2D processor grid
    fn find_optimal_2d_grid(num_procs: usize) -> (usize, usize) {
        let sqrt = (num_procs as f64).sqrt() as usize;
        for i in (1..=sqrt).rev() {
            if num_procs % i == 0 {
                return (i, num_procs / i);
            }
        }
        (1, num_procs)
    }

    /// Distribute cells among processes in one dimension
    fn distribute_dimension(total: usize, num_procs: usize, proc_id: usize) -> usize {
        let base = total / num_procs;
        let remainder = total % num_procs;

        base + if proc_id < remainder { 1 } else { 0 }
    }

    /// Compute global starting index for a process
    fn compute_global_start(total: usize, num_procs: usize, proc_id: usize) -> usize {
        let base = total / num_procs;
        let remainder = total % num_procs;

        let mut start = 0;
        for i in 0..proc_id {
            start += base + if i < remainder { 1 } else { 0 };
        }
        start
    }

    /// Compute neighbor relationships for ghost cell exchange
    fn compute_neighbors(
        local: &LocalSubdomain,
        global: &GlobalExtents,
        size: Size,
    ) -> MpiResult<HashMap<Rank, NeighborInfo>> {
        let mut neighbors = HashMap::new();

        // For 1D decomposition, only left/right neighbors
        let i_end_global = local.i_start_global + local.nx_local - 1;

        // Left neighbor
        if local.i_start_global > 0 {
            let left_rank = Self::find_rank_for_global_i(local.i_start_global - 1, global, size)?;
            neighbors.insert(
                left_rank,
                NeighborInfo {
                    direction: NeighborDirection::Left,
                    overlap: 1,
                },
            );
        }

        // Right neighbor
        if i_end_global < global.nx_global - 1 {
            let right_rank = Self::find_rank_for_global_i(i_end_global + 1, global, size)?;
            neighbors.insert(
                right_rank,
                NeighborInfo {
                    direction: NeighborDirection::Right,
                    overlap: 1,
                },
            );
        }

        // For 2D decomposition, would also check top/bottom neighbors
        // Implementation would be similar but more complex

        Ok(neighbors)
    }

    /// Find which rank owns a given global i-index
    fn find_rank_for_global_i(
        global_i: usize,
        global: &GlobalExtents,
        size: Size,
    ) -> MpiResult<Rank> {
        // For 1D decomposition
        let cells_per_proc = global.nx_global / size as usize;
        let remainder = global.nx_global % size as usize;

        let mut current_i = 0;
        for rank in 0..size {
            let extra = if rank < remainder as i32 { 1 } else { 0 };
            let proc_cells = cells_per_proc + extra;

            if global_i >= current_i && global_i < current_i + proc_cells {
                return Ok(rank);
            }
            current_i += proc_cells;
        }

        Err(MpiError::DecompositionError(format!(
            "Global index {} not found in decomposition",
            global_i
        )))
    }
}

/// Information about a neighboring process
#[derive(Debug, Clone)]
pub struct NeighborInfo {
    /// Direction of the neighbor
    pub direction: NeighborDirection,
    /// Number of overlapping cells (ghost cell layers)
    pub overlap: usize,
}

/// Direction of neighboring process
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NeighborDirection {
    /// Left neighbor (lower i indices)
    Left,
    /// Right neighbor (higher i indices)
    Right,
    /// Bottom neighbor (lower j indices)
    Bottom,
    /// Top neighbor (higher j indices)
    Top,
    /// Front neighbor (lower k indices)
    Front,
    /// Back neighbor (higher k indices)
    Back,
}

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
        let rank = self.communicator.rank() as usize;
        let size = self.communicator.size() as usize;

        // Gather workload from all processes
        let mut all_workloads = vec![0i32; size];
        let local_workload_i32 = local_workload as i32;

        // All-gather workloads
        for i in 0..size {
            if i == rank {
                all_workloads[i] = local_workload_i32;
            }
            // In real MPI, this would be MPI_Allgather
            // TODO: Replace simulated broadcast with proper MPI_Allgather implementation
            // DEPENDENCIES: Implement actual MPI communication primitives
            // BLOCKED BY: Limited MPI integration in current framework
            // PRIORITY: High - Essential for distributed memory parallelization
            // For now, simulate with broadcast
            let mut temp = local_workload_i32;
            if i == rank {
                self.communicator.broadcast(&mut temp);
                all_workloads[i] = temp;
            }
        }

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
        let neighbors =
            Self::compute_neighbors(&new_subdomain, self.current_decomp.global_extents(), size)?;

        // Create new decomposition
        let new_decomp = DomainDecomposition {
            global_extents: self.current_decomp.global_extents().clone(),
            strategy: self.current_decomp.strategy,
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
        target_work: usize,
        workloads: &[usize],
    ) -> MpiResult<LocalSubdomain> {
        // TODO: Implement a real repartitioner (e.g., recursive bisection / space-filling curve).
        let mut cumulative_work = 0;
        let mut start_idx = 0;

        for (i, &work) in workloads.iter().enumerate() {
            if cumulative_work + work > (rank as usize) * target_work {
                start_idx = i;
                break;
            }
            cumulative_work += work;
        }

        // TODO: Generalize subdomain boundary computation beyond 1D partitioning.
        let global = self.current_decomp.global_extents();
        let cells_per_work_unit = global.nx_global / workloads.len().max(1);

        let i_start = start_idx * cells_per_work_unit;
        let nx_local = (target_work * cells_per_work_unit).min(global.nx_global - i_start);

        Ok(LocalSubdomain {
            rank,
            nx_local,
            ny_local: global.ny_global,
            nz_local: global.nz_global,
            i_start_global: i_start,
            j_start_global: 0,
            k_start_global: 0,
            ghost_layers: 1,
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

/// Adaptive mesh refinement support
pub struct AdaptiveMeshRefinement {
    /// Current refinement level
    refinement_level: usize,
    /// Maximum allowed refinement level
    max_refinement_level: usize,
    /// Refinement criteria (error thresholds)
    refinement_criteria: RefinementCriteria,
    /// Load balancer for refined meshes
    load_balancer: Option<LoadBalancer>,
}

impl AdaptiveMeshRefinement {
    /// Create new AMR system
    pub fn new(
        max_refinement_level: usize,
        refinement_criteria: RefinementCriteria,
        load_balancer: Option<LoadBalancer>,
    ) -> Self {
        Self {
            refinement_level: 0,
            max_refinement_level,
            refinement_criteria,
            load_balancer,
        }
    }

    /// Check if cells need refinement based on error estimates
    pub fn needs_refinement(&self, error_estimates: &[f64]) -> Vec<bool> {
        error_estimates
            .iter()
            .map(|&error| {
                error > self.refinement_criteria.error_threshold
                    && self.refinement_level < self.max_refinement_level
            })
            .collect()
    }

    /// Check if cells need coarsening
    pub fn needs_coarsening(&self, error_estimates: &[f64]) -> Vec<bool> {
        error_estimates
            .iter()
            .map(|&error| {
                error < self.refinement_criteria.coarsening_threshold && self.refinement_level > 0
            })
            .collect()
    }

    /// Perform mesh adaptation
    pub fn adapt_mesh(
        &mut self,
        error_estimates: &[f64],
        current_decomp: &DomainDecomposition,
    ) -> MpiResult<DomainDecomposition> {
        let refine_flags = self.needs_refinement(error_estimates);
        let coarsen_flags = self.needs_coarsening(error_estimates);

        // Count cells needing adaptation
        let refine_count = refine_flags.iter().filter(|&&x| x).count();
        let coarsen_count = coarsen_flags.iter().filter(|&&x| x).count();

        // If significant adaptation needed, trigger load rebalancing
        if refine_count > error_estimates.len() / 10 || coarsen_count > error_estimates.len() / 10 {
            if let Some(ref mut balancer) = self.load_balancer {
                // Create workload estimate based on refinement
                let local_workload = error_estimates.len() + refine_count - coarsen_count;
                let metrics = balancer.assess_load_balance(local_workload)?;

                if balancer.should_repartition(&metrics) {
                    let new_workloads =
                        vec![local_workload; current_decomp.communicator.size() as usize];
                    return balancer.repartition(&new_workloads);
                }
            }
        }

        // Return current decomposition if no rebalancing needed
        Ok(current_decomp.clone())
    }

    /// Get current refinement level
    pub fn refinement_level(&self) -> usize {
        self.refinement_level
    }

    /// Increment refinement level
    pub fn increment_level(&mut self) {
        if self.refinement_level < self.max_refinement_level {
            self.refinement_level += 1;
        }
    }
}

/// Refinement criteria for adaptive mesh refinement
#[derive(Debug, Clone)]
pub struct RefinementCriteria {
    /// Error threshold above which cells are refined
    pub error_threshold: f64,
    /// Error threshold below which cells are coarsened
    pub coarsening_threshold: f64,
    /// Maximum refinement ratio between adjacent cells
    pub max_refinement_ratio: usize,
}
