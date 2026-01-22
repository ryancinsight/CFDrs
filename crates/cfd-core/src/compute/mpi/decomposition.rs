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
        let neighbors =
            Self::compute_neighbors(&local_subdomain, &global_extents, size, strategy)?;

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
    fn decompose_1d(global: &GlobalExtents, rank: Rank, size: Size) -> MpiResult<LocalSubdomain> {
        if size <= 0 {
            return Err(MpiError::DecompositionError(
                "MPI size must be positive".to_string(),
            ));
        }
        if global.nx_global == 0 {
            return Err(MpiError::DecompositionError(
                "Global domain has zero extent in x".to_string(),
            ));
        }
        let weights = vec![global.ny_global.max(1) * global.nz_global.max(1); global.nx_global];
        let partition = Self::partition_1d_by_weights(global.nx_global, &weights, size)?;
        let (i_start, nx_local) = partition[rank as usize];

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

    fn decompose_bisection(
        global: &GlobalExtents,
        rank: Rank,
        size: Size,
    ) -> MpiResult<LocalSubdomain> {
        let weights = vec![1.0; size as usize];
        let partitions = Self::partition_rcb(global, &weights)?;
        partitions
            .into_iter()
            .find(|subdomain| subdomain.rank == rank)
            .ok_or_else(|| {
                MpiError::DecompositionError("Rank not found in bisection partition".to_string())
            })
    }

    fn decompose_metis(
        global: &GlobalExtents,
        rank: Rank,
        size: Size,
    ) -> MpiResult<LocalSubdomain> {
        if size <= 0 {
            return Err(MpiError::DecompositionError(
                "MPI size must be positive".to_string(),
            ));
        }
        let (px, py, pz) = Self::find_optimal_3d_grid(
            size as usize,
            global.nx_global,
            global.ny_global,
            global.nz_global,
        );
        if px * py * pz != size as usize {
            return Err(MpiError::DecompositionError(format!(
                "Cannot create {} x {} x {} grid for {} processes",
                px, py, pz, size
            )));
        }
        let rx = rank as usize % px;
        let ry = (rank as usize / px) % py;
        let rz = rank as usize / (px * py);

        let nx_local = Self::distribute_dimension(global.nx_global, px, rx);
        let ny_local = Self::distribute_dimension(global.ny_global, py, ry);
        let nz_local = Self::distribute_dimension(global.nz_global, pz, rz);

        let i_start = Self::compute_global_start(global.nx_global, px, rx);
        let j_start = Self::compute_global_start(global.ny_global, py, ry);
        let k_start = Self::compute_global_start(global.nz_global, pz, rz);

        Ok(LocalSubdomain {
            rank,
            nx_local,
            ny_local,
            nz_local,
            i_start_global: i_start,
            j_start_global: j_start,
            k_start_global: k_start,
            ghost_layers: 1,
        })
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

    fn find_optimal_3d_grid(
        num_procs: usize,
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> (usize, usize, usize) {
        if nz <= 1 {
            let (px, py) = Self::find_optimal_2d_grid(num_procs);
            return (px, py, 1);
        }
        let mut best = (1, 1, num_procs);
        let mut best_cost = f64::INFINITY;
        for px in 1..=num_procs {
            if num_procs % px != 0 {
                continue;
            }
            let rem = num_procs / px;
            for py in 1..=rem {
                if rem % py != 0 {
                    continue;
                }
                let pz = rem / py;
                let dx = nx as f64 / px as f64;
                let dy = ny as f64 / py as f64;
                let dz = nz as f64 / pz as f64;
                let interface_area = 2.0 * (dx * dy + dx * dz + dy * dz);
                let aspect = (dx / dy - 1.0).abs()
                    + (dx / dz - 1.0).abs()
                    + (dy / dz - 1.0).abs();
                let cost = interface_area * (1.0 + aspect);
                if cost < best_cost {
                    best_cost = cost;
                    best = (px, py, pz);
                }
            }
        }
        best
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
        strategy: DecompositionStrategy,
    ) -> MpiResult<HashMap<Rank, NeighborInfo>> {
        let subdomains = Self::build_subdomains(global, size, strategy)?;
        let local_i_end = local.i_start_global + local.nx_local - 1;
        let local_j_end = local.j_start_global + local.ny_local - 1;
        let local_k_end = local.k_start_global + local.nz_local - 1;

        let mut neighbors = HashMap::new();
        for other in &subdomains {
            if other.rank == local.rank {
                continue;
            }
            let other_i_end = other.i_start_global + other.nx_local - 1;
            let other_j_end = other.j_start_global + other.ny_local - 1;
            let other_k_end = other.k_start_global + other.nz_local - 1;

            let j_overlap =
                Self::ranges_overlap(local.j_start_global, local_j_end, other.j_start_global, other_j_end);
            let k_overlap =
                Self::ranges_overlap(local.k_start_global, local_k_end, other.k_start_global, other_k_end);

            if other_i_end + 1 == local.i_start_global && j_overlap && k_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Left,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            } else if local_i_end + 1 == other.i_start_global && j_overlap && k_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Right,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            }

            let i_overlap =
                Self::ranges_overlap(local.i_start_global, local_i_end, other.i_start_global, other_i_end);

            if other_j_end + 1 == local.j_start_global && i_overlap && k_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Bottom,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            } else if local_j_end + 1 == other.j_start_global && i_overlap && k_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Top,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            }

            if other_k_end + 1 == local.k_start_global && i_overlap && j_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Front,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            } else if local_k_end + 1 == other.k_start_global && i_overlap && j_overlap {
                neighbors.insert(
                    other.rank,
                    NeighborInfo {
                        direction: NeighborDirection::Back,
                        overlap: local.ghost_layers.min(other.ghost_layers),
                    },
                );
            }
        }

        Ok(neighbors)
    }

    fn ranges_overlap(a_start: usize, a_end: usize, b_start: usize, b_end: usize) -> bool {
        a_start <= b_end && b_start <= a_end
    }

    fn build_subdomains(
        global: &GlobalExtents,
        size: Size,
        strategy: DecompositionStrategy,
    ) -> MpiResult<Vec<LocalSubdomain>> {
        let mut subdomains = Vec::with_capacity(size as usize);
        for rank in 0..size {
            let subdomain = match strategy {
                DecompositionStrategy::Simple1D => Self::decompose_1d(global, rank, size)?,
                DecompositionStrategy::Cartesian2D => Self::decompose_2d(global, rank, size)?,
                DecompositionStrategy::RecursiveBisection => {
                    Self::decompose_bisection(global, rank, size)?
                }
                DecompositionStrategy::Metis => Self::decompose_metis(global, rank, size)?,
            };
            subdomains.push(subdomain);
        }
        Ok(subdomains)
    }

    fn partition_1d_by_weights(
        nx: usize,
        weights: &[usize],
        size: Size,
    ) -> MpiResult<Vec<(usize, usize)>> {
        if weights.len() != nx {
            return Err(MpiError::DecompositionError(
                "Weight vector length must match nx".to_string(),
            ));
        }
        let total_weight: usize = weights.iter().sum();
        if total_weight == 0 {
            return Err(MpiError::DecompositionError(
                "Total weight must be positive".to_string(),
            ));
        }
        let mut partitions = Vec::with_capacity(size as usize);
        let mut start = 0usize;
        let mut acc_weight = 0usize;
        let mut target_weight = total_weight as f64 / size as f64;
        let mut current_rank = 0usize;
        for i in 0..nx {
            acc_weight += weights[i];
            let should_split = acc_weight as f64 >= target_weight && current_rank + 1 < size as usize;
            if should_split {
                let length = i + 1 - start;
                partitions.push((start, length));
                start = i + 1;
                current_rank += 1;
                target_weight =
                    (total_weight as f64 * (current_rank as f64 + 1.0)) / size as f64;
            }
        }
        if start < nx {
            partitions.push((start, nx - start));
        }
        if partitions.len() != size as usize {
            return Err(MpiError::DecompositionError(
                "Failed to partition domain across ranks".to_string(),
            ));
        }
        Ok(partitions)
    }

    fn partition_rcb(global: &GlobalExtents, weights: &[f64]) -> MpiResult<Vec<LocalSubdomain>> {
        let size = weights.len();
        if size == 0 {
            return Err(MpiError::DecompositionError(
                "Weights must not be empty".to_string(),
            ));
        }
        if global.nx_global == 0 || global.ny_global == 0 || global.nz_global == 0 {
            return Err(MpiError::DecompositionError(
                "Global extents must be non-zero in all dimensions".to_string(),
            ));
        }
        let total_weight: f64 = weights.iter().sum();
        if total_weight <= 0.0 {
            return Err(MpiError::DecompositionError(
                "Total weight must be positive".to_string(),
            ));
        }
        let mut partitions = Vec::with_capacity(size);
        let ranks: Vec<Rank> = (0..size as i32).collect();
        let region = PartitionRegion {
            i_start: 0,
            j_start: 0,
            k_start: 0,
            nx: global.nx_global,
            ny: global.ny_global,
            nz: global.nz_global,
        };
        Self::partition_rcb_recursive(region, &ranks, weights, &mut partitions)?;
        if partitions.len() != size {
            return Err(MpiError::DecompositionError(
                "Recursive bisection did not produce all partitions".to_string(),
            ));
        }
        Ok(partitions)
    }

    fn partition_rcb_recursive(
        region: PartitionRegion,
        ranks: &[Rank],
        weights: &[f64],
        partitions: &mut Vec<LocalSubdomain>,
    ) -> MpiResult<()> {
        if ranks.len() != weights.len() {
            return Err(MpiError::DecompositionError(
                "Rank list and weights must match".to_string(),
            ));
        }
        if ranks.is_empty() {
            return Err(MpiError::DecompositionError(
                "Ranks must not be empty".to_string(),
            ));
        }
        if ranks.len() == 1 {
            partitions.push(LocalSubdomain {
                rank: ranks[0],
                nx_local: region.nx,
                ny_local: region.ny,
                nz_local: region.nz,
                i_start_global: region.i_start,
                j_start_global: region.j_start,
                k_start_global: region.k_start,
                ghost_layers: 1,
            });
            return Ok(());
        }
        let total_weight: f64 = weights.iter().sum();
        if total_weight <= 0.0 {
            return Err(MpiError::DecompositionError(
                "Total weight must be positive".to_string(),
            ));
        }
        let split_dim = region.select_split_dimension().ok_or_else(|| {
            MpiError::DecompositionError(
                "Cannot split region with no available dimension".to_string(),
            )
        })?;
        let mut cumulative = 0.0;
        let mut split_index = 0usize;
        for (i, weight) in weights.iter().enumerate() {
            cumulative += *weight;
            if cumulative >= total_weight / 2.0 {
                split_index = i + 1;
                break;
            }
        }
        if split_index == 0 || split_index >= ranks.len() {
            split_index = ranks.len() / 2;
        }
        let left_weights = &weights[..split_index];
        let right_weights = &weights[split_index..];
        let left_weight: f64 = left_weights.iter().sum();
        let right_weight: f64 = right_weights.iter().sum();
        if left_weight <= 0.0 || right_weight <= 0.0 {
            return Err(MpiError::DecompositionError(
                "Partition weights must be positive".to_string(),
            ));
        }

        let (left_region, right_region) =
            region.split(split_dim, left_weight, right_weight)?;

        Self::partition_rcb_recursive(left_region, &ranks[..split_index], left_weights, partitions)?;
        Self::partition_rcb_recursive(
            right_region,
            &ranks[split_index..],
            right_weights,
            partitions,
        )?;
        Ok(())
    }


#[derive(Clone, Copy)]
struct PartitionRegion {
    i_start: usize,
    j_start: usize,
    k_start: usize,
    nx: usize,
    ny: usize,
    nz: usize,
}

#[derive(Clone, Copy)]
enum PartitionDimension {
    X,
    Y,
    Z,
}

impl PartitionRegion {
    fn select_split_dimension(&self) -> Option<PartitionDimension> {
        let mut best = None;
        let mut best_len = 0usize;
        if self.nx > 1 && self.nx >= best_len {
            best = Some(PartitionDimension::X);
            best_len = self.nx;
        }
        if self.ny > 1 && self.ny > best_len {
            best = Some(PartitionDimension::Y);
            best_len = self.ny;
        }
        if self.nz > 1 && self.nz > best_len {
            best = Some(PartitionDimension::Z);
        }
        best
    }

    fn split(
        &self,
        dim: PartitionDimension,
        left_weight: f64,
        right_weight: f64,
    ) -> MpiResult<(PartitionRegion, PartitionRegion)> {
        let total_weight = left_weight + right_weight;
        if total_weight <= 0.0 {
            return Err(MpiError::DecompositionError(
                "Partition weights must be positive".to_string(),
            ));
        }
        let (len, min_len) = match dim {
            PartitionDimension::X => (self.nx, 1usize),
            PartitionDimension::Y => (self.ny, 1usize),
            PartitionDimension::Z => (self.nz, 1usize),
        };
        if len <= 1 {
            return Err(MpiError::DecompositionError(
                "Cannot split dimension with length <= 1".to_string(),
            ));
        }
        let mut left_len =
            ((left_weight / total_weight) * len as f64).round() as usize;
        if left_len < min_len {
            left_len = min_len;
        }
        if left_len >= len {
            left_len = len - 1;
        }
        let right_len = len - left_len;
        let left = match dim {
            PartitionDimension::X => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: left_len,
                ny: self.ny,
                nz: self.nz,
            },
            PartitionDimension::Y => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: self.nx,
                ny: left_len,
                nz: self.nz,
            },
            PartitionDimension::Z => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: self.nx,
                ny: self.ny,
                nz: left_len,
            },
        };
        let right = match dim {
            PartitionDimension::X => PartitionRegion {
                i_start: self.i_start + left_len,
                j_start: self.j_start,
                k_start: self.k_start,
                nx: right_len,
                ny: self.ny,
                nz: self.nz,
            },
            PartitionDimension::Y => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start + left_len,
                k_start: self.k_start,
                nx: self.nx,
                ny: right_len,
                nz: self.nz,
            },
            PartitionDimension::Z => PartitionRegion {
                i_start: self.i_start,
                j_start: self.j_start,
                k_start: self.k_start + left_len,
                nx: self.nx,
                ny: self.ny,
                nz: right_len,
            },
        };
        Ok((left, right))
    }
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
        let neighbors = Self::compute_neighbors(
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
        let local_volume = (local_sub.nx_local as u64)
            * (local_sub.ny_local as u64)
            * (local_sub.nz_local as u64);

        let mut all_volumes = vec![0u64; size];
        self.communicator.all_gather(&local_volume, &mut all_volumes);

        // Compute weights proportional to Volume / Load (inverse density)
        // to equalize load in the new partition.
        // weight = volume / (load + epsilon)
        let weights: Vec<f64> = workloads
            .iter()
            .zip(all_volumes.iter())
            .map(|(&load, &vol)| vol as f64 / (load as f64 + 1.0))
            .collect();

        let partitions =
            DomainDecomposition::partition_rcb(self.current_decomp.global_extents(), &weights)?;
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
