//! Domain decomposition manager — struct and core impl methods.

use super::neighbors::{NeighborDirection, NeighborInfo};
use super::strategy::DecompositionStrategy;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::error::{MpiError, MpiResult};
use crate::compute::mpi::{GlobalExtents, LocalSubdomain, Rank, Size};
use std::collections::HashMap;

/// Domain decomposition manager
#[derive(Debug, Clone)]
pub struct DomainDecomposition {
    /// Global domain information
    pub(super) global_extents: GlobalExtents,
    /// Decomposition strategy
    pub(super) strategy: DecompositionStrategy,
    /// MPI communicator
    pub(super) communicator: MpiCommunicator,
    /// Local subdomain for this rank
    pub(super) local_subdomain: LocalSubdomain,
    /// Neighbor information for ghost cell exchange
    pub(super) neighbors: HashMap<Rank, NeighborInfo>,
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
        let neighbors = Self::compute_neighbors(&local_subdomain, &global_extents, size, strategy)?;

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
        let partition =
            super::partition::partition_1d_by_weights(global.nx_global, &weights, size)?;
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
        let partitions = super::partition::partition_rcb(global, &weights)?;
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
                let aspect = (dx / dy - 1.0).abs() + (dx / dz - 1.0).abs() + (dy / dz - 1.0).abs();
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
    pub(super) fn compute_neighbors(
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

            let j_overlap = Self::ranges_overlap(
                local.j_start_global,
                local_j_end,
                other.j_start_global,
                other_j_end,
            );
            let k_overlap = Self::ranges_overlap(
                local.k_start_global,
                local_k_end,
                other.k_start_global,
                other_k_end,
            );

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

            let i_overlap = Self::ranges_overlap(
                local.i_start_global,
                local_i_end,
                other.i_start_global,
                other_i_end,
            );

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
