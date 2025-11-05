//! Distributed grid management for MPI parallelization

use super::communicator::MpiCommunicator;
use super::decomposition::{DomainDecomposition, NeighborDirection};
use super::error::{MpiError, MpiResult};
use super::ghost_cells::{GhostCellManager, GhostCellUpdate};
use super::{GlobalExtents, LocalSubdomain};
use crate::error::Result;
use nalgebra::{RealField, Vector2};
use std::collections::HashMap;

/// Distributed grid that manages local subdomain with ghost cells
#[derive(Debug)]
pub struct DistributedGrid<T: RealField + Copy> {
    /// Domain decomposition information
    decomposition: DomainDecomposition,
    /// Local grid dimensions including ghost cells
    local_nx: usize,
    local_ny: usize,
    local_nz: usize,
    /// Field data storage (includes ghost cells)
    velocity_u: Vec<Vec<Vector2<T>>>,
    velocity_v: Vec<Vec<Vector2<T>>>,
    pressure: Vec<Vec<T>>,
    /// Ghost cell manager
    ghost_manager: GhostCellManager<T>,
}

impl<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> DistributedGrid<T> {
    /// Create a new distributed grid from global grid information
    pub fn new(
        global_extents: GlobalExtents,
        communicator: &MpiCommunicator,
        decomposition_strategy: super::decomposition::DecompositionStrategy,
    ) -> MpiResult<Self> {
        let decomposition = DomainDecomposition::new(
            global_extents,
            communicator,
            decomposition_strategy,
        )?;

        let local_subdomain = decomposition.local_subdomain();
        let local_nx = local_subdomain.total_nx();
        let local_ny = local_subdomain.total_ny();
        let local_nz = local_subdomain.total_nz();

        // Initialize field data with ghost cells
        let velocity_u = vec![vec![Vector2::zeros(); local_ny]; local_nx];
        let velocity_v = vec![vec![Vector2::zeros(); local_ny]; local_nx];
        let pressure = vec![vec![T::zero(); local_ny]; local_nx];

        // Create ghost cell manager
        let ghost_manager = GhostCellManager::new(
            communicator.clone(),
            decomposition.neighbors().clone(),
            local_subdomain.ghost_layers,
        );

        Ok(Self {
            decomposition,
            local_nx,
            local_ny,
            local_nz,
            velocity_u,
            velocity_v,
            pressure,
            ghost_manager,
        })
    }

    /// Get local subdomain information
    pub fn local_subdomain(&self) -> &LocalSubdomain {
        self.decomposition.local_subdomain()
    }

    /// Get global extents
    pub fn global_extents(&self) -> &GlobalExtents {
        self.decomposition.global_extents()
    }

    /// Get velocity field at local indices (including ghost cells)
    pub fn velocity_at(&self, i: usize, j: usize) -> Vector2<T> {
        Vector2::new(self.velocity_u[i][j], self.velocity_v[i][j])
    }

    /// Set velocity field at local indices
    pub fn set_velocity_at(&mut self, i: usize, j: usize, velocity: Vector2<T>) {
        self.velocity_u[i][j] = velocity;
        self.velocity_v[i][j] = velocity;
    }

    /// Get pressure at local indices
    pub fn pressure_at(&self, i: usize, j: usize) -> T {
        self.pressure[i][j]
    }

    /// Set pressure at local indices
    pub fn set_pressure_at(&mut self, i: usize, j: usize, pressure: T) {
        self.pressure[i][j] = pressure;
    }

    /// Check if local indices are within owned region (not ghost cells)
    pub fn is_owned_cell(&self, i: usize, j: usize) -> bool {
        self.local_subdomain().is_owned_cell(i, j, 0)
    }

    /// Update ghost cells by communicating with neighboring processes
    pub fn update_ghost_cells(&mut self) -> MpiResult<()> {
        // Update velocity ghost cells
        self.ghost_manager.update_ghost_cells(
            &self.velocity_u,
            &self.velocity_v,
            &mut self.pressure,
            self.local_subdomain(),
        )
    }

    /// Apply boundary conditions to owned cells only
    pub fn apply_boundary_conditions(&mut self, bc_function: impl Fn(usize, usize) -> (Vector2<T>, T)) -> Result<()> {
        let subdomain = self.local_subdomain();

        // Apply BCs only to owned cells (not ghost cells)
        for i in subdomain.ghost_layers..subdomain.ghost_layers + subdomain.nx_local {
            for j in subdomain.ghost_layers..subdomain.ghost_layers + subdomain.ny_local {
                let (velocity, pressure) = bc_function(i, j);
                self.set_velocity_at(i, j, velocity);
                self.set_pressure_at(i, j, pressure);
            }
        }

        Ok(())
    }

    /// Iterate over owned cells (excluding ghost cells)
    pub fn owned_cells_iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        let subdomain = self.local_subdomain();
        let i_start = subdomain.ghost_layers;
        let i_end = i_start + subdomain.nx_local;
        let j_start = subdomain.ghost_layers;
        let j_end = j_start + subdomain.ny_local;

        (i_start..i_end).flat_map(move |i| (j_start..j_end).map(move |j| (i, j)))
    }

    /// Get total local grid size (including ghost cells)
    pub fn local_size(&self) -> (usize, usize, usize) {
        (self.local_nx, self.local_ny, self.local_nz)
    }

    /// Global reduction operations
    pub fn global_reduce_sum(&self, local_value: T) -> T {
        let mut result = local_value;
        self.decomposition.communicator.all_reduce_sum(&mut result);
        result
    }

    pub fn global_reduce_max(&self, local_value: T) -> T {
        let mut result = local_value;
        self.decomposition.communicator.all_reduce_max(&mut result);
        result
    }

    pub fn global_reduce_min(&self, local_value: T) -> T {
        let mut result = local_value;
        self.decomposition.communicator.all_reduce_min(&mut result);
        result
    }

    /// Gather statistics across all processes
    pub fn gather_statistics(&self) -> GridStatistics<T> {
        let local_owned = (self.local_subdomain().nx_local * self.local_subdomain().ny_local) as f64;
        let local_total = (self.local_nx * self.local_ny) as f64;

        let global_owned = self.global_reduce_sum(T::from_f64(local_owned).unwrap());
        let global_total = self.global_reduce_sum(T::from_f64(local_total).unwrap());

        GridStatistics {
            local_owned_cells: local_owned as usize,
            local_total_cells: local_total as usize,
            global_owned_cells: global_owned.to_f64().unwrap() as usize,
            global_total_cells: global_total.to_f64().unwrap() as usize,
            num_processes: self.decomposition.communicator.size() as usize,
            load_imbalance: Self::compute_load_imbalance(
                local_owned as usize,
                global_owned.to_f64().unwrap() as usize / self.decomposition.communicator.size() as usize
            ),
        }
    }

    /// Compute load imbalance metric
    fn compute_load_imbalance(actual: usize, ideal: usize) -> f64 {
        if ideal == 0 {
            0.0
        } else {
            (actual as f64 - ideal as f64).abs() / ideal as f64
        }
    }
}

/// Statistics about the distributed grid
#[derive(Debug, Clone)]
pub struct GridStatistics<T: RealField> {
    /// Number of owned cells on this process
    pub local_owned_cells: usize,
    /// Total cells including ghosts on this process
    pub local_total_cells: usize,
    /// Total owned cells across all processes
    pub global_owned_cells: usize,
    /// Total cells including ghosts across all processes
    pub global_total_cells: usize,
    /// Number of MPI processes
    pub num_processes: usize,
    /// Load imbalance metric (0.0 = perfect balance)
    pub load_imbalance: T,
}

impl<T: RealField + std::fmt::Display> std::fmt::Display for GridStatistics<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Distributed Grid Statistics:")?;
        writeln!(f, "  Local: {} owned, {} total cells", self.local_owned_cells, self.local_total_cells)?;
        writeln!(f, "  Global: {} owned, {} total cells", self.global_owned_cells, self.global_total_cells)?;
        writeln!(f, "  Processes: {}", self.num_processes)?;
        writeln!(f, "  Load imbalance: {:.3}", self.load_imbalance)?;
        Ok(())
    }
}
