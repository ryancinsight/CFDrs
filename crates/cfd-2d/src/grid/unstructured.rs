//! Unstructured grid implementation for 2D CFD

use super::boundary::BoundaryType;
use super::traits::Grid2D;
use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector2};
use std::collections::HashMap;

/// 2D unstructured grid implementation
#[derive(Debug, Clone)]
pub struct UnstructuredGrid2D<T: RealField + Copy> {
    /// Cell centers
    pub centers: Vec<Vector2<T>>,
    /// Cell areas
    pub areas: Vec<T>,
    /// Cell connectivity (adjacency list)
    pub connectivity: Vec<Vec<usize>>,
    /// Boundary conditions
    pub boundaries: HashMap<usize, BoundaryType>,
    /// Grid dimensions (for compatibility)
    nx: usize,
    ny: usize,
}

impl<T: RealField + Copy> UnstructuredGrid2D<T> {
    /// Create a new unstructured grid
    pub fn new(
        centers: Vec<Vector2<T>>,
        areas: Vec<T>,
        connectivity: Vec<Vec<usize>>,
    ) -> Result<Self> {
        if centers.len() != areas.len() || centers.len() != connectivity.len() {
            return Err(Error::InvalidConfiguration(
                "Inconsistent grid data sizes".to_string(),
            ));
        }

        // Estimate nx, ny for compatibility
        let n = centers.len();
        let sqrt_n = (n as f64).sqrt() as usize;

        Ok(Self {
            centers,
            areas,
            connectivity,
            boundaries: HashMap::new(),
            nx: sqrt_n,
            ny: sqrt_n,
        })
    }

    /// Set boundary type for a cell
    pub fn set_boundary(&mut self, cell_id: usize, boundary_type: BoundaryType) {
        self.boundaries.insert(cell_id, boundary_type);
    }

    /// Convert 2D indices to cell ID
    fn cell_id(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Check if cell ID is valid
    fn check_cell_id(&self, id: usize) -> Result<()> {
        if id >= self.centers.len() {
            return Err(Error::InvalidConfiguration(format!(
                "Cell ID {} out of bounds for grid with {} cells",
                id,
                self.centers.len()
            )));
        }
        Ok(())
    }
}

impl<T: RealField + Copy> Grid2D<T> for UnstructuredGrid2D<T> {
    fn nx(&self) -> usize {
        self.nx
    }

    fn ny(&self) -> usize {
        self.ny
    }

    fn num_cells(&self) -> usize {
        self.centers.len()
    }

    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>> {
        let id = self.cell_id(i, j);
        self.check_cell_id(id)?;
        Ok(self.centers[id])
    }

    fn cell_area(&self, i: usize, j: usize) -> Result<T> {
        let id = self.cell_id(i, j);
        self.check_cell_id(id)?;
        Ok(self.areas[id])
    }

    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)> {
        let id = self.cell_id(i, j);
        if id >= self.connectivity.len() {
            return Vec::new();
        }

        self.connectivity[id]
            .iter()
            .map(|&neighbor_id| {
                let ni = neighbor_id % self.nx;
                let nj = neighbor_id / self.nx;
                (ni, nj)
            })
            .collect()
    }

    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        let id = self.cell_id(i, j);
        let nx = self.nx;

        self.connectivity
            .get(id)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .map(move |neighbor_id| {
                let ni = neighbor_id % nx;
                let nj = neighbor_id / nx;
                (ni, nj)
            })
    }

    fn is_boundary(&self, i: usize, j: usize) -> bool {
        let id = self.cell_id(i, j);
        self.boundaries.contains_key(&id)
    }

    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType> {
        let id = self.cell_id(i, j);
        self.boundaries.get(&id).copied()
    }
}
