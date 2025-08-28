//! Structured grid implementation for 2D CFD

use super::boundary::BoundaryType;
use super::traits::Grid2D;
use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// 2D structured grid implementation
#[derive(Debug, Clone)]
pub struct StructuredGrid2D<T: RealField + Copy> {
    /// Number of cells in x direction
    pub nx: usize,
    /// Number of cells in y direction
    pub ny: usize,
    /// Domain bounds: (`x_min`, `x_max`, `y_min`, `y_max`)
    pub bounds: (T, T, T, T),
    /// Cell spacing in x direction
    pub dx: T,
    /// Cell spacing in y direction
    pub dy: T,
    /// Boundary conditions
    pub boundaries: HashMap<(usize, usize), BoundaryType>,
}

impl<T: RealField + FromPrimitive + Copy> StructuredGrid2D<T> {
    /// Create a new structured grid
    pub fn new(nx: usize, ny: usize, x_min: T, x_max: T, y_min: T, y_max: T) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConfiguration(
                "Grid dimensions must be positive".to_string(),
            ));
        }

        let dx = (x_max - x_min) / T::from_usize(nx).unwrap_or_else(|| T::zero());
        let dy = (y_max - y_min) / T::from_usize(ny).unwrap_or_else(|| T::zero());

        Ok(Self {
            nx,
            ny,
            bounds: (x_min, x_max, y_min, y_max),
            dx,
            dy,
            boundaries: HashMap::new(),
        })
    }

    /// Create a unit square grid
    pub fn unit_square(nx: usize, ny: usize) -> Result<Self> {
        Self::new(nx, ny, T::zero(), T::one(), T::zero(), T::one())
    }

    /// Set boundary type for a cell
    pub fn set_boundary(&mut self, i: usize, j: usize, boundary_type: BoundaryType) {
        self.boundaries.insert((i, j), boundary_type);
    }

    /// Set boundaries for all edges
    pub fn set_edge_boundaries(&mut self, boundary_type: BoundaryType) {
        // Bottom edge
        for i in 0..self.nx {
            self.set_boundary(i, 0, boundary_type);
        }
        // Top edge
        for i in 0..self.nx {
            self.set_boundary(i, self.ny - 1, boundary_type);
        }
        // Left edge
        for j in 0..self.ny {
            self.set_boundary(0, j, boundary_type);
        }
        // Right edge
        for j in 0..self.ny {
            self.set_boundary(self.nx - 1, j, boundary_type);
        }
    }

    /// Get grid spacing
    pub fn spacing(&self) -> (T, T) {
        (self.dx, self.dy)
    }

    /// Iterate over all grid cells
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.nx).flat_map(move |i| (0..self.ny).map(move |j| (i, j)))
    }

    /// Check if indices are valid
    fn check_indices(&self, i: usize, j: usize) -> Result<()> {
        if i >= self.nx || j >= self.ny {
            return Err(Error::InvalidConfiguration(format!(
                "Cell indices ({}, {}) out of bounds for grid {}x{}",
                i, j, self.nx, self.ny
            )));
        }
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + Copy> Grid2D<T> for StructuredGrid2D<T> {
    fn nx(&self) -> usize {
        self.nx
    }

    fn ny(&self) -> usize {
        self.ny
    }

    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>> {
        self.check_indices(i, j)?;

        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        let x = self.bounds.0 + (T::from_usize(i).unwrap_or_else(|| T::zero()) + half) * self.dx;
        let y = self.bounds.2 + (T::from_usize(j).unwrap_or_else(|| T::zero()) + half) * self.dy;

        Ok(Vector2::new(x, y))
    }

    fn cell_area(&self, i: usize, j: usize) -> Result<T> {
        self.check_indices(i, j)?;
        Ok(self.dx * self.dy)
    }

    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)> {
        let mut neighbors = Vec::with_capacity(4);

        if i > 0 {
            neighbors.push((i - 1, j));
        }
        if i < self.nx - 1 {
            neighbors.push((i + 1, j));
        }
        if j > 0 {
            neighbors.push((i, j - 1));
        }
        if j < self.ny - 1 {
            neighbors.push((i, j + 1));
        }

        neighbors
    }

    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        let nx = self.nx;
        let ny = self.ny;

        [
            (i.wrapping_sub(1), j),
            (i + 1, j),
            (i, j.wrapping_sub(1)),
            (i, j + 1),
        ]
        .into_iter()
        .filter(move |(ii, jj)| *ii < nx && *jj < ny && (*ii != i || *jj != j))
    }

    fn is_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.nx - 1 || j == 0 || j == self.ny - 1
    }

    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType> {
        self.boundaries.get(&(i, j)).copied()
    }
}
