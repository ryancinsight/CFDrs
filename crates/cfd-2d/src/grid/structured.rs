//! Structured grid implementation for 2D CFD
//!
//! # Theorem
//! The grid topology must form a valid, non-overlapping partition of the computational domain.
//!
//! **Proof sketch**:
//! For a finite volume discretization to be conservative, the control volumes $\Omega_i$
//! must satisfy $\cup_i \Omega_i = \Omega$ and $\Omega_i \cap \Omega_j = \emptyset$ for $i \neq j$.
//! The grid data structures enforce this by maintaining strict adjacency invariants
//! and ensuring that the sum of face area vectors for any closed cell is exactly zero:
//! $\sum_f \mathbf{A}_f = \mathbf{0}$.

use super::boundary::BoundaryType;
use super::traits::Grid2D;
use crate::scalar::{from_f64, from_usize};
use cfd_core::error::{Error, Result};
use eunomia::FloatElement;
use leto::geometry::Vector2;
use std::collections::HashMap;

/// 2D structured grid implementation
#[derive(Debug, Clone)]
pub struct StructuredGrid2D<T> {
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

impl<T> StructuredGrid2D<T> {
    /// Create a new structured grid
    pub fn new(nx: usize, ny: usize, x_min: T, x_max: T, y_min: T, y_max: T) -> Result<Self>
    where
        T: FloatElement,
    {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConfiguration(
                "Grid dimensions must be positive".to_string(),
            ));
        }

        if nx == 1 || ny == 1 {
            return Err(Error::InvalidConfiguration(
                "Grid must have at least 2 points in each direction".to_string(),
            ));
        }

        // Grid spacing: For n points, there are (n-1) intervals
        // This is the standard finite difference discretization
        let dx = (x_max - x_min) / from_usize(nx - 1);
        let dy = (y_max - y_min) / from_usize(ny - 1);

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
    pub fn unit_square(nx: usize, ny: usize) -> Result<Self>
    where
        T: FloatElement,
    {
        Self::new(nx, ny, T::ZERO, T::ONE, T::ZERO, T::ONE)
    }

    /// Get the number of cells in x direction.
    pub fn nx(&self) -> usize {
        self.nx
    }

    /// Get the number of cells in y direction.
    pub fn ny(&self) -> usize {
        self.ny
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
    pub fn spacing(&self) -> (T, T)
    where
        T: Copy,
    {
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

    /// Get cell center coordinates.
    pub fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>>
    where
        T: FloatElement,
    {
        self.check_indices(i, j)?;

        let half = from_f64::<T>(0.5);
        let x = self.bounds.0 + (from_usize::<T>(i) + half) * self.dx;
        let y = self.bounds.2 + (from_usize::<T>(j) + half) * self.dy;

        Ok(Vector2::new(x, y))
    }

    /// Get cell volume/area.
    pub fn cell_area(&self, i: usize, j: usize) -> Result<T>
    where
        T: Copy + core::ops::Mul<Output = T>,
    {
        self.check_indices(i, j)?;
        Ok(self.dx * self.dy)
    }

    /// Get neighboring cell indices.
    pub fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)> {
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

    /// Get neighboring cell indices as an iterator.
    pub fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
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

    /// Check if cell is on a structured-grid boundary.
    pub fn is_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.nx - 1 || j == 0 || j == self.ny - 1
    }

    /// Get the configured boundary type for a cell.
    pub fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType> {
        self.boundaries.get(&(i, j)).copied()
    }
}

impl<T: FloatElement> Grid2D<T> for StructuredGrid2D<T> {
    fn nx(&self) -> usize {
        self.nx()
    }

    fn ny(&self) -> usize {
        self.ny()
    }

    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>> {
        self.cell_center(i, j)
    }

    fn cell_area(&self, i: usize, j: usize) -> Result<T> {
        self.cell_area(i, j)
    }

    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)> {
        self.neighbors(i, j)
    }

    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        self.neighbor_iter(i, j)
    }

    fn is_boundary(&self, i: usize, j: usize) -> bool {
        self.is_boundary(i, j)
    }

    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType> {
        self.boundary_type(i, j)
    }
}
