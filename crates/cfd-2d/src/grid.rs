//! 2D grid structures for CFD simulations.
//!
//! This module provides various grid types for 2D CFD simulations, including
//! structured and unstructured grids with support for adaptive refinement.

use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Trait for 2D computational grids
pub trait Grid2D<T: RealField> {
    /// Get the number of cells in x direction
    fn nx(&self) -> usize;

    /// Get the number of cells in y direction
    fn ny(&self) -> usize;

    /// Get total number of cells
    fn num_cells(&self) -> usize {
        self.nx() * self.ny()
    }

    /// Get cell center coordinates
    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>>;

    /// Get cell volume/area
    fn cell_area(&self, i: usize, j: usize) -> Result<T>;

    /// Get neighboring cell indices (prefer neighbor_iter for better performance)
    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)>;

    /// Get neighboring cell indices as iterator (more efficient)
    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)>;

    /// Check if cell is on boundary
    fn is_boundary(&self, i: usize, j: usize) -> bool;

    /// Get boundary type for cell
    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType>;
}

/// Boundary types for grid cells
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryType {
    /// Wall boundary (no-slip)
    Wall,
    /// Inlet boundary
    Inlet,
    /// Outlet boundary
    Outlet,
    /// Symmetry boundary
    Symmetry,
    /// Periodic boundary
    Periodic,
}

/// 2D structured grid implementation
#[derive(Debug, Clone)]
pub struct StructuredGrid2D<T: RealField> {
    /// Number of cells in x direction
    pub nx: usize,
    /// Number of cells in y direction
    pub ny: usize,
    /// Domain bounds: (x_min, x_max, y_min, y_max)
    pub bounds: (T, T, T, T),
    /// Cell spacing in x direction
    pub dx: T,
    /// Cell spacing in y direction
    pub dy: T,
    /// Boundary conditions
    pub boundaries: HashMap<(usize, usize), BoundaryType>,
}

impl<T: RealField + FromPrimitive> StructuredGrid2D<T> {
    /// Create a new structured grid
    pub fn new(
        nx: usize,
        ny: usize,
        x_min: T,
        x_max: T,
        y_min: T,
        y_max: T,
    ) -> Result<Self> {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConfiguration(
                "Grid dimensions must be positive".to_string(),
            ));
        }

        let dx = (x_max.clone() - x_min.clone()) / T::from_usize(nx).unwrap_or_else(|| T::zero());
        let dy = (y_max.clone() - y_min.clone()) / T::from_usize(ny).unwrap_or_else(|| T::zero());

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
        Self::new(
            nx,
            ny,
            T::zero(),
            T::one(),
            T::zero(),
            T::one(),
        )
    }

    /// Set boundary condition for a cell
    pub fn set_boundary(&mut self, i: usize, j: usize, boundary_type: BoundaryType) {
        self.boundaries.insert((i, j), boundary_type);
    }

    /// Set boundary conditions for entire edges
    pub fn set_edge_boundary(&mut self, edge: GridEdge, boundary_type: BoundaryType) {
        match edge {
            GridEdge::Left => {
                for j in 0..self.ny {
                    self.set_boundary(0, j, boundary_type);
                }
            }
            GridEdge::Right => {
                for j in 0..self.ny {
                    self.set_boundary(self.nx - 1, j, boundary_type);
                }
            }
            GridEdge::Bottom => {
                for i in 0..self.nx {
                    self.set_boundary(i, 0, boundary_type);
                }
            }
            GridEdge::Top => {
                for i in 0..self.nx {
                    self.set_boundary(i, self.ny - 1, boundary_type);
                }
            }
        }
    }

    /// Get grid spacing
    pub fn spacing(&self) -> (T, T) {
        (self.dx.clone(), self.dy.clone())
    }

    /// Get domain bounds
    pub fn bounds(&self) -> (T, T, T, T) {
        self.bounds.clone()
    }
}

/// Grid edges for boundary condition setting
#[derive(Debug, Clone, Copy)]
pub enum GridEdge {
    /// Left edge (x = x_min)
    Left,
    /// Right edge (x = x_max)
    Right,
    /// Bottom edge (y = y_min)
    Bottom,
    /// Top edge (y = y_max)
    Top,
}

impl<T: RealField + FromPrimitive> Grid2D<T> for StructuredGrid2D<T> {
    fn nx(&self) -> usize {
        self.nx
    }

    fn ny(&self) -> usize {
        self.ny
    }

    fn cell_center(&self, i: usize, j: usize) -> Result<Vector2<T>> {
        if i >= self.nx || j >= self.ny {
            return Err(Error::InvalidConfiguration(
                "Cell indices out of bounds".to_string(),
            ));
        }

        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        let x = self.bounds.0.clone() +
                (T::from_usize(i).unwrap_or_else(|| T::zero()) + half.clone()) * self.dx.clone();
        let y = self.bounds.2.clone() +
                (T::from_usize(j).unwrap_or_else(|| T::zero()) + half) * self.dy.clone();

        Ok(Vector2::new(x, y))
    }

    fn cell_area(&self, i: usize, j: usize) -> Result<T> {
        if i >= self.nx || j >= self.ny {
            return Err(Error::InvalidConfiguration(
                "Cell indices out of bounds".to_string(),
            ));
        }

        Ok(self.dx.clone() * self.dy.clone())
    }

    fn neighbors(&self, i: usize, j: usize) -> Vec<(usize, usize)> {
        // For backward compatibility, collect from the efficient iterator
        // Note: prefer neighbor_iter() for better performance
        self.neighbor_iter(i, j).collect()
    }

    fn neighbor_iter(&self, i: usize, j: usize) -> impl Iterator<Item = (usize, usize)> {
        [(0isize, 1isize), (0, -1), (1, 0), (-1, 0)]
            .iter()
            .filter_map(move |&(di, dj)| {
                let ni = i as isize + di;
                let nj = j as isize + dj;
                if ni >= 0 && ni < self.nx as isize && nj >= 0 && nj < self.ny as isize {
                    Some((ni as usize, nj as usize))
                } else {
                    None
                }
            })
    }

    fn is_boundary(&self, i: usize, j: usize) -> bool {
        i == 0 || i == self.nx - 1 || j == 0 || j == self.ny - 1
    }

    fn boundary_type(&self, i: usize, j: usize) -> Option<BoundaryType> {
        self.boundaries.get(&(i, j)).copied()
    }
}

// GridIterator struct removed in favor of standard iterators

impl<T: RealField + FromPrimitive> StructuredGrid2D<T> {
    /// Create an iterator over all cells
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.ny).flat_map(move |j| (0..self.nx).map(move |i| (i, j)))
    }

    /// Create an iterator over boundary cells only
    pub fn boundary_iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter().filter(move |(i, j)| self.is_boundary(*i, *j))
    }

    /// Create an iterator over interior cells only
    pub fn interior_iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter().filter(move |(i, j)| !self.is_boundary(*i, *j))
    }

    /// Create an iterator over cell rows for efficient processing
    pub fn row_iter(&self) -> impl Iterator<Item = impl Iterator<Item = (usize, usize)> + '_> + '_ {
        (0..self.ny).map(move |j| (0..self.nx).map(move |i| (i, j)))
    }

    /// Create an iterator over cell columns for efficient processing
    pub fn col_iter(&self) -> impl Iterator<Item = impl Iterator<Item = (usize, usize)> + '_> + '_ {
        (0..self.nx).map(move |i| (0..self.ny).map(move |j| (i, j)))
    }

    /// Create windowed iterator for stencil operations (zero-allocation)
    pub fn stencil_iter(&self, stencil_size: usize) -> impl Iterator<Item = impl Iterator<Item = (usize, usize)>> + '_ {
        let half_size = stencil_size / 2;
        self.interior_iter()
            .filter(move |(i, j)| {
                *i >= half_size && *i < self.nx - half_size &&
                *j >= half_size && *j < self.ny - half_size
            })
            .map(move |(i, j)| {
                // Return an iterator that calculates indices on the fly
                (0..stencil_size).flat_map(move |dj| {
                    (0..stencil_size).map(move |di| {
                        (i + di - half_size, j + dj - half_size)
                    })
                })
            })
    }

    // neighbor_iter is now implemented in the Grid2D trait
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_structured_grid_creation() {
        let grid = StructuredGrid2D::<f64>::new(10, 20, 0.0, 1.0, 0.0, 2.0).unwrap();

        assert_eq!(grid.nx(), 10);
        assert_eq!(grid.ny(), 20);
        assert_eq!(grid.num_cells(), 200);

        let (dx, dy) = grid.spacing();
        assert_relative_eq!(dx, 0.1, epsilon = 1e-10);
        assert_relative_eq!(dy, 0.1, epsilon = 1e-10);
    }

    #[test]
    fn test_unit_square_grid() {
        let grid = StructuredGrid2D::<f64>::unit_square(5, 5).unwrap();

        let (x_min, x_max, y_min, y_max) = grid.bounds();
        assert_eq!(x_min, 0.0);
        assert_eq!(x_max, 1.0);
        assert_eq!(y_min, 0.0);
        assert_eq!(y_max, 1.0);
    }

    #[test]
    fn test_cell_center() {
        let grid = StructuredGrid2D::<f64>::unit_square(4, 4).unwrap();

        // Test center of first cell
        let center = grid.cell_center(0, 0).unwrap();
        assert_relative_eq!(center.x, 0.125, epsilon = 1e-10);
        assert_relative_eq!(center.y, 0.125, epsilon = 1e-10);

        // Test center of last cell
        let center = grid.cell_center(3, 3).unwrap();
        assert_relative_eq!(center.x, 0.875, epsilon = 1e-10);
        assert_relative_eq!(center.y, 0.875, epsilon = 1e-10);
    }

    #[test]
    fn test_cell_area() {
        let grid = StructuredGrid2D::<f64>::unit_square(4, 4).unwrap();

        let area = grid.cell_area(0, 0).unwrap();
        assert_relative_eq!(area, 0.0625, epsilon = 1e-10); // (1/4)^2
    }

    #[test]
    fn test_neighbors() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();

        // Corner cell should have 2 neighbors
        let neighbors = grid.neighbors(0, 0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&(1, 0)));
        assert!(neighbors.contains(&(0, 1)));

        // Center cell should have 4 neighbors
        let neighbors = grid.neighbors(1, 1);
        assert_eq!(neighbors.len(), 4);
        assert!(neighbors.contains(&(0, 1)));
        assert!(neighbors.contains(&(2, 1)));
        assert!(neighbors.contains(&(1, 0)));
        assert!(neighbors.contains(&(1, 2)));
    }

    #[test]
    fn test_boundary_detection() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();

        // Corner and edge cells should be boundary
        assert!(grid.is_boundary(0, 0));
        assert!(grid.is_boundary(0, 1));
        assert!(grid.is_boundary(1, 0));
        assert!(grid.is_boundary(2, 2));

        // Center cell should not be boundary
        assert!(!grid.is_boundary(1, 1));
    }

    #[test]
    fn test_boundary_conditions() {
        let mut grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();

        // Set boundary conditions
        grid.set_edge_boundary(GridEdge::Left, BoundaryType::Inlet);
        grid.set_edge_boundary(GridEdge::Right, BoundaryType::Outlet);
        grid.set_edge_boundary(GridEdge::Bottom, BoundaryType::Wall);
        grid.set_edge_boundary(GridEdge::Top, BoundaryType::Wall);

        // Check boundary types
        assert_eq!(grid.boundary_type(0, 1), Some(BoundaryType::Inlet));
        assert_eq!(grid.boundary_type(2, 1), Some(BoundaryType::Outlet));
        assert_eq!(grid.boundary_type(1, 0), Some(BoundaryType::Wall));
        assert_eq!(grid.boundary_type(1, 2), Some(BoundaryType::Wall));
        assert_eq!(grid.boundary_type(1, 1), None);
    }

    #[test]
    fn test_grid_iterator() {
        let grid = StructuredGrid2D::<f64>::unit_square(2, 2).unwrap();

        let cells: Vec<_> = grid.iter().collect();
        assert_eq!(cells.len(), 4);
        assert_eq!(cells, vec![(0, 0), (1, 0), (0, 1), (1, 1)]);
    }

    #[test]
    fn test_boundary_iterator() {
        let grid = StructuredGrid2D::<f64>::unit_square(3, 3).unwrap();

        let boundary_cells: Vec<_> = grid.boundary_iter().collect();
        assert_eq!(boundary_cells.len(), 8); // All cells except center

        let interior_cells: Vec<_> = grid.interior_iter().collect();
        assert_eq!(interior_cells.len(), 1); // Only center cell
        assert_eq!(interior_cells[0], (1, 1));
    }
}
