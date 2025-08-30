//! Grid generation module
//!
//! This module provides structured grid generation functionality for CFD simulations.
//! The grid builder creates mesh objects with proper connectivity.

use crate::mesh::Mesh;
use crate::topology::{Cell, Vertex};
use nalgebra::{Point3, RealField};
use num_traits::FromPrimitive;

/// Builder for structured grid generation
///
/// This builder creates a structured mesh with the specified resolution and bounds.
/// It generates vertices and hexahedral cells with proper connectivity.
#[derive(Debug, Clone)]
pub struct StructuredGridBuilder<T: RealField + Copy> {
    /// Number of cells in x direction
    nx: usize,
    /// Number of cells in y direction
    ny: usize,
    /// Number of cells in z direction
    nz: usize,
    /// Domain bounds: ((x_min, x_max), (y_min, y_max), (z_min, z_max))
    bounds: ((T, T), (T, T), (T, T)),
}

impl<T: RealField + Copy + FromPrimitive> StructuredGridBuilder<T> {
    /// Create a new structured grid builder with default unit cube bounds
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            nx,
            ny,
            nz,
            bounds: (
                (T::zero(), T::one()),
                (T::zero(), T::one()),
                (T::zero(), T::one()),
            ),
        }
    }

    /// Set the domain bounds
    #[must_use]
    pub fn with_bounds(mut self, bounds: ((T, T), (T, T), (T, T))) -> Self {
        self.bounds = bounds;
        self
    }

    /// Set bounds from min/max points
    #[must_use]
    pub fn with_domain(mut self, min: Point3<T>, max: Point3<T>) -> Self {
        self.bounds = ((min.x, max.x), (min.y, max.y), (min.z, max.z));
        self
    }

    /// Generate the mesh from the grid parameters
    pub fn build(&self) -> Mesh<T> {
        let mut mesh = Mesh::new();

        // Calculate grid spacing
        let dx = (self.bounds.0 .1 - self.bounds.0 .0)
            / T::from_usize(self.nx).expect("Failed to convert nx to numeric type");
        let dy = (self.bounds.1 .1 - self.bounds.1 .0)
            / T::from_usize(self.ny).expect("Failed to convert ny to numeric type");
        let dz = (self.bounds.2 .1 - self.bounds.2 .0)
            / T::from_usize(self.nz).expect("Failed to convert nz to numeric type");

        // Generate vertices
        // We need (nx+1) x (ny+1) x (nz+1) vertices
        for k in 0..=self.nz {
            for j in 0..=self.ny {
                for i in 0..=self.nx {
                    let x = self.bounds.0 .0
                        + T::from_usize(i).expect("Failed to convert i to numeric type") * dx;
                    let y = self.bounds.1 .0
                        + T::from_usize(j).expect("Failed to convert j to numeric type") * dy;
                    let z = self.bounds.2 .0
                        + T::from_usize(k).expect("Failed to convert k to numeric type") * dz;

                    mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
                }
            }
        }

        // Generate hexahedral cells
        // We have nx x ny x nz cells
        let nx1 = self.nx + 1;
        let ny1 = self.ny + 1;

        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    // Calculate the 8 vertex indices for this hexahedron
                    // Using the standard VTK hexahedron vertex ordering
                    let v0 = k * ny1 * nx1 + j * nx1 + i;
                    let v1 = v0 + 1;
                    let v2 = v0 + nx1 + 1;
                    let v3 = v0 + nx1;
                    let v4 = v0 + ny1 * nx1;
                    let v5 = v4 + 1;
                    let v6 = v4 + nx1 + 1;
                    let v7 = v4 + nx1;

                    mesh.add_cell(Cell::hexahedron(vec![v0, v1, v2, v3, v4, v5, v6, v7]));
                }
            }
        }

        mesh
    }

    /// Get the total number of cells that will be generated
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Get the total number of vertices that will be generated
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
    }
}

/// Builder for 2D structured grid generation
///
/// This builder creates a 2D structured mesh in the XY plane.
#[derive(Debug, Clone)]
pub struct StructuredGrid2DBuilder<T: RealField + Copy> {
    /// Number of cells in x direction
    nx: usize,
    /// Number of cells in y direction
    ny: usize,
    /// Domain bounds: ((x_min, x_max), (y_min, y_max))
    bounds: ((T, T), (T, T)),
}

impl<T: RealField + Copy + FromPrimitive> StructuredGrid2DBuilder<T> {
    /// Create a new 2D structured grid builder
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            bounds: ((T::zero(), T::one()), (T::zero(), T::one())),
        }
    }

    /// Set the domain bounds
    #[must_use]
    pub fn with_bounds(mut self, x_bounds: (T, T), y_bounds: (T, T)) -> Self {
        self.bounds = (x_bounds, y_bounds);
        self
    }

    /// Generate a 2D mesh (as a thin 3D mesh with one layer)
    pub fn build(&self) -> Mesh<T> {
        // Create a thin 3D mesh with nz=1
        StructuredGridBuilder::new(self.nx, self.ny, 1)
            .with_bounds((
                self.bounds.0,
                self.bounds.1,
                (T::zero(), T::zero()), // Zero thickness in z
            ))
            .build()
    }
}
