//! Grid generation module
//!
//! This module provides structured grid generation functionality for CFD simulations.
//! The grid builder creates mesh objects with proper connectivity.

use crate::mesh::Mesh;
use crate::topology::{Cell, Face, Vertex};
use crate::{error::MeshError, error::Result};
use nalgebra::{Point3, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;

fn get_or_add_quad_face<T: RealField + Copy>(
    mesh: &mut Mesh<T>,
    face_map: &mut HashMap<[usize; 4], usize>,
    verts: [usize; 4],
) -> usize {
    let mut key = verts;
    key.sort_unstable();
    if let Some(&idx) = face_map.get(&key) {
        idx
    } else {
        let idx = mesh.add_face(Face::quad(verts[0], verts[1], verts[2], verts[3]));
        face_map.insert(key, idx);
        idx
    }
}

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
    /// Domain bounds: ((`x_min`, `x_max`), (`y_min`, `y_max`), (`z_min`, `z_max`))
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
    pub fn build(&self) -> Result<Mesh<T>> {
        if self.nx == 0 || self.ny == 0 || self.nz == 0 {
            return Err(MeshError::GridError(
                "nx, ny, and nz must be greater than zero".to_string(),
            ));
        }

        if self.bounds.0 .1 <= self.bounds.0 .0 || self.bounds.1 .1 <= self.bounds.1 .0 {
            return Err(MeshError::GridError(
                "x and y bounds must satisfy min < max".to_string(),
            ));
        }

        if self.bounds.2 .1 < self.bounds.2 .0 {
            return Err(MeshError::GridError(
                "z bounds must satisfy min <= max".to_string(),
            ));
        }

        if self.bounds.2 .1 == self.bounds.2 .0 && self.nz != 1 {
            return Err(MeshError::GridError(
                "zero-thickness z bounds require nz == 1".to_string(),
            ));
        }

        let mut mesh = Mesh::new();

        let nx_t = T::from_usize(self.nx).ok_or_else(|| {
            MeshError::GridError(format!("failed to convert nx={} to scalar type", self.nx))
        })?;
        let ny_t = T::from_usize(self.ny).ok_or_else(|| {
            MeshError::GridError(format!("failed to convert ny={} to scalar type", self.ny))
        })?;
        let nz_t = T::from_usize(self.nz).ok_or_else(|| {
            MeshError::GridError(format!("failed to convert nz={} to scalar type", self.nz))
        })?;

        let dx = (self.bounds.0 .1 - self.bounds.0 .0) / nx_t;
        let dy = (self.bounds.1 .1 - self.bounds.1 .0) / ny_t;
        let dz = (self.bounds.2 .1 - self.bounds.2 .0) / nz_t;

        // Generate vertices
        // We need (nx+1) x (ny+1) x (nz+1) vertices
        for k in 0..=self.nz {
            let k_t = T::from_usize(k).ok_or_else(|| {
                MeshError::GridError(format!("failed to convert k={k} to scalar type"))
            })?;
            for j in 0..=self.ny {
                let j_t = T::from_usize(j).ok_or_else(|| {
                    MeshError::GridError(format!("failed to convert j={j} to scalar type"))
                })?;
                for i in 0..=self.nx {
                    let i_t = T::from_usize(i).ok_or_else(|| {
                        MeshError::GridError(format!("failed to convert i={i} to scalar type"))
                    })?;

                    let x = self.bounds.0 .0 + i_t * dx;
                    let y = self.bounds.1 .0 + j_t * dy;
                    let z = self.bounds.2 .0 + k_t * dz;

                    mesh.add_vertex(Vertex::new(Point3::new(x, y, z)));
                }
            }
        }

        // Generate hexahedral cells
        // We have nx x ny x nz cells
        let nx1 = self.nx + 1;
        let ny1 = self.ny + 1;
        let mut face_map: HashMap<[usize; 4], usize> = HashMap::new();

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

                    let f_bottom = get_or_add_quad_face(&mut mesh, &mut face_map, [v0, v1, v2, v3]);
                    let f_top = get_or_add_quad_face(&mut mesh, &mut face_map, [v4, v5, v6, v7]);
                    let f_0 = get_or_add_quad_face(&mut mesh, &mut face_map, [v0, v1, v5, v4]);
                    let f_1 = get_or_add_quad_face(&mut mesh, &mut face_map, [v1, v2, v6, v5]);
                    let f_2 = get_or_add_quad_face(&mut mesh, &mut face_map, [v2, v3, v7, v6]);
                    let f_3 = get_or_add_quad_face(&mut mesh, &mut face_map, [v3, v0, v4, v7]);

                    mesh.add_cell(Cell::hexahedron(vec![f_bottom, f_top, f_0, f_1, f_2, f_3]));
                }
            }
        }

        Ok(mesh)
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
    /// Domain bounds: ((`x_min`, `x_max`), (`y_min`, `y_max`))
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
    pub fn build(&self) -> Result<Mesh<T>> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_rejects_zero_dimensions() {
        let builder = StructuredGridBuilder::<f64>::new(0, 1, 1);
        assert!(builder.build().is_err());
    }

    #[test]
    fn build_rejects_invalid_bounds() {
        let builder = StructuredGridBuilder::<f64>::new(1, 1, 1).with_bounds((
            (1.0, 0.0),
            (0.0, 1.0),
            (0.0, 1.0),
        ));
        assert!(builder.build().is_err());
    }

    #[test]
    fn build_2d_produces_mesh() {
        let mesh = StructuredGrid2DBuilder::<f64>::new(2, 3).build().unwrap();
        assert_eq!(mesh.cells().len(), 2 * 3);
        assert_eq!(mesh.vertices().len(), (2 + 1) * (3 + 1) * (1 + 1));
    }
}
