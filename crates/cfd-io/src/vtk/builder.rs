//! VTK mesh builder

use super::types::{VtkCellType, VtkMesh};
use cfd_core::error::Result;
use nalgebra::RealField;

/// Builder for constructing VTK meshes
pub struct VtkMeshBuilder<T: RealField + Copy> {
    mesh: VtkMesh<T>,
}

impl<T: RealField + Copy> VtkMeshBuilder<T> {
    /// Create a new mesh builder
    #[must_use]
    pub fn new() -> Self {
        Self {
            mesh: VtkMesh::new(),
        }
    }

    /// Add a point to the mesh
    pub fn add_point(&mut self, x: T, y: T, z: T) -> &mut Self {
        self.mesh.add_point(x, y, z);
        self
    }

    /// Add a cell to the mesh
    pub fn add_cell(&mut self, cell_type: VtkCellType, connectivity: Vec<usize>) -> &mut Self {
        self.mesh.add_cell(cell_type, connectivity);
        self
    }

    /// Add point data
    pub fn add_point_data(&mut self, name: String, data: Vec<T>) -> Result<&mut Self> {
        self.mesh.add_point_data(name, data)?;
        Ok(self)
    }

    /// Add cell data
    pub fn add_cell_data(&mut self, name: String, data: Vec<T>) -> Result<&mut Self> {
        self.mesh.add_cell_data(name, data)?;
        Ok(self)
    }

    /// Build the mesh
    #[must_use]
    pub fn build(self) -> VtkMesh<T> {
        self.mesh
    }

    /// Create a simple box mesh
    pub fn create_box(
        &mut self,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
    ) -> Result<&mut Self> {
        // Create points
        for k in 0..=nz {
            for j in 0..=ny {
                for i in 0..=nx {
                    let x = T::from_usize(i).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(format!(
                            "Cannot convert index {} to floating point",
                            i
                        ))
                    })? * dx;
                    let y = T::from_usize(j).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(format!(
                            "Cannot convert index {} to floating point",
                            j
                        ))
                    })? * dy;
                    let z = T::from_usize(k).ok_or_else(|| {
                        cfd_core::error::Error::InvalidConfiguration(format!(
                            "Cannot convert index {} to floating point",
                            k
                        ))
                    })? * dz;
                    self.add_point(x, y, z);
                }
            }
        }

        // Create hexahedral cells
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let base = k * (ny + 1) * (nx + 1) + j * (nx + 1) + i;
                    let connectivity = vec![
                        base,
                        base + 1,
                        base + (nx + 1) + 1,
                        base + (nx + 1),
                        base + (ny + 1) * (nx + 1),
                        base + (ny + 1) * (nx + 1) + 1,
                        base + (ny + 1) * (nx + 1) + (nx + 1) + 1,
                        base + (ny + 1) * (nx + 1) + (nx + 1),
                    ];
                    self.add_cell(VtkCellType::Hexahedron, connectivity);
                }
            }
        }

        Ok(self)
    }
}

impl<T: RealField + Copy> Default for VtkMeshBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}
