//! VTK data types and structures

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use std::collections::HashMap;

/// VTK data types
#[derive(Debug, Clone, Copy)]
pub enum VtkDataType {
    /// Scalar data
    Scalar,
    /// Vector data (3 components)
    Vector,
    /// Tensor data (9 components)
    Tensor,
}

/// VTK cell types
#[derive(Debug, Clone, Copy)]
pub enum VtkCellType {
    /// Vertex (1 point)
    Vertex = 1,
    /// Line (2 points)
    Line = 3,
    /// Triangle (3 points)
    Triangle = 5,
    /// Quad (4 points)
    Quad = 9,
    /// Tetrahedron (4 points)
    Tetrahedron = 10,
    /// Hexahedron (8 points)
    Hexahedron = 12,
}

impl VtkCellType {
    /// Convert from u8 representation
    pub fn from_u8(value: u8) -> Result<Self> {
        match value {
            1 => Ok(VtkCellType::Vertex),
            3 => Ok(VtkCellType::Line),
            5 => Ok(VtkCellType::Triangle),
            9 => Ok(VtkCellType::Quad),
            10 => Ok(VtkCellType::Tetrahedron),
            12 => Ok(VtkCellType::Hexahedron),
            _ => Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unknown VTK cell type: {value}"),
            ))),
        }
    }

    /// Get the number of points for this cell type
    pub fn num_points(&self) -> usize {
        match self {
            VtkCellType::Vertex => 1,
            VtkCellType::Line => 2,
            VtkCellType::Triangle => 3,
            VtkCellType::Quad => 4,
            VtkCellType::Tetrahedron => 4,
            VtkCellType::Hexahedron => 8,
        }
    }
}

/// VTK mesh data structure
#[derive(Debug, Clone)]
pub struct VtkMesh<T: RealField + Copy> {
    /// Points (vertices)
    pub points: Vec<[T; 3]>,
    /// Cells (connectivity)
    pub cells: Vec<Vec<usize>>,
    /// Cell types
    pub cell_types: Vec<VtkCellType>,
    /// Point data
    pub point_data: HashMap<String, Vec<T>>,
    /// Cell data
    pub cell_data: HashMap<String, Vec<T>>,
}

impl<T: RealField + Copy> VtkMesh<T> {
    /// Create a new empty VTK mesh
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            cells: Vec::new(),
            cell_types: Vec::new(),
            point_data: HashMap::new(),
            cell_data: HashMap::new(),
        }
    }

    /// Add a point to the mesh
    pub fn add_point(&mut self, x: T, y: T, z: T) -> usize {
        self.points.push([x, y, z]);
        self.points.len() - 1
    }

    /// Add a cell to the mesh
    pub fn add_cell(&mut self, cell_type: VtkCellType, connectivity: Vec<usize>) {
        self.cells.push(connectivity);
        self.cell_types.push(cell_type);
    }

    /// Add point data
    pub fn add_point_data(&mut self, name: String, data: Vec<T>) -> Result<()> {
        if data.len() != self.points.len() {
            return Err(Error::InvalidInput(format!(
                "Point data size mismatch: expected {}, got {}",
                self.points.len(),
                data.len()
            )));
        }
        self.point_data.insert(name, data);
        Ok(())
    }

    /// Add cell data
    pub fn add_cell_data(&mut self, name: String, data: Vec<T>) -> Result<()> {
        if data.len() != self.cells.len() {
            return Err(Error::InvalidInput(format!(
                "Cell data size mismatch: expected {}, got {}",
                self.cells.len(),
                data.len()
            )));
        }
        self.cell_data.insert(name, data);
        Ok(())
    }

    /// Get the number of points
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Get the number of cells
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }
}

impl<T: RealField + Copy> Default for VtkMesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// VTK file header information
#[derive(Debug, Clone)]
pub struct VtkHeader {
    /// Title of the dataset
    pub title: String,
    /// Data format (ASCII or BINARY)
    pub format: String,
    /// Dataset type (`STRUCTURED_GRID`, `UNSTRUCTURED_GRID`, etc.)
    pub dataset_type: String,
}
