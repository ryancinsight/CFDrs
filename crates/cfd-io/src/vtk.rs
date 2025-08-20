//! VTK file format support.
//!
//! This module provides VTK (Visualization Toolkit) file format support
//! for CFD simulation data with zero-copy operations where possible.

use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::ToPrimitive;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

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
            _ => Err(Error::InvalidInput(format!("Unknown VTK cell type: {}", value))),
        }
    }
}

/// VTK dataset types
#[derive(Debug, Clone, Copy)]
pub enum VtkDatasetType {
    /// Structured points (regular grid)
    StructuredPoints,
    /// Structured grid
    StructuredGrid,
    /// Unstructured grid
    UnstructuredGrid,
    /// Polydata
    PolyData,
    /// Rectilinear grid
    RectilinearGrid,
}

/// VTK mesh data
#[derive(Debug, Clone)]
pub struct VtkMesh<T: RealField> {
    /// Points coordinates (x, y, z) flattened
    pub points: Vec<T>,
    /// Cell connectivity
    pub cells: Vec<Vec<usize>>,
    /// Cell types
    pub cell_types: Vec<VtkCellType>,
}

impl<T: RealField> VtkMesh<T> {
    /// Create new VTK mesh
    pub fn new(points: Vec<T>, cells: Vec<Vec<usize>>, cell_types: Vec<VtkCellType>) -> Self {
        Self {
            points,
            cells,
            cell_types,
        }
    }

    /// Get number of points
    pub fn num_points(&self) -> usize {
        self.points.len() / 3
    }

    /// Get number of cells
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Get point coordinates as iterator
    pub fn points_iter(&self) -> impl Iterator<Item = [&T; 3]> {
        self.points.chunks_exact(3).map(|chunk| [&chunk[0], &chunk[1], &chunk[2]])
    }

    /// Get mutable point coordinates as iterator
    pub fn points_iter_mut(&mut self) -> impl Iterator<Item = [&mut T; 3]> {
        self.points.chunks_exact_mut(3).map(|chunk| {
            let ptr = chunk.as_mut_ptr();
            unsafe { [&mut *ptr, &mut *ptr.add(1), &mut *ptr.add(2)] }
        })
    }

    /// Get point coordinates by index
    pub fn point(&self, index: usize) -> Option<[&T; 3]> {
        if index * 3 + 2 < self.points.len() {
            Some([
                &self.points[index * 3],
                &self.points[index * 3 + 1],
                &self.points[index * 3 + 2],
            ])
        } else {
            None
        }
    }

    /// Get mutable point coordinates by index
    pub fn point_mut(&mut self, index: usize) -> Option<[&mut T; 3]> {
        if index * 3 + 2 < self.points.len() {
            let ptr = self.points[index * 3..].as_mut_ptr();
            Some(unsafe { [&mut *ptr, &mut *ptr.add(1), &mut *ptr.add(2)] })
        } else {
            None
        }
    }
}

/// VTK file header information
#[derive(Debug, Clone)]
pub struct VtkHeader {
    /// File title/description
    pub title: String,
    /// Data format (ASCII or BINARY)
    pub format: String,
    /// Dataset type
    pub dataset_type: VtkDatasetType,
}

/// VTK writer for mesh data
pub struct VtkWriter<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + ToPrimitive> VtkWriter<T> {
    /// Create a new VTK writer
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Write mesh to VTK file
    pub fn write_mesh(&self, mesh: &VtkMesh<T>, path: &Path, title: &str) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        writeln!(writer, "# vtk DataFile Version 3.0")?;
        writeln!(writer, "{}", title)?;
        writeln!(writer, "ASCII")?;
        writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;
        writeln!(writer)?;

        // Write points
        writeln!(writer, "POINTS {} float", mesh.num_points())?;
        for point in mesh.points_iter() {
            writeln!(writer, "{} {} {}", 
                point[0].to_f64().unwrap_or(0.0),
                point[1].to_f64().unwrap_or(0.0),
                point[2].to_f64().unwrap_or(0.0)
            )?;
        }
        writeln!(writer)?;

        // Write cells
        let total_cell_data: usize = mesh.cells.iter().map(|cell| cell.len() + 1).sum();
        writeln!(writer, "CELLS {} {}", mesh.num_cells(), total_cell_data)?;
        for cell in &mesh.cells {
            write!(writer, "{}", cell.len())?;
            for &vertex_id in cell {
                write!(writer, " {}", vertex_id)?;
            }
            writeln!(writer)?;
        }
        writeln!(writer)?;

        // Write cell types
        writeln!(writer, "CELL_TYPES {}", mesh.num_cells())?;
        for &cell_type in &mesh.cell_types {
            writeln!(writer, "{}", cell_type as u8)?;
        }

        Ok(())
    }

    /// Write scalar field data
    pub fn write_scalar_field(
        &self,
        mesh: &VtkMesh<T>,
        field_data: &[T],
        field_name: &str,
        path: &Path,
        title: &str,
    ) -> Result<()> {
        if field_data.len() != mesh.num_points() {
            return Err(Error::InvalidConfiguration(
                "Field data length must match number of points".to_string(),
            ));
        }

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write mesh data first
        self.write_mesh_data(&mut writer, mesh, title)?;

        // Write point data
        writeln!(writer, "POINT_DATA {}", mesh.num_points())?;
        writeln!(writer, "SCALARS {} float 1", field_name)?;
        writeln!(writer, "LOOKUP_TABLE default")?;
        for value in field_data {
            writeln!(writer, "{}", value.to_f64().unwrap_or(0.0))?;
        }

        Ok(())
    }

    /// Write vector field data
    pub fn write_vector_field(
        &self,
        mesh: &VtkMesh<T>,
        field_data: &[(T, T, T)],
        field_name: &str,
        path: &Path,
        title: &str,
    ) -> Result<()> {
        if field_data.len() != mesh.num_points() {
            return Err(Error::InvalidConfiguration(
                "Field data length must match number of points".to_string(),
            ));
        }

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write mesh data first
        self.write_mesh_data(&mut writer, mesh, title)?;

        // Write point data
        writeln!(writer, "POINT_DATA {}", mesh.num_points())?;
        writeln!(writer, "VECTORS {} float", field_name)?;
        for (x, y, z) in field_data {
            writeln!(writer, "{} {} {}", 
                x.to_f64().unwrap_or(0.0),
                y.to_f64().unwrap_or(0.0),
                z.to_f64().unwrap_or(0.0)
            )?;
        }

        Ok(())
    }

    fn write_mesh_data(&self, writer: &mut BufWriter<File>, mesh: &VtkMesh<T>, title: &str) -> Result<()> {
        // Write header
        writeln!(writer, "# vtk DataFile Version 3.0")?;
        writeln!(writer, "{}", title)?;
        writeln!(writer, "ASCII")?;
        writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;
        writeln!(writer)?;

        // Write points
        writeln!(writer, "POINTS {} float", mesh.num_points())?;
        for point in mesh.points_iter() {
            writeln!(writer, "{} {} {}", 
                point[0].to_f64().unwrap_or(0.0),
                point[1].to_f64().unwrap_or(0.0),
                point[2].to_f64().unwrap_or(0.0)
            )?;
        }
        writeln!(writer)?;

        // Write cells
        let total_cell_data: usize = mesh.cells.iter().map(|cell| cell.len() + 1).sum();
        writeln!(writer, "CELLS {} {}", mesh.num_cells(), total_cell_data)?;
        for cell in &mesh.cells {
            write!(writer, "{}", cell.len())?;
            for &vertex_id in cell {
                write!(writer, " {}", vertex_id)?;
            }
            writeln!(writer)?;
        }
        writeln!(writer)?;

        // Write cell types
        writeln!(writer, "CELL_TYPES {}", mesh.num_cells())?;
        for &cell_type in &mesh.cell_types {
            writeln!(writer, "{}", cell_type as u8)?;
        }
        writeln!(writer)?;

        Ok(())
    }
}

impl<T: RealField + ToPrimitive> Default for VtkWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// VTK reader for mesh data
pub struct VtkReader<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> VtkReader<T> {
    /// Create a new VTK reader
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Read mesh data from VTK file
    pub fn read_mesh(&self, path: &Path) -> Result<VtkMesh<T>> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        
        // Read header and get dataset type
        let header = self.read_header_from_reader(&mut reader)?;
        
        // Parse based on dataset type using the same reader
        match header.dataset_type {
            VtkDatasetType::UnstructuredGrid => self.read_unstructured_grid(reader),
            VtkDatasetType::StructuredGrid => self.read_structured_grid(reader),
            _ => Err(Error::NotImplemented(
                format!("VTK reading for {:?} not yet implemented", header.dataset_type)
            ))
        }
    }

    /// Read VTK file header from a BufReader
    fn read_header_from_reader(&self, reader: &mut BufReader<File>) -> Result<VtkHeader> {
        let mut lines = reader.lines();

        // Read header lines
        let _version = lines.next()
            .ok_or_else(|| Error::IoError(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Missing VTK version line",
            )))??;
        
        let title = lines.next()
            .ok_or_else(|| Error::IoError(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Missing title line",
            )))??;
        
        let format = lines.next()
            .ok_or_else(|| Error::IoError(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Missing format line",
            )))??;
        
        let dataset = lines.next()
            .ok_or_else(|| Error::IoError(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Missing dataset line",
            )))??;

        Ok(VtkHeader {
            title,
            format: format.trim().to_string(),
            dataset_type: self.parse_dataset_type(&dataset)?,
        })
    }
    
    /// Parse a line into a vector of values with better error handling
    fn parse_line<V: FromStr>(&self, line: &str, context: &str) -> Result<Vec<V>> {
        line.split_whitespace()
            .map(|s| s.parse::<V>().map_err(|_| 
                Error::InvalidInput(format!("Failed to parse '{}' in {}", s, context))
            ))
            .collect()
    }

    /// Read unstructured grid data
    fn read_unstructured_grid(&self, reader: BufReader<File>) -> Result<VtkMesh<T>> {
        let mut points = Vec::new();
        let mut cells = Vec::new();
        let mut cell_types = Vec::new();
        let mut lines = reader.lines();
        
        while let Some(line) = lines.next() {
            let line = line?;
            let line = line.trim();
            
            if line.starts_with("POINTS") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                let n_points = parts.get(1)
                    .ok_or_else(|| Error::InvalidInput("Missing points count".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("Invalid points count".to_string()))?;
                
                // Read point coordinates
                for i in 0..n_points {
                    if let Some(point_line) = lines.next() {
                        let coords: Vec<f64> = self.parse_line(&point_line?, 
                            &format!("point {} coordinates", i))?;
                        
                        if coords.len() < 3 {
                            return Err(Error::InvalidInput(
                                format!("Point {} has insufficient coordinates (need 3, got {})", i, coords.len())
                            ));
                        }
                        
                        // Convert to target type with proper error handling
                        points.push(T::from_f64(coords[0])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert x-coordinate {} to target type", coords[0])
                            ))?);
                        points.push(T::from_f64(coords[1])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert y-coordinate {} to target type", coords[1])
                            ))?);
                        points.push(T::from_f64(coords[2])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert z-coordinate {} to target type", coords[2])
                            ))?);
                    } else {
                        return Err(Error::InvalidInput(
                            format!("Missing point data for point {}", i)
                        ));
                    }
                }
            } else if line.starts_with("CELLS") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                let n_cells = parts.get(1)
                    .ok_or_else(|| Error::InvalidInput("Missing cells count".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("Invalid cells count".to_string()))?;
                
                // Read cell connectivity
                for i in 0..n_cells {
                    if let Some(cell_line) = lines.next() {
                        let indices: Vec<usize> = self.parse_line(&cell_line?, 
                            &format!("cell {} connectivity", i))?;
                        
                        if indices.is_empty() {
                            return Err(Error::InvalidInput(
                                format!("Cell {} has no connectivity data", i)
                            ));
                        }
                        
                        // Skip the first element (count) and take the rest as indices
                        cells.push(indices.into_iter().skip(1).collect());
                    } else {
                        return Err(Error::InvalidInput(
                            format!("Missing cell data for cell {}", i)
                        ));
                    }
                }
            } else if line.starts_with("CELL_TYPES") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                let n_types = parts.get(1)
                    .ok_or_else(|| Error::InvalidInput("Missing cell types count".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("Invalid cell types count".to_string()))?;
                
                // Read cell types
                for i in 0..n_types {
                    if let Some(type_line) = lines.next() {
                        let cell_type_val: u8 = self.parse_line(&type_line?, 
                            &format!("cell type {}", i))?
                            .into_iter()
                            .next()
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Missing cell type value for cell {}", i)
                            ))?;
                        
                        cell_types.push(VtkCellType::from_u8(cell_type_val)?);
                    } else {
                        return Err(Error::InvalidInput(
                            format!("Missing cell type data for cell {}", i)
                        ));
                    }
                }
            }
        }
        
        Ok(VtkMesh { points, cells, cell_types })
    }
    
    /// Read structured grid data
    fn read_structured_grid(&self, reader: BufReader<File>) -> Result<VtkMesh<T>> {
        let mut points = Vec::new();
        let mut lines = reader.lines();
        
        while let Some(line) = lines.next() {
            let line = line?;
            let line = line.trim();
            
            if line.starts_with("DIMENSIONS") {
                // Parse grid dimensions for structured grid
                // This determines implicit connectivity
                continue;
            } else if line.starts_with("POINTS") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                let n_points = parts.get(1)
                    .ok_or_else(|| Error::InvalidInput("Missing points count".to_string()))?
                    .parse::<usize>()
                    .map_err(|_| Error::InvalidInput("Invalid points count".to_string()))?;
                
                // Read point coordinates
                for i in 0..n_points {
                    if let Some(point_line) = lines.next() {
                        let coords: Vec<f64> = self.parse_line(&point_line?, 
                            &format!("point {} coordinates", i))?;
                        
                        if coords.len() < 3 {
                            return Err(Error::InvalidInput(
                                format!("Point {} has insufficient coordinates (need 3, got {})", i, coords.len())
                            ));
                        }
                        
                        // Convert to target type with proper error handling
                        points.push(T::from_f64(coords[0])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert x-coordinate {} to target type", coords[0])
                            ))?);
                        points.push(T::from_f64(coords[1])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert y-coordinate {} to target type", coords[1])
                            ))?);
                        points.push(T::from_f64(coords[2])
                            .ok_or_else(|| Error::InvalidInput(
                                format!("Cannot convert z-coordinate {} to target type", coords[2])
                            ))?);
                    } else {
                        return Err(Error::InvalidInput(
                            format!("Missing point data for point {}", i)
                        ));
                    }
                }
            }
        }
        
        // For structured grids, cells and cell_types are implicit
        Ok(VtkMesh { 
            points, 
            cells: Vec::new(), // Implicit connectivity
            cell_types: Vec::new() // All cells are hexahedra for structured grids
        })
    }

    fn parse_dataset_type(&self, line: &str) -> Result<VtkDatasetType> {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 2 || parts[0] != "DATASET" {
            return Err(Error::InvalidConfiguration(
                "Invalid dataset line".to_string(),
            ));
        }

        match parts[1] {
            "STRUCTURED_POINTS" => Ok(VtkDatasetType::StructuredPoints),
            "STRUCTURED_GRID" => Ok(VtkDatasetType::StructuredGrid),
            "UNSTRUCTURED_GRID" => Ok(VtkDatasetType::UnstructuredGrid),
            "POLYDATA" => Ok(VtkDatasetType::PolyData),
            "RECTILINEAR_GRID" => Ok(VtkDatasetType::RectilinearGrid),
            unknown => Err(Error::InvalidConfiguration(
                format!("Unknown dataset type: {}", unknown),
            )),
        }
    }
}

impl<T: RealField> Default for VtkReader<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// VTK mesh builder with simplified API
pub struct VtkMeshBuilder<T: RealField> {
    points: Vec<T>,
    cells: Vec<Vec<usize>>,
    cell_types: Vec<VtkCellType>,
}

impl<T: RealField> VtkMeshBuilder<T> {
    /// Create new mesh builder
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            cells: Vec::new(),
            cell_types: Vec::new(),
        }
    }

    /// Add a point
    pub fn add_point(mut self, x: T, y: T, z: T) -> Self {
        self.points.extend_from_slice(&[x, y, z]);
        self
    }

    /// Add points from iterator - unified API
    pub fn add_points<I>(mut self, points: I) -> Self
    where
        I: IntoIterator<Item = (T, T, T)>,
    {
        // Reserve capacity to optimize performance
        let iter = points.into_iter();
        let size_hint = iter.size_hint();
        if let (lower, Some(upper)) = size_hint {
            if lower == upper {
                self.points.reserve(lower * 3);
            }
        }
        
        // Use iterator combinators to flatten coordinates efficiently
        self.points.extend(
            iter.flat_map(|(x, y, z)| [x, y, z])
        );
        self
    }

    /// Add a cell
    pub fn add_cell(mut self, connectivity: Vec<usize>, cell_type: VtkCellType) -> Self {
        self.cells.push(connectivity);
        self.cell_types.push(cell_type);
        self
    }

    /// Build the mesh
    pub fn build(self) -> VtkMesh<T> {
        VtkMesh::new(self.points, self.cells, self.cell_types)
    }
}

impl<T: RealField> Default for VtkMeshBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vtk_cell_type_from_u8() {
        assert!(matches!(VtkCellType::from_u8(1), Ok(VtkCellType::Vertex)));
        assert!(matches!(VtkCellType::from_u8(3), Ok(VtkCellType::Line)));
        assert!(matches!(VtkCellType::from_u8(5), Ok(VtkCellType::Triangle)));
        assert!(matches!(VtkCellType::from_u8(9), Ok(VtkCellType::Quad)));
        assert!(matches!(VtkCellType::from_u8(10), Ok(VtkCellType::Tetrahedron)));
        assert!(matches!(VtkCellType::from_u8(12), Ok(VtkCellType::Hexahedron)));
        assert!(VtkCellType::from_u8(99).is_err());
    }

    #[test]
    fn test_vtk_mesh_builder() {
        let mesh = VtkMeshBuilder::<f64>::new()
            .add_point(0.0, 0.0, 0.0)
            .add_point(1.0, 0.0, 0.0)
            .add_point(1.0, 1.0, 0.0)
            .add_point(0.0, 1.0, 0.0)
            .add_cell(vec![0, 1, 2, 3], VtkCellType::Quad)
            .build();

        assert_eq!(mesh.num_points(), 4);
        assert_eq!(mesh.num_cells(), 1);
    }

    #[test]
    fn test_vtk_mesh_builder_unified_api() {
        let points_data = vec![(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0)];
        
        let mesh = VtkMeshBuilder::<f64>::new()
            .add_points(points_data.iter().copied()) // for Copy types
            .add_cell(vec![0, 1, 2], VtkCellType::Triangle)
            .build();

        assert_eq!(mesh.num_points(), 3);
        assert_eq!(mesh.num_cells(), 1);
    }

    #[test]
    fn test_points_iterator() {
        let points = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0];
        let mesh = VtkMesh::new(points, vec![], vec![]);

        let collected: Vec<_> = mesh.points_iter().collect();
        assert_eq!(collected.len(), 3);
        assert_eq!(collected[0], [&0.0, &0.0, &0.0]);
        assert_eq!(collected[1], [&1.0, &0.0, &0.0]);
        assert_eq!(collected[2], [&1.0, &1.0, &0.0]);
    }
}