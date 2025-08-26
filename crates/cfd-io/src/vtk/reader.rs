//! VTK file reader

use super::types::{VtkCellType, VtkHeader, VtkMesh};
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;
/// VTK file reader
pub struct VtkReader<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}
impl<T: RealField + Copy + FromStr> VtkReader<T> {
    /// Create a new VTK reader
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
    /// Read mesh from VTK file
    pub fn read_mesh(&self, path: &Path) -> Result<VtkMesh<T>> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        // Read header
        let header = self.read_header(&mut reader)?;
        // For now, we only support UNSTRUCTURED_GRID
        if !header.dataset_type.contains("UNSTRUCTURED_GRID") {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Unsupported dataset type: {}", header.dataset_type),
            )));
        let mut mesh = VtkMesh::new();
        let mut line = String::new();
        // Read the rest of the file
        while reader.read_line(&mut line)? > 0 {
            let trimmed = line.trim();
            if trimmed.starts_with("POINTS") {
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() >= 2 {
                    let num_points: usize = parts[1].parse().map_err(|_| {
                        Error::Io(std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "Invalid number of points",
                        ))
                    })?;
                    // Read points
                    for _ in 0..num_points {
                        line.clear();
                        reader.read_line(&mut line)?;
                        let coords: Vec<&str> = line.trim().split_whitespace().collect();
                        if coords.len() >= 3 {
                            let x = T::from_str(coords[0]).map_err(|_| {
                                Error::Io(std::io::Error::new(
                                    std::io::ErrorKind::InvalidData,
                                    "Invalid coordinate",
                                ))
                            })?;
                            let y = T::from_str(coords[1]).map_err(|_| {
                            let z = T::from_str(coords[2]).map_err(|_| {
                            mesh.add_point(x, y, z);
                        }
                    }
                }
            } else if trimmed.starts_with("CELLS") {
                    let num_cells: usize = parts[1].parse().map_err(|_| {
                            "Invalid number of cells",
                    // Read cells
                    for _ in 0..num_cells {
                        let indices: Vec<usize> = line
                            .trim()
                            .split_whitespace()
                            .skip(1) // First number is the count
                            .filter_map(|s| s.parse().ok())
                            .collect();
                        mesh.cells.push(indices);
            } else if trimmed.starts_with("CELL_TYPES") {
                    let num_types: usize = parts[1].parse().map_err(|_| {
                            "Invalid number of cell types",
                    // Read cell types
                    for _ in 0..num_types {
                        let type_id: u8 = line.trim().parse().map_err(|_| {
                            Error::Io(std::io::Error::new(
                                std::io::ErrorKind::InvalidData,
                                "Invalid cell type",
                            ))
                        })?;
                        mesh.cell_types.push(VtkCellType::from_u8(type_id)?);
            }
            line.clear();
        Ok(mesh)
    fn read_header(&self, reader: &mut BufReader<File>) -> Result<VtkHeader> {
        let mut lines = reader.lines();
        // Read and validate header lines
        let version = lines.next().ok_or_else(|| {
            Error::Io(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                "Missing VTK version line",
            ))
        })??;
        // Validate version format
        if !version.starts_with("# vtk DataFile Version") {
                format!("Invalid VTK version line: {}", version),
        let title = lines.next().ok_or_else(|| {
                "Missing title line",
        let format = lines.next().ok_or_else(|| {
                "Missing format line",
        let dataset_type = lines.next().ok_or_else(|| {
                "Missing dataset type line",
        Ok(VtkHeader {
            title,
            format,
            dataset_type,
        })
impl<T: RealField + Copy + FromStr> Default for VtkReader<T> {
    fn default() -> Self {
        Self::new()
