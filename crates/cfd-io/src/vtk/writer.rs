//! VTK file writer

use super::types::VtkMesh;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::cast::ToPrimitive;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// VTK file writer
pub struct VtkWriter<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + ToPrimitive> VtkWriter<T> {
    /// Create a new VTK writer
    #[must_use]
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
        writeln!(writer, "{title}")?;
        writeln!(writer, "ASCII")?;
        writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;
        writeln!(writer)?;

        // Write points (using double precision since we convert to f64)
        writeln!(writer, "POINTS {} double", mesh.num_points())?;
        for point in &mesh.points {
            let x = point[0].to_f64().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Failed to convert X coordinate to f64 for VTK output",
                )
            })?;
            let y = point[1].to_f64().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Failed to convert Y coordinate to f64 for VTK output",
                )
            })?;
            let z = point[2].to_f64().ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Failed to convert Z coordinate to f64 for VTK output",
                )
            })?;
            writeln!(writer, "{x} {y} {z}")?;
        }
        writeln!(writer)?;

        // Write cells
        let total_cell_data: usize = mesh.cells.iter().map(|cell| cell.len() + 1).sum();
        writeln!(writer, "CELLS {} {}", mesh.num_cells(), total_cell_data)?;
        for cell in &mesh.cells {
            write!(writer, "{}", cell.len())?;
            for &vertex_id in cell {
                write!(writer, " {vertex_id}")?;
            }
            writeln!(writer)?;
        }
        writeln!(writer)?;

        // Write cell types
        writeln!(writer, "CELL_TYPES {}", mesh.num_cells())?;
        for &cell_type in &mesh.cell_types {
            writeln!(writer, "{}", cell_type as u8)?;
        }

        // Write point data if available
        if !mesh.point_data.is_empty() {
            writeln!(writer)?;
            writeln!(writer, "POINT_DATA {}", mesh.num_points())?;

            for (name, data) in &mesh.point_data {
                writeln!(writer, "SCALARS {name} double 1")?;
                writeln!(writer, "LOOKUP_TABLE default")?;
                for value in data {
                    let v = value.to_f64().ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!("Failed to convert scalar '{name}' to f64 for VTK output"),
                        )
                    })?;
                    writeln!(writer, "{v}")?;
                }
            }
        }

        // Write cell data if available
        if !mesh.cell_data.is_empty() {
            writeln!(writer)?;
            writeln!(writer, "CELL_DATA {}", mesh.num_cells())?;

            for (name, data) in &mesh.cell_data {
                writeln!(writer, "SCALARS {name} double 1")?;
                writeln!(writer, "LOOKUP_TABLE default")?;
                for value in data {
                    let v = value.to_f64().ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            format!(
                                "Failed to convert cell scalar '{name}' to f64 for VTK output"
                            ),
                        )
                    })?;
                    writeln!(writer, "{v}")?;
                }
            }
        }

        Ok(())
    }
}

impl<T: RealField + Copy + ToPrimitive> Default for VtkWriter<T> {
    fn default() -> Self {
        Self::new()
    }
}
