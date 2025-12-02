//! Parallel VTK I/O for distributed MPI simulations.
//!
//! This module provides VTK file writing capabilities for distributed
//! CFD simulations running on multiple MPI processes. It handles the
//! coordination required to write complete VTK files from distributed data.
//!
//! ## Architecture
//!
//! The parallel VTK writer follows these principles:
//! - **Root Process Coordination**: One process writes the file header and metadata
//! - **Distributed Data Collection**: Each process contributes its local data
//! - **Collective Operations**: MPI collective operations for efficient data gathering
//! - **Format Compatibility**: Output files are standard VTK format readable by ParaView
//!
//! ## Usage
//!
//! ```no_run
//! use cfd_io::vtk::ParallelVtkWriter;
//! use cfd_core::compute::mpi::DistributedVector;
//! use std::collections::HashMap;
//!
//! // Create writer
//! let writer = ParallelVtkWriter::new(&communicator)?;
//!
//! // Prepare distributed data
//! let mut point_data = HashMap::new();
//! point_data.insert("velocity".to_string(), &velocity_field);
//! point_data.insert("pressure".to_string(), &pressure_field);
//!
//! // Write parallel VTK file
//! writer.write_vtk_file(
//!     "simulation_output.vtk",
//!     &mesh_points,
//!     &cell_connectivity,
//!     &cell_types,
//!     &point_data,
//!     &HashMap::new(), // cell_data
//! )?;
//! ```

#[cfg(feature = "mpi")]
use std::fs::File;
#[cfg(feature = "mpi")]
use std::io::{BufWriter, Write};
#[cfg(feature = "mpi")]
use std::path::{Path, PathBuf};

#[cfg(feature = "mpi")]
use cfd_core::compute::mpi::{DistributedVector, MpiCommunicator};
#[cfg(feature = "mpi")]
use nalgebra::RealField;
#[cfg(feature = "mpi")]
use std::collections::HashMap;

/// I/O result type
pub type IoResult<T> = Result<T, IoError>;

/// I/O error type
#[derive(Debug)]
pub enum IoError {
    /// File I/O error
    FileError(std::io::Error),
    /// Data format error
    FormatError(String),
    /// MPI communication error
    #[cfg(feature = "mpi")]
    MpiError(cfd_core::compute::mpi::MpiError),
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IoError::FileError(e) => write!(f, "File I/O error: {e}"),
            IoError::FormatError(msg) => write!(f, "Format error: {msg}"),
            #[cfg(feature = "mpi")]
            IoError::MpiError(e) => write!(f, "MPI error: {:?}", e),
        }
    }
}

impl std::error::Error for IoError {}

impl From<std::io::Error> for IoError {
    fn from(error: std::io::Error) -> Self {
        IoError::FileError(error)
    }
}

#[cfg(feature = "mpi")]
impl From<cfd_core::compute::mpi::MpiError> for IoError {
    fn from(error: cfd_core::compute::mpi::MpiError) -> Self {
        IoError::MpiError(error)
    }
}

/// Parallel VTK writer for distributed CFD data
#[cfg(feature = "mpi")]
pub struct ParallelVtkWriter<T: RealField> {
    /// MPI communicator
    communicator: MpiCommunicator,
    /// Whether this process is the root (rank 0)
    is_root: bool,
    /// Total number of processes
    num_procs: i32,
    /// This process rank
    rank: i32,
}

#[cfg(feature = "mpi")]
impl<T: RealField + std::fmt::LowerExp> ParallelVtkWriter<T> {
    /// Create a new parallel VTK writer
    pub fn new(communicator: &MpiCommunicator) -> IoResult<Self> {
        let rank = communicator.rank();
        let num_procs = communicator.size();

        Ok(Self {
            communicator: communicator.clone(),
            is_root: rank == 0i32,
            num_procs,
            rank,
        })
    }

    /// Write a distributed VTK unstructured grid file
    ///
    /// This writes a .pvtu file on the root process and .vtu files on all processes.
    pub fn write_vtk_file<P: AsRef<Path>>(
        &self,
        filename: P,
        points: &DistributedVector<T>,
        cells: &[u32],
        cell_types: &[u8],
        point_data: &HashMap<String, &DistributedVector<T>>,
        cell_data: &HashMap<String, &DistributedVector<T>>,
    ) -> IoResult<()> {
        let path = filename.as_ref();
        let file_stem = path.file_stem().unwrap_or_default().to_string_lossy();
        let parent = path.parent().unwrap_or_else(|| Path::new("."));

        // 1. Write local .vtu file
        let local_filename = format!("{}_{}.vtu", file_stem, self.rank);
        let local_path = parent.join(&local_filename);

        self.write_local_vtu(
            &local_path,
            points,
            cells,
            cell_types,
            point_data,
            cell_data,
        )?;

        // 2. Root writes .pvtu file
        if self.is_root {
            let pvtu_filename = if path.extension().map_or(false, |ext| ext == "pvtu") {
                path.to_path_buf()
            } else {
                parent.join(format!("{}.pvtu", file_stem))
            };

            self.write_pvtu(
                &pvtu_filename,
                &file_stem,
                point_data,
                cell_data
            )?;
        }

        // Barrier to ensure all files are written before returning
        self.communicator.barrier();

        Ok(())
    }

    /// Write local data to a .vtu file (XML UnstructuredGrid)
    fn write_local_vtu(
        &self,
        path: &Path,
        points: &DistributedVector<T>,
        cells: &[u32],
        cell_types: &[u8],
        point_data: &HashMap<String, &DistributedVector<T>>,
        cell_data: &HashMap<String, &DistributedVector<T>>,
    ) -> IoResult<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "<?xml version=\"1.0\"?>")?;
        writeln!(writer, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">")?;
        writeln!(writer, "  <UnstructuredGrid>")?;

        let num_points = points.local_data.len() / 3;
        let num_cells = cell_types.len();

        writeln!(writer, "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", num_points, num_cells)?;

        // Points
        writeln!(writer, "      <Points>")?;
        writeln!(writer, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")?;
        for i in 0..num_points {
            writeln!(writer, "          {:e} {:e} {:e}", points.local_data[3*i], points.local_data[3*i+1], points.local_data[3*i+2])?;
        }
        writeln!(writer, "        </DataArray>")?;
        writeln!(writer, "      </Points>")?;

        // Cells
        writeln!(writer, "      <Cells>")?;

        // Connectivity
        writeln!(writer, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")?;
        let mut idx = 0;
        let mut offsets = Vec::with_capacity(num_cells);
        let mut current_offset = 0;

        while idx < cells.len() {
            let n = cells[idx] as usize;
            idx += 1; // Skip count
            for _ in 0..n {
                write!(writer, "{} ", cells[idx])?;
                idx += 1;
            }
            writeln!(writer)?;
            current_offset += n as i32;
            offsets.push(current_offset);
        }
        writeln!(writer, "        </DataArray>")?;

        // Offsets
        writeln!(writer, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")?;
        for offset in &offsets {
            writeln!(writer, "{}", offset)?;
        }
        writeln!(writer, "        </DataArray>")?;

        // Types
        writeln!(writer, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")?;
        for &t in cell_types {
            writeln!(writer, "{}", t)?;
        }
        writeln!(writer, "        </DataArray>")?;
        writeln!(writer, "      </Cells>")?;

        // Point Data
        if !point_data.is_empty() {
             writeln!(writer, "      <PointData>")?;
             for (name, data) in point_data {
                 writeln!(writer, "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">", name)?;
                 for val in data.local_data.iter() {
                     writeln!(writer, "          {:e}", val)?;
                 }
                 writeln!(writer, "        </DataArray>")?;
             }
             writeln!(writer, "      </PointData>")?;
        }

        // Cell Data
        if !cell_data.is_empty() {
             writeln!(writer, "      <CellData>")?;
             for (name, data) in cell_data {
                 writeln!(writer, "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">", name)?;
                 for val in data.local_data.iter() {
                     writeln!(writer, "          {:e}", val)?;
                 }
                 writeln!(writer, "        </DataArray>")?;
             }
             writeln!(writer, "      </CellData>")?;
        }

        writeln!(writer, "    </Piece>")?;
        writeln!(writer, "  </UnstructuredGrid>")?;
        writeln!(writer, "</VTKFile>")?;

        Ok(())
    }

    /// Write PVTU header on root
    fn write_pvtu(
        &self,
        path: &Path,
        file_stem: &str,
        point_data: &HashMap<String, &DistributedVector<T>>,
        cell_data: &HashMap<String, &DistributedVector<T>>,
    ) -> IoResult<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        writeln!(writer, "<?xml version=\"1.0\"?>")?;
        writeln!(writer, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">")?;
        writeln!(writer, "  <PUnstructuredGrid GhostLevel=\"0\">")?;

        // PPoints
        writeln!(writer, "    <PPoints>")?;
        writeln!(writer, "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>")?;
        writeln!(writer, "    </PPoints>")?;

        // PCells
        writeln!(writer, "    <PCells>")?;
        writeln!(writer, "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>")?;
        writeln!(writer, "      <PDataArray type=\"Int32\" Name=\"offsets\"/>")?;
        writeln!(writer, "      <PDataArray type=\"UInt8\" Name=\"types\"/>")?;
        writeln!(writer, "    </PCells>")?;

        // PPointData
        if !point_data.is_empty() {
             writeln!(writer, "    <PPointData>")?;
             for (name, _) in point_data {
                 writeln!(writer, "      <PDataArray type=\"Float64\" Name=\"{}\"/>", name)?;
             }
             writeln!(writer, "    </PPointData>")?;
        }

        // PCellData
        if !cell_data.is_empty() {
             writeln!(writer, "    <PCellData>")?;
             for (name, _) in cell_data {
                 writeln!(writer, "      <PDataArray type=\"Float64\" Name=\"{}\"/>", name)?;
             }
             writeln!(writer, "    </PCellData>")?;
        }

        // Pieces
        for r in 0..self.num_procs {
            writeln!(writer, "    <Piece Source=\"{}_{}.vtu\"/>", file_stem, r)?;
        }

        writeln!(writer, "  </PUnstructuredGrid>")?;
        writeln!(writer, "</VTKFile>")?;

        Ok(())
    }
}
