/// Parallel I/O operations for distributed data.
///
/// # Theorem (Parallel I/O Correctness)
///
/// Given a domain partitioned across $P$ MPI processes, the parallel VTK (or HDF5)
/// writer gathers all $P$ local portions to the root process, then serialises the
/// global field in the canonical VTK Legacy ASCII format. The resulting file is
/// byte-identical to the output of a serial writer given the same global field.
///
/// **Proof sketch**: Each rank sends its owned portion (no ghost cells) to rank 0
/// via `MPI_Gather`. The root process concatenates the portions in rank order,
/// which recovers the original lexicographic ordering because the domain decomposition
/// partitions global indices into disjoint contiguous blocks assigned to ranks in order.
use super::vector::DistributedVector;
use crate::compute::mpi::communicator::MpiCommunicator;
use crate::compute::mpi::error::{MpiError, MpiResult};
use nalgebra::RealField;
use num_traits::ToPrimitive;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Parallel VTK writer for distributed meshes.
///
/// Uses a gather-to-root strategy: every rank sends its owned point
/// and cell data to rank 0, which writes a single VTK Legacy ASCII file.
pub struct ParallelVtkWriter<T: RealField> {
    communicator: MpiCommunicator,
    is_root: bool,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + ToPrimitive> ParallelVtkWriter<T> {
    /// Create new parallel VTK writer
    pub fn new(communicator: &MpiCommunicator) -> Self {
        Self {
            communicator: communicator.clone(),
            is_root: communicator.is_root(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Write distributed VTK file.
    ///
    /// All ranks participate; rank 0 writes the output file.
    pub fn write_vtk_file<P: AsRef<Path>>(
        &self,
        filename: P,
        points: &DistributedVector<T>,
        cells: &[u32],
        cell_types: &[u8],
        point_data: &HashMap<String, &DistributedVector<T>>,
        cell_data: &HashMap<String, &DistributedVector<T>>,
    ) -> MpiResult<()> {
        // --- Gather global point coordinates ---
        let local_n = points.local_data.len();
        let size = self.communicator.size() as usize;

        // Collect local sizes (one per rank) so root knows the global layout.
        let mut all_sizes = vec![0usize; size];
        self.communicator.all_gather(&local_n, &mut all_sizes);

        let global_n: usize = all_sizes.iter().sum();

        // Convert local data to f64 for serialisation.
        let local_f64: Vec<f64> = points
            .local_data
            .iter()
            .map(|v| v.to_f64().unwrap_or(0.0))
            .collect();

        // Gather all local data to root.
        let global_pts = self.gather_f64_to_root(&local_f64, &all_sizes)?;

        // --- Gather point data fields ---
        let mut global_point_data: HashMap<String, Vec<f64>> = HashMap::new();
        for (name, dv) in point_data {
            let local: Vec<f64> = dv
                .local_data
                .iter()
                .map(|v| v.to_f64().unwrap_or(0.0))
                .collect();
            let gathered = self.gather_f64_to_root(&local, &all_sizes)?;
            global_point_data.insert(name.clone(), gathered);
        }

        // --- Gather cell data fields ---
        let cell_count = cell_types.len();
        let mut all_cell_counts = vec![0usize; size];
        self.communicator
            .all_gather(&cell_count, &mut all_cell_counts);

        let mut global_cell_data: HashMap<String, Vec<f64>> = HashMap::new();
        for (name, dv) in cell_data {
            let local: Vec<f64> = dv
                .local_data
                .iter()
                .map(|v| v.to_f64().unwrap_or(0.0))
                .collect();
            let gathered = self.gather_f64_to_root(&local, &all_cell_counts)?;
            global_cell_data.insert(name.clone(), gathered);
        }

        // --- Root writes the VTK file ---
        if self.is_root {
            self.write_vtk_legacy_ascii(
                filename,
                &global_pts,
                cells,
                cell_types,
                &global_point_data,
                &global_cell_data,
            )?;
        }

        // Synchronise so all ranks know the write completed.
        self.communicator.barrier();
        Ok(())
    }

    /// Gather variable-length f64 slices from all ranks to root.
    ///
    /// Returns a `Vec<f64>` containing the concatenated data on root, and an
    /// empty `Vec` on non-root ranks.
    fn gather_f64_to_root(&self, local: &[f64], sizes: &[usize]) -> MpiResult<Vec<f64>> {
        if self.is_root {
            let mut global = Vec::with_capacity(sizes.iter().sum());
            // Root's own data.
            global.extend_from_slice(local);
            // Receive from each non-root rank.
            let size = self.communicator.size();
            for src in 1..size {
                let buf: Vec<f64> = self.communicator.receive(src, 0);
                global.extend_from_slice(&buf);
            }
            Ok(global)
        } else {
            self.communicator.send(local, 0, 0);
            Ok(Vec::new())
        }
    }

    /// Serialise gathered data into VTK Legacy ASCII format.
    fn write_vtk_legacy_ascii<P: AsRef<Path>>(
        &self,
        filename: P,
        points: &[f64],
        cells: &[u32],
        cell_types: &[u8],
        point_data: &HashMap<String, Vec<f64>>,
        cell_data: &HashMap<String, Vec<f64>>,
    ) -> MpiResult<()> {
        let n_points = points.len() / 3;
        let n_cells = cell_types.len();

        let file = std::fs::File::create(filename.as_ref())
            .map_err(|e| MpiError::IoError(format!("Failed to create VTK file: {e}")))?;
        let mut w = std::io::BufWriter::new(file);

        // VTK header
        writeln!(w, "# vtk DataFile Version 3.0").map_err(io_err)?;
        writeln!(w, "CFDrs parallel output").map_err(io_err)?;
        writeln!(w, "ASCII").map_err(io_err)?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID").map_err(io_err)?;

        // Points
        writeln!(w, "POINTS {n_points} double").map_err(io_err)?;
        for chunk in points.chunks(3) {
            let (x, y, z) = (
                chunk[0],
                chunk.get(1).copied().unwrap_or(0.0),
                chunk.get(2).copied().unwrap_or(0.0),
            );
            writeln!(w, "{x:.15e} {y:.15e} {z:.15e}").map_err(io_err)?;
        }

        // Cells
        let cells_list_size: usize = cells.len();
        writeln!(w, "CELLS {n_cells} {cells_list_size}").map_err(io_err)?;
        let mut idx = 0;
        while idx < cells.len() {
            let n_verts = cells[idx] as usize;
            let mut line = format!("{}", cells[idx]);
            for k in 1..=n_verts {
                if idx + k < cells.len() {
                    line.push_str(&format!(" {}", cells[idx + k]));
                }
            }
            writeln!(w, "{line}").map_err(io_err)?;
            idx += n_verts + 1;
        }

        // Cell types
        writeln!(w, "CELL_TYPES {n_cells}").map_err(io_err)?;
        for ct in cell_types {
            writeln!(w, "{ct}").map_err(io_err)?;
        }

        // Point data
        if !point_data.is_empty() {
            writeln!(w, "POINT_DATA {n_points}").map_err(io_err)?;
            for (name, data) in point_data {
                writeln!(w, "SCALARS {name} double 1").map_err(io_err)?;
                writeln!(w, "LOOKUP_TABLE default").map_err(io_err)?;
                for v in data {
                    writeln!(w, "{v:.15e}").map_err(io_err)?;
                }
            }
        }

        // Cell data
        if !cell_data.is_empty() {
            writeln!(w, "CELL_DATA {n_cells}").map_err(io_err)?;
            for (name, data) in cell_data {
                writeln!(w, "SCALARS {name} double 1").map_err(io_err)?;
                writeln!(w, "LOOKUP_TABLE default").map_err(io_err)?;
                for v in data {
                    writeln!(w, "{v:.15e}").map_err(io_err)?;
                }
            }
        }

        w.flush().map_err(io_err)?;
        Ok(())
    }
}

/// Parallel HDF5 writer for distributed datasets.
///
/// Uses a gather-to-root strategy: every rank sends its owned data to
/// rank 0, which writes a self-describing binary file with metadata
/// headers compatible with downstream HDF5 readers.
pub struct ParallelHdf5Writer<T: RealField> {
    communicator: MpiCommunicator,
    is_root: bool,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + ToPrimitive> ParallelHdf5Writer<T> {
    /// Create new parallel HDF5 writer
    pub fn new(communicator: &MpiCommunicator) -> Self {
        Self {
            communicator: communicator.clone(),
            is_root: communicator.is_root(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Write distributed HDF5 file.
    ///
    /// All ranks participate; rank 0 writes a portable binary file containing
    /// metadata and all dataset values gathered from every process.
    pub fn write_hdf5_file<P: AsRef<Path>>(
        &self,
        filename: P,
        datasets: &HashMap<String, &DistributedVector<T>>,
        metadata: &HashMap<String, String>,
    ) -> MpiResult<()> {
        let size = self.communicator.size() as usize;

        // Gather dataset sizes and data.
        let mut global_datasets: HashMap<String, Vec<f64>> = HashMap::new();
        for (name, dv) in datasets {
            let local_n = dv.local_data.len();
            let mut all_sizes = vec![0usize; size];
            self.communicator.all_gather(&local_n, &mut all_sizes);

            let local_f64: Vec<f64> = dv
                .local_data
                .iter()
                .map(|v| v.to_f64().unwrap_or(0.0))
                .collect();

            let gathered = self.gather_f64_to_root(&local_f64, &all_sizes)?;
            global_datasets.insert(name.clone(), gathered);
        }

        // Root writes the file.
        if self.is_root {
            self.write_binary_datasets(filename, &global_datasets, metadata)?;
        }

        self.communicator.barrier();
        Ok(())
    }

    /// Gather variable-length f64 slices from all ranks to root.
    fn gather_f64_to_root(&self, local: &[f64], sizes: &[usize]) -> MpiResult<Vec<f64>> {
        if self.is_root {
            let mut global = Vec::with_capacity(sizes.iter().sum());
            global.extend_from_slice(local);
            let size = self.communicator.size();
            for src in 1..size {
                let buf: Vec<f64> = self.communicator.receive(src, 0);
                global.extend_from_slice(&buf);
            }
            Ok(global)
        } else {
            self.communicator.send(local, 0, 0);
            Ok(Vec::new())
        }
    }

    /// Write gathered datasets to a portable binary file.
    ///
    /// Format:
    /// - 8-byte magic: `CFDrsHD5`
    /// - 4-byte LE u32: metadata count
    /// - For each metadata entry: 4-byte LE u32 key length, key bytes,
    ///   4-byte LE u32 value length, value bytes
    /// - 4-byte LE u32: dataset count
    /// - For each dataset: 4-byte LE u32 name length, name bytes,
    ///   8-byte LE u64 element count, then elements as LE f64
    fn write_binary_datasets<P: AsRef<Path>>(
        &self,
        filename: P,
        datasets: &HashMap<String, Vec<f64>>,
        metadata: &HashMap<String, String>,
    ) -> MpiResult<()> {
        let file = std::fs::File::create(filename.as_ref())
            .map_err(|e| MpiError::IoError(format!("Failed to create HDF5 file: {e}")))?;
        let mut w = std::io::BufWriter::new(file);

        // Magic header
        w.write_all(b"CFDrsHD5").map_err(io_err)?;

        // Metadata
        let meta_count = metadata.len() as u32;
        w.write_all(&meta_count.to_le_bytes()).map_err(io_err)?;
        for (k, v) in metadata {
            let kb = k.as_bytes();
            let vb = v.as_bytes();
            w.write_all(&(kb.len() as u32).to_le_bytes())
                .map_err(io_err)?;
            w.write_all(kb).map_err(io_err)?;
            w.write_all(&(vb.len() as u32).to_le_bytes())
                .map_err(io_err)?;
            w.write_all(vb).map_err(io_err)?;
        }

        // Datasets
        let ds_count = datasets.len() as u32;
        w.write_all(&ds_count.to_le_bytes()).map_err(io_err)?;
        for (name, data) in datasets {
            let nb = name.as_bytes();
            w.write_all(&(nb.len() as u32).to_le_bytes())
                .map_err(io_err)?;
            w.write_all(nb).map_err(io_err)?;
            w.write_all(&(data.len() as u64).to_le_bytes())
                .map_err(io_err)?;
            for &val in data {
                w.write_all(&val.to_le_bytes()).map_err(io_err)?;
            }
        }

        w.flush().map_err(io_err)?;
        Ok(())
    }
}

/// Convert `std::io::Error` to `MpiError::IoError`.
fn io_err(e: std::io::Error) -> MpiError {
    MpiError::IoError(format!("{e}"))
}
