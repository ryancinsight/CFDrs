//! I/O operations for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod binary;
pub mod csv;
pub mod formats;
pub mod json;
pub mod vtk;
pub mod checkpoint;
#[cfg(feature = "hdf5")]
pub mod hdf5;

// Re-export commonly used types
pub use csv::{CsvReader, CsvWriter, TimeSeriesData, StreamReader, StreamWriter};
pub use json::{JsonReader, JsonWriter};
pub use vtk::{VtkReader, VtkWriter, VtkMesh, VtkMeshBuilder, VtkDataType, VtkCellType, VtkDatasetType};
pub use binary::{BinaryReader, BinaryWriter};
pub use checkpoint::{Checkpoint, CheckpointManager};
#[cfg(feature = "hdf5")]
pub use hdf5::{Hdf5Reader, Hdf5Writer, DatasetMetadata, DataChunk};

// Prelude removed - use cfd_suite::prelude::* for unified SSOT interface