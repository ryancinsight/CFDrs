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

// Re-export commonly used types
pub use csv::{CsvReader, CsvWriter, TimeSeriesData, StreamReader, StreamWriter};
pub use json::{JsonReader, JsonWriter};
pub use vtk::{VtkReader, VtkWriter, VtkMesh, VtkMeshBuilder, VtkFieldData, VtkDataType, VtkCellType, VtkDatasetType};
pub use binary::{BinaryReader, BinaryWriter};
pub use checkpoint::{Checkpoint, CheckpointManager};

/// Common I/O types and traits
pub mod prelude {
    pub use crate::{
        csv::{CsvReader, CsvWriter, TimeSeriesData},
        json::{JsonReader, JsonWriter},
        vtk::{VtkReader, VtkWriter, VtkMesh, VtkMeshBuilder},
        binary::{BinaryReader, BinaryWriter},
        checkpoint::CheckpointManager,
    };
}