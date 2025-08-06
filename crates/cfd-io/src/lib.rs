//! I/O operations and file format support for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod vtk;
pub mod csv;
pub mod json;
pub mod binary;
pub mod checkpoint;

pub use vtk::{VtkWriter, VtkReader};
pub use csv::{CsvWriter, CsvReader};
pub use json::{JsonConfig, JsonWriter, JsonReader};
pub use binary::{BinaryWriter, BinaryReader};
pub use checkpoint::{Checkpoint, CheckpointManager};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        vtk::VtkWriter,
        csv::CsvWriter,
        json::{JsonConfig, JsonWriter},
        checkpoint::CheckpointManager,
    };
}