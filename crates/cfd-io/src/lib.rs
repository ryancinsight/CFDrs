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

// Re-export commonly used types
pub use csv::{CsvReader, CsvWriter};
pub use json::{JsonReader, JsonWriter};
pub use vtk::{VtkReader, VtkWriter};
pub use binary::{BinaryReader, BinaryWriter};

/// Common I/O types and traits
pub mod prelude {
    pub use crate::{
        csv::{CsvReader, CsvWriter},
        json::{JsonReader, JsonWriter},
        vtk::{VtkReader, VtkWriter},
        binary::{BinaryReader, BinaryWriter},
    };
}