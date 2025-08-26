//! VTK file format support.
//!
//! This module provides VTK (Visualization Toolkit) file format support
//! for CFD simulation data with zero-copy operations where possible.

pub mod builder;
pub mod reader;
pub mod types;
pub mod writer;
// Re-export main types for convenience
pub use builder::VtkMeshBuilder;
pub use reader::VtkReader;
pub use types::{VtkCellType, VtkDataType, VtkHeader, VtkMesh};
pub use writer::VtkWriter;
