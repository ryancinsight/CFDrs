//! VTK file format support.
//!
//! This module provides VTK (Visualization Toolkit) file format support
//! for CFD simulation data with zero-copy operations where possible.

pub mod types;
pub mod writer;
pub mod reader;
pub mod builder;

// Re-export main types for convenience
pub use types::{VtkDataType, VtkCellType, VtkMesh, VtkHeader};
pub use writer::VtkWriter;
pub use reader::VtkReader;
pub use builder::VtkMeshBuilder;