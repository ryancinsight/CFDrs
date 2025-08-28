//! Input/output functionality for CFD simulations.
//!
//! This crate provides file format support for reading and writing
//! CFD simulation data.

#![warn(missing_docs)]

pub mod binary;
pub mod checkpoint;
pub mod csv;
pub mod hdf5_module;
pub mod vtk;

pub use checkpoint::{Checkpoint, CheckpointManager, CheckpointMetadata};
pub use csv::{CsvReader, CsvWriter};
pub use hdf5_module::{ChunkingStrategy, DataChunk, DatasetMetadata, Hdf5Reader, Hdf5Writer};
pub use vtk::{VtkMesh, VtkMeshBuilder, VtkReader, VtkWriter};
