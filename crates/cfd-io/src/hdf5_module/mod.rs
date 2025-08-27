//! HDF5 support for large dataset I/O operations.
//!
//! This module provides structured I/O capabilities for CFD datasets using HDF5,
//! with support for chunking, metadata, and streaming operations.

pub mod chunking;
pub mod metadata;
pub mod reader;
pub mod writer;

// Re-export main types
pub use chunking::{ChunkingStrategy, DataChunk};
pub use metadata::DatasetMetadata;
pub use reader::Hdf5Reader;
pub use writer::Hdf5Writer;
