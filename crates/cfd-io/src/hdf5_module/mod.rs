//! HDF5 support for large dataset I/O operations.
//!
//! This module provides structured I/O capabilities for CFD datasets using HDF5,
//! with support for chunking, metadata, and streaming operations.
//!
//! ## Parallel HDF5 Support
//!
//! For distributed MPI simulations, use `ParallelHdf5Writer` for efficient
//! collective I/O operations:
//!
//! ```no_run
//! use cfd_io::hdf5_module::ParallelHdf5Writer;
//! use cfd_core::compute::mpi::*;
//!
//! let writer = ParallelHdf5Writer::new(&communicator)?;
//!
//! // Prepare metadata and datasets
//! let mut metadata = std::collections::HashMap::new();
//! metadata.insert("simulation_time".to_string(), "1.5".to_string());
//!
//! let mut datasets = std::collections::HashMap::new();
//! datasets.insert("velocity".to_string(), &velocity_field);
//! datasets.insert("pressure".to_string(), &pressure_field);
//!
//! // Write parallel HDF5 file
//! writer.write_hdf5_file(
//!     "simulation_output.h5",
//!     &datasets,
//!     &metadata,
//! )?;
//! ```

pub mod chunking;
pub mod metadata;
pub mod parallel;
pub mod reader;
pub mod writer;

// Re-export main types
pub use chunking::{ChunkingStrategy, DataChunk};
pub use metadata::DatasetMetadata;
#[cfg(feature = "mpi")]
pub use parallel::ParallelHdf5Writer;
pub use reader::Hdf5Reader;
pub use writer::Hdf5Writer;
