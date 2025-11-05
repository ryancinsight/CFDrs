//! VTK file format support.
//!
//! This module provides VTK (Visualization Toolkit) file format support
//! for CFD simulation data with zero-copy operations where possible.
//!
//! ## Parallel VTK Support
//!
//! For distributed MPI simulations, use `ParallelVtkWriter` to write
//! VTK files across multiple processes:
//!
//! ```no_run
//! use cfd_io::vtk::ParallelVtkWriter;
//! use cfd_core::compute::mpi::*;
//!
//! let writer = ParallelVtkWriter::new(&communicator)?;
//! writer.write_vtk_file(
//!     "output.vtk",
//!     &distributed_points,
//!     &cells,
//!     &cell_types,
//!     &point_data,
//!     &cell_data,
//! )?;
//! ```

pub mod builder;
pub mod parallel;
pub mod reader;
pub mod types;
pub mod writer;

// Re-export main types for convenience
pub use builder::VtkMeshBuilder;
#[cfg(feature = "mpi")]
pub use parallel::ParallelVtkWriter;
pub use reader::VtkReader;
pub use types::{VtkCellType, VtkDataType, VtkHeader, VtkMesh};
pub use writer::VtkWriter;
