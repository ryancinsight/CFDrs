//! Input/output functionality for CFD simulations.
//!
//! This crate provides file format support for reading and writing
//! CFD simulation data.

#![warn(missing_docs)]
// I/O and file format allows

mod leto_arrays;

pub mod binary;
pub mod checkpoint;
pub mod csv;
pub mod error;
#[cfg(feature = "hdf5")]
pub mod hdf5;
#[cfg(feature = "vtk")]
pub mod vtk;

// The API is now the public module hierarchy. This provides a clear,
// self-documenting structure for users.
// Example usage:
//   use cfd_io::vtk::VtkWriter;
//   use cfd_io::checkpoint::CheckpointManager;
//
// This hierarchical structure prevents namespace pollution and makes
// it clear which file format each type belongs to.
