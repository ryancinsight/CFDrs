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

// The API is now the public module hierarchy. This provides a clear,
// self-documenting structure for users.
// Example usage:
//   use cfd_io::vtk::VtkWriter;
//   use cfd_io::checkpoint::CheckpointManager;
//   use cfd_io::hdf5_module::Hdf5Reader;
//
// This hierarchical structure prevents namespace pollution and makes
// it clear which file format each type belongs to.
