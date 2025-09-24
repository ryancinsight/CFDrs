//! Input/output functionality for CFD simulations.
//!
//! This crate provides file format support for reading and writing
//! CFD simulation data.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// I/O and file format allows
#![allow(clippy::similar_names)]           // File format variables often similar
#![allow(clippy::must_use_candidate)]      // I/O utilities typically return values that are used

pub mod binary;
pub mod checkpoint;
pub mod csv;
#[cfg(feature = "hdf5")]
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
