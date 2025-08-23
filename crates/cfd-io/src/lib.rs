//! Input/output functionality for CFD simulations.
//!
//! This crate provides file format support for reading and writing
//! CFD simulation data.

#![warn(missing_docs)]

pub mod vtk;

pub use vtk::{VtkWriter, VtkReader, VtkMesh, VtkMeshBuilder};