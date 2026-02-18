//! Mesh I/O: STL, VTK, and scheme import.
//!
//! Each format is feature-gated to avoid pulling in unnecessary dependencies.

pub mod stl;

#[cfg(feature = "vtk-io")]
pub mod vtk;

#[cfg(feature = "scheme-io")]
pub mod scheme;
