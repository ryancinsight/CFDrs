//! Mesh I/O: STL, VTK, OpenFOAM, and scheme import.
//!
//! Each format is feature-gated to avoid pulling in unnecessary dependencies,
//! except for the OpenFOAM writer which has no external dependencies.

pub mod stl;

pub mod openfoam;

#[cfg(feature = "vtk-io")]
pub mod vtk;

#[cfg(feature = "scheme-io")]
pub mod scheme;
