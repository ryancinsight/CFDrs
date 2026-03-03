//! Mesh I/O: STL, OBJ, PLY, 3MF, glTF/GLB, DXF, VTK, `OpenFOAM`, and scheme import.
//!
//! Each format is feature-gated to avoid pulling in unnecessary dependencies,
//! except for formats with no external dependencies (OBJ, PLY, glTF, DXF,
//! OpenFOAM, STL).

pub mod stl;

pub mod openfoam;

pub mod obj;

pub mod ply;

pub mod gltf_export;

pub mod dxf;

#[cfg(feature = "three-mf-io")]
pub mod three_mf;

#[cfg(feature = "vtk-io")]
pub mod vtk;

pub mod scheme;
