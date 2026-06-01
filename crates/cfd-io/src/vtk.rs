//! VTK creation and I/O, provided by the `ritk-vtk` crate.
//!
//! VTK is owned by the RITK toolkit (single source of truth); CFDrs does not
//! maintain its own VTK implementation. This module re-exports `ritk-vtk`'s
//! data model (image / poly / structured / unstructured grids, attributes,
//! cell types) and its format writers/readers (legacy VTK, VTU/VTP/VTI, glTF,
//! STL/OBJ/PLY) so they are reachable as `cfd_io::vtk::*`.
//!
//! Available behind the `vtk` feature, which pulls `ritk-vtk` from its remote.

pub use ritk_vtk::*;
