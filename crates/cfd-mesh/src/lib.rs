//! Mesh handling and generation.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// Mesh generation and geometry allows
#![allow(clippy::similar_names)]           // Geometric variables (x,y,z; i,j,k) naturally similar
#![allow(clippy::cast_precision_loss)]     // Acceptable in geometric calculations
#![allow(clippy::must_use_candidate)]      // Geometric utilities often used in expressions

pub mod connectivity;
#[cfg(feature = "csg")]
pub mod csg;
pub mod error;
pub mod geometry;
pub mod grid;
pub mod mesh;
pub mod quality;
pub mod refinement;
pub mod topology;

// The public modules are the primary API.
// Users should access types through the module hierarchy:
//   use cfd_mesh::mesh::Mesh;
//   use cfd_mesh::grid::StructuredGridBuilder;
//   use cfd_mesh::quality::aspect_ratio;
// The prelude provides convenient shortcuts for common types.

/// Common mesh types and traits
pub mod prelude {
    #[cfg(feature = "csg")]
    pub use crate::csg::CsgError;
    pub use crate::{
        connectivity::Connectivity,
        geometry::Geometry,
        mesh::{Mesh, MeshStatistics},
        topology::{Cell, Edge, Face, Vertex},
    };
}
