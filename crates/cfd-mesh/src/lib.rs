//! Mesh handling and generation.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod connectivity;
#[cfg(feature = "csg")]
pub mod csg;
pub mod error;
pub mod geometry;
pub mod grid;
pub mod mesh;
pub mod quality;
pub mod refinement;

// Re-export commonly used types
pub use connectivity::Connectivity;
#[cfg(feature = "csg")]
pub use csg::CsgMeshAdapter;
pub use geometry::Geometry;
pub use mesh::{Cell, Edge, ElementType, Face, Mesh, MeshTopology, Vertex};

/// Common mesh types and traits
pub mod prelude {
    #[cfg(feature = "csg")]
    pub use crate::csg::CsgMeshAdapter;
    pub use crate::{
        connectivity::Connectivity,
        geometry::Geometry,
        mesh::{Cell, Edge, Face, Mesh, MeshTopology, Vertex},
    };
}
