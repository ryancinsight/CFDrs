//! Mesh handling and generation.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod connectivity;
pub mod geometry;
pub mod grid;
pub mod mesh;
pub mod quality;
pub mod refinement;
pub mod csg;

// Re-export commonly used types
pub use connectivity::Connectivity;
pub use geometry::Geometry;
pub use mesh::{Cell, Edge, Face, Mesh, MeshTopology, Vertex};
pub use csg::CsgMeshAdapter;

/// Common mesh types and traits
pub mod prelude {
    pub use crate::{
        connectivity::Connectivity,
        geometry::Geometry,
        mesh::{Cell, Edge, Face, Mesh, MeshTopology, Vertex},
        csg::CsgMeshAdapter,
    };
}