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

// Re-export commonly used types
pub use connectivity::Connectivity;
pub use geometry::Geometry;
pub use mesh::{Cell, Edge, Face, Mesh, MeshTopology, Vertex};
// TODO: Implement these types
// pub use grid::{StructuredGrid, UnstructuredGrid, GridGenerator};
// pub use quality::{MeshQuality, QualityMetric};
// pub use refinement::{MeshRefinement, AdaptiveRefinement};

// #[cfg(feature = "csg")]
// pub mod csg;

/// Common mesh types and traits
pub mod prelude {
    pub use crate::{
        connectivity::Connectivity,
        geometry::Geometry,
        mesh::{Cell, Edge, Face, Mesh, MeshTopology, Vertex},
        // TODO: Add these when implemented
        // grid::{StructuredGrid, UnstructuredGrid},
        // quality::MeshQuality,
        // refinement::MeshRefinement,
    };

    // #[cfg(feature = "csg")]
    // pub use crate::csg::*;
}