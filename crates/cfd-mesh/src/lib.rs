//! Mesh handling and CSGrs integration for CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod mesh;
pub mod grid;
pub mod quality;
pub mod refinement;

#[cfg(feature = "csg")]
pub mod csg_integration;

pub use mesh::{Mesh, MeshTopology, Cell, Face, Edge, Vertex};
pub use grid::{StructuredGrid, UnstructuredGrid, GridGenerator};
pub use quality::{MeshQuality, QualityMetric};
pub use refinement::{MeshRefinement, AdaptiveRefinement};

/// Prelude module for convenient imports
pub mod prelude {
    pub use crate::{
        mesh::{Mesh, Cell, Face, Vertex},
        grid::{StructuredGrid, UnstructuredGrid},
        quality::MeshQuality,
        refinement::MeshRefinement,
    };
    
    #[cfg(feature = "csg")]
    pub use crate::csg_integration::CsgMeshGenerator;
}