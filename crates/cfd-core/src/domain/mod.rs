//! Domain module consolidating spatial representations and mesh operations.
//!
//! This module provides the core abstractions for the computational domain,
//! including geometric shapes, structured and unstructured meshes, and
//! tools for mesh manipulation and quality assessment.

pub mod mesh;
pub mod shapes;

pub use mesh::{
    Connectivity, Element, ElementType, Mesh, MeshGeneration, MeshMetadata, MeshOperationsService,
    MeshQuality, MeshRefinement, MeshStatistics, QualityReport, QualityStatistics,
    RefinementCriteria,
};
pub use shapes::{AnyDomain, Domain, Domain1D, Domain2D, Domain3D, Geometry};
