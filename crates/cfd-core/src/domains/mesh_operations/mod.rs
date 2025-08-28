//! Mesh operations domain - Geometry and discretization operations.
//!
//! This module encapsulates mesh-related knowledge following DDD principles.
//! It provides abstractions for mesh generation, refinement, and quality assessment.

pub mod element;
pub mod geometry;
pub mod mesh;
pub mod quality;
pub mod refinement;
pub mod service;

pub use element::ElementType;
pub use geometry::Geometry;
pub use mesh::{Element, Mesh, MeshMetadata, MeshStatistics};
pub use quality::{MeshQuality, QualityReport, QualityStatistics};
pub use refinement::{MeshRefinement, RefinementCriteria};
pub use service::MeshOperationsService;

/// Mesh generation trait for creating meshes
pub trait MeshGeneration<T: nalgebra::RealField + Copy>: Send + Sync {
    /// Generate a structured mesh
    fn generate_structured(&self, nx: usize, ny: usize, nz: usize) -> crate::Result<Mesh<T>>;

    /// Generate an unstructured mesh
    fn generate_unstructured(&self, target_size: T) -> crate::Result<Mesh<T>>;

    /// Generate a hybrid mesh
    fn generate_hybrid(&self) -> crate::Result<Mesh<T>>;
}
