//! Mesh operations domain - Geometry and discretization operations.
//!
//! This module encapsulates mesh-related knowledge following DDD principles.
//! It provides abstractions for mesh generation, refinement, and quality assessment.

use crate::error::Result;

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
    /// 
    /// # Errors
    /// Returns an error if mesh generation fails due to invalid dimensions or memory constraints
    fn generate_structured(&self, nx: usize, ny: usize, nz: usize) -> Result<Mesh<T>>;

    /// Generate an unstructured mesh
    /// 
    /// # Errors
    /// Returns an error if mesh generation fails due to invalid target size or algorithmic constraints
    fn generate_unstructured(&self, target_size: T) -> Result<Mesh<T>>;

    /// Generate a hybrid mesh
    /// 
    /// # Errors
    /// Returns an error if hybrid mesh generation fails due to incompatible mesh types
    fn generate_hybrid(&self) -> Result<Mesh<T>>;
}
