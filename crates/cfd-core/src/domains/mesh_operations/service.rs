//! Mesh operations service for coordinating mesh-related operations

use super::{Geometry, Mesh, MeshGeneration, MeshQuality, MeshRefinement};
use crate::error::Result;
use nalgebra::RealField;
use std::sync::Arc;

/// Service for coordinating mesh operations
pub struct MeshOperationsService<T: RealField + Copy> {
    /// Mesh generator
    pub generator: Option<Arc<dyn MeshGeneration<T>>>,
    /// Mesh refiner
    pub refiner: Option<Arc<dyn MeshRefinement<T>>>,
    /// Quality assessor
    pub quality_assessor: Option<Arc<dyn MeshQuality<T>>>,
    /// Geometry handler
    pub geometry: Option<Arc<dyn Geometry<T>>>,
}

impl<T: RealField + Copy> MeshOperationsService<T> {
    /// Create a new mesh operations service
    #[must_use]
    pub fn new() -> Self {
        Self {
            generator: None,
            refiner: None,
            quality_assessor: None,
            geometry: None,
        }
    }

    /// Set the mesh generator
    #[must_use]
    pub fn with_generator(mut self, generator: Arc<dyn MeshGeneration<T>>) -> Self {
        self.generator = Some(generator);
        self
    }

    /// Set the mesh refiner
    #[must_use]
    pub fn with_refiner(mut self, refiner: Arc<dyn MeshRefinement<T>>) -> Self {
        self.refiner = Some(refiner);
        self
    }

    /// Set the quality assessor
    #[must_use]
    pub fn with_quality_assessor(mut self, assessor: Arc<dyn MeshQuality<T>>) -> Self {
        self.quality_assessor = Some(assessor);
        self
    }

    /// Set the geometry handler
    #[must_use]
    pub fn with_geometry(mut self, geometry: Arc<dyn Geometry<T>>) -> Self {
        self.geometry = Some(geometry);
        self
    }

    /// Generate and validate a mesh
    ///
    /// # Errors
    /// Returns an error if no mesh generator is set, mesh generation fails, or validation fails
    pub fn generate_validated_mesh(&self, nx: usize, ny: usize, nz: usize) -> Result<Mesh<T>> {
        let generator = self.generator.as_ref().ok_or_else(|| {
            crate::error::Error::InvalidConfiguration("No mesh generator set".into())
        })?;

        let mesh = generator.generate_structured(nx, ny, nz)?;
        mesh.validate()?;
        Ok(mesh)
    }
}

impl<T: RealField + Copy> Default for MeshOperationsService<T> {
    fn default() -> Self {
        Self::new()
    }
}
