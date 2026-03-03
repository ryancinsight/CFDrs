//! Lidinoid sphere primitive — Lidinoid TPMS clipped to a sphere.
//!
//! The Lidinoid (Lidin & Larsson 1990) is a chiral TPMS in space group I4₁32
//! with genus 3 per conventional unit cell.  It belongs to the same Bonnet
//! family as the Gyroid (association angle ≈ 38.015°).  Its chirality makes it
//! useful for photonic crystals and drug-delivery scaffolds.
//!
//! See [`crate::domain::geometry::tpms::lidinoid`] for mathematical proof of
//! minimality and full partial derivative derivations.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, lidinoid::Lidinoid, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Lidinoid sphere — Lidinoid TPMS lattice clipped to a sphere.
///
/// ## Parameters
///
/// | Field | Default | Meaning |
/// |-------|---------|---------|
/// | `radius` | 5.0 mm | Spherical clip envelope |
/// | `period` | 2.5 mm | Unit-cell period (`k = 2π/period`) |
/// | `resolution` | 64 | Voxels per axis |
/// | `iso_value` | 0.0 | Level-set threshold |
///
/// ## Example
///
/// ```rust,ignore
/// use cfd_mesh::{LidinoidSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = LidinoidSphere::default().build().expect("lidinoid sphere");
/// ```
#[derive(Clone, Debug)]
pub struct LidinoidSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for LidinoidSphere {
    fn default() -> Self {
        Self { radius: 5.0, period: 2.5, resolution: 64, iso_value: 0.0 }
    }
}

impl PrimitiveMesh for LidinoidSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(&Lidinoid, &TpmsParams {
            radius: self.radius, period: self.period,
            resolution: self.resolution, iso_value: self.iso_value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lidinoid_sphere_builds_with_faces() {
        let mesh = LidinoidSphere { resolution: 20, ..LidinoidSphere::default() }
            .build().expect("Lidinoid sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn lidinoid_sphere_invalid_radius_errors() {
        assert!(LidinoidSphere { radius: 0.0, ..LidinoidSphere::default() }.build().is_err());
    }

    #[test]
    fn lidinoid_sphere_low_resolution_errors() {
        assert!(LidinoidSphere { resolution: 1, ..LidinoidSphere::default() }.build().is_err());
    }
}
