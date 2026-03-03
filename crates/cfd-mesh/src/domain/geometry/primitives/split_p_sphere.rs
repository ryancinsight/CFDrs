//! Split P sphere primitive — Split P TPMS clipped to a sphere.
//!
//! The Split P surface (Fogden & Hyde 1992) is a tetragonally distorted variant
//! of Schwarz P in the rPD Bonnet family.  Its anisotropic channel cross-sections
//! make it useful for microfluidic membranes requiring directional permeability.
//!
//! See [`crate::domain::geometry::tpms::split_p`] for mathematical proof of
//! minimality via Bonnet isometry and full partial derivative derivations.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, split_p::SplitP, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Split P sphere — Split P TPMS lattice clipped to a sphere.
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
/// use cfd_mesh::{SplitPSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = SplitPSphere::default().build().expect("split P sphere");
/// ```
#[derive(Clone, Debug)]
pub struct SplitPSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for SplitPSphere {
    fn default() -> Self {
        Self { radius: 5.0, period: 2.5, resolution: 64, iso_value: 0.0 }
    }
}

impl PrimitiveMesh for SplitPSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(&SplitP, &TpmsParams {
            radius: self.radius, period: self.period,
            resolution: self.resolution, iso_value: self.iso_value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_p_sphere_builds_with_faces() {
        let mesh = SplitPSphere { resolution: 20, ..SplitPSphere::default() }
            .build().expect("Split P sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn split_p_sphere_invalid_radius_errors() {
        assert!(SplitPSphere { radius: -0.5, ..SplitPSphere::default() }.build().is_err());
    }

    #[test]
    fn split_p_sphere_low_resolution_errors() {
        assert!(SplitPSphere { resolution: 3, ..SplitPSphere::default() }.build().is_err());
    }
}
