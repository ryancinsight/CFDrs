//! Schwarz D sphere primitive — Diamond TPMS clipped to a sphere.
//!
//! The Schwarz D (Diamond) surface is defined by `cos(kx)·cos(ky)·cos(kz) = 0`,
//! which is equivalent to the four-term sum used in the original Schwarz (1890)
//! description.  Its diamond-like bicontinuous network features tetrahedral
//! symmetry and extremely high surface area density — making it a prime
//! candidate for heat exchangers, filtration membranes, and tissue scaffolds.
//!
//! See [`crate::domain::geometry::tpms::schwarz_d`] for the full mathematical
//! derivation and proof of zero mean curvature.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, schwarz_d::SchwarzD, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Schwarz D sphere — Diamond TPMS lattice clipped to a sphere.
///
/// ## Parameters
///
/// | Field | Default | Meaning |
/// |-------|---------|---------|
/// | `radius` | 5.0 mm | Spherical clip envelope |
/// | `period` | 2.5 mm | D-surface unit-cell period (`k = 2π/period`) |
/// | `resolution` | 64 | Voxels per axis |
/// | `iso_value` | 0.0 | Level-set threshold |
///
/// ## Example
///
/// ```rust,ignore
/// use cfd_mesh::{SchwarzDSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = SchwarzDSphere::default().build().expect("schwarz D sphere");
/// ```
#[derive(Clone, Debug)]
pub struct SchwarzDSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// D-surface unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for SchwarzDSphere {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

impl PrimitiveMesh for SchwarzDSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(
            &SchwarzD,
            &TpmsParams {
                radius: self.radius,
                period: self.period,
                resolution: self.resolution,
                iso_value: self.iso_value,
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schwarz_d_sphere_builds_with_faces() {
        let mesh = SchwarzDSphere {
            resolution: 20,
            ..SchwarzDSphere::default()
        }
        .build()
        .expect("Schwarz D sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn schwarz_d_sphere_invalid_radius_errors() {
        assert!(SchwarzDSphere { radius: -5.0, ..SchwarzDSphere::default() }.build().is_err());
    }

    #[test]
    fn schwarz_d_sphere_invalid_period_errors() {
        assert!(SchwarzDSphere { period: -1.0, ..SchwarzDSphere::default() }.build().is_err());
    }
}
