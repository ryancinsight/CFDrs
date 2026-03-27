//! Schwarz P sphere primitive — Schwarz Primitive TPMS clipped to a sphere.
//!
//! The Schwarz P surface (`cos kx + cos ky + cos kz = 0`) is the simplest
//! triply periodic minimal surface.  Its unit cell consists of bidirectionally
//! connected passages along all three Cartesian directions, giving it high
//! structural symmetry and strong isotropic stiffness.
//!
//! See [`crate::domain::geometry::tpms::schwarz_p`] for full mathematical
//! documentation and the zero-mean-curvature theorem.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, schwarz_p::SchwarzP, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Schwarz P sphere — Primitive TPMS lattice clipped to a sphere.
///
/// ## Parameters
///
/// | Field | Default | Meaning |
/// |-------|---------|---------|
/// | `radius` | 5.0 mm | Spherical clip envelope |
/// | `period` | 2.5 mm | P-surface unit-cell period (`k = 2π/period`) |
/// | `resolution` | 64 | Voxels per axis |
/// | `iso_value` | 0.0 | Level-set threshold |
///
/// ## Example
///
/// ```rust,ignore
/// use cfd_mesh::{SchwarzPSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = SchwarzPSphere::default().build().expect("schwarz P sphere");
/// ```
#[derive(Clone, Debug)]
pub struct SchwarzPSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// P-surface unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for SchwarzPSphere {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

impl PrimitiveMesh for SchwarzPSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(
            &SchwarzP,
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
    fn schwarz_p_sphere_builds_with_faces() {
        let mesh = SchwarzPSphere {
            resolution: 20,
            ..SchwarzPSphere::default()
        }
        .build()
        .expect("Schwarz P sphere");
        assert!(!mesh.faces.is_empty());
    }

    #[test]
    fn schwarz_p_sphere_invalid_radius_errors() {
        assert!(SchwarzPSphere {
            radius: 0.0,
            ..SchwarzPSphere::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn schwarz_p_sphere_low_resolution_errors() {
        assert!(SchwarzPSphere {
            resolution: 2,
            ..SchwarzPSphere::default()
        }
        .build()
        .is_err());
    }
}
