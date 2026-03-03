//! Gyroid sphere primitive — delegates to the shared TPMS marching-cubes engine.
//!
//! The gyroid TPMS (Schoen 1970) implicit field and gradient are defined in
//! [`crate::domain::geometry::tpms::gyroid::Gyroid`].  This builder is a thin
//! adapter that exposes the `PrimitiveMesh` interface used by the rest of the
//! primitive library.
//!
//! See [`crate::domain::geometry::tpms`] for full mathematical documentation.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, gyroid::Gyroid, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Gyroid sphere — TPMS mid-surface lattice clipped to a sphere.
///
/// ## Parameters
///
/// | Field | Default | Meaning |
/// |-------|---------|---------|
/// | `radius` | 5.0 mm | Spherical clip envelope |
/// | `period` | 2.5 mm | Gyroid unit-cell period (`k = 2π/period`) |
/// | `resolution` | 64 | Voxels per axis (higher → denser and more accurate) |
/// | `iso_value` | 0.0 | Level-set threshold (0 = exact minimal mid-sheet) |
///
/// ## Example
///
/// ```rust,ignore
/// use cfd_mesh::{GyroidSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = GyroidSphere { radius: 5.0, period: 2.5, resolution: 64, iso_value: 0.0 }
///     .build()
///     .expect("gyroid sphere");
/// ```
#[derive(Clone, Debug)]
pub struct GyroidSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Gyroid unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for GyroidSphere {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

impl PrimitiveMesh for GyroidSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(
            &Gyroid,
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
    fn gyroid_sphere_builds_with_faces() {
        let mesh = GyroidSphere {
            resolution: 20,
            ..GyroidSphere::default()
        }
        .build()
        .expect("gyroid sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn gyroid_sphere_invalid_radius_errors() {
        assert!(GyroidSphere {
            radius: -1.0,
            ..GyroidSphere::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn gyroid_sphere_invalid_period_errors() {
        assert!(GyroidSphere {
            period: 0.0,
            ..GyroidSphere::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn gyroid_sphere_low_resolution_errors() {
        assert!(GyroidSphere {
            resolution: 3,
            ..GyroidSphere::default()
        }
        .build()
        .is_err());
    }
}
