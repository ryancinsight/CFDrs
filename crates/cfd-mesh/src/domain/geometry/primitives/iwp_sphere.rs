//! I-WP sphere primitive — I-WP TPMS clipped to a sphere.
//!
//! The I-WP ("I-graph Wrapped Package", Schoen 1970) is a BCC-topology TPMS
//! with genus 4 per conventional unit cell.  Its two non-congruent interlocking
//! labyrinths make it attractive for bone scaffold and heat-exchanger design.
//!
//! See [`crate::domain::geometry::tpms::iwp`] for mathematical proof of
//! minimality and full partial derivative derivations.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, iwp::Iwp, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// I-WP sphere — I-WP TPMS lattice clipped to a sphere.
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
/// use cfd_mesh::{IwpSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = IwpSphere::default().build().expect("IWP sphere");
/// ```
#[derive(Clone, Debug)]
pub struct IwpSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for IwpSphere {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

impl PrimitiveMesh for IwpSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(
            &Iwp,
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
    fn iwp_sphere_builds_with_faces() {
        let mesh = IwpSphere {
            resolution: 20,
            ..IwpSphere::default()
        }
        .build()
        .expect("IWP sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn iwp_sphere_invalid_radius_errors() {
        assert!(IwpSphere {
            radius: -2.0,
            ..IwpSphere::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn iwp_sphere_low_resolution_errors() {
        assert!(IwpSphere {
            resolution: 2,
            ..IwpSphere::default()
        }
        .build()
        .is_err());
    }
}
