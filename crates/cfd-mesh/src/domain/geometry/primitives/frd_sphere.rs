//! FRD sphere primitive — Schoen FRD TPMS clipped to a sphere.
//!
//! The FRD ("Face-centred, Rhombic Dodecahedron") surface (Schoen 1970) has
//! FCC topology with two congruent interlocking labyrinths.  Its wider channel
//! junctions make it attractive for high-flow filtration membranes and fuel-cell
//! electrode scaffolds.
//!
//! See [`crate::domain::geometry::tpms::frd`] for mathematical proof of
//! minimality and full partial derivative derivations.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, frd::Frd, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// FRD sphere — Schoen F-RD TPMS lattice clipped to a sphere.
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
/// use cfd_mesh::{FrdSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = FrdSphere::default().build().expect("FRD sphere");
/// ```
#[derive(Clone, Debug)]
pub struct FrdSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for FrdSphere {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

impl PrimitiveMesh for FrdSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(
            &Frd,
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
    fn frd_sphere_builds_with_faces() {
        let mesh = FrdSphere {
            resolution: 20,
            ..FrdSphere::default()
        }
        .build()
        .expect("FRD sphere");
        assert!(!mesh.faces.is_empty());
    }

    #[test]
    fn frd_sphere_invalid_radius_errors() {
        assert!(FrdSphere {
            radius: 0.0,
            ..FrdSphere::default()
        }
        .build()
        .is_err());
    }

    #[test]
    fn frd_sphere_low_resolution_errors() {
        assert!(FrdSphere {
            resolution: 2,
            ..FrdSphere::default()
        }
        .build()
        .is_err());
    }
}
