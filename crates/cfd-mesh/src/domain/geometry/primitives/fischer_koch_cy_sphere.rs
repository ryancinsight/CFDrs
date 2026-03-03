//! Fischer-Koch C(Y) sphere primitive — C(Y) TPMS clipped to a sphere.
//!
//! The Fischer-Koch C(Y) surface (Fischer & Koch 1989) is a cubic TPMS in
//! space group *Ia*3̄*d* (#230) with genus 6 per unit cell.  Its dual
//! labyrinthine network of Y-shaped triple junctions makes it attractive for
//! triply-connected microfluidic networks and drug delivery scaffolds.
//!
//! See [`crate::domain::geometry::tpms::fischer_koch_cy`] for mathematical proof
//! of minimality via Weierstrass–Enneper construction and full partial derivatives.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, fischer_koch_cy::FischerKochCY, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Fischer-Koch C(Y) sphere — C(Y) TPMS lattice clipped to a sphere.
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
/// use cfd_mesh::{FischerKochCySphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = FischerKochCySphere::default().build().expect("Fischer-Koch C(Y) sphere");
/// ```
#[derive(Clone, Debug)]
pub struct FischerKochCySphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for FischerKochCySphere {
    fn default() -> Self {
        Self { radius: 5.0, period: 2.5, resolution: 64, iso_value: 0.0 }
    }
}

impl PrimitiveMesh for FischerKochCySphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(&FischerKochCY, &TpmsParams {
            radius: self.radius, period: self.period,
            resolution: self.resolution, iso_value: self.iso_value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fischer_koch_cy_sphere_builds_with_faces() {
        let mesh = FischerKochCySphere { resolution: 20, ..FischerKochCySphere::default() }
            .build().expect("Fischer-Koch C(Y) sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn fischer_koch_cy_sphere_invalid_radius_errors() {
        assert!(FischerKochCySphere { radius: -1.0, ..FischerKochCySphere::default() }.build().is_err());
    }

    #[test]
    fn fischer_koch_cy_sphere_low_resolution_errors() {
        assert!(FischerKochCySphere { resolution: 3, ..FischerKochCySphere::default() }.build().is_err());
    }
}
