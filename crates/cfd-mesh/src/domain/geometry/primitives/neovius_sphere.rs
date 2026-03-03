//! Neovius sphere primitive — Neovius TPMS clipped to a sphere.
//!
//! The Neovius surface (`3(cos kx+cos ky+cos kz)+4cos(kx)cos(ky)cos(kz) = 0`)
//! has cubic symmetry (Im3̄m) with genus 9 per unit cell.  Its bicontinuous
//! network of two congruent labyrinths makes it attractive for filtration
//! membranes requiring high surface-area density.
//!
//! See [`crate::domain::geometry::tpms::neovius`] for mathematical proofs and
//! partial derivatives.

use super::{PrimitiveError, PrimitiveMesh};
use crate::domain::geometry::tpms::{build_tpms_sphere, neovius::Neovius, TpmsParams};
use crate::domain::mesh::IndexedMesh;

/// Neovius sphere — Neovius TPMS lattice clipped to a sphere.
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
/// use cfd_mesh::{NeoviusSphere, domain::geometry::primitives::PrimitiveMesh};
///
/// let mesh = NeoviusSphere::default().build().expect("neovius sphere");
/// ```
#[derive(Clone, Debug)]
pub struct NeoviusSphere {
    /// Clip-sphere radius [mm].
    pub radius: f64,
    /// Unit-cell period [mm].
    pub period: f64,
    /// Voxels per axis.
    pub resolution: usize,
    /// Level-set threshold.
    pub iso_value: f64,
}

impl Default for NeoviusSphere {
    fn default() -> Self {
        Self { radius: 5.0, period: 2.5, resolution: 64, iso_value: 0.0 }
    }
}

impl PrimitiveMesh for NeoviusSphere {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build_tpms_sphere(&Neovius, &TpmsParams {
            radius: self.radius, period: self.period,
            resolution: self.resolution, iso_value: self.iso_value,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn neovius_sphere_builds_with_faces() {
        let mesh = NeoviusSphere { resolution: 20, ..NeoviusSphere::default() }
            .build().expect("Neovius sphere");
        assert!(mesh.faces.len() > 0);
    }

    #[test]
    fn neovius_sphere_invalid_radius_errors() {
        assert!(NeoviusSphere { radius: -1.0, ..NeoviusSphere::default() }.build().is_err());
    }

    #[test]
    fn neovius_sphere_low_resolution_errors() {
        assert!(NeoviusSphere { resolution: 3, ..NeoviusSphere::default() }.build().is_err());
    }
}
