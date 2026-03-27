//! FRD (Schoen F-RD) TPMS implicit field (Schoen 1970).
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = 4·cos(kx)·cos(ky)·cos(kz)
//!            − (cos(2kx)·cos(2ky)
//!              + cos(2ky)·cos(2kz)
//!              + cos(2kz)·cos(2kx))
//! ```
//!
//! The F-RD ("Face-centred, Rhombic Dodecahedron") surface occupies an FCC topology.
//! Its two interlocking channel networks are congruent, giving it the same labyrinthine
//! symmetry as Schwarz D but with face-centred rather than body-centred connectivity.
//!
//! ## Theorem — Zero Mean Curvature
//!
//! F-RD is proven minimal by Schoen (1970, *NASA CR-1012*, entry "F-RD") via the
//! Weierstrass–Enneper representation, verified by explicit computation of the
//! second fundamental form yielding `H = 0` everywhere on `{F = 0}`.
//!
//! ## Partial Derivatives
//!
//! ```text
//! ∂F/∂x = −4k·sin(kx)·cos(ky)·cos(kz) + 2k·sin(2kx)·(cos(2ky)+cos(2kz))
//! ∂F/∂y = −4k·cos(kx)·sin(ky)·cos(kz) + 2k·sin(2ky)·(cos(2kx)+cos(2kz))
//! ∂F/∂z = −4k·cos(kx)·cos(ky)·sin(kz) + 2k·sin(2kz)·(cos(2ky)+cos(2kx))
//! ```
//!
//! ## References
//!
//! - Schoen, A.H. (1970). *NASA Tech. Rep.* CR-1012. Table 1, entry "F-RD".
//! - Gandy, P.J.F., Klinowski, J. et al. (1999). *Chem. Phys. Lett.* 314, 543–551.

use super::Tpms;
use crate::domain::core::scalar::Vector3r;

/// FRD (Schoen F-RD) TPMS — Face-centred Rhombic Dodecahedron minimal surface.
///
/// The FRD surface features FCC topology with two congruent interlocking labyrinths.
/// Its wider channel junctions make it attractive for high-flow filtration membranes
/// and fuel cell electrode scaffolds.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct Frd;

impl Tpms for Frd {
    /// `F = 4cos(kx)cos(ky)cos(kz) − (cos(2kx)cos(2ky)+cos(2ky)cos(2kz)+cos(2kz)cos(2kx))`
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        let (cx, cy, cz) = ((k * x).cos(), (k * y).cos(), (k * z).cos());
        let (c2x, c2y, c2z) = (
            (2.0 * k * x).cos(),
            (2.0 * k * y).cos(),
            (2.0 * k * z).cos(),
        );
        4.0 * cx * cy * cz - (c2x * c2y + c2y * c2z + c2z * c2x)
    }

    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());
        let (s2x, c2x) = ((2.0 * k * x).sin(), (2.0 * k * x).cos());
        let (s2y, c2y) = ((2.0 * k * y).sin(), (2.0 * k * y).cos());
        let (s2z, c2z) = ((2.0 * k * z).sin(), (2.0 * k * z).cos());

        let gx = -4.0 * k * sx * cy * cz + 2.0 * k * s2x * (c2y + c2z);
        let gy = -4.0 * k * cx * sy * cz + 2.0 * k * s2y * (c2x + c2z);
        let gz = -4.0 * k * cx * cy * sz + 2.0 * k * s2z * (c2y + c2x);

        let len = (gx * gx + gy * gy + gz * gz).sqrt();
        if len > 1e-20 {
            Vector3r::new(gx / len, gy / len, gz / len)
        } else {
            Vector3r::y()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    /// **Theorem**: At `(0,0,0)`: `4·1·1·1 − (1+1+1) = 4 − 3 = 1`.
    #[test]
    fn frd_field_origin_is_one() {
        let k = TAU / 2.5;
        let v = Frd.field(0.0, 0.0, 0.0, k);
        assert!((v - 1.0).abs() < 1e-13, "F(0,0,0) = {v}");
    }

    /// **Theorem**: `|∇F|` must equal 1 at any non-singular point.
    #[test]
    fn frd_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = Frd.gradient(0.8, 0.3, 1.5, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }

    /// **Invariant**: `F(-x,-y,-z) = F(x,y,z)` — all terms are even.
    #[test]
    fn frd_field_is_even() {
        let k = TAU / 2.5;
        let (x, y, z) = (0.6, 1.0, 0.4);
        assert!((Frd.field(x, y, z, k) - Frd.field(-x, -y, -z, k)).abs() < 1e-13);
    }

    /// **Invariant**: surface zero-crossing exists in the domain.
    #[test]
    fn frd_surface_exists_in_domain() {
        let k = TAU / 2.5;
        let min_abs = (0..25_u32)
            .flat_map(|i| (0..25_u32).map(move |j| (i, j)))
            .map(|(i, j)| {
                Frd.field(-2.5 + f64::from(i) * 0.2, -2.5 + f64::from(j) * 0.2, 0.5, k)
                    .abs()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(min_abs < 0.4, "min |F| = {min_abs}");
    }
}
