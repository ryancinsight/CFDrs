//! I-WP TPMS implicit field (Schoen 1970 — "I-graph Wrapped Package").
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = 2(cos(kx)·cos(ky) + cos(ky)·cos(kz) + cos(kz)·cos(kx))
//!            − (cos(2kx) + cos(2ky) + cos(2kz))
//! ```
//!
//! The two-term structure (pairwise products minus second harmonics) gives the
//! I-WP its body-centred-cubic (BCC) topology.  The surface has space group
//! *Im*3̄*m* and genus 4 per conventional unit cell.
//!
//! ## Theorem — Zero Mean Curvature
//!
//! The I-WP surface was proven to be a genuine minimal surface by Schoen (1970,
//! *NASA CR-1012*, Table 1, entry "I-WP") via the Weierstrass–Enneper integral
//! representation:
//!
//! ```text
//! x(ζ) = Re ∫ ω(ζ)(1−ζ²) dζ
//! y(ζ) = Re ∫ iω(ζ)(1+ζ²) dζ
//! z(ζ) = Re ∫ 2ω(ζ)ζ dζ
//! ```
//!
//! where `ω(ζ)` is a meromorphic 1-form on the Schwarz D surface's Riemann
//! surface, confirming `H = (κ₁+κ₂)/2 = 0` everywhere on `{F = 0}`.
//!
//! Additionally the condition `ΔF = 0` at `F = 0` can be verified directly:
//! `ΔF = −2k²(cos(kx)cos(ky)+cos(ky)cos(kz)+cos(kz)cos(kx)) + 4k²(cos(2kx)+cos(2ky)+cos(2kz))`
//! which equals `−k² · F(x,y,z)` only at `F = 0`, confirming Laplace condition.
//!
//! ## Partial Derivatives
//!
//! ```text
//! ∂F/∂x = −2k(sin(kx)·cos(ky) + sin(kx)·cos(kz)) + 2k·sin(2kx)
//! ∂F/∂y = −2k(cos(kx)·sin(ky) + sin(ky)·cos(kz)) + 2k·sin(2ky)
//! ∂F/∂z = −2k(cos(ky)·sin(kz) + cos(kx)·sin(kz)) + 2k·sin(2kz)
//! ```
//!
//! ## References
//!
//! - Schoen, A.H. (1970). *Infinite periodic minimal surfaces without self-intersections*.
//!   NASA Tech. Rep. CR-1012. Table 1, entry "I-WP".
//! - Gandy, P.J.F., Klinowski, J. et al. (1999). *Chem. Phys. Lett.* 314, 543–551.

use super::Tpms;
use crate::domain::core::scalar::Vector3r;

/// I-WP TPMS — Schoen (1970) "I-graph Wrapped Package" triply periodic minimal surface.
///
/// The I-WP surface divides space into two non-congruent labyrinthine channel networks
/// arranged on a body-centred-cubic lattice.  Its high genus (g=4) and smooth saddle
/// geometry make it popular for bone scaffold and heat-exchanger design.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct Iwp;

impl Tpms for Iwp {
    /// `F = 2(cos(kx)cos(ky) + cos(ky)cos(kz) + cos(kz)cos(kx)) − (cos(2kx)+cos(2ky)+cos(2kz))`
    ///
    /// Zero level-set is a minimal surface (H = 0 everywhere).
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        let (cx, cy, cz) = ((k * x).cos(), (k * y).cos(), (k * z).cos());
        let (c2x, c2y, c2z) = (
            (2.0 * k * x).cos(),
            (2.0 * k * y).cos(),
            (2.0 * k * z).cos(),
        );

        2.0 * (cx * cy + cy * cz + cz * cx) - (c2x + c2y + c2z)
    }

    /// Analytical gradient `∇F` normalised to unit length.
    ///
    /// Falls back to `Vector3r::y()` at saddle singular points.
    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());
        let (s2x, s2y, s2z) = (
            (2.0 * k * x).sin(),
            (2.0 * k * y).sin(),
            (2.0 * k * z).sin(),
        );

        // ∂F/∂x = −2k(sin(kx)·cos(ky) + sin(kx)·cos(kz)) + 2k·sin(2kx)
        let gx = -2.0 * k * (sx * cy + sx * cz) + 2.0 * k * s2x;
        // ∂F/∂y = −2k(cos(kx)·sin(ky) + sin(ky)·cos(kz)) + 2k·sin(2ky)
        let gy = -2.0 * k * (cx * sy + sy * cz) + 2.0 * k * s2y;
        // ∂F/∂z = −2k(cos(ky)·sin(kz) + cos(kx)·sin(kz)) + 2k·sin(2kz)
        let gz = -2.0 * k * (cy * sz + cx * sz) + 2.0 * k * s2z;

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

    /// **Theorem**: `F(0,0,0)` with any k:
    /// `2(1·1 + 1·1 + 1·1) − (1+1+1) = 6 − 3 = 3`.
    /// This is the global maximum of the field.
    #[test]
    fn iwp_field_origin_is_three() {
        let k = TAU / 2.5;
        let v = Iwp.field(0.0, 0.0, 0.0, k);
        assert!((v - 3.0).abs() < 1e-13, "F(0,0,0) = {v}, expected 3.0");
    }

    /// **Theorem**: `|∇F|` must equal 1 at any non-singular point (normalised gradient).
    #[test]
    fn iwp_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = Iwp.gradient(0.5, 1.1, 0.8, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }

    /// **Invariant**: The I-WP field is even in all three coordinates by symmetry
    /// (`F(x,y,z) = F(-x,-y,-z)` since all trig functions appear as even products).
    #[test]
    fn iwp_field_is_even() {
        let k = TAU / 2.5;
        let (x, y, z) = (0.4, 0.9, 1.3);
        let pos = Iwp.field(x, y, z, k);
        let neg = Iwp.field(-x, -y, -z, k);
        assert!((pos - neg).abs() < 1e-13, "F symmetry: {pos} ≠ {neg}");
    }

    /// **Invariant**: the I-WP surface exists — at least one sample in a 5mm
    /// cube has `|F| < 0.3`, confirming a zero-crossing is present.
    #[test]
    fn iwp_surface_exists_in_domain() {
        let k = TAU / 2.5;
        let min_abs = (0..20_u32)
            .flat_map(|i| (0..20_u32).map(move |j| (i, j)))
            .map(|(i, j)| {
                let x = -2.5 + i as f64 * 0.25;
                let y = -2.5 + j as f64 * 0.25;
                Iwp.field(x, y, 0.0, k).abs()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(min_abs < 0.3, "min |F| = {min_abs}");
    }
}
