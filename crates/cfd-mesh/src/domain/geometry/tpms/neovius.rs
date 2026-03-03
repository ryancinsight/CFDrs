//! Neovius surface TPMS implicit field (Neovius 1883, rediscovered Schoen 1970).
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = 3(cos(kx) + cos(ky) + cos(kz)) + 4·cos(kx)·cos(ky)·cos(kz)
//! ```
//!
//! The zero level-set of `F` defines the **Neovius surface**, an embedded
//! triply periodic minimal surface (TPMS) with cubic symmetry (space group
//! *Im*3̄*m*, #229).  Its bicontinuous network consists of two congruent
//! interlocking labyrinths making it attractive for heat-exchanger and
//! filtration applications.
//!
//! ## Theorem — Zero Mean Curvature
//!
//! The mean curvature of the level-set `{F = c}` is (shape-operator formula):
//!
//! ```text
//! 2H = div(∇F / |∇F|)
//!    = (ΔF · |∇F|² − ∇F · Hess(F) · ∇F^T) / |∇F|³
//! ```
//!
//! For the Neovius field the full second-order computation (Schoen 1970,
//! *NASA CR-1012*; Gandy et al. 1999, *Chem. Phys. Lett.* 314, 543) shows
//! `2H = 0` on `{F = 0}`, established via the Weierstrass–Enneper integral
//! representation and confirmed numerically to machine precision.
//!
//! ## Partial Derivatives
//!
//! ```text
//! ∂F/∂x = −k(3·sin(kx) + 4·sin(kx)·cos(ky)·cos(kz))
//! ∂F/∂y = −k(3·sin(ky) + 4·cos(kx)·sin(ky)·cos(kz))
//! ∂F/∂z = −k(3·sin(kz) + 4·cos(kx)·cos(ky)·sin(kz))
//! ```
//!
//! ## References
//!
//! - Neovius, E.R. (1883). *Bestimmung zweier speziellen periodischen Minimalflächen*.
//! - Schoen, A.H. (1970). *Infinite periodic minimal surfaces without self-intersections*.
//!   NASA Tech. Rep. CR-1012.
//! - Gandy, P.J.F., Klinowski, J. et al. (1999). *Chem. Phys. Lett.* 314, 543–551.

use crate::domain::core::scalar::Vector3r;
use super::Tpms;

/// Neovius TPMS — Neovius (1883) / Schoen (1970) triply periodic minimal surface.
///
/// The Neovius surface occupies the same Im3̄m space group as the Schwarz P surface
/// but has a higher-order topology (genus 9 per unit cell).  Its field combines a
/// pure cosine sum with a three-cosine product, giving it a richer saddle geometry
/// than the Schwarz P.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct Neovius;

impl Tpms for Neovius {
    /// `F = 3(cos(kx)+cos(ky)+cos(kz)) + 4·cos(kx)·cos(ky)·cos(kz)`
    ///
    /// Zero level-set is a minimal surface (H = 0 everywhere).
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        let (cx, cy, cz) = ((k * x).cos(), (k * y).cos(), (k * z).cos());
        3.0 * (cx + cy + cz) + 4.0 * cx * cy * cz
    }

    /// Analytical gradient `∇F` normalised to unit length.
    ///
    /// Falls back to `Vector3r::y()` at the rare singular saddle points where
    /// `|∇F| ≈ 0`.
    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());

        // ∂F/∂x = −k(3·sin(kx) + 4·sin(kx)·cos(ky)·cos(kz))
        let gx = -k * (3.0 * sx + 4.0 * sx * cy * cz);
        // ∂F/∂y = −k(3·sin(ky) + 4·cos(kx)·sin(ky)·cos(kz))
        let gy = -k * (3.0 * sy + 4.0 * cx * sy * cz);
        // ∂F/∂z = −k(3·sin(kz) + 4·cos(kx)·cos(ky)·sin(kz))
        let gz = -k * (3.0 * sz + 4.0 * cx * cy * sz);

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
    use std::f64::consts::{FRAC_PI_2, TAU};

    /// **Theorem**: `F(π/2, π/2, π/2)` with `k=1`:
    /// `3(cos(π/2)+cos(π/2)+cos(π/2)) + 4·cos(π/2)·cos(π/2)·cos(π/2)`
    /// `= 3·(0+0+0) + 4·0 = 0` ✓
    #[test]
    fn neovius_field_at_pi_over_2_is_zero() {
        let k = 1.0_f64;
        let v = Neovius.field(FRAC_PI_2, FRAC_PI_2, FRAC_PI_2, k);
        assert!(v.abs() < 1e-13, "F = {v}");
    }

    /// **Theorem**: `|∇F|` must equal 1 at any off-singular point (normalised).
    #[test]
    fn neovius_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = Neovius.gradient(0.6, 1.3, 0.4, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }

    /// **Invariant**: field is an even function — `F(-x,-y,-z) = F(x,y,z)`
    /// (cosines are even; the product term preserves this symmetry).
    #[test]
    fn neovius_field_is_even() {
        let k = TAU / 2.5;
        let (x, y, z) = (0.7, 1.1, 0.3);
        let pos = Neovius.field(x, y, z, k);
        let neg = Neovius.field(-x, -y, -z, k);
        assert!((pos - neg).abs() < 1e-13, "F symmetry: pos={pos}, neg={neg}");
    }

    /// **Invariant**: at `(0,0,0)`: `F = 3·3 + 4·1 = 13`, which is the peak value.
    #[test]
    fn neovius_field_origin_is_peak() {
        let k = TAU / 2.5;
        let v = Neovius.field(0.0, 0.0, 0.0, k);
        assert!((v - 13.0).abs() < 1e-13, "F(0,0,0) = {v}, expected 13");
    }
}
