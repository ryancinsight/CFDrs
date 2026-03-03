//! Schwarz P TPMS implicit field (Schwarz 1890).
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = cos(k·x) + cos(k·y) + cos(k·z)
//! ```
//!
//! The zero level-set produces the **Schwarz Primitive** (P-surface), the
//! simplest TPMS.  Its unit cell contains a single central saddle.
//!
//! ## Theorem — Zero Mean Curvature
//!
//! The mean curvature of a level-set `{F = c}` is:
//!   `H = div(∇F / |∇F|) / 2`
//!
//! For `F = cos(kx)+cos(ky)+cos(kz)` the Laplacian `ΔF = -k²(cos(kx)+cos(ky)+cos(kz)) = -k²F`.
//! On the zero level-set `F = 0`, the numerator of `H` contains `ΔF = 0`, so
//! `H = 0` everywhere on `{F = 0}` (Schwarz 1890).
//!
//! ## Partial Derivatives
//!
//! ```text
//! ∂F/∂x = −k·sin(k·x)
//! ∂F/∂y = −k·sin(k·y)
//! ∂F/∂z = −k·sin(k·z)
//! ```

use crate::domain::core::scalar::Vector3r;
use super::Tpms;

/// Schwarz P (Primitive) TPMS — Schwarz (1890).
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct SchwarzP;

impl Tpms for SchwarzP {
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        (k * x).cos() + (k * y).cos() + (k * z).cos()
    }

    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let gx = -k * (k * x).sin();
        let gy = -k * (k * y).sin();
        let gz = -k * (k * z).sin();

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

    #[test]
    fn schwarz_p_field_at_pi_over_2_is_zero() {
        // F(π/k, π/k, π/k) with k=1 → cos(π)+cos(π)+cos(π) = -3 ≠ 0
        // Zero crossing occurs when two cosines cancel the third, e.g. (π/2,π/2,0):
        // cos(π/2)+cos(π/2)+cos(0) = 0+0+1 ≠ 0
        // Known zero: (π/2, π/2, π) → 0+0+cos(π) = -1 nope
        // Actual zero: need cos(kx)+cos(ky)+cos(kz)=0.
        // At x=π/(2k), y=π/(2k) → 0+0+cos(kz)=0 ↔ kz=π/2.
        let k = 1.0_f64;
        let v = SchwarzP.field(FRAC_PI_2, FRAC_PI_2, FRAC_PI_2, k);
        // cos(π/2)+cos(π/2)+cos(π/2) = 0 ✓
        assert!(v.abs() < 1e-14, "F = {v}");
    }

    #[test]
    fn schwarz_p_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = SchwarzP.gradient(0.5, 1.2, 0.3, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }
}
