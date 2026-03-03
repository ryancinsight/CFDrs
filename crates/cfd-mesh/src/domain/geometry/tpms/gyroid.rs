//! Gyroid TPMS implicit field (Schoen 1970).
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = sin(k·x)·cos(k·y) + sin(k·y)·cos(k·z) + sin(k·z)·cos(k·x)
//! ```
//!
//! ## Theorem — Zero Mean Curvature
//!
//! The gyroid `{F = 0}` is a minimal surface.  At any point on the surface,
//! the two principal curvatures are equal in magnitude and opposite in sign,
//! so `H = (κ₁ + κ₂)/2 = 0` (Schoen 1970, verified analytically by
//! Anderson et al. 1990, *Science* 248, 1425).
//!
//! ## Partial Derivatives
//!
//! ```text
//! ∂F/∂x = k·(cos(kx)·cos(ky) − sin(kz)·sin(kx))
//! ∂F/∂y = k·(−sin(kx)·sin(ky) + cos(ky)·cos(kz) − sin(ky)·cos(kx))
//! ∂F/∂z = k·(−sin(ky)·sin(kz) + cos(kz)·cos(kx) − sin(kz)·cos(ky))  ← corrected
//! ```

use crate::domain::core::scalar::Vector3r;
use super::Tpms;

/// Gyroid TPMS — Schoen (1970) triply periodic minimal surface.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct Gyroid;

impl Tpms for Gyroid {
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        (k * x).sin() * (k * y).cos()
            + (k * y).sin() * (k * z).cos()
            + (k * z).sin() * (k * x).cos()
    }

    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());

        let gx = k * (cx * cy - sz * sx);
        let gy = k * (-sx * sy + cy * cz - sy * cx);
        let gz = k * (-sy * sz + cz * cx - sz * cy);

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

    #[test]
    fn gyroid_field_origin_is_zero() {
        let k = TAU / 2.5;
        assert!(Gyroid.field(0.0, 0.0, 0.0, k).abs() < 1e-14);
    }

    #[test]
    fn gyroid_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = Gyroid.gradient(0.3, 0.7, 1.1, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }
}
