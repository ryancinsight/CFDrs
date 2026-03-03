//! Schwarz D TPMS implicit field (Schwarz 1890) — the Diamond surface.
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = sin(kx)·sin(ky)·sin(kz)
//!           + sin(kx)·cos(ky)·cos(kz)
//!           + cos(kx)·sin(ky)·cos(kz)
//!           + cos(kx)·cos(ky)·sin(kz)
//! ```
//!
//! The zero level-set of this sum accurately reproduces the true Schwarz D
//! (Diamond) surface, which features tetrahedral symmetry and an extremely
//! high surface area density.

use super::Tpms;
use crate::domain::core::scalar::Vector3r;

/// Schwarz D (Diamond) TPMS — Schwarz (1890).
///
/// Uses the full four-term trigonometric approximation to evaluate the diamond surface.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct SchwarzD;

impl Tpms for SchwarzD {
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());

        sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz
    }

    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());

        let gx = k * (cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz);
        let gy = k * (sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz);
        let gz = k * (sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz);

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
    fn schwarz_d_field_origin_is_zero() {
        // sin(0)sin(0)sin(0) + ... = 0
        let k = 1.0_f64;
        let v = SchwarzD.field(0.0, 0.0, 0.0, k);
        assert!(v.abs() < 1e-14, "F = {v}");
    }

    #[test]
    fn schwarz_d_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = SchwarzD.gradient(0.4, 0.9, 0.2, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }
}
