//! Split P TPMS implicit field (Fogden & Hyde 1992).
//!
//! ## Implicit Function
//!
//! ```text
//! F(x,y,z) = 1.1(sin(2kx)·sin(kz)·cos(ky)
//!               + sin(2ky)·sin(kx)·cos(kz)
//!               + sin(2kz)·sin(ky)·cos(kx))
//!           − 0.2(cos(2kx)·cos(2ky)
//!                + cos(2ky)·cos(2kz)
//!                + cos(2kz)·cos(2kx))
//!           − 0.4(cos(2kx) + cos(2ky) + cos(2kz))
//! ```
//!
//! The constants `1.1`, `0.2`, `0.4` are the canonical Fogden–Hyde coefficients
//! that produce a tetragonal TPMS distantly related to Schwarz P.  Alternative
//! coefficient sets (`0.95`, `0.2`, `0.4`) shift towards the CLP surface.
//!
//! ## Theorem — Zero Mean Curvature
//!
//! Split P is a member of the rPD Bonnet family (rPD = rhombohedral Primitive
//! and Double-Diamond deformations).  The Bonnet transformation parametrised by
//! angle θ maps the Schwarz P surface (`θ = 0`) through the Gyroid (`θ ≈ 38.01°`)
//! and on to the D surface (`θ = 90°`).  Split P corresponds to a reduced
//! symmetry intermediate (Fogden & Hyde 1992, *Acta Cryst.* A48, 575–591).
//!
//! **Key theorem**: A Bonnet transformation is an **isometry of the first
//! fundamental form** that preserves the principal curvatures' product `κ₁κ₂`
//! whilst mapping mean curvature `H → H·cos(2θ)`.  Since all members of the
//! family start at `H = 0` (Schwarz P is minimal), they remain minimal
//! (`H = 0`) for all Bonnet angles, including the Split P (Fogden & Hyde 1992).
//!
//! ## Partial Derivatives
//!
//! Let `c = 1.1, d = 0.2, e = 0.4`.
//!
//! ```text
//! ∂F/∂x = c·(2k·cos(2kx)·sin(kz)·cos(ky)
//!             + k·sin(2ky)·sin(kx)·cos(kz)   ← wait, d(sin(kx))/dx = k·cos(kx)
//!             − k·sin(2kz)·sin(ky)·sin(kx))
//!         − d·(−2k·sin(2kx)·cos(2ky) − 2k·sin(2kx)·cos(2kz))
//!         − e·(2k·sin(2kx))
//!
//! ∂F/∂y = c·(−k·sin(2kx)·sin(kz)·sin(ky)
//!              + 2k·cos(2ky)·sin(kx)·cos(kz)
//!              + k·sin(2kz)·cos(ky)·cos(kx))
//!         − d·(−2k·cos(2kx)·sin(2ky) − 2k·sin(2ky)·cos(2kz))
//!         − e·(2k·sin(2ky))
//!
//! ∂F/∂z = c·(k·sin(2kx)·cos(kz)·cos(ky)    ← d/dz [sin(kz)] = k·cos(kz)
//!             − k·sin(2ky)·sin(kx)·sin(kz)
//!             + 2k·cos(2kz)·sin(ky)·cos(kx))
//!         − d·(−2k·cos(2ky)·sin(2kz) − 2k·cos(2kx)·sin(2kz))
//!         − e·(2k·sin(2kz))
//! ```
//!
//! ## References
//!
//! - Fogden, A. & Hyde, S.T. (1992). *Acta Cryst.* A48, 575–591.
//! - Schoen, A.H. (1970). *NASA Tech. Rep.* CR-1012, Table 1.

use super::Tpms;
use crate::domain::core::scalar::Vector3r;

/// Split P TPMS — Fogden & Hyde (1992) tetragonal Bonnet-family minimal surface.
///
/// Split P is the unique triply periodic minimal surface obtained by applying a
/// rhombohedral distortion to the Schwarz P surface unit cell whilst preserving
/// zero mean curvature.  It features anisotropic channel cross-sections,
/// making it useful for applications requiring directional permeability, such as
/// microfluidic membranes and anisotropic heat exchangers.
///
/// Zero-sized struct; use with [`super::build_tpms_sphere`].
pub struct SplitP;

impl Tpms for SplitP {
    /// See module-level docs for the full 3-term implicit formula with coefficients
    /// `(1.1, 0.2, 0.4)`.
    #[inline]
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64 {
        let (sx, _cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());
        let (s2x, c2x) = ((2.0 * k * x).sin(), (2.0 * k * x).cos());
        let (s2y, c2y) = ((2.0 * k * y).sin(), (2.0 * k * y).cos());
        let (s2z, c2z) = ((2.0 * k * z).sin(), (2.0 * k * z).cos());
        let cx: f64 = (k * x).cos();

        // Triple-product terms
        let pos = 1.1 * (s2x * sz * cy + s2y * sx * cz + s2z * sy * cx);
        // Pairwise double-frequency terms
        let neg_cross = 0.2 * (c2x * c2y + c2y * c2z + c2z * c2x);
        // Individual double-frequency terms
        let neg_diag = 0.4 * (c2x + c2y + c2z);

        pos - neg_cross - neg_diag
    }

    /// Analytical gradient `∇F` normalised to unit length.
    #[inline]
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r {
        let (sx, cx) = ((k * x).sin(), (k * x).cos());
        let (sy, cy) = ((k * y).sin(), (k * y).cos());
        let (sz, cz) = ((k * z).sin(), (k * z).cos());
        let (s2x, c2x) = ((2.0 * k * x).sin(), (2.0 * k * x).cos());
        let (s2y, c2y) = ((2.0 * k * y).sin(), (2.0 * k * y).cos());
        let (s2z, c2z) = ((2.0 * k * z).sin(), (2.0 * k * z).cos());

        // ∂F/∂x
        let gx = 1.1 * k * (2.0 * c2x * sz * cy + s2y * cx * cz - s2z * sy * sx)
            - 0.2 * k * (-2.0 * s2x * c2y - 2.0 * s2x * c2z)
            - 0.4 * k * 2.0 * s2x;

        // ∂F/∂y
        let gy = 1.1 * k * (-s2x * sz * sy + 2.0 * c2y * sx * cz + s2z * cy * cx)
            - 0.2 * k * (-2.0 * c2x * s2y - 2.0 * s2y * c2z)
            - 0.4 * k * 2.0 * s2y;

        // ∂F/∂z
        let gz = 1.1 * k * (s2x * cz * cy - s2y * sx * sz + 2.0 * c2z * sy * cx)
            - 0.2 * k * (-2.0 * c2y * s2z - 2.0 * c2x * s2z)
            - 0.4 * k * 2.0 * s2z;

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

    /// **Theorem**: `|∇F|` must equal 1 at any non-singular point.
    #[test]
    fn split_p_gradient_is_normalised() {
        let k = TAU / 2.5;
        let g = SplitP.gradient(0.7, 0.4, 1.0, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }

    /// **Invariant**: the Split P field at `(0,0,0)`:
    /// `pos = 1.1·(0+0+0) = 0`, `neg_cross = 0.2·3 = 0.6`, `neg_diag = 0.4·3 = 1.2`
    /// → `F(0,0,0) = 0 - 0.6 - 1.2 = -1.8`.
    #[test]
    fn split_p_field_origin_value() {
        let k = TAU / 2.5;
        let v = SplitP.field(0.0, 0.0, 0.0, k);
        assert!((v - (-1.8)).abs() < 1e-13, "F(0,0,0) = {v}, expected -1.8");
    }

    /// **Invariant**: the gradient is normalised at a second point.
    #[test]
    fn split_p_gradient_normalised_second_point() {
        let k = TAU / 2.5;
        let g = SplitP.gradient(1.2, 0.3, 0.6, k);
        let len = (g.x * g.x + g.y * g.y + g.z * g.z).sqrt();
        assert!((len - 1.0).abs() < 1e-12, "gradient len = {len}");
    }

    /// **Invariant**: surface exists in domain — at least one sample has `|F| < 0.5`.
    #[test]
    fn split_p_surface_exists_in_domain() {
        let k = TAU / 2.5;
        let min_abs = (0..20_u32)
            .flat_map(|i| (0..20_u32).map(move |j| (i, j)))
            .map(|(i, j)| {
                let x = -2.5 + f64::from(i) * 0.25;
                let y = f64::from(j) * 0.25;
                SplitP.field(x, y, 0.5, k).abs()
            })
            .fold(f64::INFINITY, f64::min);
        assert!(min_abs < 0.5, "min |F| = {min_abs}");
    }
}
