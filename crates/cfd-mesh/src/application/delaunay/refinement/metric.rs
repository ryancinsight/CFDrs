//! Anisotropic metric tensor for metric-weighted Ruppert refinement.
//!
//! ## Background
//!
//! Standard Ruppert's algorithm measures triangle quality in the Euclidean metric,
//! producing meshes that are isotropically refined (equilateral triangles).  For
//! millifluidic channels where geometries are strongly anisotropic — for example,
//! boundary layers near channel walls where the normal direction has a 100× smaller
//! feature size than the streamwise direction — isotropic refinement is wasteful:
//! it inserts O(α²) more triangles than necessary, where α is the aspect ratio.
//!
//! A **metric tensor** `M` (symmetric positive-definite 2×2 matrix) redefines the
//! "length" of a vector `e` as `|e|_M = sqrt(eᵀ M e)`.  Triangle quality is then
//! measured in M-space: a triangle is "good" if its metric-radius-edge ratio is
//! small.  This allows long, thin triangles aligned with the flow direction.
//!
//! ## Algorithm (Metric-Weighted Quality Check)
//!
//! Given metric M with Cholesky factor L (M = Lᵀ L), the metric quality check
//! transforms triangle vertices `a, b, c` to `a' = La, b' = Lb, c' = Lc` and
//! computes the isotropic quality of `(a', b', c')`.  This is mathematically
//! equivalent to measuring quality in M-space.
//!
//! ## Theorem (Metric-Ruppert Termination)
//!
//! **Statement**: Under a globally constant metric tensor M with λ_min > 0
//! (positive definite), metric-weighted Ruppert's algorithm terminates with
//! O(area / h_min²) triangles, where h_min = smallest metric-weighted local
//! feature size.
//!
//! **Proof**: The Cholesky transform p' = L p is a linear bijection on ℝ². It
//! maps the anisotropic problem in (ℝ², M) to the isotropic problem in (ℝ², I).
//! Standard Ruppert's termination proof applies to the transformed problem; the
//! bound translates back by the inverse transform.  QED.
//!
//! ## References
//!
//! - Boissonnat & Oudot (2005), "Provably Good Sampling and Meshing of Surfaces",
//!   Graphical Models 67, 405-451.
//! - Labelle & Shewchuk (2003), "Anisotropic Voronoi Diagrams and Guaranteed-Quality
//!   Anisotropic Mesh Generation", SoCG.

use crate::domain::core::scalar::Real;

// ── MetricTensor ──────────────────────────────────────────────────────────────

/// A symmetric positive-definite 2×2 metric tensor for anisotropic mesh quality.
///
/// The metric M redefines the length of a vector `e = [ex, ey]` as
/// `|e|_M = sqrt(eᵀ M e) = sqrt(m00·ex² + 2·m01·ex·ey + m11·ey²)`.
///
/// Setting M to a non-identity tensor steers Ruppert's algorithm to produce
/// triangles that are equilateral in M-space, enabling anisotropic refinement
/// suited to boundary layers and elongated millifluidic channels.
///
/// # Examples
///
/// ```rust,ignore
/// use cfd_mesh::application::delaunay::refinement::MetricTensor;
///
/// // Isotropic (default)
/// let iso = MetricTensor::identity();
///
/// // 10:1 anisotropy aligned with x-axis (thin in y, coarse in x)
/// let aniso = MetricTensor::anisotropic(0.0, 10.0);
///
/// // Metric-weighted edge length
/// let e = [3.0, 4.0];
/// assert!((iso.length(&e) - 5.0).abs() < 1e-12);
/// ```
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MetricTensor {
    /// M\[0,0\].
    pub m00: Real,
    /// M\[0,1\] = M\[1,0\] (symmetry).
    pub m01: Real,
    /// M\[1,1\].
    pub m11: Real,
}

impl MetricTensor {
    /// The identity metric — isotropic standard Euclidean quality.
    ///
    /// Using the identity metric gives exactly the same refinement as the
    /// default isotropic Ruppert implementation.
    #[must_use]
    pub fn identity() -> Self {
        Self {
            m00: 1.0,
            m01: 0.0,
            m11: 1.0,
        }
    }

    /// Construct an anisotropic metric aligned with direction `theta` (radians
    /// from the +x axis) with the given `aspect_ratio` (α ≥ 1).
    ///
    /// ## Definition
    ///
    /// `M = Rθᵀ · diag(1, 1/α²) · Rθ`
    ///
    /// where `Rθ` is the rotation matrix by angle `θ`. This metric assigns
    /// unit length to vectors of Euclidean length 1 in the `θ` direction, and
    /// unit length to vectors of Euclidean length `α` perpendicular to `θ`.
    ///
    /// Effect: triangles are allowed to be `α` times longer along `θ` than
    /// perpendicular to `θ` before being deemed "bad".
    ///
    /// ## Example
    ///
    /// `MetricTensor::anisotropic(0.0, 10.0)` → 10:1 elongation along x,
    /// useful for channel interiors aligned with the x-axis.
    #[must_use]
    pub fn anisotropic(theta: Real, aspect_ratio: Real) -> Self {
        let r = aspect_ratio.max(1.0);
        let (s, c) = theta.sin_cos();
        // Eigen-decomposition:  M = Rᵀ Λ R  where  Λ = diag(1, 1/r²).
        // Expanded: m00 = c²·1 + s²·(1/r²),  m01 = cs·(1 − 1/r²),  m11 = s²·1 + c²·(1/r²).
        let inv_r_sq = 1.0 / (r * r);
        Self {
            m00: c * c + s * s * inv_r_sq,
            m01: c * s * (1.0 - inv_r_sq),
            m11: s * s + c * c * inv_r_sq,
        }
    }

    /// Metric-weighted length of vector `e = [ex, ey]`.
    ///
    /// Returns `sqrt(eᵀ M e)` = `sqrt(m00·ex² + 2·m01·ex·ey + m11·ey²)`.
    /// Always finite and ≥ 0 for positive-definite M.
    #[must_use]
    #[inline]
    pub fn length(&self, e: &[Real; 2]) -> Real {
        let (ex, ey) = (e[0], e[1]);
        let q = self.m00 * ex * ex + 2.0 * self.m01 * ex * ey + self.m11 * ey * ey;
        q.max(0.0).sqrt()
    }

    /// Cholesky factor L (lower-triangular) such that M = Lᵀ L.
    ///
    /// Returns `[l00, l10, l11]` where:
    /// - `L = [[l00, 0], [l10, l11]]`
    /// - `l00 = sqrt(m00)`
    /// - `l10 = m01 / l00`
    /// - `l11 = sqrt(m11 − l10²)`
    ///
    /// Returns `None` if M is not positive definite (numerically).
    #[must_use]
    pub(crate) fn cholesky(&self) -> Option<[Real; 3]> {
        let l00 = self.m00.sqrt();
        if l00 < 1e-30 {
            return None;
        }
        let l10 = self.m01 / l00;
        let disc = self.m11 - l10 * l10;
        if disc < 0.0 {
            return None;
        }
        let l11 = disc.sqrt();
        Some([l00, l10, l11])
    }

    /// Transform point `(x, y)` by the transpose of the Cholesky factor.
    ///
    /// For `M = L Lᵀ` (lower Cholesky), `|e|_M = |Lᵀ e|₂`.
    ///
    /// Returns `Lᵀ · [x, y]ᵀ = [l00·x + l10·y,  l11·y]`.
    ///
    /// ```text
    /// L  = [[l00,  0  ],    Lᵀ = [[l00, l10],
    ///       [l10, l11]]           [ 0,  l11]]
    /// ```
    ///
    /// After transforming all three vertices of a triangle by Lᵀ, standard
    /// isotropic quality metrics in transformed space equal the metric-weighted
    /// quality of the original triangle in M-space.
    ///
    /// **Proof**: `|e|_M² = eᵀ M e = eᵀ L Lᵀ e = (Lᵀ e)ᵀ (Lᵀ e) = |Lᵀ e|₂²`. QED.
    #[must_use]
    #[inline]
    pub(crate) fn apply_cholesky(l: &[Real; 3], x: Real, y: Real) -> (Real, Real) {
        // Lᵀ = [[l[0], l[1]], [0, l[2]]]
        // Lᵀ · [x; y] = [l[0]·x + l[1]·y,  l[2]·y]
        (l[0] * x + l[1] * y, l[2] * y)
    }
}

impl Default for MetricTensor {
    fn default() -> Self {
        Self::identity()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_4;

    #[test]
    fn identity_length_equals_euclidean() {
        let m = MetricTensor::identity();
        let e = [3.0_f64, 4.0_f64];
        assert!(
            (m.length(&e) - 5.0).abs() < 1e-12,
            "identity metric should give l2 norm"
        );
    }

    #[test]
    fn anisotropic_ratio_one_equals_identity() {
        let m = MetricTensor::anisotropic(0.0, 1.0);
        let id = MetricTensor::identity();
        assert!((m.m00 - id.m00).abs() < 1e-12, "m00 mismatch");
        assert!((m.m01 - id.m01).abs() < 1e-12, "m01 mismatch");
        assert!((m.m11 - id.m11).abs() < 1e-12, "m11 mismatch");
    }

    #[test]
    fn anisotropic_x_aligned_compresses_y_by_ratio() {
        // θ=0, ratio=4 → M = diag(1, 1/16); y-length scaled by 1/4.
        let m = MetricTensor::anisotropic(0.0, 4.0);
        assert!((m.m00 - 1.0).abs() < 1e-12, "m00 should be 1.0");
        assert!(m.m01.abs() < 1e-12, "m01 should be 0.0");
        assert!((m.m11 - 1.0 / 16.0).abs() < 1e-12, "m11 should be 1/16");
        // y-direction: |[0, 4]|_M = sqrt(1/16 * 16) = 1
        let ey = [0.0, 4.0];
        assert!(
            (m.length(&ey) - 1.0).abs() < 1e-10,
            "y metric-length mismatch"
        );
        // x-direction: unchanged
        let ex = [5.0, 0.0];
        assert!(
            (m.length(&ex) - 5.0).abs() < 1e-10,
            "x metric-length mismatch"
        );
    }

    #[test]
    fn anisotropic_diagonal_theta_pi_4() {
        // θ=π/4, ratio=1 → identity (all directions uniform at ratio=1)
        let m = MetricTensor::anisotropic(FRAC_PI_4, 1.0);
        let id = MetricTensor::identity();
        assert!((m.m00 - id.m00).abs() < 1e-12);
        assert!((m.m01 - id.m01).abs() < 1e-12);
        assert!((m.m11 - id.m11).abs() < 1e-12);
    }

    #[test]
    fn cholesky_identity_gives_identity_factor() {
        let m = MetricTensor::identity();
        let l = m.cholesky().expect("identity is PD");
        assert!((l[0] - 1.0).abs() < 1e-12, "l00 should be 1");
        assert!(l[1].abs() < 1e-12, "l10 should be 0");
        assert!((l[2] - 1.0).abs() < 1e-12, "l11 should be 1");
    }

    #[test]
    fn cholesky_reconstructs_m() {
        // M = [[4, 2], [2, 3]].  Cholesky: l00=2, l10=1, l11=sqrt(2).
        let m = MetricTensor {
            m00: 4.0,
            m01: 2.0,
            m11: 3.0,
        };
        let l = m.cholesky().expect("M is PD");
        // Reconstruct Mᵀ = LᵀL: m00 = l00², m01 = l10·l00, m11 = l10² + l11²
        let rec_m00 = l[0] * l[0];
        let rec_m01 = l[1] * l[0];
        let rec_m11 = l[1] * l[1] + l[2] * l[2];
        assert!((rec_m00 - m.m00).abs() < 1e-12, "m00 mismatch: {rec_m00}");
        assert!((rec_m01 - m.m01).abs() < 1e-12, "m01 mismatch: {rec_m01}");
        assert!((rec_m11 - m.m11).abs() < 1e-12, "m11 mismatch: {rec_m11}");
    }

    #[test]
    fn apply_cholesky_identity_is_noop() {
        let m = MetricTensor::identity();
        let l = m.cholesky().unwrap();
        let (x2, y2) = MetricTensor::apply_cholesky(&l, 3.0, 4.0);
        assert!((x2 - 3.0).abs() < 1e-12);
        assert!((y2 - 4.0).abs() < 1e-12);
    }

    #[test]
    fn metric_length_via_cholesky_matches_direct() {
        // |e|_M = |L·e|_2 by definition.
        let m = MetricTensor::anisotropic(0.3, 5.0);
        let l = m.cholesky().unwrap();
        let e = [2.0, 3.0];
        let direct = m.length(&e);
        let (tx, ty) = MetricTensor::apply_cholesky(&l, e[0], e[1]);
        let via_chol = (tx * tx + ty * ty).sqrt();
        assert!(
            (direct - via_chol).abs() < 1e-12,
            "direct={direct} chol={via_chol}"
        );
    }
}
