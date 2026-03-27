//! Cross-section shape factor (Poiseuille number) computation.
//!
//! # Theorem: Laminar Poiseuille Number (Shah-London 1978, Table 48)
//!
//! For fully-developed laminar flow in a duct of arbitrary cross-section,
//! the product of the Darcy friction factor `f` and the Reynolds number
//! `Re` is a geometric constant called the **Poiseuille number** `Po`:
//!
//! ```text
//! f · Re = Po
//! ```
//!
//! `Po` depends only on the cross-section shape, not on fluid properties
//! or flow conditions.
//!
//! ## Canonical cross-section values
//!
//! | Shape                    | Poiseuille number Po |
//! |--------------------------|---------------------|
//! | Circular                 | 64 (exact)          |
//! | Square (AR = 1)          | 56.908              |
//! | Infinite plates (AR → ∞) | 96                  |
//!
//! ## Shah-London 5-term polynomial (rectangular / trapezoidal)
//!
//! For a rectangle with aspect ratio `α = a/b ≥ 1` (a = longer, b = shorter):
//!
//! ```text
//! Po(α) = 96 · (1 − 1.3553/α + 1.9467/α² − 1.7012/α³ + 0.9564/α⁴ − 0.2537/α⁵)
//! ```
//!
//! This polynomial is accurate to within 0.05% for all `α ≥ 1`.
//!
//! **Proof sketch**: The coefficients arise from a Fourier-series solution of the
//! biharmonic equation for fully-developed flow in a rectangular duct, truncated
//! at 5 terms. The polynomial was fitted by Shah & London (1978) by minimising
//! the L∞ error against the exact infinite-series solution across `α ∈ [1, ∞)`.
//!
//! # References
//! - Shah, R. K. & London, A. L. (1978). *Laminar Flow Forced Convection in Ducts*.
//!   Academic Press, Supplement 1 to Advances in Heat Transfer. Table 48.
//! - Dryden, H. L., Murnaghan, F. D. & Bateman, H. (1932). *Hydrodynamics*.
//!   NRC Bulletin No. 84 (elliptical exact solution).
//! - Muzychka, Y. S. & Yovanovich, M. M. (2002). Laminar forced convection
//!   heat transfer in the combined entry region of non-circular ducts.
//!   *ASME J. Heat Transfer*, 124(2), 260–267.

use crate::domain::channel::cross_section::CrossSection;
use nalgebra::RealField;
use num_traits::{cast::FromPrimitive, Float};

/// Shah-London 5-term polynomial coefficients for `Po(α)`.
///
/// Source: Shah & London (1978), Table 48.
const SHAH_LONDON_C: [f64; 5] = [1.3553, 1.9467, 1.7012, 0.9564, 0.2537];

/// Evaluate the Shah-London Poiseuille number `Po(α)` for a rectangular
/// (or equivalent-rectangular) cross-section.
///
/// # Theorem
///
/// For `α = a/b ≥ 1`:
/// ```text
/// Po(α) = 96 · (1 − 1.3553/α + 1.9467/α² − 1.7012/α³ + 0.9564/α⁴ − 0.2537/α⁵)
/// ```
///
/// # Arguments
/// - `alpha` — aspect ratio `a/b ≥ 1` (longer side / shorter side)
///
/// # Returns
/// Poiseuille number `Po ∈ [56.908, 96]`.
pub fn shah_london_po<T: RealField + Copy + FromPrimitive>(alpha: T) -> T {
    let inv_a = T::one() / alpha;
    let mut power = inv_a;
    let mut correction = T::zero();
    let mut sign = -T::one();

    for &coeff in &SHAH_LONDON_C {
        let c = T::from_f64(coeff).expect("Shah-London coefficient conversion");
        correction += sign * c * power;
        power *= inv_a;
        sign = -sign;
    }

    T::from_f64(96.0).expect("Po base constant") * (T::one() + correction)
}

/// Compute the Poiseuille number `Po = f·Re` for a given cross-section shape.
///
/// Dispatches to the appropriate analytical formula based on cross-section type.
///
/// See module-level documentation for the governing theorems and references.
pub fn poiseuille_number<T: RealField + Copy + FromPrimitive + Float>(
    cross_section: &CrossSection<T>,
) -> T {
    match cross_section {
        CrossSection::Circular { .. } => {
            // Exact Hagen-Poiseuille: Po = 64
            T::from_f64(64.0).expect("Po constant")
        }
        CrossSection::Rectangular { width, height } => {
            // Shah-London polynomial for rectangular ducts
            let (a, b) = if *width >= *height {
                (*width, *height)
            } else {
                (*height, *width)
            };
            if b <= T::zero() {
                return T::from_f64(64.0).expect("Po fallback constant");
            }
            shah_london_po(a / b)
        }
        CrossSection::Elliptical {
            major_axis,
            minor_axis,
        } => {
            // Exact Stokes solution: Po = 2π(a² + b²) / (a·b)
            // (Dryden, Murnaghan & Bateman 1932)
            let two = T::one() + T::one();
            let a = *major_axis / two;
            let b = *minor_axis / two;
            if a <= T::zero() || b <= T::zero() {
                return T::from_f64(64.0).expect("Po fallback constant");
            }
            T::from_f64(2.0 * std::f64::consts::PI).expect("2π constant") * (a * a + b * b)
                / (a * b)
        }
        CrossSection::Trapezoidal {
            top_width,
            bottom_width,
            height,
        } => {
            // Muzychka & Yovanovich (2002) equivalent rectangle method:
            // use mean width as the effective rectangular width.
            let avg_width = (*top_width + *bottom_width) / (T::one() + T::one());
            if *height <= T::zero() || avg_width <= T::zero() {
                return T::from_f64(64.0).expect("Po fallback constant");
            }
            let (a, b) = if avg_width >= *height {
                (avg_width, *height)
            } else {
                (*height, avg_width)
            };
            shah_london_po(a / b)
        }
        CrossSection::Custom { .. } => {
            // For custom cross-sections, assume circular (Po = 64).
            T::from_f64(64.0).expect("Po fallback constant")
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Square channel (α = 1): Po should be ≈ 56.908 per Shah-London Table 48.
    #[test]
    fn test_shah_london_square() {
        let po: f64 = shah_london_po(1.0);
        assert!(
            (po - 56.908).abs() < 0.1,
            "Square channel Po should be ~56.908, got {po}"
        );
    }

    /// As α → ∞, Po → 96 (parallel plate limit).
    #[test]
    fn test_shah_london_wide_aspect_ratio() {
        let po: f64 = shah_london_po(100.0);
        // Shah-London 5-term polynomial approaches 96 asymptotically.
        // At α = 100, the truncation residual gives Po ≈ 94.72.
        assert!(
            (po - 96.0).abs() < 2.0,
            "Wide aspect ratio Po should approach 96, got {po}"
        );
    }

    /// Po is monotonically increasing with aspect ratio.
    #[test]
    fn test_shah_london_monotonic() {
        let mut prev: f64 = shah_london_po(1.0);
        for alpha in [1.5, 2.0, 3.0, 5.0, 10.0, 50.0] {
            let po: f64 = shah_london_po(alpha);
            assert!(
                po > prev,
                "Po should increase with α: Po({alpha}) = {po} <= Po(prev) = {prev}"
            );
            prev = po;
        }
    }

    /// Circular cross-section gives Po = 64 exactly.
    #[test]
    fn test_circular_po() {
        let cs = CrossSection::Circular { diameter: 0.001 };
        let po: f64 = poiseuille_number(&cs);
        assert!(
            (po - 64.0).abs() < 1e-12,
            "Circular Po should be exactly 64, got {po}"
        );
    }

    /// Rectangular cross-section delegates to Shah-London.
    #[test]
    fn test_rectangular_po() {
        let cs = CrossSection::Rectangular {
            width: 0.002,
            height: 0.001,
        };
        let po: f64 = poiseuille_number(&cs);
        let expected: f64 = shah_london_po(2.0);
        assert!(
            (po - expected).abs() < 1e-12,
            "Rectangular Po should match shah_london_po(2.0): got {po}, expected {expected}"
        );
    }
}
