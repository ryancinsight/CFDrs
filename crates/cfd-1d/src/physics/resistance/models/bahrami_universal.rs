//! Exact laminar resistance helpers for common microfluidic cross-sections.
//!
//! This module keeps the shape-specific formulas used by the cfd-1d resistance
//! stack:
//! - `resistance_circular`: exact Hagen-Poiseuille
//! - `resistance_rectangular`: validated rectangular-channel SSOT
//! - `resistance_elliptical`: exact closed-form ellipse
//!
//! The small geometry helpers below support those exact formulas.

use std::f64::consts::PI;

/// Compute the dimensionless polar moment $I^*_p = I_{polar}/A^2$
/// for a rectangular cross-section with width $w$ and height $h$.
///
/// ```text
/// I_{polar} = w·h·(w² + h²)/12
/// I*_p = (w² + h²) / (12·w·h)
/// ```
#[inline]
#[must_use]
pub fn rectangular_ip(width: f64, height: f64) -> f64 {
    let w = width.max(1e-15);
    let h = height.max(1e-15);
    (w * w + h * h) / (12.0 * w * h)
}

/// Compute the dimensionless polar moment $I^*_p$ for an elliptical cross-section.
///
/// ```text
/// I_{polar} = π·a·b·(a² + b²)/4
/// I*_p = (a² + b²) / (4·π·a·b)
/// ```
#[inline]
#[must_use]
pub fn elliptical_ip(semi_major: f64, semi_minor: f64) -> f64 {
    let a = semi_major.max(1e-15);
    let b = semi_minor.max(1e-15);
    (a * a + b * b) / (4.0 * PI * a * b)
}

/// Dimensionless perimeter $\Gamma = P / \sqrt{A}$.
#[inline]
#[must_use]
pub fn gamma_parameter(perimeter: f64, area: f64) -> f64 {
    perimeter / area.max(1e-30).sqrt()
}

/// Exact laminar resistance for a circular tube.
///
/// For circles the analytical HP formula is exact, so we use it directly.
/// This function serves as a comparison baseline.
///
/// ```text
/// R_circle = 128 µ L / (π D⁴)    [Hagen-Poiseuille]
/// ```
#[inline]
#[must_use]
pub fn resistance_circular(diameter_m: f64, length_m: f64, viscosity_pa_s: f64) -> f64 {
    if diameter_m <= 0.0 || length_m <= 0.0 || viscosity_pa_s <= 0.0 {
        return 0.0;
    }

    128.0 * viscosity_pa_s * length_m / (PI * diameter_m.powi(4))
}

/// Validated laminar resistance for a rectangular channel.
#[inline]
#[must_use]
pub fn resistance_rectangular(
    width_m: f64,
    height_m: f64,
    length_m: f64,
    viscosity_pa_s: f64,
) -> f64 {
    if width_m <= 0.0 || height_m <= 0.0 || length_m <= 0.0 || viscosity_pa_s <= 0.0 {
        return 0.0;
    }

    let w = width_m;
    let h = height_m;
    let area = w * h;
    let hydraulic_diameter = 2.0 * w * h / (w + h);
    let epsilon = w.min(h) / w.max(h);
    let pi = PI;
    let poiseuille_number = if epsilon == 0.0 {
        96.0
    } else {
        96.0 / ((1.0 + epsilon).powi(2)
            * (1.0 - (192.0 * epsilon / pi.powi(5)) * (pi / (2.0 * epsilon)).tanh()))
    };

    poiseuille_number * viscosity_pa_s * length_m / (2.0 * area * hydraulic_diameter.powi(2))
}

/// Exact laminar resistance for an elliptical channel.
#[inline]
#[must_use]
pub fn resistance_elliptical(
    semi_major_m: f64,
    semi_minor_m: f64,
    length_m: f64,
    viscosity_pa_s: f64,
) -> f64 {
    if semi_major_m <= 0.0 || semi_minor_m <= 0.0 || length_m <= 0.0 || viscosity_pa_s <= 0.0 {
        return 0.0;
    }

    let a = semi_major_m;
    let b = semi_minor_m;
    4.0 * viscosity_pa_s * length_m * (a * a + b * b) / (PI * a.powi(3) * b.powi(3))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::ConstantPropertyFluid;
    use crate::physics::resistance::models::{FlowConditions, RectangularChannelModel, ResistanceModel};

    const MU: f64 = 0.001; // Water viscosity [Pa·s]
    const L: f64 = 0.01;   // 1 cm channel

    fn water() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "water".to_string(),
            1000.0,
            MU,
            4186.0,
            0.598,
            1480.0,
        )
    }

    /// Circular tube: exact Hagen-Poiseuille resistance.
    #[test]
    fn circular_exact_matches_hagen_poiseuille() {
        let d = 100e-6; // 100 µm
        let r_exact = resistance_circular(d, L, MU);
        let r_hp = 128.0 * MU * L / (PI * d.powi(4));
        assert_relative_eq!(r_exact, r_hp, epsilon = 1e-15);
    }

    /// Rectangular channel: validated rectangular-channel model.
    #[test]
    fn rectangular_matches_shah_london_fit() {
        let w = 200e-6;
        let h = 100e-6; // AR = 2
        let r_rect = resistance_rectangular(w, h, L, MU);
        let model = RectangularChannelModel::new(w, h, L);
        let fluid = water();
        let conditions = FlowConditions::new(0.0);
        let (r_model, _) = model
            .calculate_coefficients(&fluid, &conditions)
            .expect("rectangular model should evaluate successfully");

        assert_relative_eq!(r_rect, r_model, epsilon = 1e-12);
    }

    /// Square channel: resistance should be positive and finite.
    #[test]
    fn square_channel_reasonable() {
        let s = 100e-6;
        let r_sq = resistance_rectangular(s, s, L, MU);
        assert!(
            r_sq > 0.0 && r_sq.is_finite(),
            "Square resistance should be finite positive, got {r_sq:.3e}"
        );
    }

    /// Resistance scales linearly with length.
    #[test]
    fn resistance_scales_with_length() {
        let r1 = resistance_circular(100e-6, 0.01, MU);
        let r2 = resistance_circular(100e-6, 0.02, MU);
        let ratio = r2 / r1;
        assert!(
            (ratio - 2.0).abs() < 0.01,
            "Resistance ratio {ratio:.4} should be 2.0 for double length"
        );
    }

    /// Resistance scales linearly with viscosity.
    #[test]
    fn resistance_scales_with_viscosity() {
        let r1 = resistance_circular(100e-6, L, 0.001);
        let r2 = resistance_circular(100e-6, L, 0.003);
        let ratio = r2 / r1;
        assert!(
            (ratio - 3.0).abs() < 0.01,
            "Resistance ratio {ratio:.4} should be 3.0 for triple viscosity"
        );
    }

    /// Elliptical channel with a=b reduces to circular.
    #[test]
    fn elliptical_circle_limit() {
        let r = 50e-6;
        let r_ellipse = resistance_elliptical(r, r, L, MU);
        let r_circle = resistance_circular(2.0 * r, L, MU);
        let rel_err = (r_ellipse - r_circle).abs() / r_circle;
        assert!(
            rel_err < 1e-12,
            "Ellipse(a=b) R={r_ellipse:.3e} should match circle R={r_circle:.3e}, err={rel_err:.4}"
        );
    }

    #[test]
    fn elliptical_exact_formula_matches_closed_form() {
        let a = 80e-6;
        let b = 40e-6;
        let expected = 4.0 * MU * L * (a * a + b * b) / (PI * a.powi(3) * b.powi(3));
        let actual = resistance_elliptical(a, b, L, MU);
        assert_relative_eq!(actual, expected, epsilon = 1e-15);
    }
}
