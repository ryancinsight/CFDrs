//! Phase-property coupling for liquid/gas mixtures.
//!
//! This module owns the phase-fraction interpolation boundary. Field-level
//! volume-of-fluid transport and interface reconstruction live under [`vof`],
//! while the [`MomentumCoupling`] selector records the pressure/velocity grid
//! arrangement used by a caller. See
//! `repos/CFDrs/docs/book/turbulence_multiphase.md` for the integration
//! contract.

use aequitas::systems::si::quantities::{Dimensionless, DynamicViscosity, MassDensity};
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement};

/// Momentum-coupling scheme selector (co-located vs staggered mesh).
///
/// Mirrors the `MomentumCoupling` enum spelled out in
/// `turbulence_multiphase.md` § "Multiphase Coupling".
#[derive(Default, Clone, Copy, Debug, PartialEq, Eq)]
pub enum MomentumCoupling {
    /// Co-located velocity-pressure grid.
    #[default]
    Colocated,
    /// Staggered (Marker-and-Cell) grid.
    Staggered,
}

/// Volume-weighted mixture properties from a phase volume fraction.
///
/// `alpha`: gas-phase volume fraction ∈ [0, 1].
/// `rho_l`, `rho_g`: liquid / gas densities (must be finite and non-negative).
/// `mu_l`, `mu_g`: liquid / gas dynamic viscosities (must be finite and
/// non-negative).
///
/// Returns `(rho_mix, mu_mix)`.
///
/// Mirrors the chapter's pseudocode verbatim:
/// ```text
/// let rho_mix = rho_l * (1.0 - alpha) + rho_g * alpha;
/// let mu_mix  = mu_l  * (1.0 - alpha) + mu_g  * alpha;
/// ```
///
/// # Errors
///
/// Returns [`Error::InvalidConfiguration`] when:
///
/// - `alpha` is non-finite or falls outside the closed physical interval
///   `[0, 1]`,
/// - any of `rho_l`, `rho_g`, `mu_l`, `mu_g` is non-finite, or
/// - any of `rho_l`, `rho_g`, `mu_l`, `mu_g` is negative.
///
/// These checks guard the **volume-fraction conservation identity**
/// `alpha + (1 - alpha) = 1`: a non-physical `alpha` produces `liquid_fraction`
/// outside `[0, 1]`, and the mixture density/viscosity stop representing a
/// state of the two-phase mixture. Negative or non-finite single-phase
/// properties are likewise rejected because the resulting mixture properties
/// violate mass and momentum conservation of the underlying transport solve.
pub fn exchange<T: FloatElement>(
    alpha: Dimensionless<T>,
    rho_l: MassDensity<T>,
    rho_g: MassDensity<T>,
    mu_l: DynamicViscosity<T>,
    mu_g: DynamicViscosity<T>,
) -> Result<(MassDensity<T>, DynamicViscosity<T>)> {
    let alpha = alpha.into_base();
    let rho_l = rho_l.into_base();
    let rho_g = rho_g.into_base();
    let mu_l = mu_l.into_base();
    let mu_g = mu_g.into_base();

    if !<T as NumericElement>::is_finite(alpha) {
        return Err(Error::InvalidConfiguration(format!(
            "multiphase::exchange: gas volume fraction must be finite, got alpha={alpha:?}"
        )));
    }
    if alpha < <T as NumericElement>::ZERO || alpha > <T as NumericElement>::ONE {
        return Err(Error::InvalidConfiguration(format!(
            "multiphase::exchange: gas volume fraction must satisfy 0 <= alpha <= 1, got alpha={alpha:?}"
        )));
    }

    for (name, value) in [
        ("rho_l", rho_l),
        ("rho_g", rho_g),
        ("mu_l", mu_l),
        ("mu_g", mu_g),
    ] {
        if !<T as NumericElement>::is_finite(value) {
            return Err(Error::InvalidConfiguration(format!(
                "multiphase::exchange: phase property {name} must be finite, got {value:?}"
            )));
        }
        if value < <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "multiphase::exchange: phase property {name} must be non-negative, got {value:?}"
            )));
        }
    }

    let liquid_fraction = <T as NumericElement>::ONE - alpha;
    let rho_mix = rho_l * liquid_fraction + rho_g * alpha;
    let mu_mix = mu_l * liquid_fraction + mu_g * alpha;
    Ok((
        MassDensity::from_base(rho_mix),
        DynamicViscosity::from_base(mu_mix),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exchange_preserves_phase_limits_and_physical_dimensions() {
        let (rho, viscosity) = exchange(
            Dimensionless::from_base(0.25_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        )
        .expect("physically valid inputs");

        assert!((rho.into_base() - 750.25).abs() < 1.0e-12);
        assert!((viscosity.into_base() - 0.0007545).abs() < 1.0e-12);
    }

    #[test]
    fn exchange_supports_single_precision_without_complex_units() {
        let (rho, viscosity) = exchange(
            Dimensionless::from_base(1.0_f32),
            MassDensity::from_base(1_000.0_f32),
            MassDensity::from_base(1.0_f32),
            DynamicViscosity::from_base(1.0e-3_f32),
            DynamicViscosity::from_base(1.8e-5_f32),
        )
        .expect("physically valid inputs");

        assert!((rho.into_base() - 1.0_f32).abs() < 1.0e-6);
        assert!((viscosity.into_base() - 1.8e-5_f32).abs() < 1.0e-8);
    }

    #[test]
    fn exchange_returns_pure_liquid_density_at_alpha_zero() {
        let (rho, viscosity) = exchange(
            Dimensionless::from_base(0.0_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        )
        .expect("alpha=0 is a valid liquid-only limit");

        assert!((rho.into_base() - 1_000.0).abs() < 1.0e-12);
        assert!((viscosity.into_base() - 1.0e-3).abs() < 1.0e-15);
    }

    #[test]
    fn exchange_volume_conservation_across_grid_alpha_range() {
        // alpha + (1 - alpha) = 1 must hold for the VOF mass-conservation
        // identity; exercise a sweep to verify the mixture returns finite
        // mixture properties without the conservation identity breaking.
        let rho_l = MassDensity::from_base(1_000.0_f64);
        let rho_g = MassDensity::from_base(1.2_f64);
        let mu_l = DynamicViscosity::from_base(1.0e-3_f64);
        let mu_g = DynamicViscosity::from_base(1.8e-5_f64);

        let alphas = [0.0_f64, 0.125, 0.25, 0.5, 0.75, 0.875, 1.0];
        let mut prev_rho = f64::INFINITY;
        for alpha in alphas {
            let (rho, viscosity) =
                exchange(Dimensionless::from_base(alpha), rho_l, rho_g, mu_l, mu_g)
                    .expect("physically valid inputs");

            let rho_value = rho.into_base();
            let mu_value = viscosity.into_base();

            // rho_mix = rho_l * (1 - alpha) + rho_g * alpha is monotone
            // decreasing in alpha when rho_l > rho_g.
            assert!(rho_value.is_finite());
            assert!(rho_value < prev_rho + 1.0e-12);
            prev_rho = rho_value;

            // Viscosity must remain bounded between the two phase values.
            assert!((1.8e-5 - 1.0e-12..=1.0e-3 + 1.0e-12).contains(&mu_value));
        }
    }

    #[test]
    fn exchange_rejects_nan_alpha() {
        let result = exchange(
            Dimensionless::from_base(f64::NAN),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("alpha")),
            "expected InvalidConfiguration error for NaN alpha, got {result:?}"
        );
    }

    #[test]
    fn exchange_rejects_negative_alpha() {
        let result = exchange(
            Dimensionless::from_base(-0.1_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("0 <= alpha")),
            "expected InvalidConfiguration error for negative alpha, got {result:?}"
        );
    }

    #[test]
    fn exchange_rejects_alpha_above_one() {
        let result = exchange(
            Dimensionless::from_base(1.5_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("0 <= alpha")),
            "expected InvalidConfiguration error for alpha > 1, got {result:?}"
        );
    }

    #[test]
    fn exchange_rejects_infinite_density() {
        let result = exchange(
            Dimensionless::from_base(0.5_f64),
            MassDensity::from_base(f64::INFINITY),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("rho_l")),
            "expected InvalidConfiguration error for infinite rho_l, got {result:?}"
        );
    }

    #[test]
    fn exchange_rejects_negative_density() {
        let result = exchange(
            Dimensionless::from_base(0.5_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(-1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("rho_g")),
            "expected InvalidConfiguration error for negative rho_g, got {result:?}"
        );
    }

    #[test]
    fn exchange_rejects_negative_viscosity() {
        let result = exchange(
            Dimensionless::from_base(0.5_f64),
            MassDensity::from_base(1_000.0),
            MassDensity::from_base(1.0),
            DynamicViscosity::from_base(1.0e-3),
            DynamicViscosity::from_base(-1.8e-5),
        );
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("mu_g")),
            "expected InvalidConfiguration error for negative mu_g, got {result:?}"
        );
    }
}
