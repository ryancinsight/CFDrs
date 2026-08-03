//! Phase-property coupling for liquid/gas mixtures.
//!
//! This module owns the phase-fraction interpolation boundary. Field-level
//! volume-of-fluid transport and interface reconstruction live under [`vof`],
//! while the [`MomentumCoupling`] selector records the pressure/velocity grid
//! arrangement used by a caller. See
//! `repos/CFDrs/docs/book/turbulence_multiphase.md` for the integration
//! contract.

use aequitas::systems::si::quantities::{Dimensionless, DynamicViscosity, MassDensity};
use eunomia::FloatElement;

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
/// `rho_l`, `rho_g`: liquid / gas densities.
/// `mu_l`, `mu_g`: liquid / gas dynamic viscosities.
///
/// Returns `(rho_mix, mu_mix)`.
///
/// Mirrors the chapter's pseudocode verbatim:
/// ```text
/// let rho_mix = rho_l * (1.0 - alpha) + rho_g * alpha;
/// let mu_mix  = mu_l  * (1.0 - alpha) + mu_g  * alpha;
/// ```
pub fn exchange<T: FloatElement>(
    alpha: Dimensionless<T>,
    rho_l: MassDensity<T>,
    rho_g: MassDensity<T>,
    mu_l: DynamicViscosity<T>,
    mu_g: DynamicViscosity<T>,
) -> (MassDensity<T>, DynamicViscosity<T>) {
    let alpha = alpha.into_base();
    let liquid_fraction = T::ONE - alpha;
    let rho_mix = rho_l.into_base() * liquid_fraction + rho_g.into_base() * alpha;
    let mu_mix = mu_l.into_base() * liquid_fraction + mu_g.into_base() * alpha;
    (
        MassDensity::from_base(rho_mix),
        DynamicViscosity::from_base(mu_mix),
    )
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
        );

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
        );

        assert!((rho.into_base() - 1.0_f32).abs() < 1.0e-6);
        assert!((viscosity.into_base() - 1.8e-5_f32).abs() < 1.0e-8);
    }
}
