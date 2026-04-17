//! Biological Cellular Cavitation Damage Models.
//!
//! # Theorem — Cellular Strain Fractional Injury
//!
//! The total cell population mass fraction is rigorously conserved across four topological injury states:
//!
//! ```math
//! f_{healthy} + f_{permeabilized} + f_{necrotic} + f_{lysed} = 1.0
//! ```
//!
//! When exposed to an external Rayleigh collapse impact pressure $P_{impact}$, the fractional distributions
//! transition sequentially based on critical mechanical strain thresholds of the cellular membrane. Unaffected
//! cells are defined dynamically based on remainder.
//!
//! **Proof sketch**: This represents a deterministic mapping from derived spatial acoustic pressure to cumulative probability
//! thresholds (bounded within `[0.0, 1.0]`). Each state represents mutually exclusive populations, ensuring
//! exact mass conservation. The applied membrane load is the local
//! cavitation-induced overpressure above ambient, so ambient-only conditions
//! produce zero induced strain and therefore preserve the fully healthy state.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Defined structural mechanics for cellular membranes.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CellularMembraneMechanics<T: RealField + Copy> {
    /// Nominal radius of the cell [m].
    pub cell_radius_m: T,
    /// Thickness of the cellular membrane [m].
    pub membrane_thickness_m: T,
    /// Elastic modulus of the cellular membrane [Pa].
    pub membrane_elastic_modulus_pa: T,
    /// Areal strain threshold inducing reversible or irreversible cell membrane pores.
    pub critical_areal_strain_permeabilization: T,
    /// Areal strain threshold inducing cell death without complete fragmentation.
    pub critical_areal_strain_necrosis: T,
    /// Areal strain threshold inducing catastrophic cell rupture / complete fragmentation.
    pub critical_areal_strain_lysis: T,
}

/// Distribution of cellular populations post-cavitation injury.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CellularInjuryProfile<T: RealField + Copy> {
    /// Fraction of unaffected robust cells in the localized domain.
    pub fraction_healthy: T,
    /// Fraction of cells surviving but acquiring membrane porosity (ideal for drug uptake).
    pub fraction_permeabilized: T,
    /// Fraction of cells undergoing mechanically-induced death without destruction.
    pub fraction_necrotic: T,
    /// Fraction of cells completely destroyed by microjet or shockwave structural collapse.
    pub fraction_lysed: T,
}

impl<T: RealField + Copy> CellularMembraneMechanics<T> {
    /// Validates the strictly monotonic ordering of cellular injury thresholds and positive dimensions.
    pub fn is_valid(&self) -> bool {
        let zero = T::zero();
        self.cell_radius_m > zero
            && self.membrane_thickness_m > zero
            && self.membrane_elastic_modulus_pa > zero
            && self.critical_areal_strain_permeabilization > zero
            && self.critical_areal_strain_permeabilization < self.critical_areal_strain_necrosis
            && self.critical_areal_strain_necrosis < self.critical_areal_strain_lysis
    }

    /// Evaluates the fractional injury profile for a bubble collapse event.
    ///
    /// Acoustic radiation follows a $1/r$ decay bounded by the microbubble minimum radius.
    /// Stress is derived using thin-walled pressure vessel approximation.
    pub fn evaluate_spatial_injury(
        &self,
        collapse_pressure_pa: T,
        bubble_radius_m: T,
        standoff_distance_m: T,
        ambient_pressure_pa: T,
    ) -> CellularInjuryProfile<T> {
        let zero = T::zero();
        let one = T::one();

        let r_impact = standoff_distance_m.max(bubble_radius_m);
        // Acoustic radiation decay P(r) = P_inf + R/r * (P_c - P_inf)
        // Bounded mathematically to ensure P(r) does not exceed collapse pressure
        let ratio = bubble_radius_m / r_impact;
        let pressure_difference = collapse_pressure_pa - ambient_pressure_pa;
        let p_r = ambient_pressure_pa + ratio * pressure_difference;

        // The membrane responds to cavitation-induced overpressure, not to the
        // ambient baseline itself. Ambient-only loading must therefore map to
        // zero induced strain.
        let two = one + one;
        let impact_overpressure = (p_r - ambient_pressure_pa).max(zero);
        let membrane_stress =
            (impact_overpressure * self.cell_radius_m) / (two * self.membrane_thickness_m);

        // Areal strain approx for Hookean membrane
        let strain = membrane_stress / self.membrane_elastic_modulus_pa;

        let evaluate_cumulative = |threshold: T| -> T {
            if strain <= zero {
                zero
            } else {
                let k = T::from_f64(5.0).unwrap_or(one); // shape factor
                let r = strain / threshold;
                // Exponent formulation representing survival reliability function 1 - exp(- (e/e_c)^k)
                let exponent = -r.powf(k);
                one - exponent.exp()
            }
        };

        let cumulative_lysed =
            evaluate_cumulative(self.critical_areal_strain_lysis).clamp(zero, one);
        let cumulative_necrotic =
            evaluate_cumulative(self.critical_areal_strain_necrosis).clamp(zero, one);
        let cumulative_permeabilized =
            evaluate_cumulative(self.critical_areal_strain_permeabilization).clamp(zero, one);

        let fraction_lysed = cumulative_lysed;
        let fraction_necrotic = (cumulative_necrotic - cumulative_lysed).max(zero);
        let fraction_permeabilized = (cumulative_permeabilized - cumulative_necrotic).max(zero);

        let subtotal = fraction_lysed + fraction_necrotic + fraction_permeabilized;
        let fraction_healthy = (one - subtotal).max(zero);

        // Ensure rigorous conservation normalization due to floating-point drift
        let total = fraction_healthy + fraction_permeabilized + fraction_necrotic + fraction_lysed;

        CellularInjuryProfile {
            fraction_healthy: fraction_healthy / total,
            fraction_permeabilized: fraction_permeabilized / total,
            fraction_necrotic: fraction_necrotic / total,
            fraction_lysed: fraction_lysed / total,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn mass_conservation_invariant(
            collapse_pressure in 100_000.0_f64..1_000_000_000.0,
            standoff_distance in 1e-6_f64..1e-3,
            radius in 1e-6_f64..50e-6,
            thickness in 10e-9_f64..100e-9,
            modulus in 10_000.0_f64..1_000_000.0,
            e_perm in 0.01_f64..0.1,
            e_necro in 0.15_f64..0.3,
            e_lysis in 0.35_f64..0.8
        ) {
            let mechanics = CellularMembraneMechanics::<f64> {
                cell_radius_m: radius,
                membrane_thickness_m: thickness,
                membrane_elastic_modulus_pa: modulus,
                critical_areal_strain_permeabilization: e_perm,
                critical_areal_strain_necrosis: e_necro,
                critical_areal_strain_lysis: e_lysis,
            };

            prop_assume!(mechanics.is_valid());

            let bubble_radius = 5e-6;
            let ambient = 101_325.0;

            let profile = mechanics.evaluate_spatial_injury(collapse_pressure, bubble_radius, standoff_distance, ambient);

            let eps = 1e-12;
            let total = profile.fraction_healthy + profile.fraction_permeabilized + profile.fraction_necrotic + profile.fraction_lysed;

            prop_assert!((total - 1.0).abs() < eps);
            prop_assert!(profile.fraction_healthy >= 0.0 && profile.fraction_healthy <= 1.0);
            prop_assert!(profile.fraction_permeabilized >= 0.0 && profile.fraction_permeabilized <= 1.0);
            prop_assert!(profile.fraction_necrotic >= 0.0 && profile.fraction_necrotic <= 1.0);
            prop_assert!(profile.fraction_lysed >= 0.0 && profile.fraction_lysed <= 1.0);
        }

        #[test]
        fn complete_healthy_at_ambient_collapse(
            e_perm in 0.01_f64..0.1,
            e_necro in 0.15_f64..0.3,
            e_lysis in 0.35_f64..0.8
        ) {
            let mechanics = CellularMembraneMechanics::<f64> {
                cell_radius_m: 10e-6,
                membrane_thickness_m: 50e-9,
                membrane_elastic_modulus_pa: 100_000.0,
                critical_areal_strain_permeabilization: e_perm,
                critical_areal_strain_necrosis: e_necro,
                critical_areal_strain_lysis: e_lysis,
            };
            prop_assume!(mechanics.is_valid());

            let ambient = 101_325.0;
            // No overpressure means no impact stress
            let profile = mechanics.evaluate_spatial_injury(ambient, 5e-6, 10e-6, ambient);
            prop_assert_eq!(profile.fraction_healthy, 1.0);
            prop_assert_eq!(profile.fraction_lysed, 0.0);
        }
    }

    #[test]
    fn injury_increases_monotonically_with_collapse_overpressure() {
        let mechanics = CellularMembraneMechanics::<f64> {
            cell_radius_m: 10e-6,
            membrane_thickness_m: 50e-9,
            membrane_elastic_modulus_pa: 100_000.0,
            critical_areal_strain_permeabilization: 0.05,
            critical_areal_strain_necrosis: 0.2,
            critical_areal_strain_lysis: 0.45,
        };
        let ambient = 101_325.0;
        let mild = mechanics.evaluate_spatial_injury(ambient + 10_000.0, 5e-6, 10e-6, ambient);
        let severe = mechanics.evaluate_spatial_injury(ambient + 2_000_000.0, 5e-6, 10e-6, ambient);

        assert!(severe.fraction_lysed >= mild.fraction_lysed);
        assert!(severe.fraction_healthy <= mild.fraction_healthy);
    }
}
