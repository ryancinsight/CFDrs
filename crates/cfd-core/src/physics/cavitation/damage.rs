//! Cavitation damage and erosion models.

use aequitas::systems::si::quantities::{Frequency, Length, MassDensity, Pressure, Time, Velocity};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Cavitation damage and erosion model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationDamage<T: FloatElement + Copy> {
    /// Material yield strength (Pa)
    pub yield_strength: Pressure<T>,
    /// Material ultimate tensile strength (Pa)
    pub ultimate_strength: Pressure<T>,
    /// Material hardness (Pa)
    pub hardness: Pressure<T>,
    /// Fatigue strength (Pa)
    pub fatigue_strength: Pressure<T>,
    /// Number of loading cycles
    pub cycles: u64,
}

impl<T: FloatElement + Copy> CavitationDamage<T> {
    /// Calculate erosion rate using Hammitt's correlation
    /// Rate ∝ (ΔP)^n where n ≈ 2-3
    pub fn erosion_rate_hammitt(&self, pressure_amplitude: Pressure<T>, exponent: T) -> T {
        let reference_pressure = self.yield_strength.into_base();
        if pressure_amplitude.into_base() > reference_pressure {
            let normalized_pressure = pressure_amplitude.into_base() / reference_pressure;
            <T as FloatElement>::powf(normalized_pressure, exponent)
        } else {
            <T as NumericElement>::ZERO
        }
    }

    /// Calculate Mean Depth of Penetration Rate (MDPR)
    /// Based on ASTM G32 standard and Plesset-Chapman erosion model
    pub fn mdpr(
        &self,
        impact_pressure: Pressure<T>,
        frequency: Frequency<T>,
        exposure_time: Time<T>,
    ) -> Velocity<T> {
        // Plesset-Chapman model: MDPR = K * (P - P_th)^n * f^m * t
        // where K is material-dependent erosion coefficient
        // n ≈ 2.0-2.5 for most materials (pressure exponent)
        // m ≈ 1.0 for frequency dependency

        let impact_pressure = impact_pressure.into_base();
        let frequency = frequency.into_base();
        let exposure_time = exposure_time.into_base();
        let threshold = self.fatigue_strength.into_base();

        if impact_pressure > threshold {
            // Material-specific erosion coefficient from ASTM G32 data
            // K varies from 1e-15 to 1e-11 m³/N²·Hz·s for different materials
            let erosion_coefficient =
                <T as FloatElement>::from_f64(super::constants::EROSION_COEFFICIENT_STEEL);
            let pressure_exponent =
                <T as FloatElement>::from_f64(super::constants::EROSION_PRESSURE_EXPONENT);

            // MDPR calculation with proper units
            let pressure_term =
                <T as FloatElement>::powf(impact_pressure - threshold, pressure_exponent);
            Velocity::from_base(erosion_coefficient * pressure_term * frequency * exposure_time)
        } else {
            Velocity::from_base(<T as NumericElement>::ZERO)
        }
    }

    /// Calculate impact pressure from bubble collapse (Rayleigh)
    pub fn collapse_impact_pressure(
        &self,
        bubble_radius: Length<T>,
        collapse_distance: Length<T>,
        liquid_density: MassDensity<T>,
        sound_speed: Velocity<T>,
    ) -> Pressure<T> {
        // Rayleigh collapse pressure: P ∝ ρc²(R/r)
        if collapse_distance.into_base() > <T as NumericElement>::ZERO {
            Pressure::from_base(
                liquid_density.into_base()
                    * sound_speed.into_base()
                    * sound_speed.into_base()
                    * bubble_radius.into_base()
                    / collapse_distance.into_base(),
            )
        } else {
            Pressure::from_base(<T as NumericElement>::ZERO)
        }
    }

    /// Calculate incubation period before damage onset using Basquin's law
    pub fn incubation_period(&self, stress_amplitude: Pressure<T>) -> u64 {
        // Maximum safe f64 value that can be precisely represented as u64
        const MAX_SAFE_F64: f64 = 9_007_199_254_740_992.0; // 2^53

        // Basquin's law for high-cycle fatigue: σ_a = σ_f' * (2N)^b
        // Rearranged: N = 0.5 * ((σ_a / σ_f')^(1/b))
        // where σ_f' is fatigue strength coefficient ≈ 0.9 * UTS
        // b is fatigue strength exponent (typically -0.085 to -0.12 for steels)

        // Fatigue strength coefficient (Morrow approximation)
        let fatigue_coeff = <T as FloatElement>::from_f64(super::constants::FATIGUE_STRENGTH_RATIO)
            * self.ultimate_strength.into_base();

        // Basquin exponent (typical for metals)
        let basquin_exponent =
            <T as FloatElement>::from_f64(super::constants::BASQUIN_EXPONENT_STEEL);

        if stress_amplitude.into_base() > self.fatigue_strength.into_base() {
            // Calculate cycles to failure
            let stress_ratio = stress_amplitude.into_base() / fatigue_coeff;
            let exponent = <T as NumericElement>::ONE / basquin_exponent;
            let two_n = <T as FloatElement>::powf(stress_ratio, exponent);
            let n = two_n / <T as FloatElement>::from_f64(2.0);

            // Safe conversion to u64 with proper bounds checking
            if n > <T as NumericElement>::ZERO && n < <T as FloatElement>::from_f64(u64::MAX as f64)
            {
                let n_f64 = <T as NumericElement>::to_f64(n);
                // Use clamp for safe conversion at the integer API boundary.
                #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
                let result = n_f64.round().clamp(0.0, MAX_SAFE_F64) as u64;
                result
            } else {
                1 // Minimum one cycle if calculation fails
            }
        } else {
            u64::MAX // Infinite life below fatigue limit
        }
    }

    /// Calculate cumulative damage using Miner's rule
    pub fn cumulative_damage(&self, stress_levels: &[(Pressure<T>, u64)]) -> T {
        let mut damage = <T as NumericElement>::ZERO;

        for &(stress, cycles) in stress_levels {
            let life_cycles = self.incubation_period(stress);
            if life_cycles < u64::MAX {
                damage += <T as FloatElement>::from_f64(cycles as f64)
                    / <T as FloatElement>::from_f64(life_cycles as f64);
            }
        }

        damage
    }

    /// Check if material has failed
    pub fn has_failed(&self, cumulative_damage: T) -> bool {
        cumulative_damage >= <T as NumericElement>::ONE
    }

    /// Calculate erosion resistance parameter
    pub fn erosion_resistance(&self) -> T {
        // Empirical correlation with material properties
        // ER ∝ (H × UTS × FS)^(1/3)
        let one_third = <T as FloatElement>::from_f64(1.0 / 3.0);
        <T as FloatElement>::powf(
            self.hardness.into_base()
                * self.ultimate_strength.into_base()
                * self.fatigue_strength.into_base(),
            one_third,
        )
    }

    /// Estimate pit depth from impact pressure
    pub fn pit_depth(&self, impact_pressure: Pressure<T>, material_constant: T) -> Length<T> {
        // Empirical model: depth ∝ (P/H)^n
        if impact_pressure.into_base() > self.hardness.into_base() {
            let ratio = impact_pressure.into_base() / self.hardness.into_base();
            Length::from_base(
                material_constant
                    * <T as FloatElement>::powf(ratio, <T as FloatElement>::from_f64(2.0)),
            )
        } else {
            Length::from_base(<T as NumericElement>::ZERO)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::CavitationDamage;
    use aequitas::systems::si::quantities::{Length, MassDensity, Pressure, Velocity};

    fn steel_like_damage() -> CavitationDamage<f64> {
        CavitationDamage {
            yield_strength: Pressure::from_base(250.0e6),
            ultimate_strength: Pressure::from_base(400.0e6),
            hardness: Pressure::from_base(1.2e9),
            fatigue_strength: Pressure::from_base(160.0e6),
            cycles: 0,
        }
    }

    #[test]
    fn hammitt_erosion_uses_pressure_ratio_power() {
        let damage = steel_like_damage();

        assert_eq!(
            damage.erosion_rate_hammitt(Pressure::from_base(200.0e6), 2.0),
            0.0
        );
        assert!(
            (damage.erosion_rate_hammitt(Pressure::from_base(500.0e6), 2.0) - 4.0).abs() <= 1.0e-12
        );
    }

    #[test]
    fn collapse_impact_pressure_matches_rayleigh_scaling() {
        let damage = steel_like_damage();
        let pressure = damage.collapse_impact_pressure(
            Length::from_base(1.0e-6),
            Length::from_base(2.0e-6),
            MassDensity::from_base(1_000.0),
            Velocity::from_base(1_500.0),
        );
        let expected = 1_000.0 * 1_500.0 * 1_500.0 * 1.0e-6 / 2.0e-6;

        assert!((pressure.into_base() - expected).abs() <= 1.0e-6);
    }

    #[test]
    fn pit_depth_uses_hardness_normalized_quadratic() {
        let damage = steel_like_damage();
        let depth = damage.pit_depth(Pressure::from_base(2.4e9), 3.0e-9);

        assert!((depth.into_base() - 12.0e-9).abs() <= 1.0e-18);
    }
}
