//! Cavitation damage and erosion models.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Cavitation damage and erosion model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationDamage<T: RealField + Copy> {
    /// Material yield strength (Pa)
    pub yield_strength: T,
    /// Material ultimate tensile strength (Pa)
    pub ultimate_strength: T,
    /// Material hardness (Pa)
    pub hardness: T,
    /// Fatigue strength (Pa)
    pub fatigue_strength: T,
    /// Number of loading cycles
    pub cycles: u64,
}

impl<T: RealField + Copy + FromPrimitive> CavitationDamage<T> {
    /// Calculate erosion rate using Hammitt's correlation
    /// Rate ∝ (ΔP)^n where n ≈ 2-3
    pub fn erosion_rate_hammitt(&self, pressure_amplitude: T, exponent: T) -> T {
        let reference_pressure = self.yield_strength;
        if pressure_amplitude > reference_pressure {
            let normalized_pressure = pressure_amplitude / reference_pressure;
            normalized_pressure.powf(exponent)
        } else {
            T::zero()
        }
    }

    /// Calculate Mean Depth of Penetration Rate (MDPR)
    /// Based on ASTM G32 standard
    pub fn mdpr(&self, impact_pressure: T, frequency: T, exposure_time: T) -> T {
        // Simplified model: MDPR ∝ (P - P_threshold)² × f × t
        let threshold = self.fatigue_strength;

        if impact_pressure > threshold {
            let coefficient = T::from_f64(1e-12).unwrap_or_else(|| T::zero()); // Material-dependent
            coefficient * (impact_pressure - threshold).powi(2) * frequency * exposure_time
        } else {
            T::zero()
        }
    }

    /// Calculate impact pressure from bubble collapse (Rayleigh)
    pub fn collapse_impact_pressure(
        &self,
        bubble_radius: T,
        collapse_distance: T,
        liquid_density: T,
        sound_speed: T,
    ) -> T {
        // Rayleigh collapse pressure: P ∝ ρc²(R/r)
        if collapse_distance > T::zero() {
            liquid_density * sound_speed * sound_speed * bubble_radius / collapse_distance
        } else {
            T::zero()
        }
    }

    /// Calculate incubation period before damage onset
    pub fn incubation_period(&self, stress_amplitude: T) -> u64 {
        // S-N curve approximation: N = C × S^(-m)
        let m = T::from_f64(3.0).unwrap_or_else(|| T::one()); // Material constant
        let c = T::from_f64(1e12).unwrap_or_else(|| T::one()); // Material constant

        if stress_amplitude > self.fatigue_strength {
            let ratio = stress_amplitude / self.fatigue_strength;
            let n = c / ratio.powf(m);
            // Convert to u64 - using a safe conversion
            // Note: This is a simplification - real implementation would need proper rounding
            let n_f64 = T::to_subset(&n).unwrap_or(0.0_f64);
            n_f64 as u64
        } else {
            u64::MAX // Infinite life below fatigue limit
        }
    }

    /// Calculate cumulative damage using Miner's rule
    pub fn cumulative_damage(&self, stress_levels: &[(T, u64)]) -> T {
        let mut damage = T::zero();

        for &(stress, cycles) in stress_levels {
            let life_cycles = self.incubation_period(stress);
            if life_cycles < u64::MAX {
                damage = damage
                    + T::from_u64(cycles).unwrap_or_else(|| T::zero())
                        / T::from_u64(life_cycles).unwrap_or_else(|| T::one());
            }
        }

        damage
    }

    /// Check if material has failed
    pub fn has_failed(&self, cumulative_damage: T) -> bool {
        cumulative_damage >= T::one()
    }

    /// Calculate erosion resistance parameter
    pub fn erosion_resistance(&self) -> T {
        // Empirical correlation with material properties
        // ER ∝ (H × UTS × FS)^(1/3)
        let one_third = T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one());
        (self.hardness * self.ultimate_strength * self.fatigue_strength).powf(one_third)
    }

    /// Estimate pit depth from impact pressure
    pub fn pit_depth(&self, impact_pressure: T, material_constant: T) -> T {
        // Empirical model: depth ∝ (P/H)^n
        if impact_pressure > self.hardness {
            let ratio = impact_pressure / self.hardness;
            material_constant * ratio.powf(T::from_f64(2.0).unwrap_or_else(|| T::one()))
        } else {
            T::zero()
        }
    }
}
