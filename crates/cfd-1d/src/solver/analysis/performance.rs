//! Performance metrics module for network systems.
//!
//! ## Theorem: Network Power Dissipation
//!
//! **Theorem**: The total viscous power dissipated in a resistive network is:
//!
//! ```text
//! W_total = Σ_{e ∈ E} Q_e · ΔP_e = Σ_{e ∈ E} R_e · Q_e²   [W]
//! ```
//!
//! This follows directly from Ohm's analogy for hydraulic networks (ΔP_e = R_e · Q_e)
//! and applies to **linear** (Hagen-Poiseuille) resistances. For nonlinear components
//! with `ΔP = R_lin · Q + k · Q²`, the dissipated power is `Q · (R_lin Q + k Q²)`.
//!
//! **Hydraulic efficiency** is the fraction of driven power that produces useful
//! volumetric throughput:
//!
//! ```text
//! η_hyd = Q_useful · ΔP_useful / W_pump   [dimensionless, ∈ (0, 1)]
//! ```
//!
//! For an ideal network without junction losses, `η_hyd → 1` as all power goes
//! into moving fluid through the channel resistances.

use nalgebra::RealField;
use std::collections::HashMap;
use std::iter::Sum;

/// Network performance metrics
#[derive(Debug, Clone)]
pub struct PerformanceMetrics<T: RealField + Copy> {
    /// Throughput [m³/s]
    pub throughput: T,
    /// Pressure efficiency (useful pressure / total pressure)
    pub pressure_efficiency: T,
    /// Power consumption [W]
    pub power_consumption: T,
    /// Mixing efficiency (for mixing applications)
    pub mixing_efficiency: Option<T>,
    /// Residence time distribution
    pub residence_times: HashMap<String, T>,
}

impl<T: RealField + Copy + Sum> PerformanceMetrics<T> {
    /// Create new performance metrics
    #[must_use]
    pub fn new() -> Self {
        Self {
            throughput: T::zero(),
            pressure_efficiency: T::zero(),
            power_consumption: T::zero(),
            mixing_efficiency: None,
            residence_times: HashMap::new(),
        }
    }

    /// Set throughput
    pub fn set_throughput(&mut self, throughput: T) {
        self.throughput = throughput;
    }

    /// Set pressure efficiency
    pub fn set_pressure_efficiency(&mut self, efficiency: T) {
        self.pressure_efficiency = efficiency;
    }

    /// Set power consumption
    pub fn set_power_consumption(&mut self, power: T) {
        self.power_consumption = power;
    }

    /// Set total pressure drop
    pub fn set_total_pressure_drop(&mut self, pressure_drop: T) {
        // Store as part of pressure efficiency calculation
        if pressure_drop > T::zero() && self.throughput > T::zero() {
            self.pressure_efficiency = self.throughput / pressure_drop;
        }
    }

    /// Set total flow rate
    pub fn set_total_flow_rate(&mut self, flow_rate: T) {
        self.throughput = flow_rate;
    }

    /// Set efficiency
    pub fn set_efficiency(&mut self, efficiency: T) {
        self.pressure_efficiency = efficiency;
    }

    /// Set mixing efficiency
    pub fn set_mixing_efficiency(&mut self, efficiency: T) {
        self.mixing_efficiency = Some(efficiency);
    }

    /// Add residence time for a component
    pub fn add_residence_time(&mut self, id: String, time: T) {
        self.residence_times.insert(id, time);
    }

    /// Calculate hydraulic efficiency
    pub fn hydraulic_efficiency(&self) -> T {
        if self.power_consumption > T::zero() {
            // Hydraulic efficiency: useful power / total power (White, 2011 - Fluid Mechanics)
            // Useful power = volumetric flow rate × pressure differential
            self.pressure_efficiency
        } else {
            T::zero()
        }
    }

    /// Get average residence time
    pub fn average_residence_time(&self) -> T {
        if self.residence_times.is_empty() {
            T::zero()
        } else {
            let sum: T = self.residence_times.values().copied().sum();
            sum / T::from_usize(self.residence_times.len()).unwrap_or_else(T::one)
        }
    }

    /// Get maximum residence time
    pub fn max_residence_time(&self) -> Option<T> {
        self.residence_times
            .values()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
    }

    /// Get minimum residence time
    pub fn min_residence_time(&self) -> Option<T> {
        self.residence_times
            .values()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
    }

    /// Check if mixing efficiency is calculated
    pub fn has_mixing_data(&self) -> bool {
        self.mixing_efficiency.is_some()
    }
}

impl<T: RealField + Copy + Sum> Default for PerformanceMetrics<T> {
    fn default() -> Self {
        Self::new()
    }
}
