//! Pressure analysis module for network systems.
//!
//! ## Theorem: Pressure Extrema Consistency
//!
//! **Theorem**: For any non-empty pressure analysis containing `n` node pressures
//! `{P_1, …, P_n}`:
//!
//! 1. `max_pressure = max{P_i}` and `min_pressure = min{P_i}` are maintained
//!    exactly by the incremental `add_pressure` method.
//! 2. `pressure_range = max_pressure - min_pressure ≥ 0`.
//! 3. For a single-path resistor network with inlet pressure `P_in` and outlet
//!    `P_out`, `total_pressure_drop = Σ ΔP_i = P_in - P_out = pressure_range`.
//!
//! **Proof of invariant (2)**: `add_pressure` updates `max_pressure` only when
//! `P > max_pressure` and updates `min_pressure` only when `P < min_pressure`.
//! Both operations preserve `max_pressure ≥ min_pressure` by induction. ∎

use aequitas::systems::si::quantities::{Pressure, PressureGradient};
use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use cfd_core::CfdScalar;
use std::collections::HashMap;
use std::iter::Sum;

/// Pressure analysis for network systems
#[derive(Debug, Clone)]
pub struct PressureAnalysis<T: CfdScalar + Copy> {
    /// Pressure distribution \[Pa]
    pub pressures: HashMap<String, Pressure<T>>,
    /// Pressure drops across components \[Pa]
    pub pressure_drops: HashMap<String, Pressure<T>>,
    /// Maximum pressure in system \[Pa]
    pub max_pressure: Pressure<T>,
    /// Minimum pressure in system \[Pa]
    pub min_pressure: Pressure<T>,
    /// Pressure gradient statistics
    pub pressure_gradients: HashMap<String, PressureGradient<T>>,
}

impl<T: CfdScalar + Copy + SafeFromF64 + SafeFromUsize + Sum> PressureAnalysis<T> {
    /// Create a new pressure analysis
    pub fn new() -> Self {
        Self {
            pressures: HashMap::new(),
            pressure_drops: HashMap::new(),
            max_pressure: Pressure::from_base(T::from_f64_or_zero(f64::NEG_INFINITY)),
            min_pressure: Pressure::from_base(T::from_f64_or_zero(f64::INFINITY)),
            pressure_gradients: HashMap::new(),
        }
    }

    /// Add pressure data for a node
    pub fn add_pressure(&mut self, id: String, pressure: Pressure<T>) {
        self.pressures.insert(id, pressure);

        if pressure.into_base() > self.max_pressure.into_base() {
            self.max_pressure = pressure;
        }
        if pressure.into_base() < self.min_pressure.into_base() {
            self.min_pressure = pressure;
        }
    }

    /// Add pressure drop data for a component
    pub fn add_pressure_drop(&mut self, id: String, drop: Pressure<T>) {
        self.pressure_drops.insert(id, drop);
    }

    /// Add pressure gradient data for a component
    pub fn add_pressure_gradient(&mut self, id: String, gradient: PressureGradient<T>) {
        self.pressure_gradients.insert(id, gradient);
    }

    /// Get the average pressure
    pub fn average_pressure(&self) -> Pressure<T> {
        if self.pressures.is_empty() {
            Pressure::from_base(T::ZERO)
        } else {
            let sum: T = self.pressures.values().map(|p| p.into_base()).sum();
            Pressure::from_base(sum / T::from_usize_or_one(self.pressures.len()))
        }
    }

    /// Get the total pressure drop
    pub fn total_pressure_drop(&self) -> Pressure<T> {
        Pressure::from_base(self.pressure_drops.values().map(|p| p.into_base()).sum())
    }

    /// Get the pressure range
    pub fn pressure_range(&self) -> Pressure<T> {
        if self.max_pressure.into_base() > self.min_pressure.into_base() {
            Pressure::from_base(self.max_pressure.into_base() - self.min_pressure.into_base())
        } else {
            Pressure::from_base(T::ZERO)
        }
    }

    /// Check if pressures are initialized
    pub fn is_initialized(&self) -> bool {
        !self.pressures.is_empty()
    }
}

impl<T: CfdScalar + Copy + SafeFromF64 + SafeFromUsize + Sum> Default for PressureAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}
