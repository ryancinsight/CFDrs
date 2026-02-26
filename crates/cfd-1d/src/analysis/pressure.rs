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


use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;
use std::iter::Sum;

/// Pressure analysis for network systems
#[derive(Debug, Clone)]
pub struct PressureAnalysis<T: RealField + Copy> {
    /// Pressure distribution [Pa]
    pub pressures: HashMap<String, T>,
    /// Pressure drops across components [Pa]
    pub pressure_drops: HashMap<String, T>,
    /// Maximum pressure in system [Pa]
    pub max_pressure: T,
    /// Minimum pressure in system [Pa]
    pub min_pressure: T,
    /// Pressure gradient statistics
    pub pressure_gradients: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive + Sum> PressureAnalysis<T> {
    /// Create a new pressure analysis
    pub fn new() -> Self {
        Self {
            pressures: HashMap::new(),
            pressure_drops: HashMap::new(),
            max_pressure: T::from_f64(f64::NEG_INFINITY).unwrap_or_else(T::zero),
            min_pressure: T::from_f64(f64::INFINITY).unwrap_or_else(T::zero),
            pressure_gradients: HashMap::new(),
        }
    }

    /// Add pressure data for a node
    pub fn add_pressure(&mut self, id: String, pressure: T) {
        self.pressures.insert(id, pressure);

        if pressure > self.max_pressure {
            self.max_pressure = pressure;
        }
        if pressure < self.min_pressure {
            self.min_pressure = pressure;
        }
    }

    /// Add pressure drop data for a component
    pub fn add_pressure_drop(&mut self, id: String, drop: T) {
        self.pressure_drops.insert(id, drop);
    }

    /// Add pressure gradient data for a component
    pub fn add_pressure_gradient(&mut self, id: String, gradient: T) {
        self.pressure_gradients.insert(id, gradient);
    }

    /// Get the average pressure
    pub fn average_pressure(&self) -> T {
        if self.pressures.is_empty() {
            T::zero()
        } else {
            let sum: T = self.pressures.values().copied().sum();
            sum / T::from_usize(self.pressures.len()).unwrap_or_else(T::one)
        }
    }

    /// Get the total pressure drop
    pub fn total_pressure_drop(&self) -> T {
        self.pressure_drops.values().copied().sum()
    }

    /// Get the pressure range
    pub fn pressure_range(&self) -> T {
        if self.max_pressure > self.min_pressure {
            self.max_pressure - self.min_pressure
        } else {
            T::zero()
        }
    }

    /// Check if pressures are initialized
    pub fn is_initialized(&self) -> bool {
        !self.pressures.is_empty()
    }
}

impl<T: RealField + Copy + FromPrimitive + Sum> Default for PressureAnalysis<T> {
    fn default() -> Self {
        Self::new()
    }
}
