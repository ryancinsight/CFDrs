//! Venturi cavitation analysis.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Venturi cavitation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiCavitation<T: RealField + Copy> {
    /// Inlet diameter (m)
    pub inlet_diameter: T,
    /// Throat diameter (m)
    pub throat_diameter: T,
    /// Outlet diameter (m)
    pub outlet_diameter: T,
    /// Convergent angle (radians)
    pub convergent_angle: T,
    /// Divergent angle (radians)
    pub divergent_angle: T,
    /// Inlet pressure (Pa)
    pub inlet_pressure: T,
    /// Inlet velocity (m/s)
    pub inlet_velocity: T,
    /// Fluid density (kg/m³)
    pub density: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
}

impl<T: RealField + Copy + FromPrimitive> VenturiCavitation<T> {
    /// Calculate throat velocity using continuity equation
    pub fn throat_velocity(&self) -> T {
        let area_ratio = (self.inlet_diameter / self.throat_diameter).powi(2);
        self.inlet_velocity * area_ratio
    }

    /// Calculate throat pressure using Bernoulli equation
    pub fn throat_pressure(&self) -> T {
        let v_inlet = self.inlet_velocity;
        let v_throat = self.throat_velocity();
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));

        self.inlet_pressure - half * self.density * (v_throat * v_throat - v_inlet * v_inlet)
    }

    /// Calculate cavitation number at throat
    pub fn cavitation_number(&self) -> T {
        let p_throat = self.throat_pressure();
        let v_throat = self.throat_velocity();
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));

        if v_throat > T::zero() {
            (p_throat - self.vapor_pressure) / (half * self.density * v_throat * v_throat)
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one())
        }
    }

    /// Check if cavitation occurs
    pub fn is_cavitating(&self) -> bool {
        self.throat_pressure() < self.vapor_pressure
    }

    /// Calculate pressure recovery in diffuser
    pub fn outlet_pressure(&self, recovery_coefficient: T) -> T {
        let v_inlet = self.inlet_velocity;
        let v_outlet = self.outlet_velocity();
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));

        let ideal_recovery = half * self.density * (v_inlet * v_inlet - v_outlet * v_outlet);
        self.inlet_pressure + recovery_coefficient * ideal_recovery
    }

    /// Calculate outlet velocity
    pub fn outlet_velocity(&self) -> T {
        let area_ratio = (self.inlet_diameter / self.outlet_diameter).powi(2);
        self.inlet_velocity * area_ratio
    }

    /// Calculate loss coefficient
    pub fn loss_coefficient(&self, measured_outlet_pressure: T) -> T {
        let ideal_outlet = self.outlet_pressure(T::one());
        let actual_recovery = measured_outlet_pressure - self.inlet_pressure;
        let ideal_recovery = ideal_outlet - self.inlet_pressure;

        if ideal_recovery.abs() > T::from_f64(1e-10).unwrap_or_else(|| T::zero()) {
            T::one() - actual_recovery / ideal_recovery
        } else {
            T::zero()
        }
    }

    /// Calculate choked flow condition
    pub fn is_choked(&self) -> bool {
        let sigma = self.cavitation_number();
        let sigma_critical = T::from_f64(1.0).unwrap_or_else(|| T::one());
        sigma < sigma_critical
    }

    /// Estimate cavity length (empirical correlation)
    pub fn cavity_length(&self, cavitation_number: T) -> T {
        // Empirical correlation: L/D = f(σ)
        // Simplified model: L/D ≈ 1/σ for σ < 1
        if cavitation_number < T::one() && cavitation_number > T::zero() {
            self.throat_diameter / cavitation_number
        } else {
            T::zero()
        }
    }
}
