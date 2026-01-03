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
        let sigma_critical =
            T::from_f64(super::constants::SIGMA_CRITICAL).unwrap_or_else(|| T::one());
        sigma < sigma_critical
    }

    /// Calculate cavity length using Nurick correlation
    /// Based on Nurick (1976) for venturi cavitation
    /// L/D = K * (1/σ - `1/σ_i)^n` where `σ_i` is incipient cavitation number
    pub fn cavity_length(&self, cavitation_number: T) -> T {
        // Nurick correlation constants from literature
        let k_coefficient =
            T::from_f64(super::constants::NURICK_K_COEFFICIENT).unwrap_or_else(|| T::one());
        let exponent = T::from_f64(super::constants::NURICK_EXPONENT)
            .unwrap_or_else(|| T::one() / (T::one() + T::one()));
        let sigma_incipient =
            T::from_f64(super::constants::SIGMA_INCIPIENT).unwrap_or_else(|| T::one());

        if cavitation_number < sigma_incipient && cavitation_number > T::zero() {
            let term = T::one() / cavitation_number - T::one() / sigma_incipient;
            if term > T::zero() {
                self.throat_diameter * k_coefficient * term.powf(exponent)
            } else {
                T::zero()
            }
        } else {
            T::zero()
        }
    }

    /// Calculate cavity closure position using Callenaere correlation
    /// Based on Callenaere et al. (2001) for cavity closure location
    pub fn cavity_closure_position(&self, cavitation_number: T) -> T {
        let cavity_len = self.cavity_length(cavitation_number);
        let divergence_factor = self.divergent_angle.tan();

        if divergence_factor > T::zero() {
            // Closure position from throat
            cavity_len
                + self.throat_diameter * divergence_factor * cavity_len / (T::one() + T::one())
        } else {
            cavity_len
        }
    }

    /// Calculate cavity volume based on conical approximation
    pub fn cavity_volume(&self, cavitation_number: T) -> T {
        let cavity_len = self.cavity_length(cavitation_number);
        let pi = T::from_f64(std::f64::consts::PI)
            .unwrap_or_else(|| T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::one()));
        let one_third = T::one() / (T::one() + T::one() + T::one());

        // Conical cavity approximation: V = (π/3) * r² * L
        let radius = self.throat_diameter / (T::one() + T::one());
        one_third * pi * radius * radius * cavity_len
    }
}
