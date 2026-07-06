//! Venturi cavitation analysis.

use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Venturi cavitation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiCavitation<T: FloatElement + Copy> {
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

impl<T: FloatElement + Copy> VenturiCavitation<T> {
    /// Calculate throat velocity using continuity equation
    pub fn throat_velocity(&self) -> T {
        let area_ratio = <T as FloatElement>::powi(self.inlet_diameter / self.throat_diameter, 2);
        self.inlet_velocity * area_ratio
    }

    /// Calculate throat pressure using Bernoulli equation
    pub fn throat_pressure(&self) -> T {
        let v_inlet = self.inlet_velocity;
        let v_throat = self.throat_velocity();
        let half = <T as FloatElement>::from_f64(0.5);

        self.inlet_pressure - half * self.density * (v_throat * v_throat - v_inlet * v_inlet)
    }

    /// Calculate cavitation number at throat
    pub fn cavitation_number(&self) -> T {
        let p_throat = self.throat_pressure();
        let v_throat = self.throat_velocity();
        let half = <T as FloatElement>::from_f64(0.5);

        if v_throat > <T as NumericElement>::ZERO {
            (p_throat - self.vapor_pressure) / (half * self.density * v_throat * v_throat)
        } else {
            <T as FloatElement>::from_f64(1e10)
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
        let half = <T as FloatElement>::from_f64(0.5);

        let ideal_recovery = half * self.density * (v_inlet * v_inlet - v_outlet * v_outlet);
        self.inlet_pressure + recovery_coefficient * ideal_recovery
    }

    /// Calculate outlet velocity
    pub fn outlet_velocity(&self) -> T {
        let area_ratio = <T as FloatElement>::powi(self.inlet_diameter / self.outlet_diameter, 2);
        self.inlet_velocity * area_ratio
    }

    /// Calculate loss coefficient
    pub fn loss_coefficient(&self, measured_outlet_pressure: T) -> T {
        let ideal_outlet = self.outlet_pressure(<T as NumericElement>::ONE);
        let actual_recovery = measured_outlet_pressure - self.inlet_pressure;
        let ideal_recovery = ideal_outlet - self.inlet_pressure;

        if <T as NumericElement>::abs(ideal_recovery) > <T as FloatElement>::from_f64(1e-10) {
            <T as NumericElement>::ONE - actual_recovery / ideal_recovery
        } else {
            <T as NumericElement>::ZERO
        }
    }

    /// Calculate choked flow condition
    pub fn is_choked(&self) -> bool {
        let sigma = self.cavitation_number();
        let sigma_critical = <T as FloatElement>::from_f64(super::constants::SIGMA_CRITICAL);
        sigma < sigma_critical
    }

    /// Calculate cavity length using Nurick correlation
    /// Based on Nurick (1976) for venturi cavitation
    /// L/D = K * (1/σ - `1/σ_i)^n` where `σ_i` is incipient cavitation number
    pub fn cavity_length(&self, cavitation_number: T) -> T {
        // Nurick correlation constants from literature
        let k_coefficient = <T as FloatElement>::from_f64(super::constants::NURICK_K_COEFFICIENT);
        let exponent = <T as FloatElement>::from_f64(super::constants::NURICK_EXPONENT);
        let sigma_incipient = <T as FloatElement>::from_f64(super::constants::SIGMA_INCIPIENT);

        if cavitation_number < sigma_incipient && cavitation_number > <T as NumericElement>::ZERO {
            let term = <T as NumericElement>::ONE / cavitation_number
                - <T as NumericElement>::ONE / sigma_incipient;
            if term > <T as NumericElement>::ZERO {
                self.throat_diameter * k_coefficient * <T as FloatElement>::powf(term, exponent)
            } else {
                <T as NumericElement>::ZERO
            }
        } else {
            <T as NumericElement>::ZERO
        }
    }

    /// Calculate cavity closure position using Callenaere correlation
    /// Based on Callenaere et al. (2001) for cavity closure location
    pub fn cavity_closure_position(&self, cavitation_number: T) -> T {
        let cavity_len = self.cavity_length(cavitation_number);
        let divergence_factor = <T as FloatElement>::tan(self.divergent_angle);

        if divergence_factor > <T as NumericElement>::ZERO {
            // Closure position from throat
            cavity_len
                + self.throat_diameter * divergence_factor * cavity_len
                    / (<T as NumericElement>::ONE + <T as NumericElement>::ONE)
        } else {
            cavity_len
        }
    }

    /// Calculate cavity volume based on conical approximation
    pub fn cavity_volume(&self, cavitation_number: T) -> T {
        let cavity_len = self.cavity_length(cavitation_number);
        let pi = <T as FloatElement>::from_f64(std::f64::consts::PI);
        let one_third = <T as NumericElement>::ONE
            / (<T as NumericElement>::ONE
                + <T as NumericElement>::ONE
                + <T as NumericElement>::ONE);

        // Conical cavity approximation: V = (π/3) * r² * L
        let radius =
            self.throat_diameter / (<T as NumericElement>::ONE + <T as NumericElement>::ONE);
        one_third * pi * radius * radius * cavity_len
    }
}

#[cfg(test)]
mod tests {
    use super::VenturiCavitation;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn throat_velocity_matches_continuity() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.05,
            throat_diameter: 0.02,
            outlet_diameter: 0.05,
            convergent_angle: 0.0,
            divergent_angle: 0.0,
            inlet_pressure: 300_000.0,
            inlet_velocity: 3.0,
            density: 998.2,
            vapor_pressure: 2339.0,
        };

        let expected = 3.0_f64 * (0.05_f64 / 0.02_f64).powi(2);
        assert!(approx_eq(venturi.throat_velocity(), expected, 1e-12));
    }

    #[test]
    fn throat_pressure_matches_bernoulli() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.05,
            throat_diameter: 0.02,
            outlet_diameter: 0.05,
            convergent_angle: 0.0,
            divergent_angle: 0.0,
            inlet_pressure: 300_000.0,
            inlet_velocity: 3.0,
            density: 998.2,
            vapor_pressure: 2339.0,
        };

        let v1 = 3.0;
        let v2 = venturi.throat_velocity();
        let expected = 300_000.0 - 0.5 * 998.2 * (v2 * v2 - v1 * v1);
        assert!(approx_eq(venturi.throat_pressure(), expected, 1e-9));
    }

    #[test]
    fn cavitation_number_matches_definition() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.05,
            throat_diameter: 0.02,
            outlet_diameter: 0.05,
            convergent_angle: 0.0,
            divergent_angle: 0.0,
            inlet_pressure: 300_000.0,
            inlet_velocity: 3.0,
            density: 998.2,
            vapor_pressure: 2339.0,
        };

        let p_throat = venturi.throat_pressure();
        let v_throat = venturi.throat_velocity();
        let expected = (p_throat - 2339.0) / (0.5 * 998.2 * v_throat * v_throat);
        assert!(approx_eq(venturi.cavitation_number(), expected, 1e-12));
        assert!(venturi.cavitation_number().is_finite());
    }

    #[test]
    fn cavity_length_is_zero_above_incipient_and_grows_as_sigma_drops() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.05,
            throat_diameter: 0.02,
            outlet_diameter: 0.05,
            convergent_angle: 0.0,
            divergent_angle: 0.0,
            inlet_pressure: 300_000.0,
            inlet_velocity: 3.0,
            density: 998.2,
            vapor_pressure: 2339.0,
        };

        let l_non = venturi.cavity_length(1.3);
        let l_weak = venturi.cavity_length(1.1);
        let l_strong = venturi.cavity_length(0.6);

        assert!(approx_eq(l_non, 0.0, 0.0));
        assert!(l_weak > 0.0);
        assert!(l_strong > l_weak);
    }

    #[test]
    fn cavity_volume_matches_conical_approximation() {
        let venturi = VenturiCavitation::<f64> {
            inlet_diameter: 0.05,
            throat_diameter: 0.02,
            outlet_diameter: 0.05,
            convergent_angle: 0.0,
            divergent_angle: 0.0,
            inlet_pressure: 300_000.0,
            inlet_velocity: 3.0,
            density: 998.2,
            vapor_pressure: 2339.0,
        };

        let sigma = 0.6;
        let l = venturi.cavity_length(sigma);
        let r = 0.02 * 0.5;
        let expected = (std::f64::consts::PI / 3.0) * r * r * l;
        assert!(approx_eq(venturi.cavity_volume(sigma), expected, 1e-15));
    }
}
