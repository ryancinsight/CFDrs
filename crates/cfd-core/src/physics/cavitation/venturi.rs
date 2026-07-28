//! Venturi cavitation analysis.

use aequitas::systems::si::quantities::{
    Angle, Dimensionless, Length, MassDensity, Pressure, Velocity, Volume,
};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Venturi cavitation parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiCavitation<T: FloatElement + Copy> {
    /// Inlet diameter (m)
    pub inlet_diameter: Length<T>,
    /// Throat diameter (m)
    pub throat_diameter: Length<T>,
    /// Outlet diameter (m)
    pub outlet_diameter: Length<T>,
    /// Convergent angle (radians)
    pub convergent_angle: Angle<T>,
    /// Divergent angle (radians)
    pub divergent_angle: Angle<T>,
    /// Inlet pressure (Pa)
    pub inlet_pressure: Pressure<T>,
    /// Inlet velocity (m/s)
    pub inlet_velocity: Velocity<T>,
    /// Fluid density (kg/m³)
    pub density: MassDensity<T>,
    /// Vapor pressure (Pa)
    pub vapor_pressure: Pressure<T>,
}

impl<T: FloatElement + Copy> VenturiCavitation<T> {
    /// Calculate throat velocity using continuity equation
    pub fn throat_velocity(&self) -> Velocity<T> {
        let area_ratio = <T as FloatElement>::powi(
            self.inlet_diameter.into_base() / self.throat_diameter.into_base(),
            2,
        );
        Velocity::from_base(self.inlet_velocity.into_base() * area_ratio)
    }

    /// Calculate throat pressure using Bernoulli equation
    pub fn throat_pressure(&self) -> Pressure<T> {
        let v_inlet = self.inlet_velocity.into_base();
        let v_throat = self.throat_velocity().into_base();
        let half = <T as FloatElement>::from_f64(0.5);

        Pressure::from_base(
            self.inlet_pressure.into_base()
                - half * self.density.into_base() * (v_throat * v_throat - v_inlet * v_inlet),
        )
    }

    /// Calculate cavitation number at throat
    pub fn cavitation_number(&self) -> Dimensionless<T> {
        let p_throat = self.throat_pressure().into_base();
        let v_throat = self.throat_velocity().into_base();
        let half = <T as FloatElement>::from_f64(0.5);

        if v_throat > <T as NumericElement>::ZERO {
            Dimensionless::from_base(
                (p_throat - self.vapor_pressure.into_base())
                    / (half * self.density.into_base() * v_throat * v_throat),
            )
        } else {
            Dimensionless::from_base(<T as FloatElement>::from_f64(1e10))
        }
    }

    /// Check if cavitation occurs
    pub fn is_cavitating(&self) -> bool {
        self.throat_pressure().into_base() < self.vapor_pressure.into_base()
    }

    /// Calculate pressure recovery in diffuser
    pub fn outlet_pressure(&self, recovery_coefficient: Dimensionless<T>) -> Pressure<T> {
        let v_inlet = self.inlet_velocity.into_base();
        let v_outlet = self.outlet_velocity().into_base();
        let half = <T as FloatElement>::from_f64(0.5);

        let ideal_recovery =
            half * self.density.into_base() * (v_inlet * v_inlet - v_outlet * v_outlet);
        Pressure::from_base(
            self.inlet_pressure.into_base() + recovery_coefficient.into_base() * ideal_recovery,
        )
    }

    /// Calculate outlet velocity
    pub fn outlet_velocity(&self) -> Velocity<T> {
        let area_ratio = <T as FloatElement>::powi(
            self.inlet_diameter.into_base() / self.outlet_diameter.into_base(),
            2,
        );
        Velocity::from_base(self.inlet_velocity.into_base() * area_ratio)
    }

    /// Calculate loss coefficient
    pub fn loss_coefficient(&self, measured_outlet_pressure: Pressure<T>) -> Dimensionless<T> {
        let ideal_outlet =
            self.outlet_pressure(Dimensionless::from_base(<T as NumericElement>::ONE));
        let actual_recovery =
            measured_outlet_pressure.into_base() - self.inlet_pressure.into_base();
        let ideal_recovery = ideal_outlet.into_base() - self.inlet_pressure.into_base();

        if <T as NumericElement>::abs(ideal_recovery) > <T as FloatElement>::from_f64(1e-10) {
            Dimensionless::from_base(<T as NumericElement>::ONE - actual_recovery / ideal_recovery)
        } else {
            Dimensionless::from_base(<T as NumericElement>::ZERO)
        }
    }

    /// Calculate choked flow condition
    pub fn is_choked(&self) -> bool {
        let sigma = self.cavitation_number().into_base();
        let sigma_critical = <T as FloatElement>::from_f64(super::constants::SIGMA_CRITICAL);
        sigma < sigma_critical
    }

    /// Calculate cavity length using Nurick correlation
    /// Based on Nurick (1976) for venturi cavitation
    /// L/D = K * (1/σ - `1/σ_i)^n` where `σ_i` is incipient cavitation number
    pub fn cavity_length(&self, cavitation_number: Dimensionless<T>) -> Length<T> {
        // Nurick correlation constants from literature
        let k_coefficient = <T as FloatElement>::from_f64(super::constants::NURICK_K_COEFFICIENT);
        let exponent = <T as FloatElement>::from_f64(super::constants::NURICK_EXPONENT);
        let sigma_incipient = <T as FloatElement>::from_f64(super::constants::SIGMA_INCIPIENT);

        let cavitation_number = cavitation_number.into_base();
        if cavitation_number < sigma_incipient && cavitation_number > <T as NumericElement>::ZERO {
            let term = <T as NumericElement>::ONE / cavitation_number
                - <T as NumericElement>::ONE / sigma_incipient;
            if term > <T as NumericElement>::ZERO {
                Length::from_base(
                    self.throat_diameter.into_base()
                        * k_coefficient
                        * <T as FloatElement>::powf(term, exponent),
                )
            } else {
                Length::from_base(<T as NumericElement>::ZERO)
            }
        } else {
            Length::from_base(<T as NumericElement>::ZERO)
        }
    }

    /// Calculate cavity closure position using Callenaere correlation
    /// Based on Callenaere et al. (2001) for cavity closure location
    pub fn cavity_closure_position(&self, cavitation_number: Dimensionless<T>) -> Length<T> {
        let cavity_len = self.cavity_length(cavitation_number);
        let divergence_factor = <T as FloatElement>::tan(self.divergent_angle.into_base());

        if divergence_factor > <T as NumericElement>::ZERO {
            // Closure position from throat
            Length::from_base(
                cavity_len.into_base()
                    + self.throat_diameter.into_base() * divergence_factor * cavity_len.into_base()
                        / (<T as NumericElement>::ONE + <T as NumericElement>::ONE),
            )
        } else {
            cavity_len
        }
    }

    /// Calculate cavity volume based on conical approximation
    pub fn cavity_volume(&self, cavitation_number: Dimensionless<T>) -> Volume<T> {
        let cavity_len = self.cavity_length(cavitation_number);
        let pi = <T as FloatElement>::from_f64(std::f64::consts::PI);
        let one_third = <T as NumericElement>::ONE
            / (<T as NumericElement>::ONE
                + <T as NumericElement>::ONE
                + <T as NumericElement>::ONE);

        // Conical cavity approximation: V = (π/3) * r² * L
        let radius = self.throat_diameter.into_base()
            / (<T as NumericElement>::ONE + <T as NumericElement>::ONE);
        Volume::from_base(one_third * pi * radius * radius * cavity_len.into_base())
    }
}

#[cfg(test)]
mod tests {
    use super::VenturiCavitation;
    use aequitas::systems::si::quantities::{
        Angle, Dimensionless, Length, MassDensity, Pressure, Velocity,
    };

    fn representative_venturi() -> VenturiCavitation<f64> {
        VenturiCavitation {
            inlet_diameter: Length::from_base(0.05),
            throat_diameter: Length::from_base(0.02),
            outlet_diameter: Length::from_base(0.05),
            convergent_angle: Angle::from_base(0.0),
            divergent_angle: Angle::from_base(0.0),
            inlet_pressure: Pressure::from_base(300_000.0),
            inlet_velocity: Velocity::from_base(3.0),
            density: MassDensity::from_base(998.2),
            vapor_pressure: Pressure::from_base(2339.0),
        }
    }

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn throat_velocity_matches_continuity() {
        let venturi = representative_venturi();

        let expected = 3.0_f64 * (0.05_f64 / 0.02_f64).powi(2);
        assert!(approx_eq(
            venturi.throat_velocity().into_base(),
            expected,
            1e-12
        ));
    }

    #[test]
    fn throat_pressure_matches_bernoulli() {
        let venturi = representative_venturi();

        let v1 = 3.0;
        let v2 = venturi.throat_velocity().into_base();
        let expected = 300_000.0 - 0.5 * 998.2 * (v2 * v2 - v1 * v1);
        assert!(approx_eq(
            venturi.throat_pressure().into_base(),
            expected,
            1e-9
        ));
    }

    #[test]
    fn cavitation_number_matches_definition() {
        let venturi = representative_venturi();

        let p_throat = venturi.throat_pressure().into_base();
        let v_throat = venturi.throat_velocity().into_base();
        let expected = (p_throat - 2339.0) / (0.5 * 998.2 * v_throat * v_throat);
        assert!(approx_eq(
            venturi.cavitation_number().into_base(),
            expected,
            1e-12
        ));
        assert!(venturi.cavitation_number().into_base().is_finite());
    }

    #[test]
    fn cavity_length_is_zero_above_incipient_and_grows_as_sigma_drops() {
        let venturi = representative_venturi();

        let l_non = venturi.cavity_length(Dimensionless::from_base(1.3));
        let l_weak = venturi.cavity_length(Dimensionless::from_base(1.1));
        let l_strong = venturi.cavity_length(Dimensionless::from_base(0.6));

        assert!(approx_eq(l_non.into_base(), 0.0, 0.0));
        assert!(l_weak.into_base() > 0.0);
        assert!(l_strong.into_base() > l_weak.into_base());
    }

    #[test]
    fn cavity_volume_matches_conical_approximation() {
        let venturi = representative_venturi();

        let sigma = 0.6;
        let l = venturi
            .cavity_length(Dimensionless::from_base(sigma))
            .into_base();
        let r = 0.02 * 0.5;
        let expected = (std::f64::consts::PI / 3.0) * r * r * l;
        assert!(approx_eq(
            venturi
                .cavity_volume(Dimensionless::from_base(sigma))
                .into_base(),
            expected,
            1e-15
        ));
    }
}
