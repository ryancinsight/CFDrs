//! Fluid dynamics service following Domain Service pattern
//!
//! This service coordinates complex fluid mechanics operations and
//! business logic that doesn't naturally fit into a single domain entity.

use crate::error::Result;
use crate::physics::fluid::ConstantPropertyFluid;
use crate::physics::fluid_dynamics::flow_regimes::{FlowClassifier, FlowRegime};
use eunomia::{FloatElement, NumericElement, RealField};

/// Service for fluid dynamics calculations
pub struct FluidDynamicsService;

impl FluidDynamicsService {
    /// Calculate Reynolds number for a given flow configuration
    pub fn reynolds_number<T: RealField + Copy>(
        fluid: &ConstantPropertyFluid<T>,
        velocity: T,
        characteristic_length: T,
    ) -> T {
        let kinematic_viscosity = fluid.viscosity / fluid.density;
        velocity * characteristic_length / kinematic_viscosity
    }

    /// Calculate Prandtl number for constant property fluid
    pub fn prandtl_number<T: RealField + Copy>(fluid: &ConstantPropertyFluid<T>) -> T {
        let cp = fluid.specific_heat;
        let k = fluid.thermal_conductivity;
        let mu = fluid.viscosity;
        mu * cp / k
    }

    /// Determine flow regime based on Reynolds number
    pub fn flow_regime<T: RealField>(reynolds: T) -> FlowRegime {
        FlowClassifier::classify_by_reynolds(reynolds)
    }

    /// Calculate pressure drop for pipe flow
    pub fn pipe_pressure_drop<T: RealField + FloatElement + Copy>(
        fluid: &ConstantPropertyFluid<T>,
        velocity: T,
        length: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let reynolds = Self::reynolds_number(fluid, velocity, diameter);
        let friction_factor = Self::friction_factor(reynolds, diameter, roughness)?;

        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;

        Ok(friction_factor * length * fluid.density * velocity * velocity / (two * diameter))
    }

    /// Calculate friction factor using appropriate correlation
    fn friction_factor<T: RealField + FloatElement + Copy>(
        reynolds: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let re_2300 = <T as FloatElement>::from_f64(
            crate::physics::constants::physics::dimensionless::reynolds::PIPE_LAMINAR_MAX,
        );
        let sixty_four = <T as FloatElement>::from_f64(64.0);

        if reynolds < re_2300 {
            // Laminar flow: f = 64/Re
            Ok(sixty_four / reynolds)
        } else if let Some(eps) = roughness {
            // Colebrook-White equation
            let relative_roughness = eps / diameter;
            Self::colebrook_white_friction_factor(reynolds, relative_roughness)
        } else {
            // Smooth pipe: Blasius (low Re) or Haaland (general explicit)
            let blasius_max = <T as FloatElement>::from_f64(
                crate::physics::constants::physics::hydraulics::BLASIUS_MAX_RE,
            );
            if reynolds < blasius_max {
                let coeff = <T as FloatElement>::from_f64(
                    crate::physics::constants::physics::hydraulics::BLASIUS_COEFFICIENT,
                );
                let exp = <T as FloatElement>::from_f64(
                    crate::physics::constants::physics::hydraulics::BLASIUS_EXPONENT,
                );
                Ok(coeff / <T as FloatElement>::powf(reynolds, exp))
            } else {
                Self::haaland_friction_factor(reynolds, <T as NumericElement>::ZERO)
            }
        }
    }

    /// Calculate friction factor using Colebrook-White equation (iterative)
    fn colebrook_white_friction_factor<T: RealField + FloatElement + Copy>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        // Initial guess using Haaland equation
        let f_init = Self::haaland_friction_factor(reynolds, relative_roughness)?;
        // Solve for x = 1/sqrt(f)
        let one = <T as NumericElement>::ONE;
        let mut x = one / <T as NumericElement>::sqrt(f_init);

        let ten = <T as FloatElement>::from_f64(10.0);
        let tolerance = <T as FloatElement>::from_f64(1e-6);
        let two = one + one;
        let three_point_seven = <T as FloatElement>::from_f64(3.7);
        let two_point_five_one = <T as FloatElement>::from_f64(2.51);
        let ln_10 = <T as FloatElement>::ln(ten);

        // Newton-Raphson iteration on x = 1/sqrt(f)
        // Equation: g(x) = x + 2 * log10(eps/(3.7*D) + 2.51*x/Re) = 0
        for _ in 0..20 {
            let argument =
                relative_roughness / three_point_seven + two_point_five_one * x / reynolds;
            let g = x + two * (<T as FloatElement>::ln(argument) / ln_10);
            let g_prime = one + (two / ln_10) * (two_point_five_one / reynolds) / argument;

            let x_new = x - g / g_prime;

            if <T as NumericElement>::abs(x_new - x) < tolerance {
                return Ok(one / (x_new * x_new));
            }
            x = x_new;
        }

        Ok(one / (x * x))
    }

    /// Calculate friction factor using Haaland explicit correlation
    fn haaland_friction_factor<T: RealField + FloatElement + Copy>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        let one = <T as NumericElement>::ONE;
        let ten = <T as FloatElement>::from_f64(10.0);
        let six_point_nine = <T as FloatElement>::from_f64(6.9);
        let three_point_seven = <T as FloatElement>::from_f64(3.7);
        let one_point_one_one = <T as FloatElement>::from_f64(1.11);
        let one_point_eight = <T as FloatElement>::from_f64(1.8);

        let term1 =
            <T as FloatElement>::powf(relative_roughness / three_point_seven, one_point_one_one);
        let term2 = six_point_nine / reynolds;
        let inv_sqrt_f = -one_point_eight
            * (<T as FloatElement>::ln(term1 + term2) / <T as FloatElement>::ln(ten));

        Ok(one / (inv_sqrt_f * inv_sqrt_f))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::fluid::ConstantPropertyFluid;

    #[test]
    fn laminar_friction_factor_matches_closed_form() -> Result<()> {
        let friction = FluidDynamicsService::friction_factor(1_000.0_f64, 0.1, None)?;

        assert!((friction - 0.064).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn test_colebrook_white_convergence() {
        // Setup fluid properties (water-like)
        let fluid = ConstantPropertyFluid {
            name: "Water".to_string(),
            density: 1000.0,
            viscosity: 0.001,
            specific_heat: 4182.0,
            thermal_conductivity: 0.6,
            speed_of_sound: 1482.0,
        };

        // Pipe parameters
        let velocity = 5.0; // m/s
        let diameter = 0.1; // m
        let length = 1.0; // m
        let roughness = Some(0.00015); // m, commercial steel

        // Calculate Reynolds number manually to verify inputs
        // Re = rho * v * D / mu = 1000 * 5 * 0.1 / 0.001 = 500,000
        let re = FluidDynamicsService::reynolds_number(&fluid, velocity, diameter);
        assert!((re - 500_000.0).abs() < 1e-1);

        // Calculate pressure drop
        let delta_p =
            FluidDynamicsService::pipe_pressure_drop(&fluid, velocity, length, diameter, roughness)
                .expect("Pressure drop calculation failed");

        // Back-calculate friction factor
        // dp = f * (L/D) * (rho * v^2 / 2)
        // f = dp * D * 2 / (L * rho * v^2)
        let dynamic_pressure = 0.5 * fluid.density * velocity * velocity;
        let f_calc = delta_p * diameter / (length * dynamic_pressure);

        println!("Calculated Friction Factor: {f_calc}");

        // Standard Colebrook-White value for Re=500,000 and eps/D=0.0015 is approx 0.0218
        // Our iterative solution yields ~0.02217, which is within 2% of the chart value.
        let expected_f = 0.0218;
        let error = (f_calc - expected_f).abs() / expected_f;

        // Verify result is within 2% of expected value (chart reading tolerance)
        assert!(
            error < 0.02,
            "Friction factor {f_calc} deviates too much from expected {expected_f}"
        );
    }
}
