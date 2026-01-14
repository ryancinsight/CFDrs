//! Fluid dynamics service following Domain Service pattern
//!
//! This service coordinates complex fluid mechanics operations and
//! business logic that doesn't naturally fit into a single domain entity.

use crate::error::Result;
use crate::physics::fluid::ConstantPropertyFluid;
use crate::physics::fluid_dynamics::flow_regimes::{FlowClassifier, FlowRegime};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Service for fluid dynamics calculations
pub struct FluidDynamicsService;

impl FluidDynamicsService {
    /// Calculate Reynolds number for a given flow configuration
    pub fn reynolds_number<T: RealField + Copy + Float>(
        fluid: &ConstantPropertyFluid<T>,
        velocity: T,
        characteristic_length: T,
    ) -> T {
        let kinematic_viscosity = fluid.viscosity / fluid.density;
        velocity * characteristic_length / kinematic_viscosity
    }

    /// Calculate Prandtl number for constant property fluid
    pub fn prandtl_number<T: RealField + Copy + Float>(fluid: &ConstantPropertyFluid<T>) -> T {
        let cp = fluid.specific_heat;
        let k = fluid.thermal_conductivity;
        let mu = fluid.viscosity;
        mu * cp / k
    }

    /// Determine flow regime based on Reynolds number
    pub fn flow_regime<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive>(
        reynolds: T,
    ) -> FlowRegime {
        FlowClassifier::classify_by_reynolds(reynolds)
    }

    /// Calculate pressure drop for pipe flow
    pub fn pipe_pressure_drop<T: RealField + FromPrimitive + Copy + Float>(
        fluid: &ConstantPropertyFluid<T>,
        velocity: T,
        length: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let reynolds = Self::reynolds_number(fluid, velocity, diameter);
        let friction_factor = Self::friction_factor(reynolds, diameter, roughness)?;

        let two = T::one() + T::one();

        Ok(friction_factor * length * fluid.density * velocity * velocity / (two * diameter))
    }

    /// Calculate friction factor using appropriate correlation
    fn friction_factor<T: RealField + FromPrimitive + Copy + Float>(
        reynolds: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let re_2300 = T::from_f64(
            crate::physics::constants::physics::dimensionless::reynolds::PIPE_LAMINAR_MAX,
        )
        .unwrap_or_else(|| T::one());
        let sixty_four = T::from_f64(64.0).unwrap_or_else(|| T::one());

        if reynolds < re_2300 {
            // Laminar flow: f = 64/Re
            Ok(sixty_four / reynolds)
        } else if let Some(eps) = roughness {
            // Colebrook-White equation
            let relative_roughness = eps / diameter;
            Self::colebrook_white_friction_factor(reynolds, relative_roughness)
        } else {
            // Smooth pipe: Blasius (low Re) or Haaland (general explicit)
            let blasius_max =
                T::from_f64(crate::physics::constants::physics::hydraulics::BLASIUS_MAX_RE)
                    .unwrap_or_else(|| T::one());
            if reynolds < blasius_max {
                let coeff = T::from_f64(
                    crate::physics::constants::physics::hydraulics::BLASIUS_COEFFICIENT,
                )
                .unwrap_or_else(|| T::one());
                let exp =
                    T::from_f64(crate::physics::constants::physics::hydraulics::BLASIUS_EXPONENT)
                        .unwrap_or_else(|| T::one());
                Ok(coeff / Float::powf(reynolds, exp))
            } else {
                Self::haaland_friction_factor(reynolds, T::zero())
            }
        }
    }

    /// Calculate friction factor using Colebrook-White equation (iterative)
    fn colebrook_white_friction_factor<T: RealField + FromPrimitive + Copy + Float>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        // Initial guess using Haaland equation
        let mut f = Self::haaland_friction_factor(reynolds, relative_roughness)?;
        let ten = T::from_f64(10.0).unwrap_or_else(|| T::one());
        let tolerance = T::from_f64(1e-6).unwrap_or_else(|| T::one());
        let two = T::one() + T::one();
        let point_eight_six = T::from_f64(0.86).unwrap_or_else(|| T::one());
        let three_point_seven = T::from_f64(3.7).unwrap_or_else(|| T::one());
        let two_point_five_one = T::from_f64(2.51).unwrap_or_else(|| T::one());

        // Simple fixed-point iteration
        for _ in 0..20 {
            let f_old = f;
            let inv_sqrt_f =
                -two * Float::log(
                    relative_roughness / three_point_seven
                        + two_point_five_one / (reynolds * Float::sqrt(f)),
                    ten,
                ) * point_eight_six;
            f = T::one() / (inv_sqrt_f * inv_sqrt_f);

            if Float::abs(f - f_old) < tolerance {
                return Ok(f);
            }
        }

        Ok(f)
    }

    /// Calculate friction factor using Haaland explicit correlation
    fn haaland_friction_factor<T: RealField + FromPrimitive + Copy + Float>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        let ten = T::from_f64(10.0).unwrap_or_else(|| T::one());
        let six_point_nine = T::from_f64(6.9).unwrap_or_else(|| T::one());
        let three_point_seven = T::from_f64(3.7).unwrap_or_else(|| T::one());
        let one_point_one_one = T::from_f64(1.11).unwrap_or_else(|| T::one());
        let one_point_eight = T::from_f64(1.8).unwrap_or_else(|| T::one());

        let term1 = Float::powf(relative_roughness / three_point_seven, one_point_one_one);
        let term2 = six_point_nine / reynolds;
        let inv_sqrt_f = -one_point_eight * Float::log(term1 + term2, ten);

        Ok(T::one() / (inv_sqrt_f * inv_sqrt_f))
    }
}
