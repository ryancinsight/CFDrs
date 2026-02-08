//! Resistance calculator with model selection and validation.

use super::geometry::ChannelGeometry;
use super::models::{
    DarcyWeisbachModel, FlowConditions, HagenPoiseuilleModel, RectangularChannelModel,
    ResistanceModel, SerpentineCrossSection, SerpentineModel, VenturiModel,
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::{Fluid, FluidTrait};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Resistance calculator with model selection and validation
pub struct ResistanceCalculator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive> ResistanceCalculator<T> {
    /// Create new resistance calculator
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Calculate resistance with automatic model selection
    pub fn calculate_auto<F: FluidTrait<T>>(
        &self,
        geometry: &ChannelGeometry<T>,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let mut local_conditions = conditions.clone();

        // Compute Reynolds number if not provided.
        if local_conditions.reynolds_number.is_none() {
            let state = fluid.properties_at(local_conditions.temperature, local_conditions.pressure)?;
            let density = state.density;
            let viscosity = state.dynamic_viscosity;

            let velocity = if let Some(v) = local_conditions.velocity {
                v
            } else if let Some(q) = local_conditions.flow_rate {
                let area = geometry.cross_sectional_area()?;
                if area <= T::zero() {
                    return Err(Error::InvalidConfiguration(
                        "Channel area must be positive to compute Reynolds number".to_string(),
                    ));
                }
                q / area
            } else {
                return Err(Error::InvalidConfiguration(
                    "Automatic resistance selection requires either velocity or flow_rate (or an explicit Reynolds number)".to_string(),
                ));
            };

            let dh = geometry.hydraulic_diameter()?;
            let re = density * velocity * dh / viscosity;

            local_conditions.velocity = Some(velocity);
            local_conditions.reynolds_number = Some(re);
        }

        match geometry {
            ChannelGeometry::Circular { diameter, length } => {
                let hp = HagenPoiseuilleModel {
                    diameter: *diameter,
                    length: *length,
                };
                if hp.is_applicable(&local_conditions) {
                    hp.validate_invariants(fluid, &local_conditions)?;
                    return hp.calculate_resistance(fluid, &local_conditions);
                }

                // Default to a smooth pipe if roughness is not specified by the geometry.
                let dw = DarcyWeisbachModel::circular(*diameter, *length, T::zero());
                if dw.is_applicable(&local_conditions) {
                    dw.validate_invariants(fluid, &local_conditions)?;
                    return dw.calculate_resistance(fluid, &local_conditions);
                }

                Err(Error::InvalidConfiguration(
                    "No applicable resistance model for circular channel at the given Reynolds number".to_string(),
                ))
            }
            ChannelGeometry::Rectangular {
                width,
                height,
                length,
            } => {
                let rect = RectangularChannelModel {
                    width: *width,
                    height: *height,
                    length: *length,
                };
                if rect.is_applicable(&local_conditions) {
                    rect.validate_invariants(fluid, &local_conditions)?;
                    rect.calculate_resistance(fluid, &local_conditions)
                } else {
                    Err(Error::InvalidConfiguration(
                        "Rectangular-channel resistance currently supports laminar flow only; provide a laminar Reynolds number or use a different model".to_string(),
                    ))
                }
            }
            _ => Err(Error::InvalidConfiguration(
                "No resistance model available for this geometry type; use Circular or Rectangular".to_string(),
            )),
        }
    }

    /// Calculate linear (R) and quadratic (k) coefficients with automatic model selection
    pub fn calculate_coefficients_auto<F: FluidTrait<T>>(
        &self,
        geometry: &ChannelGeometry<T>,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let mut local_conditions = conditions.clone();

        // Compute Reynolds number if not provided.
        if local_conditions.reynolds_number.is_none() {
            let state = fluid.properties_at(local_conditions.temperature, local_conditions.pressure)?;
            let density = state.density;
            let viscosity = state.dynamic_viscosity;

            let velocity = if let Some(v) = local_conditions.velocity {
                v
            } else if let Some(q) = local_conditions.flow_rate {
                let area = geometry.cross_sectional_area()?;
                if area <= T::zero() {
                    return Err(Error::InvalidConfiguration(
                        "Channel area must be positive to compute Reynolds number".to_string(),
                    ));
                }
                q / area
            } else {
                return Err(Error::InvalidConfiguration(
                    "Automatic resistance selection requires either velocity or flow_rate (or an explicit Reynolds number)".to_string(),
                ));
            };

            let dh = geometry.hydraulic_diameter()?;
            let re = density * velocity * dh / viscosity;

            local_conditions.velocity = Some(velocity);
            local_conditions.reynolds_number = Some(re);
        }

        match geometry {
            ChannelGeometry::Circular { diameter, length } => {
                let hp = HagenPoiseuilleModel {
                    diameter: *diameter,
                    length: *length,
                };
                if hp.is_applicable(&local_conditions) {
                    hp.validate_invariants(fluid, &local_conditions)?;
                    return hp.calculate_coefficients(fluid, &local_conditions);
                }

                // Default to a smooth pipe if roughness is not specified by the geometry.
                let dw = DarcyWeisbachModel::circular(*diameter, *length, T::zero());
                if dw.is_applicable(&local_conditions) {
                    dw.validate_invariants(fluid, &local_conditions)?;
                    return dw.calculate_coefficients(fluid, &local_conditions);
                }

                Err(Error::InvalidConfiguration(
                    "No applicable resistance model for circular channel at the given Reynolds number".to_string(),
                ))
            }
            ChannelGeometry::Rectangular {
                width,
                height,
                length,
            } => {
                let rect = RectangularChannelModel {
                    width: *width,
                    height: *height,
                    length: *length,
                };
                if rect.is_applicable(&local_conditions) {
                    rect.validate_invariants(fluid, &local_conditions)?;
                    rect.calculate_coefficients(fluid, &local_conditions)
                } else {
                    Err(Error::InvalidConfiguration(
                        "Rectangular-channel resistance currently supports laminar flow only; provide a laminar Reynolds number or use a different model".to_string(),
                    ))
                }
            }
            _ => Err(Error::InvalidConfiguration(
                "No resistance model available for this geometry type; use Circular or Rectangular".to_string(),
            )),
        }
    }

    /// Calculate linear (R) and quadratic (k) coefficients using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille_coefficients<F: FluidTrait<T>>(
        &self,
        diameter: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let model = HagenPoiseuilleModel { diameter, length };
        model.validate_invariants(fluid, conditions)?;
        model.calculate_coefficients(fluid, conditions)
    }

    /// Calculate resistance using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille<F: FluidTrait<T>>(
        &self,
        diameter: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = HagenPoiseuilleModel { diameter, length };
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate linear (R) and quadratic (k) coefficients using rectangular channel model
    pub fn calculate_rectangular_coefficients<F: FluidTrait<T>>(
        &self,
        width: T,
        height: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        let model = RectangularChannelModel {
            width,
            height,
            length,
        };
        model.validate_invariants(fluid, conditions)?;
        model.calculate_coefficients(fluid, conditions)
    }

    /// Calculate resistance using rectangular channel model
    pub fn calculate_rectangular<F: FluidTrait<T>>(
        &self,
        width: T,
        height: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = RectangularChannelModel {
            width,
            height,
            length,
        };
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate linear (R) and quadratic (k) coefficients using Darcy-Weisbach model
    pub fn calculate_darcy_weisbach_coefficients<F: FluidTrait<T>>(
        &self,
        hydraulic_diameter: T,
        length: T,
        roughness: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        // Treat the provided hydraulic diameter as the circular diameter.
        let model = DarcyWeisbachModel::circular(hydraulic_diameter, length, roughness);
        model.validate_invariants(fluid, conditions)?;
        model.calculate_coefficients(fluid, conditions)
    }

    /// Calculate resistance using Darcy-Weisbach model
    pub fn calculate_darcy_weisbach<F: FluidTrait<T>>(
        &self,
        hydraulic_diameter: T,
        length: T,
        roughness: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        // Treat the provided hydraulic diameter as the circular diameter.
        let model = DarcyWeisbachModel::circular(hydraulic_diameter, length, roughness);
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate resistance using serpentine channel model with circular cross-section.
    ///
    /// Accounts for Dean flow in curved sections and minor losses at bends.
    /// The model uses White/Ito curvature correlations and Idelchik bend
    /// loss coefficients.
    ///
    /// # Arguments
    /// - `diameter`: Channel diameter [m]
    /// - `straight_length`: Total length of all straight segments [m]
    /// - `num_segments`: Number of straight segments (bends = segments - 1)
    /// - `bend_radius`: Radius of curvature of bends [m]
    /// - `fluid`: Fluid properties (supports non-Newtonian via `FluidTrait`)
    /// - `conditions`: Flow conditions (velocity, Re, temperature, pressure)
    pub fn calculate_serpentine_circular<F: FluidTrait<T>>(
        &self,
        diameter: T,
        straight_length: T,
        num_segments: usize,
        bend_radius: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = super::factory::ResistanceModelFactory::serpentine_circular(
            diameter,
            straight_length,
            num_segments,
            bend_radius,
        );
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate resistance using serpentine channel model with rectangular cross-section.
    ///
    /// Applies Shah-London f·Re correction for rectangular ducts in addition
    /// to Dean flow curvature enhancement.
    pub fn calculate_serpentine_rectangular<F: FluidTrait<T>>(
        &self,
        width: T,
        height: T,
        straight_length: T,
        num_segments: usize,
        bend_radius: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = super::factory::ResistanceModelFactory::serpentine_rectangular(
            width,
            height,
            straight_length,
            num_segments,
            bend_radius,
        );
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate (R, k) coefficients using serpentine channel model.
    ///
    /// Returns the linear (friction-dominated) and quadratic (bend-loss-dominated)
    /// resistance coefficients: ΔP = R·Q + k·Q|Q|
    pub fn calculate_serpentine_coefficients<F: FluidTrait<T>>(
        &self,
        model: &SerpentineModel<T>,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        model.validate_invariants(fluid, conditions)?;
        model.calculate_coefficients(fluid, conditions)
    }

    /// Calculate resistance using Venturi tube model.
    ///
    /// Accounts for contraction loss (via discharge coefficient), throat
    /// friction (Darcy), and expansion recovery loss (Borda-Carnot).
    ///
    /// # Arguments
    /// - `inlet_diameter`: Upstream pipe diameter [m]
    /// - `throat_diameter`: Throat diameter [m]
    /// - `throat_length`: Length of the throat section [m]
    /// - `total_length`: Total device length [m]
    /// - `fluid`: Fluid properties
    /// - `conditions`: Flow conditions
    pub fn calculate_venturi<F: FluidTrait<T>>(
        &self,
        inlet_diameter: T,
        throat_diameter: T,
        throat_length: T,
        total_length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = VenturiModel::symmetric(
            inlet_diameter,
            throat_diameter,
            throat_length,
            total_length,
        );
        model.validate_invariants(fluid, conditions)?;
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate (R, k) coefficients using Venturi tube model.
    ///
    /// Returns the linear (friction-dominated) and quadratic (inertial)
    /// resistance coefficients: ΔP = R·Q + k·Q|Q|
    pub fn calculate_venturi_coefficients<F: FluidTrait<T>>(
        &self,
        model: &VenturiModel<T>,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        model.validate_invariants(fluid, conditions)?;
        model.calculate_coefficients(fluid, conditions)
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for ResistanceCalculator<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::ConstantFluid;

    #[test]
    fn test_hagen_poiseuille() -> Result<()> {
        // Test case: 100 μm diameter, 1 mm length pipe with water
        let diameter = 100e-6; // 100 μm
        let length = 0.001; // 1 mm
        let model = HagenPoiseuilleModel::new(diameter, length);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;
        let conditions = FlowConditions::new(0.001);

        let resistance = model.calculate_resistance(&fluid, &conditions)?;

        // Calculate expected resistance: R = 128μL/(πD⁴)
        // For water at 20°C: μ ≈ 1.002e-3 Pa·s
        let viscosity = 1.002e-3;
        let expected = 128.0_f64 * viscosity * length / (std::f64::consts::PI * diameter.powi(4));

        // Verify within 1% of theoretical value
        assert_relative_eq!(resistance, expected, epsilon = expected * 0.01);

        // Verify pressure drop for 1 μL/s flow rate
        let flow_rate = 1e-9; // 1 μL/s = 1e-9 m³/s
        let pressure_drop = resistance * flow_rate;
        let expected_pressure = expected * flow_rate;
        assert_relative_eq!(
            pressure_drop,
            expected_pressure,
            epsilon = expected_pressure * 0.01
        );

        Ok(())
    }

    #[test]
    fn test_rectangular_channel() -> Result<()> {
        let model = RectangularChannelModel::new(100e-6, 50e-6, 0.001);
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;
        let mut conditions = FlowConditions::new(0.001);
        // Compute Reynolds number based on hydraulic diameter and set it
        let dh = 2.0 * model.width * model.height / (model.width + model.height);
        let density = fluid.density;
        let viscosity = fluid.dynamic_viscosity();
        let velocity = conditions.velocity.unwrap();
        let re = density * velocity * dh / viscosity;
        conditions.reynolds_number = Some(re);

        let resistance = model.calculate_resistance(&fluid, &conditions)?;

        // Verify resistance is positive
        assert!(resistance > 0.0);

        // Verify hydraulic diameter calculation
        // Calculate hydraulic diameter for rectangular channel
        // D_h = 2 * w * h / (w + h)
        let hydraulic_diameter = 2.0 * model.width * model.height / (model.width + model.height);
        assert_relative_eq!(
            hydraulic_diameter,
            2.0 * 100e-6 * 50e-6 / (100e-6 + 50e-6),
            epsilon = 1e-10
        );

        Ok(())
    }

    #[test]
    fn test_darcy_weisbach() -> Result<()> {
        use cfd_core::physics::constants::physics::fluid::WATER_DENSITY;

        let diameter = 10e-3; // 10mm pipe for turbulent flow
        let model = DarcyWeisbachModel::circular(diameter, 0.001, 1e-6);

        // Test with turbulent flow
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;
        let velocity = 1.0;
        let density = WATER_DENSITY;
        let viscosity = fluid.viscosity;
        let re = density * velocity * diameter / viscosity;

        // Create flow conditions with Reynolds number
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(re);

        assert!(re > cfd_core::physics::constants::physics::dimensionless::reynolds::PIPE_TURBULENT_MIN);

        let resistance = model.calculate_resistance(&fluid, &conditions)?;
        assert!(resistance > 0.0);

        Ok(())
    }

    #[test]
    fn test_model_applicability() {
        let conditions_laminar = FlowConditions {
            reynolds_number: Some(1000.0),
            velocity: Some(0.001),
            flow_rate: None,
            shear_rate: None,
            temperature: 293.15,
            pressure: 101_325.0,
        };

        let conditions_turbulent = FlowConditions {
            reynolds_number: Some(10_000.0),
            velocity: Some(1.0),
            flow_rate: None,
            shear_rate: None,
            temperature: 293.15,
            pressure: 101_325.0,
        };

        let hp_model = HagenPoiseuilleModel::new(0.001f64, 0.01);
        assert!(hp_model.is_applicable(&conditions_laminar));
        assert!(!hp_model.is_applicable(&conditions_turbulent));

        let dw_model = DarcyWeisbachModel::circular(0.001f64, 0.01, 1e-6);
        assert!(!dw_model.is_applicable(&conditions_laminar));
        assert!(dw_model.is_applicable(&conditions_turbulent));
    }

    #[test]
    fn test_calculator_auto_selection() -> Result<()> {
        let calculator = ResistanceCalculator::<f64>::new();
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;
        let mut conditions = FlowConditions::new(0.001);

        // Test circular geometry
        let circular = ChannelGeometry::Circular {
            diameter: 100e-6,
            length: 0.001,
        };
        let resistance_circular = calculator.calculate_auto(&circular, &fluid, &conditions)?;
        assert!(resistance_circular > 0.0);

        // Test rectangular geometry
        let rectangular = ChannelGeometry::Rectangular {
            width: 100e-6,
            height: 50e-6,
            length: 0.001,
        };
        // Provide Reynolds number for rectangular geometry
        let dh = 2.0 * 100e-6 * 50e-6 / (100e-6 + 50e-6);
        let density = fluid.density;
        let viscosity = fluid.dynamic_viscosity();
        let velocity = conditions.velocity.unwrap();
        let re = density * velocity * dh / viscosity;
        conditions.reynolds_number = Some(re);

        let resistance_rectangular =
            calculator.calculate_auto(&rectangular, &fluid, &conditions)?;
        assert!(resistance_rectangular > 0.0);

        Ok(())
    }
}
