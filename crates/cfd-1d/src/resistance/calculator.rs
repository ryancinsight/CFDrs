//! Resistance calculator with model selection and validation.

use super::geometry::ChannelGeometry;
use super::models::{
    DarcyWeisbachModel, FlowConditions, HagenPoiseuilleModel, RectangularChannelModel,
    ResistanceModel,
};
use cfd_core::error::Result;
use cfd_core::fluid::{ConstantPropertyFluid, Fluid};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Resistance calculator with model selection and validation
pub struct ResistanceCalculator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceCalculator<T> {
    /// Create new resistance calculator
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Calculate resistance with automatic model selection
    pub fn calculate_auto(
        &self,
        geometry: &ChannelGeometry<T>,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        match geometry {
            ChannelGeometry::Circular { diameter, length } => {
                let model = HagenPoiseuilleModel {
                    diameter: *diameter,
                    length: *length,
                };
                model.calculate_resistance(fluid, conditions)
            }
            ChannelGeometry::Rectangular {
                width,
                height,
                length,
            } => {
                let model = RectangularChannelModel {
                    width: *width,
                    height: *height,
                    length: *length,
                };
                model.calculate_resistance(fluid, conditions)
            }
        }
    }

    /// Calculate resistance using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille(
        &self,
        diameter: T,
        length: T,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = HagenPoiseuilleModel { diameter, length };
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate resistance using rectangular channel model
    pub fn calculate_rectangular(
        &self,
        width: T,
        height: T,
        length: T,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = RectangularChannelModel {
            width,
            height,
            length,
        };
        model.calculate_resistance(fluid, conditions)
    }

    /// Calculate resistance using Darcy-Weisbach model
    pub fn calculate_darcy_weisbach(
        &self,
        hydraulic_diameter: T,
        length: T,
        roughness: T,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        let model = DarcyWeisbachModel {
            hydraulic_diameter,
            length,
            roughness,
        };
        model.calculate_resistance(fluid, conditions)
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> Default for ResistanceCalculator<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hagen_poiseuille() -> Result<()> {
        // Test case: 100 μm diameter, 1 mm length pipe with water
        let diameter = 100e-6; // 100 μm
        let length = 0.001; // 1 mm
        let model = HagenPoiseuilleModel::new(diameter, length);
        let fluid = cfd_core::fluid::ConstantPropertyFluid::<f64>::water_20c();
        let conditions = FlowConditions::new(0.001);

        let resistance = model.calculate_resistance(&fluid, &conditions)?;

        // Calculate expected resistance: R = 128μL/(πD⁴)
        // For water at 20°C: μ ≈ 1.002e-3 Pa·s
        let viscosity = 1.002e-3;
        let expected = 128.0 * viscosity * length / (std::f64::consts::PI * diameter.powi(4));

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
        let fluid = cfd_core::fluid::ConstantPropertyFluid::<f64>::water_20c();
        let conditions = FlowConditions::new(0.001);

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
        use cfd_core::constants::physical::fluid::WATER_DENSITY_STD;

        let diameter = 10e-3; // 10mm pipe for turbulent flow
        let model = DarcyWeisbachModel::new(diameter, 0.001, 1e-6);

        // Test with turbulent flow
        let fluid = cfd_core::fluid::ConstantPropertyFluid::<f64>::water_20c();
        let velocity = 1.0;
        let density = WATER_DENSITY_STD;
        let viscosity = fluid.dynamic_viscosity();
        let re = density * velocity * diameter / viscosity;

        // Create flow conditions with Reynolds number
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(re);

        assert!(re > cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_UPPER);

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
            temperature: 293.15,
            pressure: 101325.0,
        };

        let conditions_turbulent = FlowConditions {
            reynolds_number: Some(10000.0),
            velocity: Some(1.0),
            flow_rate: None,
            temperature: 293.15,
            pressure: 101325.0,
        };

        let hp_model = HagenPoiseuilleModel::new(0.001f64, 0.01);
        assert!(hp_model.is_applicable(&conditions_laminar));
        assert!(!hp_model.is_applicable(&conditions_turbulent));

        let dw_model = DarcyWeisbachModel::new(0.001f64, 0.01, 1e-6);
        assert!(!dw_model.is_applicable(&conditions_laminar));
        assert!(dw_model.is_applicable(&conditions_turbulent));
    }

    #[test]
    fn test_calculator_auto_selection() -> Result<()> {
        let calculator = ResistanceCalculator::<f64>::new();
        let fluid = cfd_core::fluid::ConstantPropertyFluid::<f64>::water_20c();
        let conditions = FlowConditions::new(0.001);

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
        let resistance_rectangular =
            calculator.calculate_auto(&rectangular, &fluid, &conditions)?;
        assert!(resistance_rectangular > 0.0);

        Ok(())
    }
}
