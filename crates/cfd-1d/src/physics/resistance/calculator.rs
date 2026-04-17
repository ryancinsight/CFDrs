//! Resistance calculator with model selection and validation.
//!
//! ## Dispatch Pattern
//!
//! `ResistanceCalculator` is a thin facade delegating all physics to the
//! `dispatch` and `coefficients` sub-modules. The dispatch logic selects the
//! appropriate resistance model based on channel geometry:
//!
//! | Geometry variant | Delegated model |
//! |------------------|-----------------|
//! | `Circular` | `HagenPoiseuilleModel` |
//! | `Rectangular` | `RectangularChannelModel` (Shah-London) |
//! | `Serpentine` | `SerpentineModel` (Dean-enhanced friction) |
//! | `Venturi` | `VenturiModel` (contraction + throat + expansion) |
//! | `DarcyWeisbach` | `DarcyWeisbachModel` (general turbulent) |
//! | `Membrane` | `MembranePoreModel` (Hagen-Poiseuille per pore) |
//! | `JunctionLoss` | `JunctionLossModel` (K-factor minor losses) |
//!
//! Each model implements the `ResistanceModel` trait, returning a (R, k) tuple
//! where the total pressure drop is ΔP = R·Q + k·Q² (linear + quadratic terms).
//! The quadratic term is non-zero only for turbulent or geometry-dependent models
//! (Darcy-Weisbach, Venturi expansion losses).

/// Resistance and quadratic coefficient calculation dispatch.
pub mod coefficients;
pub mod dispatch;

use super::geometry::ChannelGeometry;
use super::models::{FlowConditions, SerpentineModel, VenturiModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

pub(crate) fn populate_shear_aware_conditions<T, F>(
    geometry: &ChannelGeometry<T>,
    fluid: &F,
    local_conditions: &mut FlowConditions<T>,
) -> Result<()>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T>,
{
    if local_conditions.reynolds_number.is_some() {
        return Ok(());
    }

    let density = fluid
        .properties_at(local_conditions.temperature, local_conditions.pressure)?
        .density;

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

    let velocity_abs = if velocity >= T::zero() {
        velocity
    } else {
        -velocity
    };
    let dh = geometry.hydraulic_diameter()?;
    if dh <= T::zero() {
        return Err(Error::InvalidConfiguration(
            "Hydraulic diameter must be positive to compute Reynolds number".to_string(),
        ));
    }

    let shear_rate = if let Some(sr) = local_conditions.shear_rate {
        sr
    } else {
        T::from_f64(8.0).expect("Mathematical constant conversion compromised") * velocity_abs / dh
    };
    let apparent_viscosity = fluid.viscosity_at_shear(
        shear_rate,
        local_conditions.temperature,
        local_conditions.pressure,
    )?;
    if apparent_viscosity <= T::zero() {
        return Err(Error::InvalidConfiguration(
            "Viscosity must be positive to compute Reynolds number".to_string(),
        ));
    }

    local_conditions.velocity = Some(velocity);
    local_conditions.shear_rate = Some(shear_rate);
    local_conditions.reynolds_number = Some(density * velocity_abs * dh / apparent_viscosity);
    Ok(())
}

pub(crate) fn rectangular_auto_selection_error<T: RealField + Copy>(
    conditions: &FlowConditions<T>,
) -> Error {
    match conditions.reynolds_number {
        Some(reynolds) => Error::InvalidConfiguration(format!(
            "Automatic rectangular resistance selection is only validated for the laminar rectangular-channel model (Re < 2300). Got Re = {reynolds}. Select RectangularChannelModel explicitly for validated laminar use, or select Darcy-Weisbach explicitly if you intend a hydraulic-diameter surrogate."
        )),
        None => Error::InvalidConfiguration(
            "Automatic rectangular resistance selection requires an explicit Reynolds number or enough flow information to compute one".to_string(),
        ),
    }
}

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
        dispatch::calculate_auto(geometry, fluid, conditions)
    }

    /// Calculate linear (R) and quadratic (k) coefficients with automatic model selection
    pub fn calculate_coefficients_auto<F: FluidTrait<T>>(
        &self,
        geometry: &ChannelGeometry<T>,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        coefficients::calculate_coefficients_auto(geometry, fluid, conditions)
    }

    /// Calculate linear (R) and quadratic (k) coefficients using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille_coefficients<F: FluidTrait<T>>(
        &self,
        diameter: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        coefficients::calculate_hagen_poiseuille_coefficients(diameter, length, fluid, conditions)
    }

    /// Calculate resistance using Hagen-Poiseuille model
    pub fn calculate_hagen_poiseuille<F: FluidTrait<T>>(
        &self,
        diameter: T,
        length: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        dispatch::calculate_hagen_poiseuille(diameter, length, fluid, conditions)
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
        coefficients::calculate_rectangular_coefficients(width, height, length, fluid, conditions)
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
        dispatch::calculate_rectangular(width, height, length, fluid, conditions)
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
        coefficients::calculate_darcy_weisbach_coefficients(
            hydraulic_diameter,
            length,
            roughness,
            fluid,
            conditions,
        )
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
        dispatch::calculate_darcy_weisbach(hydraulic_diameter, length, roughness, fluid, conditions)
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
        dispatch::calculate_serpentine_circular(
            diameter,
            straight_length,
            num_segments,
            bend_radius,
            fluid,
            conditions,
        )
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
        dispatch::calculate_serpentine_rectangular(
            width,
            height,
            straight_length,
            num_segments,
            bend_radius,
            fluid,
            conditions,
        )
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
        coefficients::calculate_serpentine_coefficients(model, fluid, conditions)
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
        dispatch::calculate_venturi(
            inlet_diameter,
            throat_diameter,
            throat_length,
            total_length,
            fluid,
            conditions,
        )
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
        coefficients::calculate_venturi_coefficients(model, fluid, conditions)
    }

    /// Calculate hydraulic resistance of a porous membrane.
    ///
    /// Uses the parallel-pore model from `MembranePoreModel`.
    pub fn calculate_membrane_porous<F: FluidTrait<T>>(
        &self,
        thickness: T,
        width: T,
        height: T,
        pore_radius: T,
        porosity: T,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        dispatch::calculate_membrane_porous(
            thickness,
            width,
            height,
            pore_radius,
            porosity,
            fluid,
            conditions,
        )
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
    use crate::physics::resistance::models::{
        DarcyWeisbachModel, HagenPoiseuilleModel, RectangularChannelModel, ResistanceModel,
    };
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
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

        assert!(
            re > cfd_core::physics::constants::physics::dimensionless::reynolds::PIPE_TURBULENT_MIN
        );

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

    #[test]
    fn test_calculator_auto_uses_shear_aware_reynolds_for_non_newtonian_circular() -> Result<()> {
        let calculator = ResistanceCalculator::<f64>::new();
        let fluid = CarreauYasudaBlood::<f64>::normal_blood();
        let circular = ChannelGeometry::Circular {
            diameter: 0.04,
            length: 2.5,
        };
        let fluid_state = fluid.properties_at(293.15, 101_325.0)?;
        let target_reference_re = 2400.0;
        let flow_rate =
            target_reference_re * std::f64::consts::PI * 0.04 * fluid_state.dynamic_viscosity
                / (4.0 * fluid_state.density);
        let area = circular.cross_sectional_area()?;
        let velocity = flow_rate / area;
        let shear_rate = 8.0 * velocity.abs() / 0.04;
        let apparent_viscosity = fluid.viscosity_at_shear(shear_rate, 293.15, 101_325.0)?;
        let shear_aware_re = fluid_state.density * velocity.abs() * 0.04 / apparent_viscosity;

        assert!(target_reference_re > 2300.0);
        assert!(
            shear_aware_re < 2300.0,
            "expected laminar shear-aware Reynolds, got {shear_aware_re}"
        );

        let mut conditions = FlowConditions::new(0.0);
        conditions.velocity = None;
        conditions.flow_rate = Some(flow_rate);

        let resistance = calculator.calculate_auto(&circular, &fluid, &conditions)?;
        let (coeff_resistance, coeff_k) =
            calculator.calculate_coefficients_auto(&circular, &fluid, &conditions)?;

        let mut expected_conditions = conditions.clone();
        expected_conditions.shear_rate = Some(shear_rate);
        expected_conditions.reynolds_number = Some(shear_aware_re);

        let expected =
            calculator.calculate_hagen_poiseuille(0.04, 2.5, &fluid, &expected_conditions)?;
        let (expected_r, expected_k) = calculator.calculate_hagen_poiseuille_coefficients(
            0.04,
            2.5,
            &fluid,
            &expected_conditions,
        )?;

        assert_relative_eq!(resistance, expected, max_relative = 1e-12);
        assert_relative_eq!(coeff_resistance, expected_r, max_relative = 1e-12);
        assert_relative_eq!(coeff_k, expected_k, epsilon = 1e-15);
        assert_relative_eq!(coeff_k, 0.0, epsilon = 1e-15);

        Ok(())
    }

    #[test]
    fn test_calculator_auto_rejects_non_laminar_rectangular_auto_selection() -> Result<()> {
        let calculator = ResistanceCalculator::<f64>::new();
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>()?;
        let rectangular = ChannelGeometry::Rectangular {
            width: 1.0e-3,
            height: 5.0e-4,
            length: 2.0e-2,
        };
        let conditions = FlowConditions::new(5.0);

        let auto_error = calculator
            .calculate_auto(&rectangular, &fluid, &conditions)
            .expect_err("rectangular auto-selection must not silently fall back above Re=2300");
        assert!(
            auto_error
                .to_string()
                .contains("select Darcy-Weisbach explicitly"),
            "unexpected rectangular auto-selection error: {auto_error}"
        );

        let coeff_error = calculator
            .calculate_coefficients_auto(&rectangular, &fluid, &conditions)
            .expect_err(
                "rectangular coefficient auto-selection must not silently fall back above Re=2300",
            );
        assert!(
            coeff_error
                .to_string()
                .contains("select Darcy-Weisbach explicitly"),
            "unexpected rectangular coefficient auto-selection error: {coeff_error}"
        );

        let dh = 2.0 * 1.0e-3 * 5.0e-4 / (1.0e-3 + 5.0e-4);
        let explicit = calculator.calculate_darcy_weisbach(dh, 2.0e-2, 0.0, &fluid, &conditions)?;
        assert!(explicit > 0.0);

        Ok(())
    }
}
