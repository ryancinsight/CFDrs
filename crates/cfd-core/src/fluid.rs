//! Fluid properties and models.
//!
//! This module provides representations for various fluid types including
//! Newtonian and non-Newtonian fluids with different viscosity models.

use crate::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use num_traits::Float;
use serde::{Deserialize, Serialize};

/// Represents fluid properties
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Fluid<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: Option<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: Option<T>,
    /// Viscosity model (contains viscosity data)
    pub viscosity_model: ViscosityModel<T>,
}

/// Viscosity models for different fluid types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ViscosityModel<T: RealField + Copy> {
    /// Constant viscosity model
    Constant {
        /// Dynamic viscosity [Pa·s]
        viscosity: T,
    },
    /// Power-law model: μ = K * γ^(n-1)
    PowerLaw {
        /// Consistency index [Pa·s^n]
        consistency_index: T,
        /// Flow behavior index [-]
        flow_index: T,
    },
    /// Carreau model for shear-thinning fluids
    Carreau {
        /// Zero-shear viscosity [Pa·s]
        mu_zero: T,
        /// Infinite-shear viscosity [Pa·s]
        mu_inf: T,
        /// Relaxation time [s]
        lambda: T,
        /// Power-law index [-]
        n: T,
    },
    /// Cross model for polymer solutions
    Cross {
        /// Zero-shear viscosity [Pa·s]
        mu_zero: T,
        /// Infinite-shear viscosity [Pa·s]
        mu_inf: T,
        /// Natural time [s]
        lambda: T,
        /// Rate constant [-]
        m: T,
    },
}

impl<T: RealField + Copy + FromPrimitive + Float> Fluid<T> {
    /// Create a fluid with constant viscosity
    pub fn constant_viscosity(name: impl Into<String>, density: T, viscosity: T) -> Self {
        Self {
            name: name.into(),
            density,
            specific_heat: None,
            thermal_conductivity: None,
            viscosity_model: ViscosityModel::Constant { viscosity },
        }
    }

    /// Create a new power-law fluid
    pub fn current_power_law(
        name: impl Into<String>,
        density: T,
        consistency_index: T,
        flow_index: T,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            specific_heat: None,
            thermal_conductivity: None,
            viscosity_model: ViscosityModel::PowerLaw {
                consistency_index,
                flow_index,
            },
        }
    }

    /// Create a new Carreau fluid
    pub fn current_carreau(
        name: impl Into<String>,
        density: T,
        mu_zero: T,
        mu_inf: T,
        lambda: T,
        n: T,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            specific_heat: None,
            thermal_conductivity: None,
            viscosity_model: ViscosityModel::Carreau {
                mu_zero,
                mu_inf,
                lambda,
                n,
            },
        }
    }

    /// Create a new Cross fluid
    pub fn current_cross(
        name: impl Into<String>,
        density: T,
        mu_zero: T,
        mu_inf: T,
        lambda: T,
        m: T,
    ) -> Self {
        Self {
            name: name.into(),
            density,
            specific_heat: None,
            thermal_conductivity: None,
            viscosity_model: ViscosityModel::Cross {
                mu_zero,
                mu_inf,
                lambda,
                m,
            },
        }
    }

    /// Predefined water properties at 20°C, 1 atm
    pub fn water() -> Result<Self> {
        Ok(Self {
            name: "Water (20°C)".to_string(),
            density: T::from_f64(998.2)
                .ok_or_else(|| Error::InvalidInput("Cannot convert density value".into()))?,
            specific_heat: Some(
                T::from_f64(4182.0)
                    .ok_or_else(|| Error::InvalidInput("Cannot convert specific heat".into()))?,
            ),
            thermal_conductivity: Some(T::from_f64(0.598).ok_or_else(|| {
                Error::InvalidInput("Cannot convert thermal conductivity".into())
            })?),
            viscosity_model: ViscosityModel::Constant {
                viscosity: T::from_f64(1.002e-3)
                    .ok_or_else(|| Error::InvalidInput("Cannot convert viscosity value".into()))?,
            },
        })
    }

    /// Predefined air properties at 20°C, 1 atm
    pub fn air() -> Result<Self> {
        Ok(Self {
            name: "Air (20°C, 1 atm)".to_string(),
            density: T::from_f64(1.204)
                .ok_or_else(|| Error::InvalidInput("Cannot convert density value".into()))?,
            specific_heat: Some(
                T::from_f64(718.0)
                    .ok_or_else(|| Error::InvalidInput("Cannot convert specific heat cv".into()))?,
            ),
            thermal_conductivity: Some(T::from_f64(0.0257).ok_or_else(|| {
                Error::InvalidInput("Cannot convert thermal conductivity".into())
            })?),
            viscosity_model: ViscosityModel::Constant {
                viscosity: T::from_f64(1.825e-5)
                    .ok_or_else(|| Error::InvalidInput("Cannot convert viscosity value".into()))?,
            },
        })
    }

    /// Calculate kinematic viscosity [m²/s]
    ///
    /// For variable viscosity fluids, this uses the characteristic viscosity.
    pub fn kinematic_viscosity(&self) -> T {
        self.characteristic_viscosity() / self.density
    }

    /// Get the characteristic viscosity for the fluid
    ///
    /// For constant viscosity fluids, this returns the constant value.
    /// For variable viscosity fluids, this returns the zero-shear viscosity
    /// or consistency index as appropriate for Reynolds number calculations.
    pub fn characteristic_viscosity(&self) -> T {
        match &self.viscosity_model {
            ViscosityModel::Constant { viscosity } => *viscosity,
            ViscosityModel::PowerLaw {
                consistency_index, ..
            } => *consistency_index,
            ViscosityModel::Carreau { mu_zero, .. } => *mu_zero,
            ViscosityModel::Cross { mu_zero, .. } => *mu_zero,
        }
    }

    /// Get dynamic viscosity for a given shear rate
    ///
    /// # Errors
    ///
    /// Returns an error if numeric conversions fail during calculation.
    pub fn dynamic_viscosity(&self, shear_rate: T) -> Result<T> {
        match &self.viscosity_model {
            ViscosityModel::Constant { viscosity } => Ok(*viscosity),
            ViscosityModel::PowerLaw {
                consistency_index,
                flow_index,
            } => Ok(
                *consistency_index * num_traits::Float::powf(shear_rate, *flow_index - T::one())
            ),
            ViscosityModel::Carreau {
                mu_zero,
                mu_inf,
                lambda,
                n,
            } => {
                let two = T::from_f64(2.0).ok_or_else(|| {
                    Error::InvalidConfiguration("Cannot convert 2.0 to target type".into())
                })?;
                let exponent = (*n - T::one()) / two;
                let factor = T::one() + num_traits::Float::powf(*lambda * shear_rate, two);

                Ok(*mu_inf + (*mu_zero - *mu_inf) * num_traits::Float::powf(factor, exponent))
            }
            ViscosityModel::Cross {
                mu_zero,
                mu_inf,
                lambda,
                m,
            } => {
                let factor = T::one() + num_traits::Float::powf(*lambda * shear_rate, *m);
                Ok(*mu_inf + (*mu_zero - *mu_inf) / factor)
            }
        }
    }

    /// Calculate Reynolds number
    ///
    /// Uses the characteristic viscosity for variable viscosity fluids.
    pub fn reynolds_number(&self, velocity: T, length: T) -> T {
        velocity * length / self.kinematic_viscosity()
    }

    /// Calculate Prandtl number
    ///
    /// Returns None if thermal properties are not defined.
    pub fn prandtl_number(&self) -> Option<T> {
        match (self.specific_heat, self.thermal_conductivity) {
            (Some(cp), Some(k)) => Some(self.characteristic_viscosity() * cp / k),
            _ => None,
        }
    }

    /// Set specific heat capacity
    pub fn with_specific_heat(mut self, cp: T) -> Self {
        self.specific_heat = Some(cp);
        self
    }

    /// Set thermal conductivity
    pub fn with_thermal_conductivity(mut self, k: T) -> Self {
        self.thermal_conductivity = Some(k);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_water_properties() -> Result<()> {
        let water = Fluid::<f64>::water()?;
        assert_eq!(water.name, "Water (20°C)");
        assert_relative_eq!(water.density, 998.2, epsilon = 0.1);

        // Check that viscosity is stored in the model
        let ViscosityModel::Constant { viscosity } = water.viscosity_model else {
            return Err(Error::InvalidInput(
                "Water should have constant viscosity model".into(),
            ));
        };
        assert_relative_eq!(viscosity, 1.002e-3, epsilon = 1e-6);

        // Check kinematic viscosity calculation
        assert_relative_eq!(
            water.kinematic_viscosity(),
            1.002e-3 / 998.2,
            epsilon = 1e-9
        );
        Ok(())
    }

    #[test]
    fn test_air_properties() -> Result<()> {
        let air = Fluid::<f64>::air()?;
        assert_eq!(air.name, "Air (20°C, 1 atm)");
        assert_relative_eq!(air.density, 1.204, epsilon = 0.001);

        let ViscosityModel::Constant { viscosity } = air.viscosity_model else {
            return Err(Error::InvalidInput(
                "Air should have constant viscosity model".into(),
            ));
        };
        assert_relative_eq!(viscosity, 1.825e-5, epsilon = 1e-8);
        Ok(())
    }

    #[test]
    fn test_reynolds_number() -> Result<()> {
        let water = Fluid::<f64>::water()?;
        let re = water.reynolds_number(1.0, 0.1);
        // Re = ρVL/μ = 998.2 * 1.0 * 0.1 / 1.002e-3 ≈ 99,620
        assert_relative_eq!(re, 99_620.76, epsilon = 1.0);
        Ok(())
    }

    #[test]
    fn test_power_law_viscosity() -> Result<()> {
        // Example: Ketchup-like fluid
        let fluid = Fluid::<f64>::current_power_law(
            "Power-law fluid",
            1000.0, // density
            0.5,    // consistency index K [Pa·s^n]
            0.8,    // flow index n [-]
        );

        let shear_rate = 100.0; // [s^-1]
                                // μ = K * γ^(n-1) = 0.5 * 100^(0.8-1) = 0.5 * 100^(-0.2)
        let expected_viscosity = 0.5 * 100.0_f64.powf(-0.2);
        let calculated = fluid.dynamic_viscosity(shear_rate)?;

        assert_relative_eq!(calculated, expected_viscosity, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn test_carreau_model() -> Result<()> {
        // Example: Polymer solution
        let fluid = Fluid::<f64>::current_carreau(
            "Carreau fluid",
            1050.0, // density [kg/m³]
            1.0,    // zero-shear viscosity [Pa·s]
            0.001,  // infinite-shear viscosity [Pa·s]
            0.1,    // relaxation time [s]
            0.5,    // power-law index [-]
        );

        // Test at zero shear (should approach mu_zero)
        let visc_low = fluid.dynamic_viscosity(0.001)?;
        assert_relative_eq!(visc_low, 1.0, epsilon = 0.01);

        // Test at high shear (should approach mu_inf)
        // At shear_rate=1000: viscosity ≈ mu_inf + (mu_zero - mu_inf) * (1 + (lambda*gamma)^2)^((n-1)/2)
        // ≈ 0.001 + 0.999 * (1 + 10000)^(-0.25) ≈ 0.001 + 0.999 * 0.1 ≈ 0.1
        let visc_high = fluid.dynamic_viscosity(1000.0)?;
        assert_relative_eq!(visc_high, 0.1, epsilon = 0.01);

        // Test characteristic viscosity (should be mu_zero)
        assert_relative_eq!(fluid.characteristic_viscosity(), 1.0, epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn test_cross_model() -> Result<()> {
        // Example: Polymer melt
        let fluid = Fluid::<f64>::current_cross(
            "Cross fluid",
            950.0, // density [kg/m³]
            10.0,  // zero-shear viscosity [Pa·s]
            0.01,  // infinite-shear viscosity [Pa·s]
            1.0,   // natural time [s]
            0.75,  // rate constant [-]
        );

        // Test at low shear rate
        let visc_low = fluid.dynamic_viscosity(0.001)?;
        assert!(visc_low > 9.0); // Should be close to mu_zero

        // Test at high shear rate
        let visc_high = fluid.dynamic_viscosity(1000.0)?;
        assert!(visc_high < 1.0); // Should approach mu_inf

        // Test intermediate shear rate
        let visc_mid = fluid.dynamic_viscosity(1.0)?;
        assert!(visc_mid > 0.01 && visc_mid < 10.0);
        Ok(())
    }

    #[test]
    fn test_kinematic_viscosity_calculation() -> Result<()> {
        // Test that kinematic viscosity is calculated correctly
        let water = Fluid::<f64>::water()?;
        let kinematic = water.kinematic_viscosity();
        let expected = 1.002e-3 / 998.2;
        assert_relative_eq!(kinematic, expected, epsilon = 1e-9);

        // Test for variable viscosity fluid
        let power_law = Fluid::<f64>::current_power_law(
            "Test fluid",
            1200.0,
            0.8, // This becomes the characteristic viscosity
            0.9,
        );
        assert_relative_eq!(
            power_law.kinematic_viscosity(),
            0.8 / 1200.0,
            epsilon = 1e-9
        );
        Ok(())
    }

    #[test]
    fn test_thermal_properties() -> Result<()> {
        // Water now includes thermal properties by default
        let water = Fluid::<f64>::water()?;

        // Check that thermal properties are present
        assert!(water.specific_heat.is_some());
        assert!(water.thermal_conductivity.is_some());

        // Calculate Prandtl number
        let pr = water
            .prandtl_number()
            .expect("Water should have thermal properties");
        // Pr = μ·cp/k = 1.002e-3 * 4182 / 0.598 ≈ 7.0
        assert_relative_eq!(pr, 7.0, epsilon = 0.1);

        // Test fluid without thermal properties
        let fluid = Fluid::<f64>::constant_viscosity("Test", 1000.0, 0.001);
        assert!(fluid.prandtl_number().is_none());

        Ok(())
    }

    #[test]
    fn test_prandtl_number() -> Result<()> {
        let water = Fluid::<f64>::water()?;
        let pr = water
            .prandtl_number()
            .ok_or_else(|| Error::InvalidInput("Cannot calculate Prandtl number".into()))?;

        // Water at 20°C has Pr ≈ 7
        assert!(pr > 6.0 && pr < 8.0);

        Ok(())
    }

    #[test]
    fn test_error_handling() {
        // This would fail if we had a type that can't represent these values
        // For f64, it should always succeed
        assert!(Fluid::<f64>::water().is_ok());
        assert!(Fluid::<f64>::air().is_ok());
    }
}
