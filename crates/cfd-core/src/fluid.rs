//! Fluid properties and models.

use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use num_traits::Float;
use serde::{Deserialize, Serialize};

/// Represents fluid properties
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Fluid<T: RealField> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Kinematic viscosity [m²/s]
    pub kinematic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: Option<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: Option<T>,
    /// Viscosity model
    pub viscosity_model: ViscosityModel<T>,
}

/// Viscosity models for non-Newtonian fluids
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum ViscosityModel<T: RealField> {
    /// Newtonian fluid with constant viscosity
    Newtonian,
    /// Power-law model: μ = K * γ^(n-1)
    PowerLaw {
        /// Consistency index [Pa·s^n]
        consistency_index: T,
        /// Flow behavior index [-]
        flow_index: T,
    },
    /// Carreau model
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
    /// Cross model
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

impl<T: RealField + FromPrimitive + Float> Fluid<T> {
    /// Create a new Newtonian fluid
    pub fn new_newtonian(name: impl Into<String>, density: T, viscosity: T) -> Self {
        let kinematic_viscosity = viscosity.clone() / density.clone();
        Self {
            name: name.into(),
            density: density.clone(),
            viscosity: viscosity.clone(),
            kinematic_viscosity,
            specific_heat: None,
            thermal_conductivity: None,
            viscosity_model: ViscosityModel::Newtonian,
        }
    }

    /// Create water at 20°C
    pub fn water() -> Self {
        Self::new_newtonian(
            "Water (20°C)",
            T::from_f64(998.2).unwrap(),
            T::from_f64(1.002e-3).unwrap(),
        )
    }

    /// Create air at 20°C, 1 atm
    pub fn air() -> Self {
        Self::new_newtonian(
            "Air (20°C, 1 atm)",
            T::from_f64(1.204).unwrap(),
            T::from_f64(1.825e-5).unwrap(),
        )
    }

    /// Get dynamic viscosity for a given shear rate
    pub fn dynamic_viscosity(&self, shear_rate: T) -> T {
        match &self.viscosity_model {
            ViscosityModel::Newtonian => self.viscosity.clone(),
            ViscosityModel::PowerLaw {
                consistency_index,
                flow_index,
            } => {
                *consistency_index * num_traits::Float::powf(shear_rate, *flow_index - T::one())
            }
            ViscosityModel::Carreau {
                mu_zero,
                mu_inf,
                lambda,
                n,
            } => {
                let exponent = (*n - T::one()) / T::from_f64(2.0).unwrap();
                let factor = T::one() + num_traits::Float::powf(*lambda * shear_rate, T::from_f64(2.0).unwrap());
                
                *mu_inf
                    + (*mu_zero - *mu_inf)
                        * num_traits::Float::powf(factor, exponent)
            }
            ViscosityModel::Cross {
                mu_zero,
                mu_inf,
                lambda,
                m,
            } => {
                let factor = T::one() + num_traits::Float::powf(*lambda * shear_rate, *m);
                *mu_inf
                    + (*mu_zero - *mu_inf) / factor
            }
        }
    }

    /// Calculate Reynolds number
    pub fn reynolds_number(&self, velocity: T, length: T) -> T {
        velocity * length / self.kinematic_viscosity.clone()
    }

    /// Calculate Prandtl number  
    pub fn prandtl_number(&self) -> Option<T> {
        match (self.specific_heat.clone(), self.thermal_conductivity.clone()) {
            (Some(cp), Some(k)) => Some(self.viscosity.clone() * cp / k),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_water_properties() {
        let water = Fluid::<f64>::water();
        assert_eq!(water.name, "Water (20°C)");
        assert_relative_eq!(water.density, 998.2, epsilon = 0.1);
        assert_relative_eq!(water.viscosity, 1.002e-3, epsilon = 1e-6);
    }

    #[test]
    fn test_reynolds_number() {
        let water = Fluid::<f64>::water();
        let re = water.reynolds_number(1.0, 0.1);
        assert_relative_eq!(re, 99_820.0, epsilon = 10.0);
    }
}