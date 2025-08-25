//! Value objects for CFD domain modeling.
//!
//! This module provides immutable value objects that represent important
//! domain concepts with built-in validation and behavior.

use crate::error::{Error, Result};
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Flow geometry type for Reynolds number interpretation
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum FlowGeometry {
    /// Pipe or channel flow
    Pipe,
    /// Flow over flat plate
    FlatPlate,
    /// Flow around sphere
    Sphere,
    /// Flow around cylinder
    Cylinder,
    /// General internal flow
    Internal,
    /// General external flow
    External,
}

/// Reynolds number with geometry-aware flow regime detection
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ReynoldsNumber<T: RealField + Copy> {
    value: T,
    geometry: FlowGeometry,
}

impl<T: RealField + Copy + FromPrimitive> ReynoldsNumber<T> {
    /// Create a new Reynolds number with specified geometry
    pub fn new(value: T, geometry: FlowGeometry) -> crate::error::Result<Self> {
        if value < T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Reynolds number must be non-negative".into(),
            ));
        }
        Ok(Self { value, geometry })
    }

    /// Create Reynolds number for pipe flow
    pub fn pipe(value: T) -> crate::error::Result<Self> {
        Self::new(value, FlowGeometry::Pipe)
    }

    /// Get the raw value
    pub fn value(&self) -> T {
        self.value
    }

    /// Check if flow is laminar based on geometry
    pub fn is_laminar(&self) -> bool {
        let threshold = match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => {
                T::from_f64(crate::constants::physics::dimensionless::reynolds::PIPE_CRITICAL_LOWER)
            }
            FlowGeometry::FlatPlate | FlowGeometry::External => {
                T::from_f64(crate::constants::physics::dimensionless::reynolds::PLATE_CRITICAL)
            }
            FlowGeometry::Sphere => {
                T::from_f64(crate::constants::physics::dimensionless::reynolds::SPHERE_CRITICAL)
            }
            FlowGeometry::Cylinder => {
                T::from_f64(crate::constants::physics::dimensionless::reynolds::CYLINDER_CRITICAL)
            }
        }
        .unwrap_or_else(|| T::zero());

        self.value < threshold
    }

    /// Check if flow is turbulent based on geometry
    pub fn is_turbulent(&self) -> bool {
        let threshold = match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => T::from_f64(4000.0),
            FlowGeometry::FlatPlate | FlowGeometry::External => T::from_f64(1e6),
            FlowGeometry::Sphere | FlowGeometry::Cylinder => T::from_f64(3e5),
        }
        .unwrap_or_else(|| T::zero());

        self.value >= threshold
    }

    /// Check if flow is in transition region
    pub fn is_transitional(&self) -> bool {
        !self.is_laminar() && !self.is_turbulent()
    }

    /// Get transition probability (0 = fully laminar, 1 = fully turbulent)
    /// Uses smooth transition function
    pub fn transition_probability(&self) -> T {
        let (lower, upper) = match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => {
                (T::from_f64(2300.0), T::from_f64(4000.0))
            }
            _ => (T::from_f64(5e5), T::from_f64(1e6)),
        };

        let lower = lower.unwrap_or_else(|| T::zero());
        let upper = upper.unwrap_or_else(|| T::one());

        if self.value <= lower {
            T::zero()
        } else if self.value >= upper {
            T::one()
        } else {
            // Smooth transition using tanh function
            let mid = (lower + upper) / T::from_f64(2.0).unwrap_or_else(|| T::one());
            let width = upper - lower;
            let normalized =
                (self.value - mid) / (width * T::from_f64(0.25).unwrap_or_else(|| T::one()));
            (T::one() + normalized.tanh()) / T::from_f64(2.0).unwrap_or_else(|| T::one())
        }
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for ReynoldsNumber<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Re = {} ({:?})", self.value, self.geometry)
    }
}

/// Pressure value object with unit awareness
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Pressure<T: RealField + Copy> {
    value: T,
    unit: PressureUnit,
}

impl<T: RealField + FromPrimitive + Copy> Pressure<T> {
    /// Create pressure in Pascals
    pub fn pascals(value: T) -> Self {
        Self {
            value,
            unit: PressureUnit::Pascal,
        }
    }

    /// Create pressure in atmospheres
    pub fn atmospheres(value: T) -> Self {
        Self {
            value,
            unit: PressureUnit::Atmosphere,
        }
    }

    /// Create pressure in bar
    pub fn bar(value: T) -> Self {
        Self {
            value,
            unit: PressureUnit::Bar,
        }
    }

    /// Convert to Pascals
    pub fn to_pascals(&self) -> T {
        match self.unit {
            PressureUnit::Pascal => self.value,
            PressureUnit::Atmosphere => {
                self.value * T::from_f64(101_325.0).unwrap_or_else(|| T::one())
            }
            PressureUnit::Bar => self.value * T::from_f64(100_000.0).unwrap_or_else(|| T::one()),
        }
    }

    /// Get value in current unit
    pub fn value(&self) -> T {
        self.value
    }

    /// Get unit
    pub fn unit(&self) -> PressureUnit {
        self.unit
    }
}

/// Pressure units
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PressureUnit {
    /// Pascal (N/m²)
    Pascal,
    /// Atmosphere (atm)
    Atmosphere,
    /// Bar
    Bar,
}

impl fmt::Display for PressureUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pascal => write!(f, "Pa"),
            Self::Atmosphere => write!(f, "atm"),
            Self::Bar => write!(f, "bar"),
        }
    }
}

/// Velocity value object with magnitude and direction
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Velocity<T: RealField + Copy> {
    vector: Vector3<T>,
}

impl<T: RealField + Copy> Velocity<T> {
    /// Create velocity from components
    pub fn new(x: T, y: T, z: T) -> Self {
        Self {
            vector: Vector3::new(x, y, z),
        }
    }

    /// Create velocity from vector
    pub fn from_vector(vector: Vector3<T>) -> Self {
        Self { vector }
    }

    /// Get magnitude
    pub fn magnitude(&self) -> T {
        self.vector.norm()
    }

    /// Get unit vector (direction)
    pub fn direction(&self) -> Vector3<T> {
        let mag = self.magnitude();
        if mag > T::zero() {
            self.vector / mag
        } else {
            Vector3::zeros()
        }
    }

    /// Get velocity vector
    pub fn vector(&self) -> &Vector3<T> {
        &self.vector
    }

    /// Get x component
    pub fn x(&self) -> T {
        self.vector.x
    }

    /// Get y component
    pub fn y(&self) -> T {
        self.vector.y
    }

    /// Get z component
    pub fn z(&self) -> T {
        self.vector.z
    }
}

/// Temperature value object with unit conversion
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Temperature<T: RealField + Copy> {
    value: T,
    unit: TemperatureUnit,
}

impl<T: RealField + FromPrimitive + Copy> Temperature<T> {
    /// Create temperature in Kelvin
    pub fn kelvin(value: T) -> Result<Self> {
        if value < T::zero() {
            return Err(Error::InvalidConfiguration(
                "Temperature cannot be below absolute zero".to_string(),
            ));
        }
        Ok(Self {
            value,
            unit: TemperatureUnit::Kelvin,
        })
    }

    /// Create temperature in Celsius
    pub fn celsius(value: T) -> Result<Self> {
        let kelvin_value = value
            + T::from_f64(crate::constants::thermo::CELSIUS_TO_KELVIN).unwrap_or_else(|| T::one());
        Self::kelvin(kelvin_value)
    }

    /// Create temperature in Fahrenheit
    pub fn fahrenheit(value: T) -> Result<Self> {
        use crate::constants::physical::temperature::{
            FAHRENHEIT_TO_CELSIUS_FACTOR, FAHRENHEIT_ZERO_OFFSET,
        };
        let celsius_value = (value
            - T::from_f64(FAHRENHEIT_ZERO_OFFSET).unwrap_or_else(|| T::one()))
            * T::from_f64(FAHRENHEIT_TO_CELSIUS_FACTOR).unwrap_or_else(|| T::one());
        Self::celsius(celsius_value)
    }

    /// Convert to Kelvin
    pub fn to_kelvin(&self) -> T {
        match self.unit {
            TemperatureUnit::Kelvin => self.value,
            TemperatureUnit::Celsius => {
                self.value
                    + T::from_f64(crate::constants::thermo::CELSIUS_TO_KELVIN)
                        .unwrap_or_else(|| T::zero())
            }
            TemperatureUnit::Fahrenheit => {
                use crate::constants::physical::temperature::{
                    FAHRENHEIT_TO_CELSIUS_FACTOR, FAHRENHEIT_ZERO_OFFSET,
                };
                let celsius = (self.value
                    - T::from_f64(FAHRENHEIT_ZERO_OFFSET).unwrap_or_else(|| T::zero()))
                    * T::from_f64(FAHRENHEIT_TO_CELSIUS_FACTOR).unwrap_or_else(|| T::zero());
                celsius
                    + T::from_f64(crate::constants::thermo::CELSIUS_TO_KELVIN)
                        .unwrap_or_else(|| T::zero())
            }
        }
    }

    /// Get value in current unit
    pub fn value(&self) -> T {
        self.value
    }

    /// Get unit
    pub fn unit(&self) -> TemperatureUnit {
        self.unit
    }
}

/// Temperature units
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TemperatureUnit {
    /// Kelvin
    Kelvin,
    /// Celsius
    Celsius,
    /// Fahrenheit
    Fahrenheit,
}

impl fmt::Display for TemperatureUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Kelvin => write!(f, "K"),
            Self::Celsius => write!(f, "°C"),
            Self::Fahrenheit => write!(f, "°F"),
        }
    }
}

/// Dimensionless number value object
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DimensionlessNumber<T: RealField + Copy> {
    value: T,
    name: &'static str,
}

impl<T: RealField + Copy> DimensionlessNumber<T> {
    /// Create a new dimensionless number
    pub fn new(value: T, name: &'static str) -> Self {
        Self { value, name }
    }

    /// Create Prandtl number
    pub fn prandtl(value: T) -> Self {
        Self::new(value, "Prandtl")
    }

    /// Create Nusselt number
    pub fn nusselt(value: T) -> Self {
        Self::new(value, "Nusselt")
    }

    /// Create Grashof number
    pub fn grashof(value: T) -> Self {
        Self::new(value, "Grashof")
    }

    /// Get value
    pub fn value(&self) -> T {
        self.value
    }

    /// Get name
    pub fn name(&self) -> &'static str {
        self.name
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for DimensionlessNumber<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} = {}", self.name, self.value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reynolds_number() {
        let re = ReynoldsNumber::new(1500.0, FlowGeometry::Pipe).expect("Valid Reynolds number");
        assert!(re.is_laminar());
        assert!(!re.is_turbulent());
        assert!(!re.is_transitional());

        let re_turb =
            ReynoldsNumber::new(5000.0, FlowGeometry::Pipe).expect("Valid Reynolds number");
        assert!(!re_turb.is_laminar());
        assert!(re_turb.is_turbulent());
        assert!(!re_turb.is_transitional());

        let re_trans =
            ReynoldsNumber::new(3000.0, FlowGeometry::Pipe).expect("Valid Reynolds number");
        assert!(!re_trans.is_laminar());
        assert!(!re_trans.is_turbulent());
        assert!(re_trans.is_transitional());
    }

    #[test]
    fn test_pressure_conversion() {
        let p_atm = Pressure::atmospheres(1.0f64);
        let p_pa = p_atm.to_pascals();
        assert!((p_pa - 101_325.0).abs() < 1e-6);

        let p_bar = Pressure::bar(1.0f64);
        let p_pa_bar = p_bar.to_pascals();
        assert!((p_pa_bar - 100_000.0).abs() < 1e-6);
    }

    #[test]
    fn test_temperature_conversion() {
        let t_c = Temperature::celsius(0.0f64).expect("CRITICAL: Add proper error handling");
        let t_k = t_c.to_kelvin();
        assert!((t_k - 273.15).abs() < 1e-6);

        let t_f = Temperature::fahrenheit(32.0f64).expect("CRITICAL: Add proper error handling");
        let t_k_f = t_f.to_kelvin();
        assert!((t_k_f - 273.15).abs() < 1e-6);
    }

    #[test]
    fn test_velocity() {
        let vel = Velocity::new(3.0f64, 4.0f64, 0.0f64);
        assert!((vel.magnitude() - 5.0).abs() < 1e-10);

        let dir = vel.direction();
        assert!((dir.x - 0.6).abs() < 1e-10);
        assert!((dir.y - 0.8).abs() < 1e-10);
    }
}
