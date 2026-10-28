//! Flow-related value objects

use super::dimensionless::DimensionlessNumber;
use crate::error::Result;
use crate::physics::constants::physics::dimensionless::reynolds::{
    PIPE_LAMINAR_MAX, PIPE_TURBULENT_MIN, PLATE_CRITICAL,
};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

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
pub struct ReynoldsNumber<T: FloatElement + Copy> {
    value: T,
    geometry: FlowGeometry,
}

impl<T: FloatElement + Copy> ReynoldsNumber<T> {
    /// Create a new Reynolds number with specified geometry
    ///
    /// # Errors
    /// Returns an error if the value is negative
    pub fn new(value: T, geometry: FlowGeometry) -> Result<Self> {
        if value < <T as NumericElement>::ZERO {
            return Err(crate::error::Error::InvalidConfiguration(
                "Reynolds number must be non-negative".into(),
            ));
        }
        Ok(Self { value, geometry })
    }

    /// Create Reynolds number for pipe flow
    ///
    /// # Errors
    /// Returns an error if the value is negative
    pub fn pipe(value: T) -> Result<Self> {
        Self::new(value, FlowGeometry::Pipe)
    }

    /// Create Reynolds number for flat plate flow
    ///
    /// # Errors
    ///
    /// Returns an error if the value is not finite or is negative.
    pub fn flat_plate(value: T) -> Result<Self> {
        Self::new(value, FlowGeometry::FlatPlate)
    }

    /// Get the Reynolds number value
    pub fn value(&self) -> T {
        self.value
    }

    /// Get the flow geometry
    pub fn geometry(&self) -> FlowGeometry {
        self.geometry
    }

    /// Determine if flow is laminar
    pub fn is_laminar(&self) -> bool {
        match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => {
                self.value < <T as FloatElement>::from_f64(PIPE_LAMINAR_MAX)
            }
            FlowGeometry::FlatPlate | FlowGeometry::External => {
                self.value < <T as FloatElement>::from_f64(PLATE_CRITICAL)
            }
            _ => self.value < <T as FloatElement>::from_f64(PIPE_LAMINAR_MAX),
        }
    }

    /// Determine if flow is turbulent
    pub fn is_turbulent(&self) -> bool {
        match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => {
                self.value > <T as FloatElement>::from_f64(PIPE_TURBULENT_MIN)
            }
            FlowGeometry::FlatPlate | FlowGeometry::External => {
                self.value > <T as FloatElement>::from_f64(PLATE_CRITICAL)
            }
            _ => self.value > <T as FloatElement>::from_f64(PIPE_TURBULENT_MIN),
        }
    }

    /// Determine if flow is in transition
    pub fn is_transitional(&self) -> bool {
        !self.is_laminar() && !self.is_turbulent()
    }

    /// Calculate from flow parameters
    ///
    /// The Reynolds value uses the single [`DimensionlessNumber::reynolds`]
    /// formula; this constructor only adds the geometry-aware regime
    /// interpretation on top.
    ///
    /// # Errors
    ///
    /// Returns an error if any parameter is not finite or if density, velocity, length, or viscosity are negative or zero.
    pub fn from_flow_parameters(
        density: T,
        velocity: T,
        length: T,
        viscosity: T,
        geometry: FlowGeometry,
    ) -> Result<Self> {
        let re = DimensionlessNumber::reynolds(density, velocity, length, viscosity)?;
        Self::new(re.value(), geometry)
    }
}

#[cfg(test)]
mod tests {
    use super::super::dimensionless::DimensionlessType;
    use super::*;

    #[test]
    fn from_flow_parameters_matches_taxonomy_reynolds() {
        let density = 1_000.0_f64;
        let velocity = 2.0_f64;
        let length = 0.01_f64;
        let viscosity = 0.001_f64;
        let expected = density * velocity * length / viscosity;

        let specialized = ReynoldsNumber::from_flow_parameters(
            density,
            velocity,
            length,
            viscosity,
            FlowGeometry::Pipe,
        )
        .expect("positive viscosity produces a Reynolds number");
        let taxonomy = DimensionlessNumber::reynolds(density, velocity, length, viscosity)
            .expect("positive viscosity produces a Reynolds number");

        assert_eq!(specialized.value().to_bits(), expected.to_bits());
        assert_eq!(taxonomy.value().to_bits(), expected.to_bits());
        assert_eq!(taxonomy.number_type(), DimensionlessType::Reynolds);
        assert!(specialized.is_turbulent());
        assert!(!specialized.is_laminar());
        assert!(!specialized.is_transitional());
    }

    #[test]
    fn pipe_regime_thresholds_follow_central_constants() {
        let laminar = ReynoldsNumber::pipe(1_000.0_f64).expect("non-negative value is valid");
        let transitional = ReynoldsNumber::pipe(3_000.0_f64).expect("non-negative value is valid");
        assert!(laminar.is_laminar());
        assert!(transitional.is_transitional());
    }

    #[test]
    fn non_positive_viscosity_is_rejected_once() {
        assert!(
            ReynoldsNumber::from_flow_parameters(1_000.0, 2.0, 0.01, 0.0, FlowGeometry::Pipe)
                .is_err()
        );
        assert!(DimensionlessNumber::reynolds(1_000.0, 2.0, 0.01, 0.0).is_err());
    }
}
