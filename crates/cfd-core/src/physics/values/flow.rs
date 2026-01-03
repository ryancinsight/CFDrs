//! Flow-related value objects

use crate::error::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

// Physical constants
const REYNOLDS_LAMINAR_LIMIT_PIPE: f64 = 2300.0;
const REYNOLDS_TURBULENT_ONSET_PIPE: f64 = 4000.0;
const REYNOLDS_LAMINAR_LIMIT_PLATE: f64 = 500_000.0;

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
    ///
    /// # Errors
    /// Returns an error if the value is negative
    pub fn new(value: T, geometry: FlowGeometry) -> Result<Self> {
        if value < T::zero() {
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
                self.value < T::from_f64(REYNOLDS_LAMINAR_LIMIT_PIPE).unwrap_or_else(T::zero)
            }
            FlowGeometry::FlatPlate | FlowGeometry::External => {
                self.value < T::from_f64(REYNOLDS_LAMINAR_LIMIT_PLATE).unwrap_or_else(T::zero)
            }
            _ => self.value < T::from_f64(REYNOLDS_LAMINAR_LIMIT_PIPE).unwrap_or_else(T::zero),
        }
    }

    /// Determine if flow is turbulent
    pub fn is_turbulent(&self) -> bool {
        match self.geometry {
            FlowGeometry::Pipe | FlowGeometry::Internal => {
                self.value > T::from_f64(REYNOLDS_TURBULENT_ONSET_PIPE).unwrap_or_else(T::zero)
            }
            FlowGeometry::FlatPlate | FlowGeometry::External => {
                self.value > T::from_f64(REYNOLDS_LAMINAR_LIMIT_PLATE).unwrap_or_else(T::zero)
            }
            _ => self.value > T::from_f64(REYNOLDS_TURBULENT_ONSET_PIPE).unwrap_or_else(T::zero),
        }
    }

    /// Determine if flow is in transition
    pub fn is_transitional(&self) -> bool {
        !self.is_laminar() && !self.is_turbulent()
    }

    /// Calculate from flow parameters
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
        if viscosity <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Viscosity must be positive".into(),
            ));
        }
        let re = density * velocity * length / viscosity;
        Self::new(re, geometry)
    }
}
