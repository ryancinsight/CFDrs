//! Inlet and outlet boundary conditions
//!
//! Reference: Patankar (1980), Numerical Heat Transfer and Fluid Flow

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Inlet boundary condition types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum InletCondition<T: RealField + Copy> {
    /// Velocity inlet with prescribed velocity
    Velocity {
        /// Velocity vector [m/s]
        velocity: Vector3<T>,
    },

    /// Pressure inlet with total pressure
    Pressure {
        /// Total pressure [Pa]
        pressure: T,
        /// Optional velocity direction (normalized)
        direction: Option<Vector3<T>>,
    },

    /// Mass flow inlet
    MassFlow {
        /// Mass flow rate [kg/s]
        rate: T,
        /// Optional temperature [K]
        temperature: Option<T>,
    },

    /// Volume flow inlet
    VolumeFlow {
        /// Volume flow rate [m³/s]
        rate: T,
    },
}

/// Outlet boundary condition types
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum OutletCondition<T: RealField + Copy> {
    /// Pressure outlet with static pressure
    Pressure {
        /// Static pressure [Pa]
        pressure: T,
    },

    /// Outflow with zero gradient
    Outflow,
}

impl<T: RealField + Copy> InletCondition<T> {
    /// Create velocity inlet
    pub fn velocity(velocity: Vector3<T>) -> Self {
        Self::Velocity { velocity }
    }

    /// Create pressure inlet
    pub fn pressure(pressure: T, direction: Option<Vector3<T>>) -> Self {
        Self::Pressure {
            pressure,
            direction,
        }
    }

    /// Create mass flow inlet
    pub fn mass_flow(rate: T, temperature: Option<T>) -> Self {
        Self::MassFlow { rate, temperature }
    }

    /// Create volume flow inlet
    pub fn volume_flow(rate: T) -> Self {
        Self::VolumeFlow { rate }
    }
}

impl<T: RealField + Copy> OutletCondition<T> {
    /// Create pressure outlet
    pub fn pressure(pressure: T) -> Self {
        Self::Pressure { pressure }
    }

    /// Create outflow boundary
    #[must_use]
    pub const fn outflow() -> Self {
        Self::Outflow
    }
}
