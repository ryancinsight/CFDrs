//! Fundamental boundary condition types
//!
//! Mathematical classification of boundary conditions following
//! standard PDE theory.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Fundamental physical type of a boundary condition
///
/// Reference: LeVeque (2002). Finite Volume Methods for Hyperbolic Problems
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FundamentalBCType {
    /// Fixed value (first type)
    Dirichlet,
    /// Fixed gradient (second type)
    Neumann,
    /// Mixed (third type)
    Robin,
    /// Other specialized types
    Other,
}

/// Core boundary condition types for CFD
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum BoundaryCondition<T: RealField + Copy> {
    /// Dirichlet: u = g
    Dirichlet {
        /// The fixed value to apply at the boundary
        value: T,
    },

    /// Neumann: ∂u/∂n = g
    Neumann {
        /// The gradient normal to the boundary
        gradient: T,
    },

    /// Robin: αu + β∂u/∂n = γ
    Robin {
        /// Coefficient for the value term
        alpha: T,
        /// Coefficient for the gradient term
        beta: T,
        /// Right-hand side constant
        gamma: T,
    },

    /// Periodic boundary
    Periodic {
        /// Name of the partner boundary for periodicity
        partner: String,
    },

    /// Velocity inlet
    VelocityInlet {
        /// Velocity vector at the inlet
        velocity: Vector3<T>,
    },

    /// Pressure inlet
    PressureInlet {
        /// Pressure at the inlet
        pressure: T,
        /// Optional velocity direction vector
        velocity_direction: Option<Vector3<T>>,
    },

    /// Pressure outlet
    PressureOutlet {
        /// Pressure at the outlet
        pressure: T,
    },

    /// Mass flow inlet
    MassFlowInlet {
        /// Mass flow rate at the inlet
        mass_flow_rate: T,
        /// Optional temperature specification
        temperature: Option<T>,
    },

    /// Volume flow inlet
    VolumeFlowInlet {
        /// Volume flow rate at the inlet
        volume_flow_rate: T,
    },

    /// Wall boundary
    Wall {
        /// Type of wall (no-slip, slip, etc.)
        wall_type: super::WallType<T>,
    },

    /// Symmetry plane
    Symmetry,

    /// Outflow (zero gradient)
    Outflow,
}

impl<T: RealField + Copy> BoundaryCondition<T> {
    /// Create velocity inlet boundary condition
    pub fn velocity_inlet(velocity: Vector3<T>) -> Self {
        Self::VelocityInlet { velocity }
    }

    /// Create pressure outlet boundary condition
    pub fn pressure_outlet(pressure: T) -> Self {
        Self::PressureOutlet { pressure }
    }

    /// Create no-slip wall boundary condition
    #[must_use]
    pub fn wall_no_slip() -> Self {
        Self::Wall {
            wall_type: super::WallType::NoSlip,
        }
    }

    /// Create slip wall boundary condition
    #[must_use]
    pub fn wall_slip() -> Self {
        Self::Wall {
            wall_type: super::WallType::Slip,
        }
    }

    /// Get fundamental type classification
    pub const fn fundamental_type(&self) -> FundamentalBCType {
        match self {
            Self::Dirichlet { .. }
            | Self::VelocityInlet { .. }
            | Self::PressureInlet { .. }
            | Self::Wall { .. } => FundamentalBCType::Dirichlet,

            Self::Neumann { .. } | Self::Outflow | Self::Symmetry => FundamentalBCType::Neumann,

            Self::Robin { .. } => FundamentalBCType::Robin,

            Self::Periodic { .. }
            | Self::MassFlowInlet { .. }
            | Self::VolumeFlowInlet { .. }
            | Self::PressureOutlet { .. } => FundamentalBCType::Other,
        }
    }

    /// Check if Dirichlet-type
    pub const fn is_dirichlet(&self) -> bool {
        matches!(self.fundamental_type(), FundamentalBCType::Dirichlet)
    }

    /// Check if Neumann-type
    pub const fn is_neumann(&self) -> bool {
        matches!(self.fundamental_type(), FundamentalBCType::Neumann)
    }
}
