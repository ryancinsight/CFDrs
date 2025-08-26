//! Boundary condition types and implementations.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};
/// Boundary condition types for CFD simulations
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundaryCondition<T: RealField + Copy> {
    /// Dirichlet boundary condition (fixed value)
    Dirichlet {
        /// Fixed value at the boundary
        value: T,
    },
    /// Neumann boundary condition (fixed gradient)
    Neumann {
        /// Fixed gradient at the boundary
        gradient: T,
    },
    /// Robin boundary condition (mixed type)
    Robin {
        /// Coefficient for the value term
        alpha: T,
        /// Coefficient for the gradient term
        beta: T,
        /// Right-hand side value
        gamma: T,
    },
    /// Periodic boundary condition
    Periodic {
        /// Name of the paired boundary (can also store indices as string)
        partner: String,
    },
    /// Inlet boundary with specified velocity
    VelocityInlet {
        /// Velocity vector at inlet
        velocity: Vector3<T>,
    },
    /// Pressure inlet
    PressureInlet {
        /// Total pressure at inlet
        pressure: T,
    },
    /// Pressure outlet
    PressureOutlet {
        /// Static pressure at outlet
        pressure: T,
    },
    /// Mass flow inlet
    MassFlowInlet {
        /// Mass flow rate [kg/s]
        mass_flow: T,
    },
    /// Volume flow rate inlet (for 1D networks)
    VolumeFlowInlet {
        /// Volume flow rate [m³/s]
        flow_rate: T,
    },
    /// Wall boundary
    Wall {
        /// Wall type
        wall_type: WallType<T>,
    },
    /// Symmetry boundary
    Symmetry,
    /// Outflow boundary (zero gradient)
    Outflow,
}
/// Wall boundary types
pub enum WallType<T: RealField + Copy> {
    /// No-slip wall (zero velocity)
    NoSlip,
    /// Slip wall (zero normal velocity)
    Slip,
    /// Moving wall with specified velocity
    Moving {
        /// Wall velocity
        velocity: Vector3<T>,
    },
    /// Rotating wall
    Rotating {
        /// Angular velocity [rad/s]
        omega: Vector3<T>,
        /// Center of rotation
        center: Vector3<T>,
    },
}

impl<T: RealField + Copy> BoundaryCondition<T> {
    /// Create a pressure inlet boundary condition
    pub const fn pressure_inlet(pressure: T) -> Self {
        Self::PressureInlet { pressure }
    }
    /// Create a pressure outlet boundary condition
    pub const fn pressure_outlet(pressure: T) -> Self {
        Self::PressureOutlet { pressure }
    }
    
    /// Create a velocity inlet boundary condition
    pub const fn velocity_inlet(velocity: Vector3<T>) -> Self {
        Self::VelocityInlet { velocity }
    }
    
    /// Create a no-slip wall boundary condition
    #[must_use]
    pub const fn wall_no_slip() -> Self {
        Self::Wall {
            wall_type: WallType::NoSlip,
        }
    }
    
    /// Create a slip wall boundary condition
    pub const fn wall_slip() -> Self {
        Self::Wall {
            wall_type: WallType::Slip,
        }
    }
    
    /// Create a moving wall boundary condition
    pub const fn wall_moving(velocity: Vector3<T>) -> Self {
        Self::Wall {
            wall_type: WallType::Moving { velocity },
        }
    }
    
    /// Create a volume flow rate inlet boundary condition (for 1D networks)
    pub const fn flow_rate_inlet(flow_rate: T) -> Self {
        Self::VolumeFlowInlet { flow_rate }
    }
    
    /// Create a 1D velocity inlet (scalar velocity)
    pub fn velocity_inlet_1d(velocity: T) -> Self {
        Self::VelocityInlet {
            velocity: Vector3::new(velocity, T::zero(), T::zero()),
        }
    }
    
    /// Create a zero gradient outflow boundary condition
    pub const fn zero_gradient() -> Self {
        Self::Outflow
    }
    
    /// Check if this is a Dirichlet-type boundary condition
    pub const fn is_dirichlet(&self) -> bool {
        matches!(
            self,
            Self::Dirichlet { .. }
                | Self::VelocityInlet { .. }
                | Self::PressureInlet { .. }
                | Self::Wall { .. }
        )
    }
    
    /// Check if this is a Neumann-type boundary condition
    pub const fn is_neumann(&self) -> bool {
        matches!(self, Self::Neumann { .. } | Self::Outflow | Self::Symmetry)
    }
    
    /// Check if this is a wall boundary condition
    pub const fn is_wall(&self) -> bool {
        matches!(self, Self::Wall { .. })
    }
    
    /// Extract pressure value if this is a pressure boundary condition
    pub fn pressure_value(&self) -> Option<T>
    where
        T: Copy,
    {
        match self {
            Self::PressureInlet { pressure } | Self::PressureOutlet { pressure } => Some(*pressure),
            _ => None,
        }
    }
    
    /// Extract flow rate value if this is a flow rate boundary condition
    pub fn flow_rate_value(&self) -> Option<T> {
        match self {
            Self::VolumeFlowInlet { flow_rate } => Some(*flow_rate),
            _ => None,
        }
    }
    /// Extract 1D velocity value if this is a velocity boundary condition
    pub fn velocity_1d_value(&self) -> Option<T> {
        match self {
            Self::VelocityInlet { velocity } => Some(velocity.x),
            _ => None,
        }
    }
    /// Check if this is a zero gradient boundary condition
    pub const fn is_zero_gradient(&self) -> bool {
        matches!(self, Self::Outflow | Self::Symmetry)
    }
}

/// Boundary condition set for a complete problem
#[derive(Debug, Clone)]
pub struct BoundaryConditionSet<T: RealField + Copy> {
    /// Map of boundary names to conditions
    pub conditions: indexmap::IndexMap<String, BoundaryCondition<T>>,
}

impl<T: RealField + Copy> BoundaryConditionSet<T> {
    /// Create a new empty boundary condition set
    pub fn new() -> Self {
        Self {
            conditions: indexmap::IndexMap::new(),
        }
    }
    
    /// Add a boundary condition
    pub fn add(&mut self, name: impl Into<String>, condition: BoundaryCondition<T>) -> &mut Self {
        self.conditions.insert(name.into(), condition);
        self
    }
    
    /// Get a boundary condition by name
    pub fn get(&self, name: &str) -> Option<&BoundaryCondition<T>> {
        self.conditions.get(name)
    }
    
    /// Validate that periodic boundary conditions are properly paired
    ///
    /// # Errors
    /// Returns an error if:
    /// - A periodic boundary references a non-existent partner
    /// - A periodic boundary's partner is not also periodic
    /// - Periodic boundaries do not reference each other mutually
    pub fn validate_periodic(&self) -> crate::error::Result<()> {
        use std::collections::HashSet;
        let mut validated_partners = HashSet::new();
        for (name, bc) in &self.conditions {
            if let BoundaryCondition::Periodic { partner } = bc {
                // Skip if we already validated this boundary as part of checking its partner
                if validated_partners.contains(name.as_str()) {
                    continue;
                }
                // Check partner exists
                let partner_bc = self.conditions.get(partner).ok_or_else(|| {
                    crate::error::Error::InvalidConfiguration(format!(
                        "Periodic boundary '{name}' references non-existent partner '{partner}'"
                    ))
                })?;
                // Check partner is also periodic and points back to this boundary
                if let BoundaryCondition::Periodic {
                    partner: partner_partner,
                } = partner_bc
                {
                    if partner_partner != name {
                        return Err(crate::error::Error::InvalidConfiguration(format!(
                            "Periodic boundaries '{name}' and '{partner}' do not reference each other mutually"
                        )));
                    }
                } else {
                    return Err(crate::error::Error::InvalidConfiguration(format!(
                        "Periodic boundary '{name}' has a partner '{partner}' which is not periodic"
                    )));
                }
                // Mark the partner as validated to avoid redundant checking
                validated_partners.insert(partner.as_str());
            }
        }
        Ok(())
    }
}

impl<T: RealField + Copy> Default for BoundaryConditionSet<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::vector;
    #[test]
    fn test_boundary_condition_creation() {
        let bc = BoundaryCondition::pressure_inlet(101325.0);
        assert!(bc.is_dirichlet());
        assert!(!bc.is_neumann());
        let bc: BoundaryCondition<f64> = BoundaryCondition::wall_no_slip();
        assert!(bc.is_wall());
    }
    
    #[test]
    fn test_boundary_condition_set() {
        let mut bc_set = BoundaryConditionSet::new();
        bc_set
            .add(
                "inlet",
                BoundaryCondition::velocity_inlet(vector![1.0, 0.0, 0.0]),
            )
            .add("outlet", BoundaryCondition::pressure_outlet(0.0))
            .add("wall", BoundaryCondition::wall_no_slip());
        assert_eq!(bc_set.conditions.len(), 3);
        assert!(bc_set.get("inlet").is_some());
    }
}
