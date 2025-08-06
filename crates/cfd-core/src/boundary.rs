//! Boundary condition types and implementations.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Boundary condition types for CFD simulations
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum BoundaryCondition<T: RealField> {
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
        /// Paired boundary identifier
        paired_with: String,
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
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum WallType<T: RealField> {
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

impl<T: RealField> BoundaryCondition<T> {
    /// Create a pressure inlet boundary condition
    pub fn pressure_inlet(pressure: T) -> Self {
        Self::PressureInlet { pressure }
    }

    /// Create a pressure outlet boundary condition
    pub fn pressure_outlet(pressure: T) -> Self {
        Self::PressureOutlet { pressure }
    }

    /// Create a velocity inlet boundary condition
    pub fn velocity_inlet(velocity: Vector3<T>) -> Self {
        Self::VelocityInlet { velocity }
    }

    /// Create a no-slip wall boundary condition
    pub fn wall_no_slip() -> Self {
        Self::Wall {
            wall_type: WallType::NoSlip,
        }
    }

    /// Create a slip wall boundary condition
    pub fn wall_slip() -> Self {
        Self::Wall {
            wall_type: WallType::Slip,
        }
    }

    /// Create a moving wall boundary condition
    pub fn wall_moving(velocity: Vector3<T>) -> Self {
        Self::Wall {
            wall_type: WallType::Moving { velocity },
        }
    }

    /// Check if this is a Dirichlet-type boundary condition
    pub fn is_dirichlet(&self) -> bool {
        matches!(
            self,
            Self::Dirichlet { .. }
                | Self::VelocityInlet { .. }
                | Self::PressureInlet { .. }
                | Self::Wall { .. }
        )
    }

    /// Check if this is a Neumann-type boundary condition
    pub fn is_neumann(&self) -> bool {
        matches!(
            self,
            Self::Neumann { .. } | Self::Outflow | Self::Symmetry
        )
    }

    /// Check if this is a wall boundary condition
    pub fn is_wall(&self) -> bool {
        matches!(self, Self::Wall { .. })
    }
}

/// Boundary condition set for a complete problem
#[derive(Debug, Clone)]
pub struct BoundaryConditionSet<T: RealField> {
    /// Map of boundary names to conditions
    pub conditions: indexmap::IndexMap<String, BoundaryCondition<T>>,
}

impl<T: RealField> BoundaryConditionSet<T> {
    /// Create a new empty boundary condition set
    pub fn new() -> Self {
        Self {
            conditions: indexmap::IndexMap::new(),
        }
    }

    /// Add a boundary condition
    pub fn add(
        &mut self,
        name: impl Into<String>,
        condition: BoundaryCondition<T>,
    ) -> &mut Self {
        self.conditions.insert(name.into(), condition);
        self
    }

    /// Get a boundary condition by name
    pub fn get(&self, name: &str) -> Option<&BoundaryCondition<T>> {
        self.conditions.get(name)
    }

    /// Check if all periodic boundaries are properly paired
    pub fn validate_periodic(&self) -> Result<(), String> {
        for (name, bc) in &self.conditions {
            if let BoundaryCondition::Periodic { paired_with } = bc {
                match self.conditions.get(paired_with) {
                    Some(BoundaryCondition::Periodic {
                        paired_with: other,
                    }) if other == name => {}
                    Some(_) => {
                        return Err(format!(
                            "Boundary '{}' is paired with '{}' which is not periodic",
                            name, paired_with
                        ))
                    }
                    None => {
                        return Err(format!(
                            "Boundary '{}' is paired with non-existent boundary '{}'",
                            name, paired_with
                        ))
                    }
                }
            }
        }
        Ok(())
    }
}

impl<T: RealField> Default for BoundaryConditionSet<T> {
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

        let bc = BoundaryCondition::wall_no_slip();
        assert!(bc.is_wall());
        assert!(bc.is_dirichlet());
    }

    #[test]
    fn test_boundary_condition_set() {
        let mut bc_set = BoundaryConditionSet::new();
        bc_set
            .add("inlet", BoundaryCondition::velocity_inlet(vector![1.0, 0.0, 0.0]))
            .add("outlet", BoundaryCondition::pressure_outlet(0.0))
            .add("wall", BoundaryCondition::wall_no_slip());

        assert_eq!(bc_set.conditions.len(), 3);
        assert!(bc_set.get("inlet").is_some());
    }
}