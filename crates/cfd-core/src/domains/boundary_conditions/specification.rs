//! Boundary condition specification

use super::time_dependent::TimeDependentSpec;
use crate::boundary::BoundaryCondition;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Boundary condition specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryConditionSpec<T: RealField + Copy> {
    /// Boundary condition
    pub condition: BoundaryCondition<T>,
    /// Boundary region identifier
    pub region_id: String,
    /// Time-dependent specification
    pub time_dependent: Option<TimeDependentSpec<T>>,
}

impl<T: RealField + Copy> BoundaryConditionSpec<T> {
    /// Create a new boundary condition specification
    pub fn new(condition: BoundaryCondition<T>, region_id: String) -> Self {
        Self {
            condition,
            region_id,
            time_dependent: None,
        }
    }

    /// Add time-dependent behavior
    pub fn with_time_dependence(mut self, spec: TimeDependentSpec<T>) -> Self {
        self.time_dependent = Some(spec);
        self
    }

    /// Evaluate the boundary condition at a given time
    pub fn evaluate_at_time(&self, time: T) -> BoundaryCondition<T> {
        if let Some(ref time_spec) = self.time_dependent {
            time_spec.apply_to_condition(&self.condition, time)
        } else {
            self.condition.clone()
        }
    }

    /// Get the boundary condition (returns reference when not time-dependent)
    pub fn get_condition(&self) -> &BoundaryCondition<T> {
        &self.condition
    }

    /// Check if this specification is time-dependent
    pub fn is_time_dependent(&self) -> bool {
        self.time_dependent.is_some()
    }

    /// Get the base condition type
    pub fn condition_type(&self) -> &str {
        match self.condition {
            BoundaryCondition::Dirichlet { .. } => "Dirichlet",
            BoundaryCondition::Neumann { .. } => "Neumann",
            BoundaryCondition::Robin { .. } => "Robin",
            BoundaryCondition::Periodic { .. } => "Periodic",
            BoundaryCondition::Outflow => "Outflow",
            BoundaryCondition::VelocityInlet { .. } => "VelocityInlet",
            BoundaryCondition::PressureInlet { .. } => "PressureInlet",
            BoundaryCondition::PressureOutlet { .. } => "PressureOutlet",
            BoundaryCondition::MassFlowInlet { .. } => "MassFlowInlet",
            BoundaryCondition::VolumeFlowInlet { .. } => "VolumeFlowInlet",
            BoundaryCondition::Wall { .. } => "Wall",
            BoundaryCondition::Symmetry => "Symmetry",
        }
    }
}
