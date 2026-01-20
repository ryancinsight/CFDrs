use super::solver::MomentumSolver;
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use std::collections::HashMap;

/// Builder for configuring boundary conditions for the MomentumSolver
pub struct BoundarySetup<T: RealField + Copy> {
    conditions: HashMap<String, BoundaryCondition<T>>,
}

impl<T: RealField + Copy> BoundarySetup<T> {
    /// Create a new BoundarySetup
    pub fn new() -> Self {
        Self {
            conditions: HashMap::new(),
        }
    }

    /// Set the north boundary condition
    pub fn north(mut self, bc: BoundaryCondition<T>) -> Self {
        self.conditions.insert("north".to_string(), bc);
        self
    }

    /// Set the south boundary condition
    pub fn south(mut self, bc: BoundaryCondition<T>) -> Self {
        self.conditions.insert("south".to_string(), bc);
        self
    }

    /// Set the east boundary condition
    pub fn east(mut self, bc: BoundaryCondition<T>) -> Self {
        self.conditions.insert("east".to_string(), bc);
        self
    }

    /// Set the west boundary condition
    pub fn west(mut self, bc: BoundaryCondition<T>) -> Self {
        self.conditions.insert("west".to_string(), bc);
        self
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> BoundarySetup<T> {
    /// Apply the configured boundaries to the solver
    pub fn apply(self, solver: &mut MomentumSolver<T>) {
        for (name, bc) in self.conditions {
            solver.set_boundary(name, bc);
        }
    }
}

impl<T: RealField + Copy> Default for BoundarySetup<T> {
    fn default() -> Self {
        Self::new()
    }
}
