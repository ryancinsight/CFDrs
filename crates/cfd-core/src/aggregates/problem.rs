//! Problem aggregate for CFD problem definitions

use crate::boundary::BoundaryCondition;
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::fluid::ConstantPropertyFluid;
use crate::values::{Pressure, Velocity};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Problem aggregate that defines a CFD problem setup
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProblemAggregate<T: RealField + Copy, D: Domain<T>> {
    /// Problem name
    pub name: String,
    /// Problem description
    pub description: String,
    /// Computational domain
    pub domain: D,
    /// Fluid properties (constant property model)
    pub fluid: ConstantPropertyFluid<T>,
    /// Initial conditions
    pub initial_conditions: InitialConditions<T>,
    /// Boundary conditions
    pub boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Problem type
    pub problem_type: ProblemType,
}

/// Initial conditions for the problem
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InitialConditions<T: RealField + Copy> {
    /// Initial velocity field
    pub velocity: Velocity<T>,
    /// Initial pressure field
    pub pressure: Pressure<T>,
    /// Initial temperature (if thermal)
    pub temperature: Option<T>,
}

/// Problem type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ProblemType {
    /// Steady-state flow
    SteadyState,
    /// Transient flow
    Transient,
    /// Thermal flow
    Thermal,
    /// Multiphase flow
    Multiphase,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float, D: Domain<T>> ProblemAggregate<T, D> {
    /// Create a new problem aggregate
    pub fn new(name: String, domain: D, fluid: ConstantPropertyFluid<T>) -> Self {
        Self {
            name,
            description: String::new(),
            domain,
            fluid,
            initial_conditions: InitialConditions::default(),
            boundary_conditions: HashMap::new(),
            problem_type: ProblemType::Transient,
        }
    }

    /// Set problem description
    #[must_use]
    pub fn with_description(mut self, description: String) -> Self {
        self.description = description;
        self
    }

    /// Set initial conditions
    #[must_use]
    pub fn with_initial_conditions(mut self, conditions: InitialConditions<T>) -> Self {
        self.initial_conditions = conditions;
        self
    }

    /// Add boundary condition
    ///
    /// # Errors
    /// Returns an error if:
    /// - A boundary condition with the same name already exists
    pub fn add_boundary_condition(
        &mut self,
        name: String,
        condition: BoundaryCondition<T>,
    ) -> Result<()> {
        if self.boundary_conditions.contains_key(&name) {
            return Err(Error::InvalidConfiguration(format!(
                "Boundary condition '{name}' already exists"
            )));
        }
        self.boundary_conditions.insert(name, condition);
        Ok(())
    }

    /// Set problem type
    #[must_use]
    pub fn with_type(mut self, problem_type: ProblemType) -> Self {
        self.problem_type = problem_type;
        self
    }

    /// Validate problem setup
    ///
    /// # Errors
    /// Returns an error if:
    /// - Domain has zero or negative volume
    /// - Fluid density is zero or negative
    /// - Fluid viscosity is zero or negative
    /// - No boundary conditions are defined
    /// - Initial conditions are invalid
    pub fn validate(&self) -> Result<()> {
        // Check domain validity
        if self.domain.volume() <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Domain has zero volume".to_string(),
            ));
        }

        // Check fluid properties
        if self.fluid.density <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Fluid density must be positive".to_string(),
            ));
        }

        // Check boundary conditions
        if self.boundary_conditions.is_empty() {
            return Err(Error::InvalidConfiguration(
                "No boundary conditions specified".to_string(),
            ));
        }

        Ok(())
    }

    /// Check if problem is thermal
    pub fn is_thermal(&self) -> bool {
        matches!(self.problem_type, ProblemType::Thermal)
    }

    /// Check if problem is multiphase
    pub fn is_multiphase(&self) -> bool {
        matches!(self.problem_type, ProblemType::Multiphase)
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for InitialConditions<T> {
    fn default() -> Self {
        Self {
            velocity: Velocity::zero(),
            pressure: Pressure::zero(),
            temperature: None,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> InitialConditions<T> {
    /// Create with velocity
    pub fn with_velocity(velocity: Velocity<T>) -> Self {
        Self {
            velocity,
            ..Default::default()
        }
    }

    /// Create with pressure
    pub fn with_pressure(pressure: Pressure<T>) -> Self {
        Self {
            pressure,
            ..Default::default()
        }
    }

    /// Set temperature
    #[must_use]
    pub fn with_temperature(mut self, temperature: T) -> Self {
        self.temperature = Some(temperature);
        self
    }
}
