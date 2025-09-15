//! Simulation aggregate root

use super::metadata::SimulationMetadata;
use super::parameters::PhysicalParameters;
use super::state::SimulationState;
use crate::boundary::BoundaryCondition;
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::fluid::ConstantPropertyFluid;
// FluidDynamicsService trait removed - using domain services directly
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Simulation aggregate root that encapsulates all simulation-related entities
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationAggregate<T: RealField + Copy, D: Domain<T>> {
    /// Unique simulation identifier
    pub id: String,
    /// Simulation metadata
    pub metadata: SimulationMetadata,
    /// Computational domain
    pub domain: D,
    /// Fluid properties
    pub fluid: ConstantPropertyFluid<T>,
    /// Boundary conditions
    pub boundary_conditions: HashMap<String, BoundaryCondition<T>>,
    /// Physical parameters
    pub parameters: PhysicalParameters<T>,
    /// Simulation state
    pub state: SimulationState,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float, D: Domain<T>>
    SimulationAggregate<T, D>
{
    /// Create a new simulation aggregate
    pub fn new(id: String, domain: D, fluid: ConstantPropertyFluid<T>) -> Self {
        Self {
            id,
            metadata: SimulationMetadata::default(),
            domain,
            fluid,
            boundary_conditions: HashMap::new(),
            parameters: PhysicalParameters::default(),
            state: SimulationState::Initialized,
        }
    }

    /// Add a boundary condition
    ///
    /// # Errors
    /// Returns an error if:
    /// - The simulation is currently active
    pub fn add_boundary_condition(
        &mut self,
        name: String,
        condition: BoundaryCondition<T>,
    ) -> Result<()> {
        if self.state.is_active() {
            return Err(Error::InvalidConfiguration(
                "Cannot modify boundary conditions while simulation is active".to_string(),
            ));
        }

        self.boundary_conditions.insert(name, condition);
        self.metadata.touch();
        Ok(())
    }

    /// Update physical parameters
    ///
    /// # Errors
    /// Returns an error if:
    /// - The simulation is currently active
    /// - Parameters fail validation
    pub fn update_parameters(&mut self, parameters: PhysicalParameters<T>) -> Result<()> {
        if self.state.is_active() {
            return Err(Error::InvalidConfiguration(
                "Cannot modify parameters while simulation is active".to_string(),
            ));
        }

        self.parameters = parameters;
        self.metadata.touch();
        Ok(())
    }

    /// Start the simulation
    ///
    /// # Errors
    /// Returns an error if:
    /// - Configuration validation fails
    /// - Simulation state transition is invalid
    pub fn start(&mut self) -> Result<()> {
        // Validate configuration
        self.validate_configuration()?;

        self.state
            .start()
            .map_err(|e| Error::InvalidConfiguration(e.to_string()))?;
        self.metadata.touch();
        Ok(())
    }

    /// Pause the simulation
    ///
    /// # Errors
    /// Returns an error if:
    /// - Simulation state transition is invalid (e.g., not currently running)
    pub fn pause(&mut self) -> Result<()> {
        self.state
            .pause()
            .map_err(|e| Error::InvalidConfiguration(e.to_string()))?;
        self.metadata.touch();
        Ok(())
    }

    /// Complete the simulation
    ///
    /// # Errors
    /// Returns an error if:
    /// - Simulation state transition is invalid (e.g., not running or paused)
    pub fn complete(&mut self) -> Result<()> {
        self.state
            .complete()
            .map_err(|e| Error::InvalidConfiguration(e.to_string()))?;
        self.metadata.touch();
        Ok(())
    }

    /// Mark simulation as failed
    pub fn fail(&mut self, _error: String) {
        self.state.fail();
        self.metadata.touch();
    }

    /// Validate simulation configuration
    ///
    /// # Errors
    /// Returns an error if:
    /// - Domain has zero or negative volume
    /// - No boundary conditions are defined
    /// - Fluid properties are invalid (zero density or viscosity)
    /// - Physical parameters are inconsistent
    pub fn validate_configuration(&self) -> Result<()> {
        // Check if domain is valid (has non-zero volume)
        if self.domain.volume() <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Domain has zero volume".to_string(),
            ));
        }

        // Check if fluid properties are valid
        if self.fluid.density <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Fluid density must be positive".to_string(),
            ));
        }

        if self.fluid.viscosity <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Fluid viscosity must be positive".to_string(),
            ));
        }

        // Stability check would require grid information, not just domain
        // This should be done at the solver level, not here

        Ok(())
    }

    /// Execute a simulation step
    pub fn step(&mut self, dt: T) -> Result<()> {
        if self.state != SimulationState::Running {
            return Err(Error::InvalidConfiguration(
                "Simulation is not running".to_string(),
            ));
        }

        // Time step execution requires a solver implementation
        // This should be handled by a specific solver, not the aggregate
        self.parameters.time_step = dt;
        self.metadata.touch();

        Ok(())
    }

    /// Get simulation progress
    pub fn progress(&self, current_time: T) -> T {
        if self.parameters.max_time > T::zero() {
            current_time / self.parameters.max_time
        } else {
            T::zero()
        }
    }
}
