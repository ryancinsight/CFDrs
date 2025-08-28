//! Simulation aggregate root

use super::metadata::SimulationMetadata;
use super::parameters::PhysicalParameters;
use super::state::SimulationState;
use crate::boundary::BoundaryCondition;
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::fluid::Fluid;
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
    pub fluid: Fluid<T>,
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
    pub fn new(id: String, domain: D, fluid: Fluid<T>) -> Self {
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
    pub fn pause(&mut self) -> Result<()> {
        self.state
            .pause()
            .map_err(|e| Error::InvalidConfiguration(e.to_string()))?;
        self.metadata.touch();
        Ok(())
    }

    /// Complete the simulation
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
    pub fn validate_configuration(&self) -> Result<()> {
        // Check if domain is valid
        if self.domain.num_cells() == 0 {
            return Err(Error::InvalidConfiguration(
                "Domain has no cells".to_string(),
            ));
        }

        // Check if fluid properties are valid
        if self.fluid.density() <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Fluid density must be positive".to_string(),
            ));
        }

        if self.fluid.dynamic_viscosity() <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Fluid viscosity must be positive".to_string(),
            ));
        }

        // Check stability
        let dx = self.domain.cell_size();
        if !self.parameters.is_stable(dx) {
            return Err(Error::InvalidConfiguration(format!(
                "CFL condition violated: CFL = {}",
                self.parameters.cfl_number(dx)
            )));
        }

        Ok(())
    }

    /// Execute a simulation step
    pub fn step(&mut self, dt: T) -> Result<()> {
        if self.state != SimulationState::Running {
            return Err(Error::InvalidConfiguration(
                "Simulation is not running".to_string(),
            ));
        }

        // Simulation step logic would go here
        // This is a placeholder for the actual implementation
        let _ = dt;

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
