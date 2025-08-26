//! Domain aggregates for CFD simulations.
//!
//! This module provides aggregate roots that encapsulate related entities
//! and enforce business rules and invariants.

use crate::boundary::BoundaryCondition;
use crate::domain::Domain;
use crate::error::{Error, Result};
use crate::fluid::Fluid;
use crate::services::FluidDynamicsService;
use crate::values::{Pressure, ReynoldsNumber, Velocity};
use nalgebra::{RealField, Vector3};
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

    /// Add a boundary condition with validation
    pub fn add_boundary_condition(
        &mut self,
        name: String,
        condition: BoundaryCondition<T>,
    ) -> Result<()> {
        Self::validate_boundary_condition(&condition)?;
        self.boundary_conditions.insert(name, condition);
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339();
        Ok(())
    }

    /// Set reference conditions
    pub fn set_reference_conditions(
        &mut self,
        pressure: Pressure<T>,
        velocity: &Velocity<T>,
        length: T,
    ) -> Result<()> {
        self.parameters.reference_pressure = pressure;
        self.parameters.reference_velocity = velocity.magnitude();
        self.parameters.reference_length = length;
        
        let reynolds_value = FluidDynamicsService::reynolds_number(
            &self.fluid,
            self.parameters.reference_velocity,
            self.parameters.reference_length,
        );
        
        self.parameters.reynolds_number = Some(ReynoldsNumber::new(
            reynolds_value,
            crate::values::FlowGeometry::Pipe,
        )?);
        
        Ok(())
    }

    /// Get characteristic Reynolds number
    pub fn reynolds_number(&self) -> Option<&ReynoldsNumber<T>> {
        self.parameters.reynolds_number.as_ref()
    }

    /// Check if simulation is ready to run
    pub fn is_ready(&self) -> bool {
        !self.boundary_conditions.is_empty()
            && self.parameters.reynolds_number.is_some()
            && matches!(
                self.state,
                SimulationState::Initialized | SimulationState::Configured
            )
    }

    /// Transition to configured state
    pub fn configure(&mut self) -> Result<()> {
        if !self.is_ready() {
            return Err(Error::InvalidConfiguration(
                "Simulation is not properly configured".to_string(),
            ));
        }
        self.state = SimulationState::Configured;
        Ok(())
    }

    /// Start the simulation
    pub fn start(&mut self) -> Result<()> {
        if self.state != SimulationState::Configured {
            return Err(Error::InvalidState(
                "Simulation must be configured before starting".to_string(),
            ));
        }
        self.state = SimulationState::Running;
        Ok(())
    }

    /// Complete the simulation
    pub fn complete(&mut self) -> Result<()> {
        if self.state != SimulationState::Running {
            return Err(Error::InvalidState(
                "Simulation must be running to complete".to_string(),
            ));
        }
        self.state = SimulationState::Completed;
        Ok(())
    }

    /// Validate boundary condition
    fn validate_boundary_condition(condition: &BoundaryCondition<T>) -> Result<()> {
        // Validation logic here
        Ok(())
    }
}

/// Simulation metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationMetadata {
    /// Creation timestamp
    pub created_at: String,
    /// Last modification timestamp
    pub modified_at: String,
    /// Simulation description
    pub description: String,
    /// Version
    pub version: String,
}

impl Default for SimulationMetadata {
    fn default() -> Self {
        Self {
            created_at: chrono::Utc::now().to_rfc3339(),
            modified_at: chrono::Utc::now().to_rfc3339(),
            description: String::new(),
            version: "0.1.0".to_string(),
        }
}

/// Physical parameters for simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicalParameters<T: RealField + Copy> {
    /// Reference pressure
    pub reference_pressure: Pressure<T>,
    /// Reference velocity magnitude
    pub reference_velocity: T,
    /// Reference length scale
    pub reference_length: T,
    /// Reynolds number
    pub reynolds_number: Option<ReynoldsNumber<T>>,
}

impl<T: RealField + FromPrimitive + Copy> Default for PhysicalParameters<T> {
    fn default() -> Self {
        Self {
            reference_pressure: Pressure::pascals(T::from_f64(101325.0).unwrap_or_else(T::one)),
            reference_velocity: T::zero(),
            reference_length: T::one(),
            reynolds_number: None,
        }
}

/// Simulation state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimulationState {
    /// Initial state
    Initialized,
    /// Configured and ready to run
    Configured,
    /// Currently running
    Running,
    /// Paused
    Paused,
    /// Completed successfully
    Completed,
    /// Failed with error
    Failed,
}

/// Mesh aggregate for managing computational meshes
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshAggregate<T: RealField + Copy> {
    /// Unique mesh identifier
    pub id: String,
    /// Mesh metadata
    pub metadata: MeshMetadata,
    /// Number of nodes
    pub nodes: usize,
    /// Number of elements
    pub elements: usize,
    /// Mesh quality metrics
    pub quality: MeshQualityMetrics<T>,
    /// Mesh state
    pub state: MeshState,
}

impl<T: RealField + Copy> MeshAggregate<T> {
    /// Create a new mesh aggregate
    pub fn new(id: String, nodes: usize, elements: usize) -> Self {
        Self {
            id,
            metadata: MeshMetadata::default(),
            nodes,
            elements,
            quality: MeshQualityMetrics::default(),
            state: MeshState::Created,
        }
}

/// Mesh metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshMetadata {
    /// Mesh type
    pub mesh_type: MeshType,
    /// Dimension (1D, 2D, 3D)
    pub dimension: usize,
    /// Generator used
    pub generator: String,
}

impl Default for MeshMetadata {
    fn default() -> Self {
        Self {
            mesh_type: MeshType::Structured,
            dimension: 3,
            generator: "internal".to_string(),
        }
}

/// Mesh type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeshType {
    /// Structured mesh
    Structured,
    /// Unstructured mesh
    Unstructured,
    /// Hybrid mesh
    Hybrid,
}

/// Mesh quality metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshQualityMetrics<T: RealField + Copy> {
    /// Minimum quality
    pub min_quality: T,
    /// Maximum quality
    pub max_quality: T,
    /// Average quality
    pub avg_quality: T,
    /// Minimum orthogonality
    pub min_orthogonality: T,
    /// Maximum aspect ratio
    pub max_aspect_ratio: T,
}

impl<T: RealField + Copy> Default for MeshQualityMetrics<T> {
    fn default() -> Self {
        Self {
            min_quality: T::zero(),
            max_quality: T::one(),
            avg_quality: T::from_f64(0.5).unwrap_or_else(T::one),
            min_orthogonality: T::zero(),
            max_aspect_ratio: T::one(),
        }
}

/// Mesh state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MeshState {
    /// Just created
    Created,
    /// Validated
    Validated,
    /// Refined
    Refined,
    /// Invalid
    Invalid,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::Domain1D;

    #[test]
    fn test_simulation_aggregate_lifecycle() {
        let domain = Domain1D::new(0.0, 1.0);
        let fluid = Fluid::water().expect("Failed to create water fluid");
        let mut sim = SimulationAggregate::new("test-sim".to_string(), domain, fluid);
        
        assert!(!sim.is_ready());
        
        sim.add_boundary_condition(
            "inlet".to_string(),
            BoundaryCondition::inlet(1.0, 0.0, 0.0),
        ).expect("Failed to add boundary condition");
        
        sim.set_reference_conditions(
            Pressure::pascals(101_325.0),
            &Velocity::new(1.0, 0.0, 0.0),
            1.0,
        ).expect("Failed to set reference conditions");
        
        assert!(sim.is_ready());
        
        sim.configure().expect("Failed to configure");
        assert_eq!(sim.state, SimulationState::Configured);
        
        sim.start().expect("Failed to start");
        assert_eq!(sim.state, SimulationState::Running);
        
        sim.complete().expect("Failed to complete");
        assert_eq!(sim.state, SimulationState::Completed);

    }


}

}
}
}
}
}
