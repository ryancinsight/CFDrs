//! Domain aggregates for CFD simulations.
//!
//! This module provides aggregate roots that encapsulate related entities
//! and enforce business rules and invariants.

use crate::{
    Error, Result, Fluid, BoundaryCondition, Domain, 
    FluidDynamicsService, ReynoldsNumber, Pressure, Velocity
};
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Simulation aggregate root that encapsulates all simulation-related entities
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationAggregate<T: RealField, D: Domain<T>> {
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

impl<T: RealField + FromPrimitive, D: Domain<T>> SimulationAggregate<T, D> {
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
        // Validate boundary condition compatibility
        self.validate_boundary_condition(&condition)?;
        
        self.boundary_conditions.insert(name, condition);
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339(); // Update timestamp
        Ok(())
    }

    /// Set reference conditions
    pub fn set_reference_conditions(
        &mut self,
        pressure: Pressure<T>,
        velocity: Velocity<T>,
        length: T,
    ) -> Result<()> {
        self.parameters.reference_pressure = pressure;
        self.parameters.reference_velocity = velocity.magnitude();
        self.parameters.reference_length = length;
        
        // Calculate Reynolds number
        let reynolds_value = FluidDynamicsService::reynolds_number(
            &self.fluid,
            self.parameters.reference_velocity.clone(),
            self.parameters.reference_length.clone(),
        );
        self.parameters.reynolds_number = Some(ReynoldsNumber::new(reynolds_value)?);
        
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339(); // Update timestamp
        Ok(())
    }

    /// Get characteristic Reynolds number
    pub fn reynolds_number(&self) -> Option<&ReynoldsNumber<T>> {
        self.parameters.reynolds_number.as_ref()
    }

    /// Check if simulation is ready to run
    pub fn is_ready(&self) -> bool {
        !self.boundary_conditions.is_empty() && 
        self.parameters.reynolds_number.is_some() &&
        matches!(self.state, SimulationState::Initialized | SimulationState::Configured)
    }

    /// Transition to configured state
    pub fn configure(&mut self) -> Result<()> {
        if !self.is_ready() {
            return Err(Error::InvalidConfiguration(
                "Simulation is not properly configured".to_string()
            ));
        }
        
        self.state = SimulationState::Configured;
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339(); // Update timestamp
        Ok(())
    }

    /// Start simulation
    pub fn start(&mut self) -> Result<()> {
        if !matches!(self.state, SimulationState::Configured) {
            return Err(Error::InvalidConfiguration(
                "Simulation must be configured before starting".to_string()
            ));
        }
        
        self.state = SimulationState::Running;
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339(); // Update timestamp
        Ok(())
    }

    /// Complete simulation
    pub fn complete(&mut self) -> Result<()> {
        if !matches!(self.state, SimulationState::Running) {
            return Err(Error::InvalidConfiguration(
                "Simulation must be running to complete".to_string()
            ));
        }
        
        self.state = SimulationState::Completed;
        self.metadata.modified_at = chrono::Utc::now().to_rfc3339(); // Update timestamp
        Ok(())
    }

    /// Validate boundary condition compatibility
    fn validate_boundary_condition(&self, condition: &BoundaryCondition<T>) -> Result<()> {
        // Add domain-specific validation logic here
        match condition {
            BoundaryCondition::VelocityInlet { velocity: _ } => {
                // Validate inlet conditions
                Ok(())
            }
            BoundaryCondition::PressureOutlet { pressure: _ } => {
                // Validate outlet conditions
                Ok(())
            }
            BoundaryCondition::Wall { wall_type: _ } => {
                // Validate wall conditions
                Ok(())
            }
            _ => Ok(()),
        }
    }
}

/// Simulation metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationMetadata {
    /// Simulation name
    pub name: String,
    /// Description
    pub description: String,
    /// Creation timestamp
    pub created_at: String,
    /// Last modified timestamp
    pub modified_at: String,
    /// Version
    pub version: String,
    /// Tags for categorization
    pub tags: Vec<String>,
}

impl Default for SimulationMetadata {
    fn default() -> Self {
        Self {
            name: "Untitled Simulation".to_string(),
            description: String::new(),
            created_at: chrono::Utc::now().to_rfc3339(),
            modified_at: chrono::Utc::now().to_rfc3339(),
            version: "1.0.0".to_string(),
            tags: Vec::new(),
        }
    }
}

/// Physical parameters for the simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicalParameters<T: RealField> {
    /// Reference pressure
    pub reference_pressure: Pressure<T>,
    /// Reference velocity magnitude
    pub reference_velocity: T,
    /// Reference length
    pub reference_length: T,
    /// Reynolds number
    pub reynolds_number: Option<ReynoldsNumber<T>>,
    /// Gravitational acceleration
    pub gravity: Option<Vector3<T>>,
    /// Time step for transient simulations
    pub time_step: Option<T>,
}

impl<T: RealField + FromPrimitive> Default for PhysicalParameters<T> {
    fn default() -> Self {
        Self {
            reference_pressure: Pressure::pascals(T::from_f64(101_325.0).unwrap()),
            reference_velocity: T::one(),
            reference_length: T::one(),
            reynolds_number: None,
            gravity: Some(Vector3::new(T::zero(), T::from_f64(-9.81).unwrap(), T::zero())),
            time_step: None,
        }
    }
}

/// Simulation state enumeration
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SimulationState {
    /// Simulation is initialized but not configured
    Initialized,
    /// Simulation is configured and ready to run
    Configured,
    /// Simulation is currently running
    Running,
    /// Simulation is paused
    Paused,
    /// Simulation completed successfully
    Completed,
    /// Simulation failed with error
    Failed(String),
}

/// Mesh aggregate for managing mesh-related entities
#[derive(Debug, Clone)]
pub struct MeshAggregate<T: RealField> {
    /// Mesh identifier
    pub id: String,
    /// Mesh metadata
    pub metadata: MeshMetadata,
    /// Mesh quality metrics
    pub quality_metrics: Option<MeshQualityMetrics<T>>,
    /// Mesh state
    pub state: MeshState,
}

impl<T: RealField> MeshAggregate<T> {
    /// Create a new mesh aggregate
    pub fn new(id: String) -> Self {
        Self {
            id,
            metadata: MeshMetadata::default(),
            quality_metrics: None,
            state: MeshState::Created,
        }
    }

    /// Set quality metrics
    pub fn set_quality_metrics(&mut self, metrics: MeshQualityMetrics<T>) {
        self.quality_metrics = Some(metrics);
        self.state = MeshState::Analyzed;
    }

    /// Check if mesh is suitable for simulation
    pub fn is_suitable_for_simulation(&self) -> bool {
        if let Some(ref metrics) = self.quality_metrics {
            metrics.overall_quality_score > T::from_f64(0.7).unwrap()
        } else {
            false
        }
    }
}

/// Mesh metadata
#[derive(Debug, Clone)]
pub struct MeshMetadata {
    /// Mesh name
    pub name: String,
    /// Number of vertices
    pub num_vertices: usize,
    /// Number of cells
    pub num_cells: usize,
    /// Number of faces
    pub num_faces: usize,
    /// Mesh type
    pub mesh_type: MeshType,
}

impl Default for MeshMetadata {
    fn default() -> Self {
        Self {
            name: "Untitled Mesh".to_string(),
            num_vertices: 0,
            num_cells: 0,
            num_faces: 0,
            mesh_type: MeshType::Unstructured,
        }
    }
}

/// Mesh type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshType {
    /// Structured mesh
    Structured,
    /// Unstructured mesh
    Unstructured,
    /// Hybrid mesh
    Hybrid,
}

/// Mesh quality metrics
#[derive(Debug, Clone)]
pub struct MeshQualityMetrics<T: RealField> {
    /// Overall quality score (0-1)
    pub overall_quality_score: T,
    /// Minimum aspect ratio
    pub min_aspect_ratio: T,
    /// Maximum aspect ratio
    pub max_aspect_ratio: T,
    /// Average skewness
    pub average_skewness: T,
    /// Maximum skewness
    pub max_skewness: T,
    /// Minimum orthogonality
    pub min_orthogonality: T,
}

/// Mesh state enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshState {
    /// Mesh created but not analyzed
    Created,
    /// Mesh analyzed for quality
    Analyzed,
    /// Mesh validated for simulation
    Validated,
    /// Mesh rejected due to poor quality
    Rejected,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::Domain1D;

    #[test]
    fn test_simulation_aggregate_lifecycle() {
        let domain = Domain1D::new(0.0, 1.0);
        let fluid = Fluid::water();
        let mut sim = SimulationAggregate::new("test-sim".to_string(), domain, fluid);

        // Initially not ready
        assert!(!sim.is_ready());

        // Add boundary condition
        sim.add_boundary_condition(
            "inlet".to_string(),
            BoundaryCondition::velocity_inlet(Vector3::new(1.0, 0.0, 0.0))
        ).unwrap();

        // Set reference conditions
        sim.set_reference_conditions(
            Pressure::pascals(101_325.0),
            Velocity::new(1.0, 0.0, 0.0),
            1.0
        ).unwrap();

        // Now should be ready
        assert!(sim.is_ready());

        // Configure and start
        sim.configure().unwrap();
        assert_eq!(sim.state, SimulationState::Configured);

        sim.start().unwrap();
        assert_eq!(sim.state, SimulationState::Running);

        sim.complete().unwrap();
        assert_eq!(sim.state, SimulationState::Completed);
    }
}
