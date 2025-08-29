//! Problem definition and configuration.

use crate::boundary::BoundaryConditionSet;
use crate::domain::Domain;
use crate::error::Result;
use crate::fluid::ConstantPropertyFluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use std::sync::Arc;

/// Trait for defining CFD problems
pub trait Problem<T: RealField + Copy>: Send + Sync {
    /// Domain type for this problem
    type Domain: Domain<T>;
    /// State type for this problem
    type State;

    /// Get the computational domain
    fn domain(&self) -> &Self::Domain;

    /// Get the fluid properties
    fn fluid(&self) -> &Fluid<T>;

    /// Get the boundary conditions
    fn boundary_conditions(&self) -> &BoundaryConditionSet<T>;

    /// Get the initial state
    fn initial_state(&self) -> Result<Self::State>;

    /// Validate the problem setup
    fn validate(&self) -> Result<()> {
        // Validate boundary conditions
        self.boundary_conditions().validate_periodic()?;
        Ok(())
    }
}

/// Generic problem configuration
#[derive(Debug, Clone)]
pub struct ProblemConfig<T: RealField + Copy, D: Domain<T>> {
    /// Computational domain
    pub domain: Arc<D>,
    /// Fluid properties
    pub fluid: ConstantPropertyFluid<T>,
    /// Boundary conditions
    pub boundary_conditions: BoundaryConditionSet<T>,
    /// Problem parameters
    pub parameters: ProblemParameters<T>,
}

/// Common problem parameters
#[derive(Debug, Clone)]
pub struct ProblemParameters<T: RealField + Copy> {
    /// Reference pressure [Pa]
    pub reference_pressure: T,
    /// Reference temperature [K]
    pub reference_temperature: Option<T>,
    /// Gravity vector [m/s²]
    pub gravity: Option<nalgebra::Vector3<T>>,
    /// Time-dependent problem
    pub transient: bool,
    /// Include energy equation
    pub energy: bool,
}

impl<T: RealField + FromPrimitive + Copy> Default for ProblemParameters<T> {
    fn default() -> Self {
        Self {
            reference_pressure: T::from_f64(101_325.0).unwrap_or_else(|| T::one()), // 1 atm
            reference_temperature: Some(T::from_f64(293.15).unwrap_or_else(|| T::one())), // 20°C
            gravity: None,
            transient: false,
            energy: false,
        }
    }
}

/// Problem builder for convenient construction
pub struct ProblemBuilder<T: RealField + Copy, D: Domain<T>> {
    domain: Option<Arc<D>>,
    fluid: Option<Fluid<T>>,
    boundary_conditions: BoundaryConditionSet<T>,
    parameters: ProblemParameters<T>,
}

impl<T: RealField + Copy, D: Domain<T>> ProblemBuilder<T, D> {
    /// Create a new problem builder
    #[must_use]
    pub fn new() -> Self {
        Self {
            domain: None,
            fluid: None,
            boundary_conditions: BoundaryConditionSet::new(),
            parameters: ProblemParameters::default(),
        }
    }

    /// Set the domain
    pub fn domain(mut self, domain: D) -> Self {
        self.domain = Some(Arc::new(domain));
        self
    }

    /// Set the fluid
    pub fn fluid(mut self, fluid: ConstantPropertyFluid<T>) -> Self {
        self.fluid = Some(fluid);
        self
    }

    /// Add a boundary condition
    pub fn boundary_condition(
        mut self,
        name: impl Into<String>,
        condition: crate::boundary::BoundaryCondition<T>,
    ) -> Self {
        self.boundary_conditions.add(name, condition);
        self
    }

    /// Set reference pressure
    pub fn reference_pressure(mut self, pressure: T) -> Self {
        self.parameters.reference_pressure = pressure;
        self
    }

    /// Set gravity
    pub fn gravity(mut self, gravity: nalgebra::Vector3<T>) -> Self {
        self.parameters.gravity = Some(gravity);
        self
    }

    /// Enable transient simulation
    pub fn transient(mut self, transient: bool) -> Self {
        self.parameters.transient = transient;
        self
    }

    /// Build the problem configuration
    pub fn build(self) -> Result<ProblemConfig<T, D>> {
        let domain = self.domain.ok_or_else(|| {
            crate::error::Error::InvalidConfiguration("Domain not set".to_string())
        })?;

        let fluid = self.fluid.ok_or_else(|| {
            crate::error::Error::InvalidConfiguration("Fluid not set".to_string())
        })?;

        let config = ProblemConfig {
            domain,
            fluid,
            boundary_conditions: self.boundary_conditions,
            parameters: self.parameters,
        };

        Ok(config)
    }
}

impl<T: RealField + Copy, D: Domain<T>> Default for ProblemBuilder<T, D> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{boundary::BoundaryCondition, domain::Domain2D};
    use nalgebra::vector;

    #[test]
    fn test_problem_builder() {
        let problem = ProblemBuilder::new()
            .domain(Domain2D::from_scalars(0.0, 0.0, 1.0, 1.0))
            .fluid(ConstantPropertyFluid::water_20c())
            .boundary_condition(
                "inlet",
                BoundaryCondition::velocity_inlet(vector![1.0, 0.0, 0.0]),
            )
            .boundary_condition("outlet", BoundaryCondition::pressure_outlet(0.0))
            .boundary_condition("walls", BoundaryCondition::wall_no_slip())
            .reference_pressure(101325.0)
            .transient(true)
            .build()
            .expect("CRITICAL: Add proper error handling");

        assert_eq!(problem.fluid.name, "Water at 20°C");
        assert_eq!(problem.boundary_conditions.conditions.len(), 3);
        assert!(problem.parameters.transient);
    }
}
