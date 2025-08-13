//! Factory patterns implementing GRASP Creator principle.
//!
//! This module provides factory abstractions that assign creation responsibility
//! to classes that have the initializing data, following GRASP principles.

use crate::{Error, Result, Solver, SolverConfiguration};
use nalgebra::RealField;
use std::collections::HashMap;
use std::any::Any;
use std::sync::Arc;

/// Concrete factory trait for type-safe creation
/// Follows Open/Closed Principle - open for extension, closed for modification
pub trait ConcreteSolverFactory<T: RealField>: Send + Sync {
    /// Solver type this factory creates
    type Solver: Solver<T>;
    /// Configuration type for the solver
    type Config: SolverConfiguration<T>;

    /// Create a solver with the given configuration
    fn create(&self, config: Self::Config) -> Result<Self::Solver>;

    /// Get factory name for identification
    fn name(&self) -> &str;

    /// Check if this factory can create a solver for the given problem type
    fn can_handle(&self, problem_type: &str) -> bool;
}

/// Factory capability trait following Interface Segregation Principle
/// Separates capability checking from creation
pub trait FactoryCapability {
    /// Check if factory supports the given capability
    fn supports_capability(&self, capability: &str) -> bool;

    /// Get all supported capabilities
    fn capabilities(&self) -> Vec<&str>;
}

/// Type-erased solver trait for dynamic dispatch
pub trait DynamicSolver<T: RealField>: Send + Sync {
    /// Solve the problem
    fn solve(&mut self, problem: &dyn Any) -> Result<Box<dyn Any>>;
    
    /// Get solver name
    fn name(&self) -> &str;
}

/// Wrapper to convert concrete solvers to dynamic solvers
pub struct DynamicSolverWrapper<T, S>
where
    T: RealField,
    S: Solver<T>,
{
    solver: S,
    _phantom: std::marker::PhantomData<T>,
}

impl<T, S> DynamicSolverWrapper<T, S>
where
    T: RealField + 'static,
    S: Solver<T> + 'static,
    S::Problem: 'static,
    S::Solution: 'static,
{
    pub fn new(solver: S) -> Self {
        Self {
            solver,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T, S> DynamicSolver<T> for DynamicSolverWrapper<T, S>
where
    T: RealField + 'static,
    S: Solver<T> + 'static,
    S::Problem: 'static,
    S::Solution: 'static,
{
    fn solve(&mut self, problem: &dyn Any) -> Result<Box<dyn Any>> {
        let typed_problem = problem
            .downcast_ref::<S::Problem>()
            .ok_or_else(|| Error::InvalidConfiguration("Invalid problem type".into()))?;
        
        let solution = self.solver.solve(typed_problem)?;
        Ok(Box::new(solution))
    }
    
    fn name(&self) -> &str {
        self.solver.name()
    }
}

/// Type-erased factory wrapper for heterogeneous storage
pub trait DynamicFactory<T: RealField>: Send + Sync {
    /// Create a solver from a configuration
    fn create_solver(&self, config: &dyn Any) -> Result<Box<dyn DynamicSolver<T>>>;
    
    /// Get factory name
    fn name(&self) -> &str;
    
    /// Check if this factory can handle a problem type
    fn can_handle(&self, problem_type: &str) -> bool;
}

/// Wrapper to convert concrete factories to dynamic factories
pub struct DynamicFactoryWrapper<T, F>
where
    T: RealField,
    F: ConcreteSolverFactory<T>,
{
    factory: F,
    _phantom: std::marker::PhantomData<T>,
}

impl<T, F> DynamicFactoryWrapper<T, F>
where
    T: RealField + 'static,
    F: ConcreteSolverFactory<T> + 'static,
    F::Solver: 'static,
    F::Config: 'static,
    <F::Solver as Solver<T>>::Problem: 'static,
    <F::Solver as Solver<T>>::Solution: 'static,
{
    pub fn new(factory: F) -> Self {
        Self {
            factory,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T, F> DynamicFactory<T> for DynamicFactoryWrapper<T, F>
where
    T: RealField + 'static,
    F: ConcreteSolverFactory<T> + 'static,
    F::Solver: 'static,
    F::Config: 'static,
    <F::Solver as Solver<T>>::Problem: 'static,
    <F::Solver as Solver<T>>::Solution: 'static,
{
    fn create_solver(&self, config: &dyn Any) -> Result<Box<dyn DynamicSolver<T>>> {
        let typed_config = config
            .downcast_ref::<F::Config>()
            .ok_or_else(|| Error::InvalidConfiguration("Invalid configuration type".into()))?;
        
        let solver = self.factory.create(typed_config.clone())?;
        Ok(Box::new(DynamicSolverWrapper::new(solver)))
    }
    
    fn name(&self) -> &str {
        self.factory.name()
    }
    
    fn can_handle(&self, problem_type: &str) -> bool {
        self.factory.can_handle(problem_type)
    }
}

/// Complete factory registry with proper factory pattern implementation
pub struct SolverFactoryRegistry<T: RealField> {
    factories: HashMap<String, Arc<dyn DynamicFactory<T>>>,
    factory_metadata: HashMap<String, FactoryMetadata>,
}

/// Metadata for registered factories
#[derive(Debug, Clone)]
pub struct FactoryMetadata {
    pub name: String,
    pub factory_type: String,
    pub version: String,
    pub capabilities: Vec<String>,
}

impl<T: RealField> Default for SolverFactoryRegistry<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField> SolverFactoryRegistry<T> {
    /// Create new factory registry
    pub fn new() -> Self {
        Self {
            factories: HashMap::new(),
            factory_metadata: HashMap::new(),
        }
    }

    /// Register a factory with metadata
    pub fn register_factory(
        &mut self, 
        name: String, 
        factory: Arc<dyn DynamicFactory<T>>,
        metadata: FactoryMetadata
    ) -> Result<()> {
        if self.factories.contains_key(&name) {
            return Err(Error::InvalidInput(format!("Factory '{}' already registered", name)));
        }

        self.factories.insert(name.clone(), factory);
        self.factory_metadata.insert(name, metadata);
        Ok(())
    }

    /// Get a factory by name
    pub fn get_factory(&self, name: &str) -> Option<&dyn DynamicFactory<T>> {
        self.factories.get(name).map(|f| f.as_ref())
    }

    /// List available factories
    pub fn list_factories(&self) -> Vec<&str> {
        self.factories.keys().map(|s| s.as_str()).collect()
    }

    /// Get factory metadata
    pub fn get_factory_metadata(&self, name: &str) -> Option<&FactoryMetadata> {
        self.factory_metadata.get(name)
    }
    
    /// Create a solver using the specified factory with a configuration
    pub fn create_solver<C: Any>(&self, factory_name: &str, config: C) -> Result<Box<dyn DynamicSolver<T>>> {
        self.factories
            .get(factory_name)
            .ok_or_else(|| Error::InvalidInput(format!("Factory '{}' not found", factory_name)))?
            .create_solver(&config)
    }
    
    /// Get all factories of a specific type
    pub fn get_factories_by_type(&self, factory_type: &str) -> Vec<&str> {
        self.factory_metadata
            .iter()
            .filter(|(_, metadata)| metadata.factory_type == factory_type)
            .map(|(name, _)| name.as_str())
            .collect()
    }
    
    /// Check if a factory supports a specific capability
    pub fn factory_supports(&self, factory_name: &str, capability: &str) -> bool {
        self.factory_metadata
            .get(factory_name)
            .map(|metadata| metadata.capabilities.contains(&capability.to_string()))
            .unwrap_or(false)
    }
}

/// Builder pattern for complex object creation (Creator principle)
pub trait Builder<T> {
    /// Build the final object
    fn build(self) -> Result<T>;

    /// Validate the builder state before building
    fn validate(&self) -> Result<()>;
}

/// Configuration builder implementing Builder pattern
pub struct ConfigurationBuilder<T: RealField> {
    tolerance: Option<T>,
    max_iterations: Option<usize>,
    verbosity: Option<u8>,
    parallel: Option<bool>,
    custom_params: HashMap<String, String>,
}

impl<T: RealField> Default for ConfigurationBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField> ConfigurationBuilder<T> {
    /// Create new configuration builder
    pub fn new() -> Self {
        Self {
            tolerance: None,
            max_iterations: None,
            verbosity: None,
            parallel: None,
            custom_params: HashMap::new(),
        }
    }

    /// Set tolerance (fluent interface)
    pub fn tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = Some(tolerance);
        self
    }

    /// Set maximum iterations
    pub fn max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = Some(max_iterations);
        self
    }

    /// Set verbosity level
    pub fn verbosity(mut self, verbosity: u8) -> Self {
        self.verbosity = Some(verbosity);
        self
    }

    /// Enable/disable parallel execution
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.parallel = Some(parallel);
        self
    }

    /// Add custom parameter
    pub fn custom_param(mut self, key: String, value: String) -> Self {
        self.custom_params.insert(key, value);
        self
    }
}

impl<T: RealField + num_traits::FromPrimitive> Builder<crate::SolverConfig<T>> for ConfigurationBuilder<T> {
    fn build(self) -> Result<crate::SolverConfig<T>> {
        self.validate()?;

        Ok(crate::SolverConfig::builder()
            .tolerance(self.tolerance.unwrap_or_else(|| T::from_f64(1e-6).unwrap()))
            .max_iterations(self.max_iterations.unwrap_or(1000))
            .relaxation_factor(T::one())
            .verbosity(self.verbosity.unwrap_or(1))
            .parallel(self.parallel.unwrap_or(true))
            .num_threads(None)
            .build())
    }

    fn validate(&self) -> Result<()> {
        if let Some(tolerance) = &self.tolerance {
            if *tolerance <= T::zero() {
                return Err(Error::InvalidInput("Tolerance must be positive".to_string()));
            }
        }

        if let Some(max_iterations) = self.max_iterations {
            if max_iterations == 0 {
                return Err(Error::InvalidInput("Max iterations must be positive".to_string()));
            }
        }

        if let Some(verbosity) = self.verbosity {
            if verbosity > 3 {
                return Err(Error::InvalidInput("Verbosity level must be 0-3".to_string()));
            }
        }

        Ok(())
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_configuration_builder() {
        let config = ConfigurationBuilder::<f64>::new()
            .tolerance(1e-8)
            .max_iterations(500)
            .verbosity(2)
            .parallel(true)
            .build()
            .unwrap();

        assert_eq!(config.tolerance(), 1e-8);
        assert_eq!(config.max_iterations(), 500);
        assert_eq!(config.verbosity(), 2);
        assert!(config.parallel());
    }


}
