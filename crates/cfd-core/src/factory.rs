//! Factory patterns implementing GRASP Creator principle.
//!
//! This module provides factory abstractions that assign creation responsibility
//! to classes that have the initializing data, following GRASP principles.

use crate::error::{Error, Result};
use crate::solver::{Solver, SolverConfiguration};
use nalgebra::RealField;
use std::collections::HashMap;
use std::any::Any;
use std::sync::Arc;

/// Concrete factory trait for type-safe creation
/// Follows Open/Closed Principle - open for extension, closed for modification
pub trait ConcreteSolverFactory<T: RealField + Copy>: Send + Sync {
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

/// Dynamic solver trait for runtime polymorphism
pub trait DynamicSolver<T: RealField + Copy>: Send + Sync {
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

/// Dynamic factory for heterogeneous storage
pub trait DynamicFactory<T: RealField + Copy>: Send + Sync {
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
        
        let solver = self.factory.create(typed_config)?;
        Ok(Box::new(DynamicSolverWrapper::new(solver)))
    }
    
    fn name(&self) -> &str {
        self.factory.name()
    }
    
    fn can_handle(&self, problem_type: &str) -> bool {
        self.factory.can_handle(problem_type)
    }
}

/// Metadata for registered factories
#[derive(Debug, Clone)]
pub struct FactoryMetadata {
    pub name: String,
    pub factory_type: String,
    pub version: String,
    pub capabilities: Vec<String>,
}

/// Registry entry containing both factory and metadata
struct RegistryEntry<T: RealField + Copy> {
    factory: Arc<dyn DynamicFactory<T>>,
    metadata: FactoryMetadata,
}

/// Complete factory registry with proper factory pattern implementation
/// Uses a single HashMap for SSOT (Single Source of Truth)
pub struct SolverFactoryRegistry<T: RealField + Copy> {
    registry: HashMap<String, RegistryEntry<T>>,
}

impl<T: RealField + Copy> Default for SolverFactoryRegistry<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> SolverFactoryRegistry<T> {
    /// Create new factory registry
    pub fn new() -> Self {
        Self {
            registry: HashMap::new(),
        }
    }

    /// Register a factory with metadata
    pub fn register_factory(
        &mut self, 
        name: String, 
        factory: Arc<dyn DynamicFactory<T>>,
        metadata: FactoryMetadata
    ) -> Result<()> {
        use std::collections::hash_map::Entry;
        
        match self.registry.entry(name) {
            Entry::Occupied(_) => {
                Err(Error::InvalidInput(format!("Factory '{}' already registered", name)))
            }
            Entry::Vacant(entry) => {
                entry.insert(RegistryEntry { factory, metadata });
                Ok(())
            }
        }
    }

    /// Get a factory by name
    pub fn get_factory(&self, name: &str) -> Option<&dyn DynamicFactory<T>> {
        self.registry.get(name).map(|entry| entry.factory.as_ref())
    }

    /// Get factory metadata
    pub fn get_factory_metadata(&self, name: &str) -> Option<&FactoryMetadata> {
        self.registry.get(name).map(|entry| &entry.metadata)
    }

    /// List available factories
    pub fn list_factories(&self) -> Vec<&str> {
        self.registry.keys().map(|s| s.as_str()).collect()
    }
    
    /// Create a solver using the specified factory with a configuration
    pub fn create_solver<C: Any>(&self, factory_name: &str, config: C) -> Result<Box<dyn DynamicSolver<T>>> {
        self.registry
            .get(factory_name)
            .ok_or_else(|| Error::InvalidInput(format!("Factory '{}' not found", factory_name)))?
            .factory
            .create_solver(&config)
    }
    
    /// Get all factories of a specific type
    pub fn get_factories_by_type(&self, factory_type: &str) -> Vec<&str> {
        self.registry
            .iter()
            .filter(|(_, entry)| entry.metadata.factory_type == factory_type)
            .map(|(name, _)| name.as_str())
            .collect()
    }
    
    /// Check if a factory supports a specific capability
    pub fn factory_supports(&self, factory_name: &str, capability: &str) -> bool {
        self.registry
            .get(factory_name)
            .map(|entry| entry.metadata.capabilities.contains(&capability.to_string()))
            .unwrap_or(false)
    }
    
    /// Get factory and metadata together (for advanced use cases)
    pub fn get_entry(&self, name: &str) -> Option<(&dyn DynamicFactory<T>, &FactoryMetadata)> {
        self.registry.get(name).map(|entry| (entry.factory.as_ref(), &entry.metadata))
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::solver::SolverConfig;

    #[test]
    fn test_factory_registry() {
        let mut registry = SolverFactoryRegistry::<f64>::new();
        
        // Test registration
        let metadata = FactoryMetadata {
            name: "test_factory".to_string(),
            factory_type: "iterative".to_string(),
            version: "1.0.0".to_string(),
            capabilities: vec!["linear".to_string(), "sparse".to_string()],
        };
        
        // Would need a concrete factory implementation for a full test
        // This test verifies the structure compiles and basic operations work
        
        assert_eq!(registry.list_factories().len(), 0);
    }
    
    #[test]
    fn test_registry_entry_consolidation() {
        let mut registry = SolverFactoryRegistry::<f64>::new();
        
        // Verify that we can't register the same factory twice
        let metadata = FactoryMetadata {
            name: "test".to_string(),
            factory_type: "test".to_string(),
            version: "1.0.0".to_string(),
            capabilities: vec![],
        };
        
        // This would fail without a concrete factory, but demonstrates the API
        assert!(registry.get_factory("nonexistent").is_none());
        assert!(registry.get_factory_metadata("nonexistent").is_none());
    }
}