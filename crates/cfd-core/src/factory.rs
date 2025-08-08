//! Factory patterns implementing GRASP Creator principle.
//!
//! This module provides factory abstractions that assign creation responsibility
//! to classes that have the initializing data, following GRASP principles.

use crate::{Error, Result, Solver, SolverConfiguration};
use nalgebra::RealField;
use std::collections::HashMap;

/// Abstract factory trait following GRASP Creator principle
/// Assigns creation responsibility to classes with initializing data
/// Uses type erasure to avoid dyn compatibility issues
pub trait AbstractSolverFactory<T: RealField>: Send + Sync {
    /// Create a solver with the given configuration (simplified)
    fn create_solver_simple(&self, name: &str) -> Result<String>;

    /// Get factory name for identification
    fn name(&self) -> &str;

    /// Get supported problem types
    fn supported_types(&self) -> Vec<&str>;
}

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

/// Simplified factory registry avoiding trait object issues
pub struct SolverFactoryRegistry {
    factory_names: HashMap<String, String>,
}

impl Default for SolverFactoryRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl SolverFactoryRegistry {
    /// Create new factory registry
    pub fn new() -> Self {
        Self {
            factory_names: HashMap::new(),
        }
    }

    /// Register a factory name (simplified approach)
    pub fn register(&mut self, name: String, factory_type: String) -> Result<()> {
        if self.factory_names.contains_key(&name) {
            return Err(Error::InvalidInput(format!("Factory '{}' already registered", name)));
        }

        self.factory_names.insert(name, factory_type);
        Ok(())
    }

    /// List available factories
    pub fn list_factories(&self) -> Vec<&str> {
        self.factory_names.keys().map(|s| s.as_str()).collect()
    }

    /// Get factory type for a given name
    pub fn get_factory_type(&self, name: &str) -> Option<&str> {
        self.factory_names.get(name).map(|s| s.as_str())
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

/// Resource manager implementing ACID principles for data operations
pub struct ResourceManager<T> {
    resources: HashMap<String, T>,
    transaction_log: Vec<TransactionEntry>,
}

/// Transaction entry for audit logging
#[derive(Debug, Clone)]
pub struct TransactionEntry {
    /// Operation type (ADD, REMOVE, etc.)
    pub operation: String,
    /// Resource identifier
    pub resource_id: String,
    /// Timestamp of the operation
    pub timestamp: std::time::SystemTime,
}

impl<T: Clone> ResourceManager<T> {
    /// Create new resource manager
    pub fn new() -> Self {
        Self {
            resources: HashMap::new(),
            transaction_log: Vec::new(),
        }
    }

    /// Atomic operation to add resource
    pub fn add_resource(&mut self, id: String, resource: T) -> Result<()> {
        // Atomicity: Either succeeds completely or fails completely
        if self.resources.contains_key(&id) {
            return Err(Error::InvalidInput(format!("Resource '{}' already exists", id)));
        }

        // Log transaction for durability
        self.transaction_log.push(TransactionEntry {
            operation: "ADD".to_string(),
            resource_id: id.clone(),
            timestamp: std::time::SystemTime::now(),
        });

        self.resources.insert(id, resource);
        Ok(())
    }

    /// Atomic operation to remove resource
    pub fn remove_resource(&mut self, id: &str) -> Result<T> {
        let resource = self.resources.remove(id)
            .ok_or_else(|| Error::InvalidInput(format!("Resource '{}' not found", id)))?;

        // Log transaction
        self.transaction_log.push(TransactionEntry {
            operation: "REMOVE".to_string(),
            resource_id: id.to_string(),
            timestamp: std::time::SystemTime::now(),
        });

        Ok(resource)
    }

    /// Get resource (Isolation: read-only access)
    pub fn get_resource(&self, id: &str) -> Option<&T> {
        self.resources.get(id)
    }

    /// Get transaction log for auditing (Durability)
    pub fn get_transaction_log(&self) -> &[TransactionEntry] {
        &self.transaction_log
    }

    /// Consistency check
    pub fn validate_consistency(&self) -> Result<()> {
        // Implement consistency checks based on business rules
        // For now, just check that all logged resources exist
        for entry in &self.transaction_log {
            if entry.operation == "ADD" && !self.resources.contains_key(&entry.resource_id) {
                return Err(Error::InvalidState(
                    format!("Inconsistency: Resource '{}' in log but not in storage", entry.resource_id)
                ));
            }
        }
        Ok(())
    }
}

impl<T: Clone> Default for ResourceManager<T> {
    fn default() -> Self {
        Self::new()
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

    #[test]
    fn test_resource_manager() {
        let mut manager = ResourceManager::new();
        
        // Test atomic add
        manager.add_resource("test".to_string(), 42).unwrap();
        assert_eq!(manager.get_resource("test"), Some(&42));
        
        // Test consistency
        manager.validate_consistency().unwrap();
        
        // Test atomic remove
        let removed = manager.remove_resource("test").unwrap();
        assert_eq!(removed, 42);
        assert_eq!(manager.get_resource("test"), None);
    }
}
