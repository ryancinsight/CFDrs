//! Simulation orchestration implementing Separation of Concerns (SOC) principle.
//!
//! This module separates the concerns of simulation execution, coordination,
//! and resource management, following CLEAN architecture principles.

use crate::{Error, Result};
use nalgebra::RealField;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Execution coordinator following Single Responsibility Principle
/// Responsible only for coordinating execution flow
pub struct ExecutionCoordinator<T: RealField> {
    /// Execution context
    context: ExecutionContext<T>,
}

/// Resource manager following Single Responsibility Principle
/// Responsible only for managing computational resources
pub struct ResourceManager {
    /// Performance metrics
    metrics: Arc<Mutex<PerformanceMetrics>>,
    /// Resource limits
    memory_limit: Option<usize>,
    /// Thread pool size
    thread_pool_size: Option<usize>,
}

/// Solver registry following Single Responsibility Principle
/// Responsible only for managing solver registration and lookup
pub struct SolverRegistry {
    /// Registered solver names and types
    solver_types: HashMap<String, String>,
    /// Solver capabilities
    capabilities: HashMap<String, Vec<String>>,
}

/// Simulation orchestrator implementing Controller pattern from GRASP
/// Coordinates between different components following Low Coupling principle
pub struct SimulationOrchestrator<T: RealField> {
    /// Execution coordinator
    coordinator: ExecutionCoordinator<T>,
    /// Resource manager
    resource_manager: ResourceManager,
    /// Solver registry
    solver_registry: SolverRegistry,
}

/// Execution context for simulation runs
#[derive(Debug, Clone)]
pub struct ExecutionContext<T: RealField> {
    /// Simulation name
    pub name: String,
    /// Start time
    pub start_time: Option<Instant>,
    /// Maximum execution time
    pub max_duration: Option<Duration>,
    /// Resource limits
    pub resource_limits: ResourceLimits,
    /// Execution mode
    pub mode: ExecutionMode,
    /// Custom parameters
    pub parameters: HashMap<String, T>,
}

/// Resource limits for simulation execution
#[derive(Debug, Clone)]
pub struct ResourceLimits {
    /// Maximum memory usage in bytes
    pub max_memory: Option<usize>,
    /// Maximum CPU time
    pub max_cpu_time: Option<Duration>,
    /// Maximum number of iterations
    pub max_iterations: Option<usize>,
}

/// Execution mode for simulations
#[derive(Debug, Clone, PartialEq)]
pub enum ExecutionMode {
    /// Sequential execution
    Sequential,
    /// Parallel execution
    Parallel,
    /// Adaptive execution (switches based on problem size)
    Adaptive,
}

/// Performance metrics collection
#[derive(Debug, Default, Clone)]
pub struct PerformanceMetrics {
    /// Total execution time
    pub total_time: Duration,
    /// Solver-specific timings
    pub solver_times: HashMap<String, Duration>,
    /// Memory usage statistics
    pub memory_usage: MemoryUsage,
    /// Iteration counts
    pub iterations: HashMap<String, usize>,
    /// Convergence information
    pub convergence_info: Vec<ConvergenceRecord>,
}

/// Memory usage statistics
#[derive(Debug, Default, Clone)]
pub struct MemoryUsage {
    /// Peak memory usage
    pub peak_bytes: usize,
    /// Current memory usage
    pub current_bytes: usize,
    /// Memory allocations count
    pub allocations: usize,
}

/// Convergence record for analysis
#[derive(Debug, Clone)]
pub struct ConvergenceRecord {
    /// Solver name
    pub solver_name: String,
    /// Iteration number
    pub iteration: usize,
    /// Residual value
    pub residual: f64,
    /// Timestamp
    pub timestamp: Instant,
}

impl<T: RealField> Default for ExecutionContext<T> {
    fn default() -> Self {
        Self {
            name: "default".to_string(),
            start_time: None,
            max_duration: None,
            resource_limits: ResourceLimits::default(),
            mode: ExecutionMode::Sequential,
            parameters: HashMap::new(),
        }
    }
}

impl Default for ResourceLimits {
    fn default() -> Self {
        Self {
            max_memory: None,
            max_cpu_time: None,
            max_iterations: Some(10000),
        }
    }
}

impl<T: RealField> ExecutionCoordinator<T> {
    /// Create new execution coordinator
    pub fn new(context: ExecutionContext<T>) -> Self {
        Self { context }
    }

    /// Get execution context
    pub fn context(&self) -> &ExecutionContext<T> {
        &self.context
    }

    /// Set execution context
    pub fn set_context(&mut self, context: ExecutionContext<T>) {
        self.context = context;
    }
}

impl ResourceManager {
    /// Create new resource manager
    pub fn new() -> Self {
        Self {
            metrics: Arc::new(Mutex::new(PerformanceMetrics::default())),
            memory_limit: None,
            thread_pool_size: None,
        }
    }

    /// Set memory limit
    pub fn set_memory_limit(&mut self, limit: usize) {
        self.memory_limit = Some(limit);
    }

    /// Set thread pool size
    pub fn set_thread_pool_size(&mut self, size: usize) {
        self.thread_pool_size = Some(size);
    }

    /// Get performance metrics
    pub fn metrics(&self) -> Arc<Mutex<PerformanceMetrics>> {
        Arc::clone(&self.metrics)
    }
}

impl SolverRegistry {
    /// Create new solver registry
    pub fn new() -> Self {
        Self {
            solver_types: HashMap::new(),
            capabilities: HashMap::new(),
        }
    }

    /// Register a solver type
    pub fn register_solver(&mut self, name: String, solver_type: String, capabilities: Vec<String>) -> Result<()> {
        if self.solver_types.contains_key(&name) {
            return Err(Error::InvalidInput(format!("Solver '{}' already registered", name)));
        }

        self.solver_types.insert(name.clone(), solver_type);
        self.capabilities.insert(name, capabilities);
        Ok(())
    }

    /// Get solver type
    pub fn get_solver_type(&self, name: &str) -> Option<&str> {
        self.solver_types.get(name).map(|s| s.as_str())
    }

    /// Get solver capabilities
    pub fn get_capabilities(&self, name: &str) -> Option<&Vec<String>> {
        self.capabilities.get(name)
    }
}

impl<T: RealField> SimulationOrchestrator<T> {
    /// Create new orchestrator following Dependency Injection
    pub fn new() -> Self {
        Self {
            coordinator: ExecutionCoordinator::new(ExecutionContext::default()),
            resource_manager: ResourceManager::new(),
            solver_registry: SolverRegistry::new(),
        }
    }

    /// Create orchestrator with custom components (Dependency Injection)
    pub fn with_components(
        coordinator: ExecutionCoordinator<T>,
        resource_manager: ResourceManager,
        solver_registry: SolverRegistry,
    ) -> Self {
        Self {
            coordinator,
            resource_manager,
            solver_registry,
        }
    }

    /// Register a solver type with capabilities
    pub fn register_solver_type(&mut self, name: String, solver_type: String, capabilities: Vec<String>) -> Result<()> {
        self.solver_registry.register_solver(name, solver_type, capabilities)
    }

    /// Set execution context
    pub fn set_context(&mut self, context: ExecutionContext<T>) {
        self.coordinator.set_context(context);
    }

    /// Get resource manager
    pub fn resource_manager(&self) -> &ResourceManager {
        &self.resource_manager
    }

    /// Get solver registry
    pub fn solver_registry(&self) -> &SolverRegistry {
        &self.solver_registry
    }

    /// Execute simulation with proper separation of concerns
    pub fn execute(&mut self, solver_name: &str) -> Result<String> {
        // Preparation phase (separated concern)
        self.prepare_execution()?;

        // Execution phase (separated concern)
        let result = self.execute_solver_by_name(solver_name)?;

        // Cleanup phase (separated concern)
        self.finalize_execution()?;

        Ok(result)
    }

    /// Prepare execution environment (SOC: preparation logic)
    fn prepare_execution(&mut self) -> Result<()> {
        self.coordinator.context.start_time = Some(Instant::now());

        // Validate resource limits
        if let Some(max_memory) = self.coordinator.context.resource_limits.max_memory {
            if max_memory == 0 {
                return Err(Error::InvalidInput("Max memory must be positive".to_string()));
            }
        }

        // Initialize metrics
        let metrics_arc = self.resource_manager.metrics();
        let mut metrics = metrics_arc.lock().unwrap();
        *metrics = PerformanceMetrics::default();

        tracing::info!("Simulation '{}' prepared for execution", self.coordinator.context.name);
        Ok(())
    }

    /// Execute solver with monitoring (SOC: execution logic)
    fn execute_solver_by_name(&mut self, solver_name: &str) -> Result<String> {
        let solver_type = self.solver_registry.get_solver_type(solver_name)
            .ok_or_else(|| Error::InvalidInput(format!("Solver '{}' not found", solver_name)))?;

        let start_time = Instant::now();

        // Check resource limits before execution
        self.check_resource_limits()?;

        // For now, we simulate execution since actual solver dispatch would require
        // the concrete solver implementations to be instantiated
        let result = format!("Executed {} solver of type {}", solver_name, solver_type);

        // Record metrics
        let execution_time = start_time.elapsed();
        let metrics_arc = self.resource_manager.metrics();
        let mut metrics = metrics_arc.lock().unwrap();
        metrics.solver_times.insert(solver_name.to_string(), execution_time);
        metrics.total_time += execution_time;

        tracing::info!("Solver '{}' completed in {:?}", solver_name, execution_time);
        Ok(result)
    }

    /// Finalize execution (SOC: cleanup logic)
    fn finalize_execution(&mut self) -> Result<()> {
        let total_time = self.coordinator.context.start_time
            .map(|start| start.elapsed())
            .unwrap_or_default();

        let metrics_arc = self.resource_manager.metrics();
        let mut metrics = metrics_arc.lock().unwrap();
        metrics.total_time = total_time;

        tracing::info!("Simulation '{}' completed in {:?}", self.coordinator.context.name, total_time);
        Ok(())
    }

    /// Check resource limits (SOC: resource management)
    fn check_resource_limits(&self) -> Result<()> {
        if let Some(max_duration) = self.coordinator.context.max_duration {
            if let Some(start_time) = self.coordinator.context.start_time {
                if start_time.elapsed() > max_duration {
                    return Err(Error::TimeoutError("Execution time limit exceeded".to_string()));
                }
            }
        }

        // Additional resource checks would go here
        Ok(())
    }

    /// Get registered solver types
    pub fn get_solver_types(&self) -> Vec<(&str, &str)> {
        self.solver_registry.solver_types.iter()
            .map(|(name, solver_type)| (name.as_str(), solver_type.as_str()))
            .collect()
    }

    /// Get performance metrics
    pub fn get_metrics(&self) -> PerformanceMetrics {
        let metrics_arc = self.resource_manager.metrics();
        let result = metrics_arc.lock().unwrap().clone();
        result
    }

    /// Reset orchestrator state
    pub fn reset(&mut self) {
        self.coordinator.set_context(ExecutionContext::default());
        let metrics_arc = self.resource_manager.metrics();
        let mut metrics = metrics_arc.lock().unwrap();
        *metrics = PerformanceMetrics::default();
    }
}

impl<T: RealField> Default for SimulationOrchestrator<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for ResourceManager {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for SolverRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Builder for execution context (Builder pattern from GRASP)
pub struct ExecutionContextBuilder<T: RealField> {
    context: ExecutionContext<T>,
}

impl<T: RealField> ExecutionContextBuilder<T> {
    /// Create new builder
    pub fn new(name: String) -> Self {
        Self {
            context: ExecutionContext {
                name,
                ..Default::default()
            },
        }
    }

    /// Set maximum duration
    pub fn max_duration(mut self, duration: Duration) -> Self {
        self.context.max_duration = Some(duration);
        self
    }

    /// Set execution mode
    pub fn mode(mut self, mode: ExecutionMode) -> Self {
        self.context.mode = mode;
        self
    }

    /// Set resource limits
    pub fn resource_limits(mut self, limits: ResourceLimits) -> Self {
        self.context.resource_limits = limits;
        self
    }

    /// Add parameter
    pub fn parameter(mut self, key: String, value: T) -> Self {
        self.context.parameters.insert(key, value);
        self
    }

    /// Build the context
    pub fn build(self) -> ExecutionContext<T> {
        self.context
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_execution_context_builder() {
        let context = ExecutionContextBuilder::<f64>::new("test".to_string())
            .max_duration(Duration::from_secs(60))
            .mode(ExecutionMode::Parallel)
            .parameter("tolerance".to_string(), 1e-6)
            .build();

        assert_eq!(context.name, "test");
        assert_eq!(context.max_duration, Some(Duration::from_secs(60)));
        assert_eq!(context.mode, ExecutionMode::Parallel);
        assert_eq!(context.parameters.get("tolerance"), Some(&1e-6));
    }

    #[test]
    fn test_orchestrator_creation() {
        let orchestrator = SimulationOrchestrator::<f64>::new();
        assert_eq!(orchestrator.solver_registry.solver_types.len(), 0);
        assert_eq!(orchestrator.coordinator.context.name, "default");
    }
}
