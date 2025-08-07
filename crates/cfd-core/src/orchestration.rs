//! Simulation orchestration implementing Separation of Concerns (SOC) principle.
//!
//! This module separates the concerns of simulation execution, coordination,
//! and resource management, following CLEAN architecture principles.

use crate::{Error, Result};
use nalgebra::RealField;
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Simulation orchestrator implementing Controller pattern from GRASP
pub struct SimulationOrchestrator<T: RealField> {
    /// Registered solver names and types
    solver_registry: HashMap<String, String>,
    /// Execution context
    context: ExecutionContext<T>,
    /// Performance metrics
    metrics: Arc<Mutex<PerformanceMetrics>>,
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

impl<T: RealField> SimulationOrchestrator<T> {
    /// Create new orchestrator
    pub fn new() -> Self {
        Self {
            solver_registry: HashMap::new(),
            context: ExecutionContext::default(),
            metrics: Arc::new(Mutex::new(PerformanceMetrics::default())),
        }
    }

    /// Register a solver type (simplified approach)
    pub fn register_solver_type(&mut self, name: String, solver_type: String) -> Result<()> {
        if self.solver_registry.contains_key(&name) {
            return Err(Error::InvalidInput(format!("Solver '{}' already registered", name)));
        }

        self.solver_registry.insert(name, solver_type);
        Ok(())
    }

    /// Set execution context
    pub fn set_context(&mut self, context: ExecutionContext<T>) {
        self.context = context;
    }

    /// Execute simulation with proper separation of concerns (simplified)
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
        self.context.start_time = Some(Instant::now());

        // Validate resource limits
        if let Some(max_memory) = self.context.resource_limits.max_memory {
            if max_memory == 0 {
                return Err(Error::InvalidInput("Max memory must be positive".to_string()));
            }
        }

        // Initialize metrics
        let mut metrics = self.metrics.lock().unwrap();
        *metrics = PerformanceMetrics::default();

        tracing::info!("Simulation '{}' prepared for execution", self.context.name);
        Ok(())
    }

    /// Execute solver with monitoring (SOC: execution logic, simplified)
    fn execute_solver_by_name(&mut self, solver_name: &str) -> Result<String> {
        let solver_type = self.solver_registry.get(solver_name)
            .ok_or_else(|| Error::InvalidInput(format!("Solver '{}' not found", solver_name)))?;

        let start_time = Instant::now();

        // Check resource limits before execution
        self.check_resource_limits()?;

        // Simplified execution - in practice, this would dispatch to actual solver
        let result = format!("Executed {} solver of type {}", solver_name, solver_type);

        // Record metrics
        let execution_time = start_time.elapsed();
        let mut metrics = self.metrics.lock().unwrap();
        metrics.solver_times.insert(solver_name.to_string(), execution_time);
        metrics.total_time += execution_time;

        tracing::info!("Solver '{}' completed in {:?}", solver_name, execution_time);
        Ok(result)
    }

    /// Finalize execution (SOC: cleanup logic)
    fn finalize_execution(&mut self) -> Result<()> {
        let total_time = self.context.start_time
            .map(|start| start.elapsed())
            .unwrap_or_default();

        let mut metrics = self.metrics.lock().unwrap();
        metrics.total_time = total_time;

        tracing::info!("Simulation '{}' completed in {:?}", self.context.name, total_time);
        Ok(())
    }

    /// Check resource limits (SOC: resource management)
    fn check_resource_limits(&self) -> Result<()> {
        if let Some(max_duration) = self.context.max_duration {
            if let Some(start_time) = self.context.start_time {
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
        self.solver_registry.iter()
            .map(|(name, solver_type)| (name.as_str(), solver_type.as_str()))
            .collect()
    }

    /// Get performance metrics
    pub fn get_metrics(&self) -> PerformanceMetrics {
        self.metrics.lock().unwrap().clone()
    }

    /// Reset orchestrator state
    pub fn reset(&mut self) {
        self.context = ExecutionContext::default();
        let mut metrics = self.metrics.lock().unwrap();
        *metrics = PerformanceMetrics::default();
    }
}

impl<T: RealField> Default for SimulationOrchestrator<T> {
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
        assert_eq!(orchestrator.solver_registry.len(), 0);
        assert_eq!(orchestrator.context.name, "default");
    }
}
