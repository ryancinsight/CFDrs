//! Plugin system for extensible CFD solvers.
//!
//! This module provides a flexible plugin architecture following SOLID principles,
//! enabling runtime extensibility while maintaining type safety and zero-cost abstractions.

use crate::{Error, Result};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

/// Trait for simulation plugins following the Interface Segregation Principle
pub trait SimulationPlugin: Send + Sync + 'static {
    /// Configuration type for this plugin
    type Config: Send + Sync + 'static;
    /// State type maintained by this plugin
    type State: Send + Sync + 'static;
    /// Output type produced by this plugin
    type Output: Send + Sync + 'static;

    /// Initialize the plugin with given configuration
    fn initialize(&self, config: Self::Config) -> Result<Self::State>;

    /// Perform one simulation step
    fn step(&self, state: &mut Self::State, dt: f64) -> Result<()>;

    /// Generate output from current state
    fn output(&self, state: &Self::State) -> Self::Output;

    /// Get plugin metadata
    fn metadata(&self) -> PluginMetadata {
        PluginMetadata::default()
    }
}

/// Plugin lifecycle hooks for simulation control
pub trait PluginLifecycle: SimulationPlugin {
    /// Called before simulation starts
    fn on_start(&self, _state: &mut Self::State) -> Result<()> {
        Ok(())
    }

    /// Called after simulation ends
    fn on_stop(&self, _state: &Self::State) -> Result<()> {
        Ok(())
    }

    /// Called on simulation reset
    fn on_reset(&self, state: &mut Self::State, config: Self::Config) -> Result<()> {
        *state = self.initialize(config)?;
        Ok(())
    }
}

/// Plugin capability interface following Interface Segregation Principle
/// Allows fine-grained capability queries without forcing implementation
pub trait PluginCapabilities {
    /// Check if plugin supports parallel execution
    fn supports_parallel(&self) -> bool {
        false
    }
    
    /// Check if plugin supports adaptive time stepping
    fn supports_adaptive_timestep(&self) -> bool {
        false
    }
    
    /// Check if plugin supports state persistence
    fn supports_persistence(&self) -> bool {
        false
    }
    
    /// Get supported input/output formats
    fn supported_formats(&self) -> Vec<&'static str> {
        vec![]
    }
}

/// Plugin validation interface for ensuring physics correctness
pub trait PluginValidation: SimulationPlugin {
    /// Validate plugin state for physics consistency
    fn validate_state(&self, state: &Self::State) -> Result<()> {
        let _ = state; // Avoid unused parameter warning
        Ok(())
    }
    
    /// Validate configuration parameters
    fn validate_config(&self, config: &Self::Config) -> Result<()> {
        let _ = config; // Avoid unused parameter warning
        Ok(())
    }
    
    /// Get conservation properties that must be maintained
    fn conservation_properties(&self) -> Vec<&'static str> {
        vec![]
    }
}

/// Plugin health monitoring for production environments
pub trait PluginHealth: SimulationPlugin {
    /// Check plugin health status
    fn health_check(&self) -> PluginHealthStatus {
        PluginHealthStatus::Healthy
    }
    
    /// Get plugin performance metrics
    fn performance_metrics(&self) -> PluginMetrics {
        PluginMetrics::default()
    }
    
    /// Check if plugin can handle increased load
    fn can_scale(&self, load_factor: f64) -> bool {
        let _ = load_factor;
        true
    }
}

/// Plugin health status
#[derive(Debug, Clone, PartialEq)]
pub enum PluginHealthStatus {
    /// Plugin is operating normally
    Healthy,
    /// Plugin has minor issues but is functional
    Degraded(String),
    /// Plugin has serious issues
    Unhealthy(String),
    /// Plugin is not responding
    Unresponsive,
}

/// Plugin performance metrics
#[derive(Debug, Clone, Default)]
pub struct PluginMetrics {
    /// CPU usage percentage (0.0 to 100.0)
    pub cpu_usage: f64,
    /// Memory usage in bytes
    pub memory_usage: u64,
    /// Operations per second
    pub operations_per_second: f64,
    /// Average response time in milliseconds
    pub avg_response_time_ms: f64,
    /// Error rate (0.0 to 1.0)
    pub error_rate: f64,
}



/// Generic plugin trait for type erasure
pub trait Plugin: Any + Send + Sync {
    /// Get plugin name
    fn name(&self) -> &str;

    /// Get plugin version
    fn version(&self) -> &str;

    /// Get plugin description
    fn description(&self) -> &str;

    /// Get type ID for downcasting
    fn type_id(&self) -> TypeId;

    /// As any for downcasting
    fn as_any(&self) -> &dyn Any;

    /// As any mut for downcasting
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Plugin metadata with dependency management
#[derive(Debug, Clone, Default)]
pub struct PluginMetadata {
    /// Plugin name
    pub name: String,
    /// Plugin version following semver
    pub version: String,
    /// Plugin description
    pub description: String,
    /// Plugin author
    pub author: Option<String>,
    /// Plugin license
    pub license: Option<String>,
    /// Plugin capabilities
    pub capabilities: Vec<String>,
    /// Required dependencies
    pub dependencies: Vec<String>,
    /// Minimum API version required
    pub min_api_version: String,
    /// Plugin priority for loading order
    pub priority: i32,
}

/// Plugin storage following Single Responsibility Principle
#[derive(Clone)]
pub struct PluginStorage {
    plugins: Arc<RwLock<HashMap<String, Arc<dyn Plugin>>>>,
    metadata: Arc<RwLock<HashMap<String, PluginMetadata>>>,
}

/// Dependency resolver following SRP and Open/Closed Principle
#[derive(Clone)]
pub struct DependencyResolver {
    dependencies: Arc<RwLock<HashMap<String, Vec<String>>>>,
    load_order: Arc<RwLock<Vec<String>>>,
}

/// Factory registry following Interface Segregation Principle
#[derive(Clone)]
pub struct FactoryRegistry {
    factories: Arc<RwLock<HashMap<TypeId, Arc<dyn Any + Send + Sync>>>>,
}

/// Main plugin registry coordinating storage, dependencies, and factories
#[derive(Clone)]
pub struct PluginRegistry {
    storage: PluginStorage,
    dependency_resolver: DependencyResolver,
    factory_registry: FactoryRegistry,
    monitoring: PluginMonitoring,
}

/// Plugin monitoring system for production environments
#[derive(Clone)]
pub struct PluginMonitoring {
    health_status: Arc<RwLock<HashMap<String, PluginHealthStatus>>>,
    metrics: Arc<RwLock<HashMap<String, PluginMetrics>>>,
    last_check: Arc<RwLock<HashMap<String, std::time::Instant>>>,
}



impl PluginMetadata {
    /// Create default metadata for a plugin
    pub fn default_for_plugin(name: &str) -> Self {
        Self {
            name: name.to_string(),
            version: "1.0.0".to_string(),
            description: format!("Plugin: {}", name),
            author: None,
            license: None,
            capabilities: Vec::new(),
            dependencies: Vec::new(),
            min_api_version: "1.0.0".to_string(),
            priority: 0,
        }
    }
}

impl PluginStorage {
    /// Create new plugin storage
    pub fn new() -> Self {
        Self {
            plugins: Arc::new(RwLock::new(HashMap::new())),
            metadata: Arc::new(RwLock::new(HashMap::new())),
        }
    }
    
    /// Store a plugin with its metadata
    pub fn store<P: Plugin + 'static>(&self, plugin: P, metadata: PluginMetadata) -> Result<()> {
        let name = plugin.name().to_string();
        
        // Store plugin
        {
            let mut plugins = self.plugins.write().map_err(|_| {
                Error::PluginError("Failed to acquire write lock on plugin storage".to_string())
            })?;
            
            if plugins.contains_key(&name) {
                return Err(Error::PluginError(format!(
                    "Plugin '{}' is already stored",
                    name
                )));
            }
            
            plugins.insert(name.clone(), Arc::new(plugin));
        }
        
        // Store metadata
        {
            let mut plugin_metadata = self.metadata.write().map_err(|_| {
                Error::PluginError("Failed to acquire write lock on metadata".to_string())
            })?;
            plugin_metadata.insert(name, metadata);
        }
        
        Ok(())
    }
    
    /// Get a plugin by name
    pub fn get(&self, name: &str) -> Result<Arc<dyn Plugin>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin storage".to_string())
        })?;

        plugins
            .get(name)
            .cloned()
            .ok_or_else(|| Error::PluginError(format!("Plugin '{}' not found", name)))
    }
    
    /// List all stored plugins
    pub fn list(&self) -> Result<Vec<String>> {
        let plugins = self.plugins.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on plugin storage".to_string())
        })?;

        Ok(plugins.keys().cloned().collect())
    }
    
    /// Get metadata for a plugin
    pub fn get_metadata(&self, name: &str) -> Result<PluginMetadata> {
        let metadata = self.metadata.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on metadata".to_string())
        })?;

        metadata
            .get(name)
            .cloned()
            .ok_or_else(|| Error::PluginError(format!("Metadata for '{}' not found", name)))
    }
}

impl DependencyResolver {
    /// Create new dependency resolver
    pub fn new() -> Self {
        Self {
            dependencies: Arc::new(RwLock::new(HashMap::new())),
            load_order: Arc::new(RwLock::new(Vec::new())),
        }
    }
    
    /// Add dependencies for a plugin
    pub fn add_dependencies(&self, plugin_name: String, deps: Vec<String>) -> Result<()> {
        {
            let mut dependencies = self.dependencies.write().map_err(|_| {
                Error::PluginError("Failed to acquire write lock on dependencies".to_string())
            })?;
            dependencies.insert(plugin_name.clone(), deps);
        }
        
        self.update_load_order()
    }
    
    /// Validate dependencies for a plugin
    pub fn validate_dependencies(&self, deps: &[String], storage: &PluginStorage) -> Result<()> {
        for dep in deps {
            storage.get(dep).map_err(|_| {
                Error::InvalidConfiguration(format!("Missing dependency: {}", dep))
            })?;
        }
        Ok(())
    }
    
    /// Get the computed load order
    pub fn get_load_order(&self) -> Result<Vec<String>> {
        let load_order = self.load_order.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on load order".to_string())
        })?;
        Ok(load_order.clone())
    }
    
    /// Update plugin load order based on dependencies (topological sort)
    fn update_load_order(&self) -> Result<()> {
        // Topological sort implementation
        let dependencies = self.dependencies.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on dependencies".to_string())
        })?;
        
        let mut visited = std::collections::HashSet::new();
        let mut result = Vec::new();
        
        fn visit(
            plugin: &str,
            deps: &HashMap<String, Vec<String>>,
            visited: &mut std::collections::HashSet<String>,
            result: &mut Vec<String>,
        ) -> Result<()> {
            if visited.contains(plugin) {
                return Ok(());
            }
            
            visited.insert(plugin.to_string());
            
            if let Some(plugin_deps) = deps.get(plugin) {
                for dep in plugin_deps {
                    visit(dep, deps, visited, result)?;
                }
            }
            
            result.push(plugin.to_string());
            Ok(())
        }
        
        for plugin in dependencies.keys() {
            visit(plugin, &dependencies, &mut visited, &mut result)?;
        }
        
        {
            let mut load_order = self.load_order.write().map_err(|_| {
                Error::PluginError("Failed to acquire write lock on load order".to_string())
            })?;
            *load_order = result;
        }
        
        Ok(())
    }
}

impl FactoryRegistry {
    /// Create new factory registry
    pub fn new() -> Self {
        Self {
            factories: Arc::new(RwLock::new(HashMap::new())),
        }
    }
    
    /// Register a solver factory
    pub fn register_factory<F: SolverFactory + 'static>(&self, factory: F) -> Result<()> {
        let type_id = TypeId::of::<F::Solver>();
        let mut factories = self.factories.write().map_err(|_| {
            Error::PluginError("Failed to acquire write lock on factory registry".to_string())
        })?;

        factories.insert(type_id, Arc::new(factory));
        Ok(())
    }
    
    /// Get a solver factory by type
    pub fn get_factory<S: SimulationPlugin + 'static>(&self) -> Result<Arc<dyn Any + Send + Sync>> {
        let type_id = TypeId::of::<S>();
        let factories = self.factories.read().map_err(|_| {
            Error::PluginError("Failed to acquire read lock on factory registry".to_string())
        })?;

        factories
            .get(&type_id)
            .cloned()
            .ok_or_else(|| Error::PluginError("Factory not found for solver type".to_string()))
    }
}

impl PluginMonitoring {
    /// Create new plugin monitoring system
    pub fn new() -> Self {
        Self {
            health_status: Arc::new(RwLock::new(HashMap::new())),
            metrics: Arc::new(RwLock::new(HashMap::new())),
            last_check: Arc::new(RwLock::new(HashMap::new())),
        }
    }
    
    /// Update plugin health status
    pub fn update_health(&self, plugin_name: &str, status: PluginHealthStatus) {
        let mut health = self.health_status.write().unwrap();
        let mut last_check = self.last_check.write().unwrap();
        
        health.insert(plugin_name.to_string(), status);
        last_check.insert(plugin_name.to_string(), std::time::Instant::now());
    }
    
    /// Update plugin metrics
    pub fn update_metrics(&self, plugin_name: &str, metrics: PluginMetrics) {
        let mut plugin_metrics = self.metrics.write().unwrap();
        plugin_metrics.insert(plugin_name.to_string(), metrics);
    }
    
    /// Get overall system health
    pub fn system_health(&self) -> Result<SystemHealthReport> {
        let health = self.health_status.read().unwrap();
        let metrics = self.metrics.read().unwrap();
        
        let total_plugins = health.len();
        let healthy_plugins = health.values()
            .filter(|&status| matches!(status, PluginHealthStatus::Healthy))
            .count();
        
        let average_cpu = if !metrics.is_empty() {
            metrics.values().map(|m| m.cpu_usage).sum::<f64>() / metrics.len() as f64
        } else {
            0.0
        };
        
        let total_memory = metrics.values().map(|m| m.memory_usage).sum::<u64>();
        
        Ok(SystemHealthReport {
            total_plugins,
            healthy_plugins,
            degraded_plugins: health.values()
                .filter(|status| matches!(status, PluginHealthStatus::Degraded(_)))
                .count(),
            unhealthy_plugins: health.values()
                .filter(|status| matches!(status, PluginHealthStatus::Unhealthy(_)))
                .count(),
            average_cpu_usage: average_cpu,
            total_memory_usage: total_memory,
            system_status: if healthy_plugins == total_plugins {
                SystemStatus::Healthy
            } else if healthy_plugins as f64 / total_plugins as f64 > 0.8 {
                SystemStatus::Degraded
            } else {
                SystemStatus::Critical
            },
        })
    }
}

/// System health report
#[derive(Debug, Clone)]
pub struct SystemHealthReport {
    /// Total number of plugins
    pub total_plugins: usize,
    /// Number of healthy plugins
    pub healthy_plugins: usize,
    /// Number of degraded plugins
    pub degraded_plugins: usize,
    /// Number of unhealthy plugins
    pub unhealthy_plugins: usize,
    /// Average CPU usage across all plugins
    pub average_cpu_usage: f64,
    /// Total memory usage across all plugins
    pub total_memory_usage: u64,
    /// Overall system status
    pub system_status: SystemStatus,
}

/// Overall system status
#[derive(Debug, Clone, PartialEq)]
pub enum SystemStatus {
    /// All plugins healthy
    Healthy,
    /// Some plugins degraded but system functional
    Degraded,
    /// Critical issues requiring attention
    Critical,
}

impl PluginRegistry {
    /// Create a new plugin registry using composition
    pub fn new() -> Self {
        Self {
            storage: PluginStorage::new(),
            dependency_resolver: DependencyResolver::new(),
            factory_registry: FactoryRegistry::new(),
            monitoring: PluginMonitoring::new(),
        }
    }

    /// Register a plugin
    pub fn register<P: Plugin + 'static>(&self, plugin: P) -> Result<()> {
        let name = plugin.name().to_string();
        let metadata = PluginMetadata::default_for_plugin(&name);
        self.register_with_metadata(plugin, metadata)
    }

    /// Register a plugin with metadata using composition pattern
    pub fn register_with_metadata<P: Plugin + 'static>(
        &self,
        plugin: P,
        metadata: PluginMetadata,
    ) -> Result<()> {
        // Validate dependencies using dependency resolver
        self.dependency_resolver.validate_dependencies(&metadata.dependencies, &self.storage)?;

        // Store plugin and metadata using storage component
        self.storage.store(plugin, metadata.clone())?;
        
        // Add dependencies using dependency resolver
        self.dependency_resolver.add_dependencies(
            metadata.name.clone(), 
            metadata.dependencies.clone()
        )?;

        Ok(())
    }

    /// Register a solver factory using factory registry
    pub fn register_factory<F: SolverFactory + 'static>(&self, factory: F) -> Result<()> {
        self.factory_registry.register_factory(factory)
    }

    /// Get a plugin by name using storage
    pub fn get(&self, name: &str) -> Result<Arc<dyn Plugin>> {
        self.storage.get(name)
    }

    /// Get a solver factory by type using factory registry
    pub fn get_factory<S: SimulationPlugin + 'static>(&self) -> Result<Arc<dyn Any + Send + Sync>> {
        self.factory_registry.get_factory::<S>()
    }

    /// List all registered plugins using storage
    pub fn list(&self) -> Result<Vec<String>> {
        self.storage.list()
    }
    
    /// Get plugin metadata using storage
    pub fn get_metadata(&self, name: &str) -> Result<PluginMetadata> {
        self.storage.get_metadata(name)
    }
    
    /// Get load order using dependency resolver
    pub fn get_load_order(&self) -> Result<Vec<String>> {
        self.dependency_resolver.get_load_order()
    }

    /// Filter plugins by capability using metadata
    pub fn filter_by_capability(&self, capability: &str) -> Result<Vec<String>> {
        let plugins = self.storage.list()?;
        
        Ok(plugins
            .into_iter()
            .filter_map(|name| {
                self.storage.get_metadata(&name).ok().and_then(|meta| {
                    if meta.capabilities.contains(&capability.to_string()) {
                        Some(name)
                    } else {
                        None
                    }
                })
            })
            .collect())
    }
    
    /// Update plugin health status
    pub fn update_plugin_health(&self, plugin_name: &str, status: PluginHealthStatus) {
        self.monitoring.update_health(plugin_name, status);
    }
    
    /// Update plugin metrics
    pub fn update_plugin_metrics(&self, plugin_name: &str, metrics: PluginMetrics) {
        self.monitoring.update_metrics(plugin_name, metrics);
    }
    
    /// Get overall system health report
    pub fn system_health(&self) -> Result<SystemHealthReport> {
        self.monitoring.system_health()
    }
    
    /// Perform comprehensive system health check across all plugins
    pub fn comprehensive_health_check(&self) -> Result<ComprehensiveHealthReport> {
        let plugins = self.storage.list()?;
        let mut plugin_reports = Vec::new();
        
        for plugin_name in plugins {
            if let Ok(plugin) = self.storage.get(&plugin_name) {
                let health_status = if plugin.name() == plugin_name {
                    PluginHealthStatus::Healthy
                } else {
                    PluginHealthStatus::Degraded("Name mismatch".to_string())
                };
                
                let metrics = PluginMetrics::default();
                
                plugin_reports.push(PluginHealthReport {
                    name: plugin_name.clone(),
                    status: health_status.clone(),
                    metrics: metrics.clone(),
                    dependencies_satisfied: true,
                    last_check: std::time::Instant::now(),
                });
                
                self.monitoring.update_health(&plugin_name, health_status);
                self.monitoring.update_metrics(&plugin_name, metrics);
            }
        }
        
        Ok(ComprehensiveHealthReport {
            system_health: self.monitoring.system_health()?,
            plugin_reports,
            check_timestamp: std::time::Instant::now(),
        })
    }
}

/// Individual plugin health report
#[derive(Debug, Clone)]
pub struct PluginHealthReport {
    /// Plugin name
    pub name: String,
    /// Current health status
    pub status: PluginHealthStatus,
    /// Performance metrics
    pub metrics: PluginMetrics,
    /// Whether all dependencies are satisfied
    pub dependencies_satisfied: bool,
    /// Timestamp of last health check
    pub last_check: std::time::Instant,
}

/// Comprehensive health report for the entire system
#[derive(Debug, Clone)]
pub struct ComprehensiveHealthReport {
    /// Overall system health summary
    pub system_health: SystemHealthReport,
    /// Individual plugin health reports
    pub plugin_reports: Vec<PluginHealthReport>,
    /// Timestamp when this comprehensive check was performed
    pub check_timestamp: std::time::Instant,
}

impl Default for PluginRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// Solver factory trait for creating solver instances
pub trait SolverFactory: Send + Sync + 'static {
    /// Solver type produced by this factory
    type Solver: SimulationPlugin;

    /// Create a new solver instance
    fn create_solver(&self, params: SolverParams) -> Result<Self::Solver>;

    /// Get factory metadata
    fn metadata(&self) -> FactoryMetadata {
        FactoryMetadata::default()
    }

    /// Validate parameters before solver creation
    fn validate_params(&self, _params: &SolverParams) -> Result<()> {
        // Default implementation does no validation
        Ok(())
    }
}

/// Builder pattern for solver parameters
#[derive(Debug, Clone, Default)]
pub struct SolverParams {
    /// Solver-specific parameters
    params: HashMap<String, serde_json::Value>,
}

impl SolverParams {
    /// Create new solver parameters
    pub fn new() -> Self {
        Self {
            params: HashMap::new(),
        }
    }

    /// Builder method to add a parameter
    pub fn with<T: serde::Serialize>(mut self, key: impl Into<String>, value: T) -> Result<Self> {
        let value = serde_json::to_value(value)
            .map_err(|e| Error::SerializationError(e.to_string()))?;
        self.params.insert(key.into(), value);
        Ok(self)
    }

    /// Get a parameter
    pub fn get<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<T> {
        let value = self
            .params
            .get(key)
            .ok_or_else(|| Error::InvalidConfiguration(format!("Parameter '{}' not found", key)))?;

        serde_json::from_value(value.clone())
            .map_err(|e| Error::SerializationError(e.to_string()))
    }

    /// Get optional parameter
    pub fn get_optional<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.params.get(key) {
            Some(value) => serde_json::from_value(value.clone())
                .map(Some)
                .map_err(|e| Error::SerializationError(e.to_string())),
            None => Ok(None),
        }
    }

    /// Check if parameter exists
    pub fn contains(&self, key: &str) -> bool {
        self.params.contains_key(key)
    }

    /// Iterate over all parameters
    pub fn iter(&self) -> impl Iterator<Item = (&String, &serde_json::Value)> {
        self.params.iter()
    }

    /// Merge with another set of parameters
    pub fn merge(mut self, other: Self) -> Self {
        self.params.extend(other.params);
        self
    }
}

/// Factory metadata
#[derive(Debug, Clone, Default)]
pub struct FactoryMetadata {
    /// Factory name
    pub name: String,
    /// Supported solver types
    pub solver_types: Vec<String>,
    /// Factory description
    pub description: String,
    /// Required parameters
    pub required_params: Vec<String>,
    /// Optional parameters with defaults
    pub optional_params: HashMap<String, serde_json::Value>,
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestPlugin {
        name: String,
    }

    impl Plugin for TestPlugin {
        fn name(&self) -> &str {
            &self.name
        }

        fn version(&self) -> &str {
            "1.0.0"
        }

        fn description(&self) -> &str {
            "Test plugin"
        }

        fn type_id(&self) -> TypeId {
            TypeId::of::<Self>()
        }

        fn as_any(&self) -> &dyn Any {
            self
        }

        fn as_any_mut(&mut self) -> &mut dyn Any {
            self
        }
    }

    #[test]
    fn test_plugin_registry() {
        let registry = PluginRegistry::new();
        let plugin = TestPlugin {
            name: "test".to_string(),
        };

        // Register plugin
        registry.register(plugin).unwrap();

        // Get plugin
        let retrieved = registry.get("test").unwrap();
        assert_eq!(retrieved.name(), "test");

        // List plugins
        let list = registry.list().unwrap();
        assert_eq!(list.len(), 1);
        assert!(list.contains(&"test".to_string()));
    }

    #[test]
    fn test_solver_params() {
        let params = SolverParams::new()
            .with("tolerance", 1e-6)
            .unwrap()
            .with("max_iterations", 1000)
            .unwrap();

        let tolerance: f64 = params.get("tolerance").unwrap();
        let max_iter: i32 = params.get("max_iterations").unwrap();

        assert_eq!(tolerance, 1e-6);
        assert_eq!(max_iter, 1000);

        // Test optional parameter
        let optional: Option<String> = params.get_optional("missing").unwrap();
        assert!(optional.is_none());
    }

    #[test]
    fn test_params_builder_pattern() {
        let params = SolverParams::new()
            .with("a", 1)
            .unwrap()
            .with("b", 2.0)
            .unwrap()
            .with("c", "test")
            .unwrap();

        assert_eq!(params.get::<i32>("a").unwrap(), 1);
        assert_eq!(params.get::<f64>("b").unwrap(), 2.0);
        assert_eq!(params.get::<String>("c").unwrap(), "test");
    }
}