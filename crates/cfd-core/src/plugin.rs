//! Plugin system for extensible CFD solvers.
//!
//! This module provides a simplified plugin architecture following SOLID principles,
//! with coarse-grained locking for safety and strongly-typed configurations.

use crate::error::{Error, Result, PluginErrorKind};

use std::collections::{HashMap, HashSet};
use std::sync::Arc;


/// Core plugin trait - the single abstraction for all plugin functionality
pub trait Plugin: Send + Sync {
    /// Plugin name (must be unique)
    fn name(&self) -> &str;
    
    /// Plugin version
    fn version(&self) -> &str;
    
    /// Plugin description
    fn description(&self) -> &str;
    
    /// Plugin author (optional)
    fn author(&self) -> Option<&str> {
        None
    }
    
    /// Plugin license (optional)
    fn license(&self) -> Option<&str> {
        None
    }
    
    /// Plugin dependencies (other plugin names)
    fn dependencies(&self) -> Vec<&str> {
        vec![]
    }
    
    /// Plugin capabilities
    fn capabilities(&self) -> Vec<&str> {
        vec![]
    }
    
    
    /// Check if plugin supports parallel execution
    fn supports_parallel(&self) -> bool {
        false
    }
    
    /// Check if plugin supports adaptive time stepping
    fn supports_adaptive_timestep(&self) -> bool {
        false
    }
    
    /// Initialize the plugin (called once when registered)
    fn initialize(&self) -> Result<()> {
        Ok(())
    }
    
    /// Shutdown the plugin (called when unregistered)
    fn shutdown(&self) -> Result<()> {
        Ok(())
    }
}

/// Trait for simulation plugins with typed configuration and state
pub trait SimulationPlugin: Plugin {
    /// Configuration type for this plugin
    type Config: Send + Sync + 'static;
    /// State type maintained by this plugin
    type State: Send + Sync + 'static;
    /// Output type produced by this plugin
    type Output: Send + Sync + 'static;

    /// Create initial state from configuration
    fn create_state(&self, config: Self::Config) -> Result<Self::State>;

    /// Perform one simulation step
    fn step(&self, state: &mut Self::State, dt: f64) -> Result<()>;

    /// Generate output from current state
    fn output(&self, state: &Self::State) -> Self::Output;
    
    /// Validate configuration parameters
    fn validate_config(&self, config: &Self::Config) -> Result<()> {
        let _ = config; // Avoid unused parameter warning
        Ok(())
    }
    
    /// Validate plugin state for physics consistency
    fn validate_state(&self, state: &Self::State) -> Result<()> {
        let _ = state; // Avoid unused parameter warning
        Ok(())
    }
}

/// Plugin health status for monitoring
#[derive(Debug, Clone)]
pub enum PluginHealthStatus {
    /// Plugin is functioning normally
    Healthy,
    /// Plugin is functional but with issues
    Degraded(String),
    /// Plugin is non-functional
    Unhealthy(String),
}

/// Plugin metrics for monitoring
#[derive(Debug, Clone, Default)]
pub struct PluginMetrics {
    /// CPU usage percentage
    pub cpu_usage: f64,
    /// Memory usage in bytes
    pub memory_usage: u64,
    /// Number of successful operations
    pub success_count: u64,
    /// Number of failed operations
    pub error_count: u64,
    /// Average operation time in milliseconds
    pub avg_operation_time_ms: f64,
}

/// Standard data storage for plugins (no internal locking)
struct PluginStorage {
    plugins: HashMap<String, Arc<dyn Plugin>>,
}

impl PluginStorage {
    fn new() -> Self {
        Self {
            plugins: HashMap::new(),
        }
    }
    
    fn insert(&mut self, plugin: Arc<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        if self.plugins.contains_key(&name) {
            return Err(Error::Plugin(PluginErrorKind::AlreadyRegistered { name }));
        }
        self.plugins.insert(name, plugin);
        Ok(())
    }
    
    fn get(&self, name: &str) -> Option<&Arc<dyn Plugin>> {
        self.plugins.get(name)
    }
    
    fn remove(&mut self, name: &str) -> Option<Arc<dyn Plugin>> {
        self.plugins.remove(name)
    }
    
    fn list(&self) -> Vec<String> {
        self.plugins.keys().cloned().collect()
    }
}

/// Standard dependency resolver (no internal locking)
struct DependencyResolver {
    dependencies: HashMap<String, Vec<String>>,
    load_order: Vec<String>,
}

impl DependencyResolver {
    fn new() -> Self {
        Self {
            dependencies: HashMap::new(),
            load_order: Vec::new(),
        }
    }
    
    fn add_plugin(&mut self, name: String, deps: Vec<String>) -> Result<()> {
        self.dependencies.insert(name, deps);
        self.update_load_order()
    }
    
    fn remove_plugin(&mut self, name: &str) {
        self.dependencies.remove(name);
        let _ = self.update_load_order(); // Ignore errors during removal
    }
    
    fn validate_dependencies(&self, deps: &[String], available_plugins: &[String]) -> Result<()> {
        for dep in deps {
            if !available_plugins.contains(dep) {
                return Err(Error::Plugin(PluginErrorKind::DependencyNotSatisfied {
                    plugin: "unknown".to_string(),
                    dependency: dep.to_string(),
                }));
            }
        }
        Ok(())
    }
    
    fn get_load_order(&self) -> &[String] {
        &self.load_order
    }
    
    /// Update plugin load order using topological sort with cycle detection
    fn update_load_order(&mut self) -> Result<()> {
        let mut visited = HashSet::new();
        let mut recursion_stack = HashSet::new();
        let mut result = Vec::new();
        
        fn visit(
            plugin: &str,
            deps: &HashMap<String, Vec<String>>,
            recursion_stack: &mut HashSet<String>,
            visited: &mut HashSet<String>,
            result: &mut Vec<String>,
        ) -> Result<()> {
            // Check for cycles
            if recursion_stack.contains(plugin) {
                return Err(Error::Plugin(PluginErrorKind::CircularDependency {
                    chain: vec![plugin.to_string()]
                }));
            }
            
            // Already processed
            if visited.contains(plugin) {
                return Ok(());
            }
            
            // Mark as being processed
            recursion_stack.insert(plugin.to_string());
            visited.insert(plugin.to_string());
            
            // Process dependencies
            if let Some(plugin_deps) = deps.get(plugin) {
                for dep in plugin_deps {
                    visit(dep, deps, recursion_stack, visited, result)?;
                }
            }
            
            // Remove from recursion stack and add to result
            recursion_stack.remove(plugin);
            result.push(plugin.to_string());
            Ok(())
        }
        
        // Process all plugins
        for plugin in self.dependencies.keys() {
            visit(plugin, &self.dependencies, &mut recursion_stack, &mut visited, &mut result)?;
        }
        
        self.load_order = result;
        Ok(())
    }
}

/// Standard monitoring system (no internal locking)
struct PluginMonitoring {
    health_status: HashMap<String, PluginHealthStatus>,
    metrics: HashMap<String, PluginMetrics>,
    last_check: HashMap<String, std::time::Instant>,
}

impl PluginMonitoring {
    fn new() -> Self {
        Self {
            health_status: HashMap::new(),
            metrics: HashMap::new(),
            last_check: HashMap::new(),
        }
    }
    
    fn update_health(&mut self, plugin_name: &str, status: PluginHealthStatus) {
        self.health_status.insert(plugin_name.to_string(), status);
        self.last_check.insert(plugin_name.to_string(), std::time::Instant::now());
    }
    
    fn update_metrics(&mut self, plugin_name: &str, metrics: PluginMetrics) {
        self.metrics.insert(plugin_name.to_string(), metrics);
    }
    
    fn get_health(&self, plugin_name: &str) -> Option<&PluginHealthStatus> {
        self.health_status.get(plugin_name)
    }
    
    fn get_metrics(&self, plugin_name: &str) -> Option<&PluginMetrics> {
        self.metrics.get(plugin_name)
    }
}

/// Main plugin registry - single owner of all plugin data
/// 
/// To share between threads, wrap the entire registry in Arc<`RwLock`<PluginRegistry>>
pub struct PluginRegistry {
    storage: PluginStorage,
    resolver: DependencyResolver,
    monitoring: PluginMonitoring,
}

impl PluginRegistry {
    /// Create a new plugin registry
    #[must_use] pub fn new() -> Self {
        Self {
            storage: PluginStorage::new(),
            resolver: DependencyResolver::new(),
            monitoring: PluginMonitoring::new(),
        }
    }
    
    /// Register a plugin
    pub fn register(&mut self, plugin: Arc<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        let deps: Vec<String> = plugin.dependencies().iter().map(|s| (*s).to_string()).collect();
        
        // Validate dependencies exist
        let available: Vec<String> = self.storage.list();
        self.resolver.validate_dependencies(&deps, &available)?;
        
        // Initialize the plugin
        plugin.initialize()?;
        
        // Store the plugin
        self.storage.insert(plugin)?;
        
        // Update dependency graph
        self.resolver.add_plugin(name.clone(), deps)?;
        
        // Set initial health status
        self.monitoring.update_health(&name, PluginHealthStatus::Healthy);
        
        Ok(())
    }
    
    /// Unregister a plugin
    pub fn unregister(&mut self, name: &str) -> Result<()> {
        // Check if other plugins depend on this one
        for (plugin_name, deps) in &self.resolver.dependencies {
            if plugin_name != name && deps.contains(&name.to_string()) {
                return Err(Error::Plugin(PluginErrorKind::DependencyNotSatisfied {
                    plugin: plugin_name.to_string(),
                    dependency: name.to_string(),
                }));
            }
        }
        
        // Remove and shutdown the plugin
        if let Some(plugin) = self.storage.remove(name) {
            plugin.shutdown()?;
        }
        
        // Update dependency graph
        self.resolver.remove_plugin(name);
        
        Ok(())
    }
    
    /// Get a plugin by name
    #[must_use] pub fn get(&self, name: &str) -> Option<&Arc<dyn Plugin>> {
        self.storage.get(name)
    }
    

    
    /// List all registered plugins
    #[must_use] pub fn list(&self) -> Vec<String> {
        self.storage.list()
    }
    
    /// Get plugins in dependency order
    #[must_use] pub fn get_load_order(&self) -> &[String] {
        self.resolver.get_load_order()
    }
    
    /// Filter plugins by capability
    #[must_use] pub fn filter_by_capability(&self, capability: &str) -> Vec<String> {
        self.storage
            .plugins
            .iter()
            .filter_map(|(name, plugin)| {
                if plugin.capabilities().contains(&capability) {
                    Some(name.clone())
                } else {
                    None
                }
            })
            .collect()
    }
    
    /// Update plugin health status
    pub fn update_health(&mut self, plugin_name: &str, status: PluginHealthStatus) {
        self.monitoring.update_health(plugin_name, status);
    }
    
    /// Update plugin metrics
    pub fn update_metrics(&mut self, plugin_name: &str, metrics: PluginMetrics) {
        self.monitoring.update_metrics(plugin_name, metrics);
    }
    
    /// Get plugin health status
    #[must_use] pub fn get_health(&self, plugin_name: &str) -> Option<&PluginHealthStatus> {
        self.monitoring.get_health(plugin_name)
    }
    
    /// Get plugin metrics
    #[must_use] pub fn get_metrics(&self, plugin_name: &str) -> Option<&PluginMetrics> {
        self.monitoring.get_metrics(plugin_name)
    }
    
    /// Get a summary of system health
    #[must_use] pub fn system_health_summary(&self) -> SystemHealthSummary {
        let mut healthy = 0;
        let mut degraded = 0;
        let mut unhealthy = 0;
        
        for status in self.monitoring.health_status.values() {
            match status {
                PluginHealthStatus::Healthy => healthy += 1,
                PluginHealthStatus::Degraded(_) => degraded += 1,
                PluginHealthStatus::Unhealthy(_) => unhealthy += 1,
            }
        }
        
        let total = healthy + degraded + unhealthy;
        let status = if unhealthy > 0 {
            SystemStatus::Critical
        } else if degraded > 0 {
            SystemStatus::Degraded
        } else {
            SystemStatus::Healthy
        };
        
        SystemHealthSummary {
            total_plugins: total,
            healthy_plugins: healthy,
            degraded_plugins: degraded,
            unhealthy_plugins: unhealthy,
            system_status: status,
        }
    }
}

impl Default for PluginRegistry {
    fn default() -> Self {
        Self::new()
    }
}

/// System health summary
#[derive(Debug, Clone)]
pub struct SystemHealthSummary {
    /// Total number of plugins
    pub total_plugins: usize,
    /// Number of healthy plugins
    pub healthy_plugins: usize,
    /// Number of degraded plugins
    pub degraded_plugins: usize,
    /// Number of unhealthy plugins
    pub unhealthy_plugins: usize,
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

#[cfg(test)]
mod tests {
    use super::*;
    
    /// Test plugin implementation
    struct TestPlugin {
        name: String,
        dependencies: Vec<String>,
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
        
        fn dependencies(&self) -> Vec<&str> {
            self.dependencies.iter().map(|s| s.as_str()).collect()
        }
    }
    
    #[test]
    fn test_plugin_registration() {
        let mut registry = PluginRegistry::new();
        
        let plugin = Arc::new(TestPlugin {
            name: "test".to_string(),
            dependencies: vec![],
        });
        
        assert!(registry.register(plugin).is_ok());
        assert!(registry.get("test").is_some());
        assert_eq!(registry.list(), vec!["test"]);
    }
    
    #[test]
    fn test_dependency_resolution() {
        let mut registry = PluginRegistry::new();
        
        // Register plugins with dependencies
        let plugin_a = Arc::new(TestPlugin {
            name: "A".to_string(),
            dependencies: vec![],
        });
        
        let plugin_b = Arc::new(TestPlugin {
            name: "B".to_string(),
            dependencies: vec!["A".to_string()],
        });
        
        let plugin_c = Arc::new(TestPlugin {
            name: "C".to_string(),
            dependencies: vec!["A".to_string(), "B".to_string()],
        });
        
        // Register in any order - should work
        assert!(registry.register(plugin_a).is_ok());
        assert!(registry.register(plugin_b).is_ok());
        assert!(registry.register(plugin_c).is_ok());
        
        // Check load order
        let load_order = registry.get_load_order();
        let a_idx = load_order.iter().position(|x| x == "A").expect("CRITICAL: Add proper error handling");
        let b_idx = load_order.iter().position(|x| x == "B").expect("CRITICAL: Add proper error handling");
        let c_idx = load_order.iter().position(|x| x == "C").expect("CRITICAL: Add proper error handling");
        
        assert!(a_idx < b_idx);
        assert!(b_idx < c_idx);
    }
    
    #[test]
    fn test_circular_dependency_detection() {
        let mut resolver = DependencyResolver::new();
        
        // Create circular dependency: A -> B -> C -> A
        resolver.dependencies.insert("A".to_string(), vec!["B".to_string()]);
        resolver.dependencies.insert("B".to_string(), vec!["C".to_string()]);
        resolver.dependencies.insert("C".to_string(), vec!["A".to_string()]);
        
        // Should detect the cycle
        assert!(resolver.update_load_order().is_err());
    }
    
    #[test]
    fn test_missing_dependency() {
        let mut registry = PluginRegistry::new();
        
        let plugin = Arc::new(TestPlugin {
            name: "test".to_string(),
            dependencies: vec!["nonexistent".to_string()],
        });
        
        // Should fail due to missing dependency
        assert!(registry.register(plugin).is_err());
    }
    
    #[test]
    fn test_health_monitoring() {
        let mut registry = PluginRegistry::new();
        
        let plugin = Arc::new(TestPlugin {
            name: "test".to_string(),
            dependencies: vec![],
        });
        
        registry.register(plugin).expect("CRITICAL: Add proper error handling");
        
        // Check initial health
        assert!(matches!(
            registry.get_health("test"),
            Some(PluginHealthStatus::Healthy)
        ));
        
        // Update health
        registry.update_health("test", PluginHealthStatus::Degraded("Test issue".to_string()));
        
        assert!(matches!(
            registry.get_health("test"),
            Some(PluginHealthStatus::Degraded(_))
        ));
        
        // Check system summary
        let summary = registry.system_health_summary();
        assert_eq!(summary.degraded_plugins, 1);
        assert_eq!(summary.system_status, SystemStatus::Degraded);
    }
    
    #[test]
    fn test_plugin_removal() {
        let mut registry = PluginRegistry::new();
        
        let plugin_a = Arc::new(TestPlugin {
            name: "A".to_string(),
            dependencies: vec![],
        });
        
        let plugin_b = Arc::new(TestPlugin {
            name: "B".to_string(),
            dependencies: vec!["A".to_string()],
        });
        
        registry.register(plugin_a).expect("CRITICAL: Add proper error handling");
        registry.register(plugin_b).expect("CRITICAL: Add proper error handling");
        
        // Cannot remove A because B depends on it
        assert!(registry.unregister("A").is_err());
        
        // Can remove B
        assert!(registry.unregister("B").is_ok());
        
        // Now can remove A
        assert!(registry.unregister("A").is_ok());
        
        assert_eq!(registry.list().len(), 0);
    }
}