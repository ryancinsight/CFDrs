//! Plugin registry for managing plugin lifecycle.

use crate::error::{Error, PluginErrorKind, Result};
use crate::plugin::{
    dependency::DependencyResolver,
    health::{PluginHealthStatus, PluginMetrics, PluginMonitoring, SystemHealthSummary},
    storage::PluginStorage,
    traits::Plugin,
};
use std::sync::Arc;

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
    #[must_use]
    pub fn new() -> Self {
        Self {
            storage: PluginStorage::new(),
            resolver: DependencyResolver::new(),
            monitoring: PluginMonitoring::new(),
        }
    }

    /// Register a plugin with the registry
    pub fn register(&mut self, plugin: Arc<dyn Plugin>) -> Result<()> {
        let name = plugin.name().to_string();
        let deps: Vec<String> = plugin
            .dependencies()
            .iter()
            .map(|s| s.to_string())
            .collect();

        // Validate dependencies exist
        let available: Vec<String> = self.storage.list();
        self.resolver.validate_dependencies(&deps, &available)?;

        // Add to resolver first (can fail)
        self.resolver.add_plugin(name.clone(), deps)?;

        // Then add to storage
        if let Err(e) = self.storage.insert(plugin) {
            // Rollback resolver change
            self.resolver.remove_plugin(&name);
            return Err(e);
        }

        // Initialize health status
        self.monitoring
            .update_health(&name, PluginHealthStatus::Unknown);

        Ok(())
    }

    /// Unregister a plugin from the registry
    pub fn unregister(&mut self, name: &str) -> Result<()> {
        // Check if other plugins depend on this one
        for (plugin_name, plugin) in self.storage.iter() {
            if plugin_name != name && plugin.dependencies().contains(&name) {
                return Err(Error::Plugin(PluginErrorKind::HasDependents {
                    plugin: name.to_string(),
                    dependents: vec![plugin_name.clone()],
                }));
            }
        }

        // Remove from all subsystems
        self.storage.remove(name);
        self.resolver.remove_plugin(name);
        // Monitoring data can stay for historical purposes

        Ok(())
    }

    /// Get a plugin by name
    pub fn get(&self, name: &str) -> Option<Arc<dyn Plugin>> {
        self.storage.get(name).cloned()
    }

    /// List all registered plugins
    pub fn list(&self) -> Vec<String> {
        self.storage.list()
    }

    /// Get plugin load order based on dependencies
    pub fn get_load_order(&self) -> Vec<String> {
        self.resolver.get_load_order().to_vec()
    }

    /// Update health status for a plugin
    pub fn update_health(&mut self, plugin_name: &str, status: PluginHealthStatus) {
        self.monitoring.update_health(plugin_name, status);
    }

    /// Update metrics for a plugin
    pub fn update_metrics(&mut self, plugin_name: &str, metrics: PluginMetrics) {
        self.monitoring.update_metrics(plugin_name, metrics);
    }

    /// Get health status for a plugin
    pub fn get_health(&self, plugin_name: &str) -> Option<PluginHealthStatus> {
        self.monitoring.get_health(plugin_name).cloned()
    }

    /// Get metrics for a plugin
    pub fn get_metrics(&self, plugin_name: &str) -> Option<PluginMetrics> {
        self.monitoring.get_metrics(plugin_name).cloned()
    }

    /// Get overall system health
    pub fn get_system_health(&self) -> SystemHealthSummary {
        let plugin_names = self.storage.list();
        self.monitoring.get_system_health(&plugin_names)
    }

    /// Check if a plugin is registered
    pub fn contains(&self, name: &str) -> bool {
        self.storage.contains(name)
    }

    /// Get the number of registered plugins
    pub fn len(&self) -> usize {
        self.storage.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.storage.is_empty()
    }

    /// Clear all plugins from the registry
    pub fn clear(&mut self) {
        self.storage.clear();
        self.resolver = DependencyResolver::new();
        self.monitoring = PluginMonitoring::new();
    }

    /// Validate plugin configuration
    pub fn validate_all(&self) -> Result<()> {
        for name in self.storage.list() {
            if let Some(plugin) = self.storage.get(&name) {
                // Check dependencies
                let deps: Vec<String> = plugin
                    .dependencies()
                    .iter()
                    .map(|s| s.to_string())
                    .collect();
                let available = self.storage.list();
                self.resolver.validate_dependencies(&deps, &available)?;
            }
        }
        Ok(())
    }
}

impl Default for PluginRegistry {
    fn default() -> Self {
        Self::new()
    }
}
