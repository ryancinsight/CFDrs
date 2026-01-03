//! Health monitoring and metrics for plugins.

use std::collections::HashMap;
use std::time::Instant;

/// Plugin health status
#[derive(Debug, Clone)]
pub enum PluginHealthStatus {
    /// Plugin is healthy and operating normally
    Healthy,
    /// Plugin is degraded but still functional
    Degraded(String),
    /// Plugin has failed
    Failed(String),
    /// Plugin health is unknown
    Unknown,
}

/// Plugin performance metrics
#[derive(Debug, Clone)]
pub struct PluginMetrics {
    /// Number of successful executions
    pub executions: u64,
    /// Number of failed executions
    pub failures: u64,
    /// Average execution time in milliseconds
    pub avg_execution_time_ms: f64,
    /// Peak memory usage in bytes
    pub peak_memory_bytes: u64,
}

/// System health summary
#[derive(Debug, Clone)]
pub struct SystemHealthSummary {
    /// Overall system status
    pub status: SystemStatus,
    /// Total number of plugins
    pub total_plugins: usize,
    /// Number of healthy plugins
    pub healthy_plugins: usize,
    /// Number of degraded plugins
    pub degraded_plugins: usize,
    /// Number of failed plugins
    pub failed_plugins: usize,
    /// Plugin-specific health details
    pub plugin_details: HashMap<String, PluginHealthStatus>,
}

/// Overall system status
#[derive(Debug, Clone)]
pub enum SystemStatus {
    /// All plugins healthy
    Healthy,
    /// Some plugins degraded but system functional
    Degraded,
    /// Critical failures affecting system operation
    Critical,
}

/// Standard monitoring system (no internal locking)
pub(crate) struct PluginMonitoring {
    health_status: HashMap<String, PluginHealthStatus>,
    metrics: HashMap<String, PluginMetrics>,
    last_check: HashMap<String, Instant>,
}

impl PluginMonitoring {
    pub fn new() -> Self {
        Self {
            health_status: HashMap::new(),
            metrics: HashMap::new(),
            last_check: HashMap::new(),
        }
    }

    pub fn update_health(&mut self, plugin_name: &str, status: PluginHealthStatus) {
        self.health_status.insert(plugin_name.to_string(), status);
        self.last_check
            .insert(plugin_name.to_string(), Instant::now());
    }

    pub fn update_metrics(&mut self, plugin_name: &str, metrics: PluginMetrics) {
        self.metrics.insert(plugin_name.to_string(), metrics);
    }

    pub fn get_health(&self, plugin_name: &str) -> Option<&PluginHealthStatus> {
        self.health_status.get(plugin_name)
    }

    pub fn get_metrics(&self, plugin_name: &str) -> Option<&PluginMetrics> {
        self.metrics.get(plugin_name)
    }

    pub fn get_system_health(&self, plugin_names: &[String]) -> SystemHealthSummary {
        let mut healthy = 0;
        let mut degraded = 0;
        let mut failed = 0;
        let mut details = HashMap::new();

        for name in plugin_names {
            let status = self
                .health_status
                .get(name)
                .cloned()
                .unwrap_or(PluginHealthStatus::Unknown);

            match &status {
                PluginHealthStatus::Healthy => healthy += 1,
                PluginHealthStatus::Degraded(_) => degraded += 1,
                PluginHealthStatus::Failed(_) => failed += 1,
                PluginHealthStatus::Unknown => {}
            }

            details.insert(name.clone(), status);
        }

        let system_status = if failed > 0 {
            SystemStatus::Critical
        } else if degraded > 0 {
            SystemStatus::Degraded
        } else {
            SystemStatus::Healthy
        };

        SystemHealthSummary {
            status: system_status,
            total_plugins: plugin_names.len(),
            healthy_plugins: healthy,
            degraded_plugins: degraded,
            failed_plugins: failed,
            plugin_details: details,
        }
    }
}
