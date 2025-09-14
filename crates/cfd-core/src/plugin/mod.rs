//! Plugin system for extensible CFD solvers.
//!
//! This module provides a production-grade plugin architecture following SOLID principles,
//! with coarse-grained locking for safety and strongly-typed configurations.

pub mod dependency;
pub mod health;
pub mod registry;
pub mod storage;
pub mod traits;

// Re-export main types for convenience
pub use health::{PluginHealthStatus, PluginMetrics, SystemHealthSummary, SystemStatus};
pub use registry::PluginRegistry;
pub use traits::{Plugin, SimulationPlugin};

#[cfg(test)]
mod tests;
