//! Network analysis and performance metrics for 1D CFD simulations.
//!
//! This module provides comprehensive analysis tools for network systems,
//! organized into domain-specific submodules.

pub mod flow;
pub mod pressure;
pub mod resistance;
pub mod performance;
pub mod analyzer;

// Re-export main types for convenience
pub use flow::FlowAnalysis;
pub use pressure::PressureAnalysis;
pub use resistance::ResistanceAnalysis;
pub use performance::PerformanceMetrics;
pub use analyzer::{NetworkAnalyzer, NetworkAnalysisResult};