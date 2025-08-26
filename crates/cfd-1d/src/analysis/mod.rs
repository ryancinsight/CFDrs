//! Network analysis and performance metrics for 1D CFD simulations.
//!
//! This module provides comprehensive analysis tools for network systems,
//! organized into domain-specific submodules.

pub mod analyzer;
pub mod flow;
pub mod performance;
pub mod pressure;
pub mod resistance;

// Re-export main types for convenience
pub use analyzer::{NetworkAnalysisResult, NetworkAnalyzer};
pub use flow::FlowAnalysis;
pub use performance::PerformanceMetrics;
pub use pressure::PressureAnalysis;
pub use resistance::ResistanceAnalysis;
