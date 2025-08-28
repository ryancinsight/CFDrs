//! Network analysis and performance metrics for 1D CFD simulations.
//!
//! This module provides comprehensive analysis tools for network systems,
//! organized into domain-specific submodules.

pub mod analyzers;
pub mod flow;
pub mod performance;
pub mod pressure;
pub mod resistance;

// Re-export main types for convenience
pub use analyzers::{NetworkAnalysisResult, NetworkAnalyzerOrchestrator};
pub use flow::FlowAnalysis;
pub use performance::PerformanceMetrics;
pub use pressure::PressureAnalysis;
pub use resistance::ResistanceAnalysis;
