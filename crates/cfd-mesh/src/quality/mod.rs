//! Mesh quality assessment module
//!
//! Provides comprehensive mesh quality metrics following established CFD literature.

pub mod analyzer;
pub mod criteria;
pub mod metrics;
pub mod statistics;

pub use analyzer::QualityAnalyzer;
pub use criteria::{QualityCriteria, QualityThresholds};
pub use metrics::QualityMetrics;
pub use statistics::QualityStatistics;
