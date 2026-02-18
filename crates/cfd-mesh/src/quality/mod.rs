//! Mesh quality assessment.
//!
//! Triangle quality metrics relevant to CFD meshing:
//! aspect ratio, skewness, minimum angle, etc.

pub mod metrics;
pub mod triangle;
pub mod validation;

pub use metrics::QualityMetric;
pub use validation::MeshValidator;
