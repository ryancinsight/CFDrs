//! Mesh quality assessment.
//!
//! Triangle quality metrics relevant to CFD meshing:
//! aspect ratio, skewness, minimum angle, mean curvature.
//!
//! ## Module Overview
//!
//! | Module | Purpose |
//! |--------|---------|
//! | `metrics` | `QualityMetric` — min/max/mean per-scalar |
//! | `triangle` | Per-triangle measurement functions |
//! | `validation` | `MeshValidator`, `QualityReport`, `QualityThresholds` |
//! | `normals` | Normal consistency analysis |
//! | `histograms` | `Histogram` — fixed-width distribution |
//! | `report` | `FullQualityReport` — base + histograms + counts |
//! | `analyzer` | `QualityAnalyzer` trait + `StandardQualityAnalyzer` |
//! | `curvature` | `vertex_mean_curvature` via cotangent Laplacian |

pub mod analyzer;
pub mod curvature;
pub mod histograms;
pub mod metrics;
pub mod normals;
pub mod report;
pub mod triangle;
pub mod validation;

pub use analyzer::{QualityAnalyzer, StandardQualityAnalyzer};
pub use curvature::vertex_mean_curvature;
pub use histograms::Histogram;
pub use metrics::QualityMetric;
pub use normals::{analyze_normals, NormalAnalysis};
pub use report::FullQualityReport;
pub use validation::MeshValidator;
