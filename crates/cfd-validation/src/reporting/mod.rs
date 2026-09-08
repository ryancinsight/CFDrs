//! Automated validation report generation
//!
//! Generates comprehensive validation reports including:
//! - Test results and coverage analysis
//! - Performance benchmarks and regressions
//! - Convergence studies and error analysis
//! - Code quality metrics

mod automation;
mod builder;
mod data;
#[cfg(test)]
mod tests;

pub mod fidelity_limits;
pub mod html;
pub mod json;
pub mod markdown;
pub mod summary;

pub use automation::AutomatedReporter;
pub use builder::ReportBuilder;
pub use data::{
    CodeQualityReport, PerformanceReport, Reporter, TestCategory, TestResult, TestStatus,
    ValidationReport,
};
pub use html::HtmlReporter;
pub use json::JsonReporter;
pub use markdown::MarkdownReporter;
pub use summary::{PerformanceMetrics, ValidationSummary};
