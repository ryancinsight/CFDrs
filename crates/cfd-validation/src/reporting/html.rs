//! HTML report generation for validation results

use super::{Reporter, ValidationReport};
use cfd_core::error::Result;

/// HTML reporter for validation results (stub implementation)
pub struct HtmlReporter;

impl Reporter for HtmlReporter {
    fn generate_report(&self, _report: &ValidationReport) -> Result<String> {
        // Stub implementation - would generate full HTML report
        Ok("<html><body><h1>CFD Validation Report</h1><p>HTML report generation not yet implemented</p></body></html>".to_string())
    }
}
