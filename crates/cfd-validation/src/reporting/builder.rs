use super::{
    CodeQualityReport, PerformanceReport, TestCategory, ValidationReport, ValidationSummary,
};
use std::collections::HashMap;
use std::time::SystemTime;

/// Validation report builder
pub struct ReportBuilder {
    title: String,
    summary: ValidationSummary,
    test_results: HashMap<String, TestCategory>,
    performance: Vec<PerformanceReport>,
    code_quality: CodeQualityReport,
    recommendations: Vec<String>,
}

impl ReportBuilder {
    /// Create a new report builder with the specified title
    pub fn new(title: String) -> Self {
        Self {
            title,
            summary: ValidationSummary::default(),
            test_results: HashMap::new(),
            performance: Vec::new(),
            code_quality: CodeQualityReport::default(),
            recommendations: Vec::new(),
        }
    }

    /// Set the validation summary
    pub fn with_summary(mut self, summary: ValidationSummary) -> Self {
        self.summary = summary;
        self
    }

    /// Add a test category result
    pub fn add_test_category(mut self, category: TestCategory) -> Self {
        self.test_results.insert(category.name.clone(), category);
        self
    }

    /// Add a performance benchmark report
    pub fn add_performance_report(mut self, report: PerformanceReport) -> Self {
        self.performance.push(report);
        self
    }

    /// Set code quality metrics
    pub fn with_code_quality(mut self, quality: CodeQualityReport) -> Self {
        self.code_quality = quality;
        self
    }

    /// Add a recommendation for improvement
    pub fn add_recommendation(mut self, recommendation: String) -> Self {
        self.recommendations.push(recommendation);
        self
    }

    /// Build the final validation report
    pub fn build(self) -> ValidationReport {
        ValidationReport {
            timestamp: SystemTime::now(),
            title: self.title,
            summary: self.summary,
            test_results: self.test_results,
            performance: self.performance,
            code_quality: self.code_quality,
            recommendations: self.recommendations,
        }
    }
}
