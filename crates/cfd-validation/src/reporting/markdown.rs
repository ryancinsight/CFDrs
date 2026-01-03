//! Markdown report generation for validation results

use super::{Reporter, ValidationReport};
use cfd_core::error::Result;

/// Markdown reporter for validation results
pub struct MarkdownReporter {
    include_details: bool,
    include_performance: bool,
}

impl MarkdownReporter {
    /// Create a new Markdown reporter with default settings
    pub fn new() -> Self {
        Self {
            include_details: true,
            include_performance: true,
        }
    }

    /// Configure whether to include detailed test results
    pub fn with_details(mut self, include: bool) -> Self {
        self.include_details = include;
        self
    }

    /// Configure whether to include performance benchmarks
    pub fn with_performance(mut self, include: bool) -> Self {
        self.include_performance = include;
        self
    }
}

impl Default for MarkdownReporter {
    fn default() -> Self {
        Self::new()
    }
}

impl Reporter for MarkdownReporter {
    fn generate_report(&self, report: &ValidationReport) -> Result<String> {
        let mut markdown = String::new();

        // Header
        markdown.push_str("# CFD Validation Report\n\n");
        markdown.push_str(&format!(
            "**Generated:** {}\n\n",
            report
                .timestamp
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_secs()
        ));

        // Executive Summary
        markdown.push_str("## Executive Summary\n\n");
        markdown.push_str(&format!(
            "- **Total Tests:** {}\n",
            report.summary.total_tests
        ));
        markdown.push_str(&format!(
            "- **Passed:** {} ({:.1}%)\n",
            report.summary.passed_tests,
            (report.summary.passed_tests as f64 / report.summary.total_tests.max(1) as f64) * 100.0
        ));
        markdown.push_str(&format!("- **Failed:** {}\n", report.summary.failed_tests));
        markdown.push_str(&format!(
            "- **Skipped:** {}\n",
            report.summary.skipped_tests
        ));
        markdown.push_str(&format!(
            "- **Test Coverage:** {:.1}%\n",
            report.summary.coverage_percentage
        ));
        markdown.push_str(&format!(
            "- **Duration:** {:.1}s\n",
            report.summary.total_duration.as_secs_f64()
        ));
        markdown.push_str(&format!(
            "- **Health Score:** {:.3}\n\n",
            report.health_score()
        ));

        // Critical Issues
        let critical_issues = report.critical_issues();
        if !critical_issues.is_empty() {
            markdown.push_str("### 🚨 Critical Issues\n\n");
            for issue in &critical_issues {
                markdown.push_str(&format!("- {issue}\n"));
            }
            markdown.push('\n');
        }

        // Test Results by Category
        if !report.test_results.is_empty() {
            markdown.push_str("## Test Results\n\n");

            for (category_name, category) in &report.test_results {
                markdown.push_str(&format!("### {category_name}\n\n"));
                markdown.push_str(&format!("- **Total:** {}\n", category.total));
                markdown.push_str(&format!(
                    "- **Passed:** {} ({:.1}%)\n",
                    category.passed,
                    (category.passed as f64 / category.total.max(1) as f64) * 100.0
                ));
                markdown.push_str(&format!("- **Failed:** {}\n", category.failed));
                markdown.push_str(&format!(
                    "- **Coverage:** {:.1}%\n",
                    category.coverage_percentage
                ));

                if self.include_details && !category.details.is_empty() {
                    markdown.push_str("\n**Test Details:**\n\n");
                    markdown.push_str("| Test | Status | Duration |\n");
                    markdown.push_str("|------|--------|----------|\n");

                    for test in &category.details {
                        let status_emoji = match test.status {
                            super::TestStatus::Passed => "✅",
                            super::TestStatus::Failed => "❌",
                            super::TestStatus::Skipped => "⏭️",
                            super::TestStatus::Timeout => "⏰",
                        };

                        markdown.push_str(&format!(
                            "| {} | {} | {:.3}ms |\n",
                            test.name, status_emoji, test.duration_ms
                        ));
                    }
                }
                markdown.push('\n');
            }
        }

        // Code Quality
        markdown.push_str("## Code Quality Metrics\n\n");
        markdown.push_str(&format!(
            "- **Lines of Code:** {}\n",
            report.code_quality.lines_of_code
        ));
        markdown.push_str(&format!(
            "- **Test Coverage:** {:.1}%\n",
            report.code_quality.test_coverage
        ));
        markdown.push_str(&format!(
            "- **Documentation Coverage:** {:.1}%\n",
            report.code_quality.documentation_coverage
        ));
        markdown.push_str(&format!(
            "- **Clippy Warnings:** {}\n",
            report.code_quality.clippy_warnings
        ));
        markdown.push_str(&format!(
            "- **Compiler Errors:** {}\n",
            report.code_quality.compiler_errors
        ));
        markdown.push_str(&format!(
            "- **Cyclomatic Complexity:** {:.2}\n",
            report.code_quality.cyclomatic_complexity
        ));
        markdown.push_str(&format!(
            "- **Maintainability Index:** {:.2}\n\n",
            report.code_quality.maintainability_index
        ));

        // Performance Benchmarks
        if self.include_performance && !report.performance.is_empty() {
            markdown.push_str("## Performance Benchmarks\n\n");

            for benchmark in &report.performance {
                markdown.push_str(&format!("### {}\n\n", benchmark.benchmark_name));
                markdown.push_str(&format!(
                    "- **Mean:** {:.3}ms\n",
                    benchmark.metrics.mean * 1000.0
                ));
                markdown.push_str(&format!(
                    "- **Std Dev:** {:.3}ms\n",
                    benchmark.metrics.std_dev * 1000.0
                ));
                markdown.push_str(&format!(
                    "- **Min:** {:.3}ms\n",
                    benchmark.metrics.min * 1000.0
                ));
                markdown.push_str(&format!(
                    "- **Max:** {:.3}ms\n",
                    benchmark.metrics.max * 1000.0
                ));
                markdown.push_str(&format!("- **Samples:** {}\n", benchmark.metrics.samples));
                markdown.push_str(&format!(
                    "- **Stable:** {}\n",
                    if benchmark.metrics.is_stable(0.05) {
                        "✅"
                    } else {
                        "❌"
                    }
                ));

                if let Some(regression) = benchmark.regression_detected {
                    let direction = if regression > 0.0 { "slower" } else { "faster" };
                    markdown.push_str(&format!(
                        "- **Regression:** {:.2}% {}\n",
                        regression.abs(),
                        direction
                    ));
                }

                markdown.push('\n');
            }
        }

        // Recommendations
        if !report.recommendations.is_empty() {
            markdown.push_str("## Recommendations\n\n");
            for recommendation in &report.recommendations {
                markdown.push_str(&format!("- {recommendation}\n"));
            }
            markdown.push('\n');
        }

        // Footer
        markdown.push_str("---\n");
        markdown.push_str("*Generated by CFD Validation Framework*\n");

        Ok(markdown)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::benchmarking::PerformanceReport;
    use crate::reporting::{CodeQualityReport, TestCategory, ValidationSummary};
    use std::collections::HashMap;

    #[test]
    fn test_markdown_report_generation() {
        let report = ValidationReport {
            timestamp: std::time::SystemTime::UNIX_EPOCH,
            title: "Test CFD Report".to_string(),
            summary: ValidationSummary {
                total_tests: 100,
                passed_tests: 95,
                failed_tests: 5,
                skipped_tests: 0,
                total_duration: std::time::Duration::from_secs(30),
                coverage_percentage: 87.5,
            },
            test_results: {
                let mut results = HashMap::new();
                results.insert(
                    "unit_tests".to_string(),
                    TestCategory {
                        name: "Unit Tests".to_string(),
                        passed: 95,
                        failed: 5,
                        skipped: 0,
                        total: 100,
                        coverage_percentage: 87.5,
                        details: Vec::new(),
                    },
                );
                results
            },
            performance: vec![crate::reporting::PerformanceReport {
                benchmark_name: "Matrix Multiplication".to_string(),
                metrics: crate::reporting::PerformanceMetrics {
                    mean: 0.001,
                    std_dev: 0.0001,
                    min: 0.0009,
                    max: 0.0012,
                    median: 0.001,
                    samples: 10,
                },
                regression_detected: Some(-2.5),
                baseline_comparison: None,
            }],
            code_quality: CodeQualityReport {
                lines_of_code: 15420,
                test_coverage: 87.5,
                documentation_coverage: 73.2,
                clippy_warnings: 3,
                compiler_errors: 0,
                cyclomatic_complexity: 2.1,
                maintainability_index: 78.5,
            },
            recommendations: vec![
                "Increase test coverage to >90%".to_string(),
                "Add performance regression monitoring".to_string(),
            ],
        };

        let reporter = MarkdownReporter::new();
        let markdown = reporter.generate_report(&report).unwrap();

        // Verify key sections are present
        assert!(markdown.contains("# CFD Validation Report"));
        assert!(markdown.contains("Passed:** 95"));
        assert!(markdown.contains("Failed:** 5"));
        assert!(markdown.contains("Matrix Multiplication"));
        assert!(markdown.contains("Recommendations"));
        assert!(markdown.contains("87.5%"));
    }
}
