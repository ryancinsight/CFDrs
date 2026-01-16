//! HTML report generation for validation results

use super::{Reporter, ValidationReport};
use cfd_core::error::Result;
use std::collections::HashMap;
use std::time::SystemTime;

/// HTML reporter for validation results
pub struct HtmlReporter;

impl HtmlReporter {
    /// Create a new HTML reporter
    pub fn new() -> Self {
        Self
    }

    /// Generate HTML header with CSS styling
    fn generate_header(&self, report: &ValidationReport) -> String {
        format!(
            r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            margin-bottom: 40px;
            padding-bottom: 20px;
            border-bottom: 2px solid #e9ecef;
        }}
        .header h1 {{
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 2.5em;
        }}
        .timestamp {{
            color: #6c757d;
            font-size: 0.9em;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 40px;
        }}
        .summary-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}
        .summary-card.success {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }}
        .summary-card.warning {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        .summary-card.danger {{
            background: linear-gradient(135deg, #eb3349 0%, #f45c43 100%);
        }}
        .summary-number {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 5px;
        }}
        .summary-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        .section {{
            margin-bottom: 40px;
        }}
        .section h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .test-category {{
            margin-bottom: 30px;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            overflow: hidden;
        }}
        .category-header {{
            background: #f8f9fa;
            padding: 15px 20px;
            border-bottom: 1px solid #dee2e6;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        .category-name {{
            font-weight: bold;
            color: #495057;
        }}
        .category-stats {{
            display: flex;
            gap: 15px;
            font-size: 0.9em;
        }}
        .stat {{
            padding: 4px 8px;
            border-radius: 4px;
            color: white;
        }}
        .stat.passed {{ background: #28a745; }}
        .stat.failed {{ background: #dc3545; }}
        .stat.skipped {{ background: #6c757d; }}
        .test-details {{
            padding: 0;
        }}
        .test-row {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 12px 20px;
            border-bottom: 1px solid #f1f3f4;
        }}
        .test-row:last-child {{
            border-bottom: none;
        }}
        .test-name {{
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            flex: 1;
        }}
        .test-status {{
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.8em;
            font-weight: bold;
            text-transform: uppercase;
        }}
        .test-status.passed {{
            background: #d4edda;
            color: #155724;
        }}
        .test-status.failed {{
            background: #f8d7da;
            color: #721c24;
        }}
        .test-status.skipped {{
            background: #e2e3e5;
            color: #383d41;
        }}
        .test-duration {{
            color: #6c757d;
            font-size: 0.9em;
            margin-left: 15px;
        }}
        .performance-table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }}
        .performance-table th,
        .performance-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }}
        .performance-table th {{
            background: #f8f9fa;
            font-weight: bold;
            color: #495057;
        }}
        .regression {{
            color: #dc3545;
            font-weight: bold;
        }}
        .improvement {{
            color: #28a745;
            font-weight: bold;
        }}
        .quality-metrics {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
        }}
        .metric-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #3498db;
        }}
        .metric-name {{
            font-weight: bold;
            color: #495057;
            margin-bottom: 8px;
        }}
        .metric-value {{
            font-size: 1.5em;
            color: #2c3e50;
        }}
        .recommendations {{
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 8px;
            padding: 20px;
        }}
        .recommendations h3 {{
            color: #856404;
            margin-top: 0;
        }}
        .recommendations ul {{
            margin-bottom: 0;
        }}
        .recommendations li {{
            margin-bottom: 8px;
        }}
        .health-score {{
            text-align: center;
            margin-bottom: 30px;
        }}
        .health-score-value {{
            font-size: 4em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        .health-score.excellent {{ color: #28a745; }}
        .health-score.good {{ color: #17a2b8; }}
        .health-score.fair {{ color: #ffc107; }}
        .health-score.poor {{ color: #dc3545; }}
        @media (max-width: 768px) {{
            .summary-grid {{
                grid-template-columns: 1fr;
            }}
            .category-header {{
                flex-direction: column;
                align-items: flex-start;
                gap: 10px;
            }}
            .test-row {{
                flex-direction: column;
                align-items: flex-start;
                gap: 8px;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="timestamp">Generated on {timestamp}</div>
        </div>"#,
            title = report.title,
            timestamp = format_timestamp(report.timestamp)
        )
    }

    /// Generate summary section
    fn generate_summary(&self, report: &ValidationReport) -> String {
        let health_score = report.health_score();
        let health_class = match health_score {
            score if score >= 0.9 => "excellent",
            score if score >= 0.7 => "good",
            score if score >= 0.5 => "fair",
            _ => "poor",
        };

        format!(
            r#"
        <div class="health-score">
            <div class="health-score-value {health_class}">{health_score:.0%}</div>
            <div>Overall Health Score</div>
        </div>

        <div class="summary-grid">
            <div class="summary-card {passed_class}">
                <div class="summary-number">{passed}</div>
                <div class="summary-label">Tests Passed</div>
            </div>
            <div class="summary-card {failed_class}">
                <div class="summary-number">{failed}</div>
                <div class="summary-label">Tests Failed</div>
            </div>
            <div class="summary-card {skipped_class}">
                <div class="summary-number">{skipped}</div>
                <div class="summary-label">Tests Skipped</div>
            </div>
            <div class="summary-card {coverage_class}">
                <div class="summary-number">{coverage:.1}%</div>
                <div class="summary-label">Coverage</div>
            </div>
        </div>"#,
            health_score = health_score,
            health_class = health_class,
            passed = report.summary.passed_tests,
            failed = report.summary.failed_tests,
            skipped = report.summary.skipped_tests,
            coverage = report.summary.coverage_percentage,
            passed_class = if report.summary.failed_tests == 0 { "success" } else { "warning" },
            failed_class = if report.summary.failed_tests > 0 { "danger" } else { "success" },
            skipped_class = "warning",
            coverage_class = if report.summary.coverage_percentage >= 80.0 { "success" } else { "warning" }
        )
    }

    /// Generate test results section
    fn generate_test_results(&self, report: &ValidationReport) -> String {
        let mut html = String::from(
            r#"
        <div class="section">
            <h2>Test Results</h2>"#,
        );

        if report.test_results.is_empty() {
            html.push_str(
                r#"
            <p>No test results available.</p>"#,
            );
        } else {
            for (category_name, category) in &report.test_results {
                html.push_str(&format!(
                    r#"
            <div class="test-category">
                <div class="category-header">
                    <div class="category-name">{name}</div>
                    <div class="category-stats">
                        <span class="stat passed">{passed} passed</span>
                        <span class="stat failed">{failed} failed</span>
                        <span class="stat skipped">{skipped} skipped</span>
                    </div>
                </div>
                <div class="test-details">"#,
                    name = category.name,
                    passed = category.passed,
                    failed = category.failed,
                    skipped = category.skipped
                ));

                for test in &category.details {
                    let status_class = match test.status {
                        TestStatus::Passed => "passed",
                        TestStatus::Failed => "failed",
                        TestStatus::Skipped | TestStatus::Timeout => "skipped",
                    };

                    let status_text = match test.status {
                        TestStatus::Passed => "passed",
                        TestStatus::Failed => "failed",
                        TestStatus::Skipped => "skipped",
                        TestStatus::Timeout => "timeout",
                    };

                    html.push_str(&format!(
                        r#"
                    <div class="test-row">
                        <div class="test-name">{name}</div>
                        <div>
                            <span class="test-status {status_class}">{status_text}</span>
                            <span class="test-duration">{duration:.1}ms</span>
                        </div>
                    </div>"#,
                        name = test.name,
                        status_class = status_class,
                        status_text = status_text,
                        duration = test.duration_ms
                    ));

                    if let Some(error) = &test.error_message {
                        html.push_str(&format!(
                            r#"
                    <div class="test-row" style="background: #f8d7da;">
                        <div style="color: #721c24; font-family: monospace; font-size: 0.9em;">{error}</div>
                    </div>"#,
                            error = error
                        ));
                    }
                }

                html.push_str(
                    r#"
                </div>
            </div>"#,
                );
            }
        }

        html.push_str(
            r#"
        </div>"#,
        );

        html
    }

    /// Generate performance section
    fn generate_performance(&self, report: &ValidationReport) -> String {
        let mut html = String::from(
            r#"
        <div class="section">
            <h2>Performance Benchmarks</h2>"#,
        );

        if report.performance.is_empty() {
            html.push_str(
                r#"
            <p>No performance data available.</p>"#,
            );
        } else {
            html.push_str(
                r#"
            <table class="performance-table">
                <thead>
                    <tr>
                        <th>Benchmark</th>
                        <th>Mean (ms)</th>
                        <th>Std Dev</th>
                        <th>Min (ms)</th>
                        <th>Max (ms)</th>
                        <th>Regression</th>
                    </tr>
                </thead>
                <tbody>"#,
            );

            for perf in &report.performance {
                let regression = if let Some(reg) = perf.regression_detected {
                    if reg > 0.0 {
                        format!(r#"<span class="regression">+{reg:.1}%</span>"#)
                    } else {
                        format!(r#"<span class="improvement">{reg:.1}%</span>"#)
                    }
                } else {
                    "N/A".to_string()
                };

                html.push_str(&format!(
                    r#"
                    <tr>
                        <td>{name}</td>
                        <td>{mean:.2}</td>
                        <td>{std_dev:.2}</td>
                        <td>{min:.2}</td>
                        <td>{max:.2}</td>
                        <td>{regression}</td>
                    </tr>"#,
                    name = perf.benchmark_name,
                    mean = perf.metrics.mean,
                    std_dev = perf.metrics.std_dev,
                    min = perf.metrics.min,
                    max = perf.metrics.max,
                    regression = regression
                ));
            }

            html.push_str(
                r#"
                </tbody>
            </table>"#,
            );
        }

        html.push_str(
            r#"
        </div>"#,
        );

        html
    }

    /// Generate code quality section
    fn generate_code_quality(&self, report: &ValidationReport) -> String {
        format!(
            r#"
        <div class="section">
            <h2>Code Quality Metrics</h2>
            <div class="quality-metrics">
                <div class="metric-card">
                    <div class="metric-name">Lines of Code</div>
                    <div class="metric-value">{loc}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-name">Test Coverage</div>
                    <div class="metric-value">{coverage:.1}%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-name">Documentation Coverage</div>
                    <div class="metric-value">{doc_coverage:.1}%</div>
                </div>
                <div class="metric-card">
                    <div class="metric-name">Clippy Warnings</div>
                    <div class="metric-value">{warnings}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-name">Cyclomatic Complexity</div>
                    <div class="metric-value">{complexity:.1}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-name">Maintainability Index</div>
                    <div class="metric-value">{maintainability:.1}</div>
                </div>
            </div>
        </div>"#,
            loc = report.code_quality.lines_of_code,
            coverage = report.code_quality.test_coverage,
            doc_coverage = report.code_quality.documentation_coverage,
            warnings = report.code_quality.clippy_warnings,
            complexity = report.code_quality.cyclomatic_complexity,
            maintainability = report.code_quality.maintainability_index
        )
    }

    /// Generate recommendations section
    fn generate_recommendations(&self, report: &ValidationReport) -> String {
        let mut html = String::from(
            r#"
        <div class="section">
            <h2>Recommendations</h2>"#,
        );

        let critical_issues = report.critical_issues();
        if !critical_issues.is_empty() {
            html.push_str(
                r#"
            <div class="recommendations">
                <h3>Critical Issues</h3>
                <ul>"#,
            );

            for issue in critical_issues {
                html.push_str(&format!(
                    r#"
                    <li>{}</li>"#,
                    issue
                ));
            }

            html.push_str(
                r#"
                </ul>
            </div>"#,
            );
        }

        if !report.recommendations.is_empty() {
            html.push_str(
                r#"
            <div class="recommendations">
                <h3>Improvement Suggestions</h3>
                <ul>"#,
            );

            for rec in &report.recommendations {
                html.push_str(&format!(
                    r#"
                    <li>{}</li>"#,
                    rec
                ));
            }

            html.push_str(
                r#"
                </ul>
            </div>"#,
            );
        }

        html.push_str(
            r#"
        </div>"#,
        );

        html
    }

    /// Generate HTML footer
    fn generate_footer(&self) -> String {
        r#"
    </div>
</body>
</html>"#
            .to_string()
    }
}

impl Reporter for HtmlReporter {
    fn generate_report(&self, report: &ValidationReport) -> Result<String> {
        let html = format!(
            "{}{}{}{}{}{}{}",
            self.generate_header(report),
            self.generate_summary(report),
            self.generate_test_results(report),
            self.generate_performance(report),
            self.generate_code_quality(report),
            self.generate_recommendations(report),
            self.generate_footer()
        );

        Ok(html)
    }
}

/// Format timestamp for display
fn format_timestamp(timestamp: SystemTime) -> String {
    match timestamp.duration_since(std::time::UNIX_EPOCH) {
        Ok(duration) => {
            let datetime = chrono::DateTime::from_timestamp(duration.as_secs() as i64, 0);
            match datetime {
                Some(dt) => dt.format("%Y-%m-%d %H:%M:%S UTC").to_string(),
                None => "Unknown timestamp".to_string(),
            }
        }
        Err(_) => "Invalid timestamp".to_string(),
    }
}

impl Default for HtmlReporter {
    fn default() -> Self {
        Self::new()
    }
}
