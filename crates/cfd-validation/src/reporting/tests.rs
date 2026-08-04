use super::*;
use std::collections::HashMap;
use std::time::SystemTime;

#[test]
fn test_report_builder() {
    let report = ReportBuilder::new("Test Report".to_string())
        .with_summary(ValidationSummary {
            total_tests: 10,
            passed_tests: 9,
            failed_tests: 1,
            skipped_tests: 0,
            total_duration: std::time::Duration::from_secs(5),
            coverage_percentage: 85.0,
        })
        .add_recommendation("Fix failing test".to_string())
        .build();

    assert_eq!(report.title, "Test Report");
    assert_eq!(report.summary.passed_tests, 9);
    assert_eq!(report.summary.failed_tests, 1);
    assert_eq!(report.recommendations.len(), 1);
}

#[test]
fn test_health_score() {
    let report = ValidationReport {
        timestamp: SystemTime::now(),
        title: "Test".to_string(),
        summary: ValidationSummary {
            total_tests: 10,
            passed_tests: 8,
            failed_tests: 2,
            skipped_tests: 0,
            total_duration: std::time::Duration::from_secs(1),
            coverage_percentage: 80.0,
        },
        test_results: HashMap::new(),
        performance: Vec::new(),
        code_quality: CodeQualityReport {
            test_coverage: 80.0,
            clippy_warnings: 5,
            ..Default::default()
        },
        recommendations: Vec::new(),
    };

    let score = report.health_score();
    assert!((0.0..=1.0).contains(&score));

    let issues = report.critical_issues();
    assert!(!issues.is_empty());
}
