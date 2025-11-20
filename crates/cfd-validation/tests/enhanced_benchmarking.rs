//! Enhanced benchmarking and performance analysis tests
//!
//! Tests advanced benchmarking features including:
//! - Performance trend analysis
//! - Regression detection
//! - Statistical analysis
//! - Automated performance reporting

use cfd_validation::benchmarking::{
    analysis::{
        PerformanceAnalyzer, PerformanceReport, PerformanceTrend, RegressionAlert,
        RegressionConfig, TrendType,
    },
    suite::{BenchmarkConfig, BenchmarkSuite},
    PerformanceMetrics,
};
use std::time::Duration;

/// Test performance trend analysis
#[test]
fn test_performance_trend_analysis() {
    println!("Testing Performance Trend Analysis...");

    let mut analyzer = PerformanceAnalyzer::with_default_config();

    // Simulate improving performance over time
    let improving_data = vec![
        PerformanceMetrics {
            mean: 1.0,
            std_dev: 0.05,
            min: 0.95,
            max: 1.05,
            median: 1.0,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.95,
            std_dev: 0.05,
            min: 0.90,
            max: 1.00,
            median: 0.95,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.90,
            std_dev: 0.05,
            min: 0.85,
            max: 0.95,
            median: 0.90,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.85,
            std_dev: 0.05,
            min: 0.80,
            max: 0.90,
            median: 0.85,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.80,
            std_dev: 0.05,
            min: 0.75,
            max: 0.85,
            median: 0.80,
            samples: 10,
        },
    ];

    for metric in improving_data {
        analyzer.add_result("improving_benchmark", metric);
    }

    let trend = analyzer.analyze_trend("improving_benchmark").unwrap();

    assert_eq!(trend.trend_type, TrendType::Improving);
    assert!(
        trend.slope < 0.0,
        "Slope should be negative for improvement"
    );
    assert!(
        trend.r_squared > 0.8,
        "Should have good correlation: R² = {}",
        trend.r_squared
    );
    assert!(
        trend.p_value < 0.05,
        "Should be statistically significant: p = {}",
        trend.p_value
    );

    println!(
        "✓ Performance trend analysis passed (slope: {:.4}, R²: {:.3})",
        trend.slope, trend.r_squared
    );
}

/// Test performance regression detection
#[test]
fn test_regression_detection() {
    println!("Testing Performance Regression Detection...");

    let config = RegressionConfig {
        degradation_threshold: 2.0, // 2% threshold
        min_samples: 5,
        confidence_level: 0.95,
        lookback_window: 10,
    };

    let mut analyzer = PerformanceAnalyzer::new(config);

    // Simulate degrading performance
    let degrading_data = vec![
        PerformanceMetrics {
            mean: 1.0,
            std_dev: 0.01,
            min: 0.99,
            max: 1.01,
            median: 1.0,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.025,
            std_dev: 0.01,
            min: 1.015,
            max: 1.035,
            median: 1.025,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.050,
            std_dev: 0.01,
            min: 1.040,
            max: 1.060,
            median: 1.050,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.075,
            std_dev: 0.01,
            min: 1.065,
            max: 1.085,
            median: 1.075,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.100,
            std_dev: 0.01,
            min: 1.090,
            max: 1.110,
            median: 1.100,
            samples: 10,
        },
    ];

    for metric in degrading_data {
        analyzer.add_result("degrading_benchmark", metric);
    }

    let regression = analyzer.detect_regression("degrading_benchmark").unwrap();

    assert!(regression.is_some(), "Should detect performance regression");

    let alert = regression.unwrap();
    assert!(
        alert.degradation_rate > 2.0,
        "Degradation rate should exceed threshold: {:.2}%",
        alert.degradation_rate
    );
    assert!(
        alert.confidence > 0.9,
        "Should have high confidence: {:.2}",
        alert.confidence
    );
    assert_eq!(alert.trend.trend_type, TrendType::Degrading);

    println!(
        "✓ Regression detection passed (degradation: {:.2}%, confidence: {:.1}%)",
        alert.degradation_rate,
        alert.confidence * 100.0
    );
}

/// Test stable performance analysis
#[test]
fn test_stable_performance_analysis() {
    println!("Testing Stable Performance Analysis...");

    let mut analyzer = PerformanceAnalyzer::with_default_config();

    // Simulate stable performance with some noise
    let stable_data = vec![
        PerformanceMetrics {
            mean: 1.00,
            std_dev: 0.02,
            min: 0.98,
            max: 1.02,
            median: 1.00,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.01,
            std_dev: 0.02,
            min: 0.99,
            max: 1.03,
            median: 1.01,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.99,
            std_dev: 0.02,
            min: 0.97,
            max: 1.01,
            median: 0.99,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.02,
            std_dev: 0.02,
            min: 1.00,
            max: 1.04,
            median: 1.02,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.00,
            std_dev: 0.02,
            min: 0.98,
            max: 1.02,
            median: 1.00,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.01,
            std_dev: 0.02,
            min: 0.99,
            max: 1.03,
            median: 1.01,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 0.99,
            std_dev: 0.02,
            min: 0.97,
            max: 1.01,
            median: 0.99,
            samples: 10,
        },
    ];

    for metric in stable_data {
        analyzer.add_result("stable_benchmark", metric);
    }

    let trend = analyzer.analyze_trend("stable_benchmark").unwrap();

    assert_eq!(trend.trend_type, TrendType::Stable);
    assert!(
        trend.slope.abs() < 0.01,
        "Slope should be near zero for stable performance: {}",
        trend.slope
    );

    // Should not detect regression for stable performance
    let regression = analyzer.detect_regression("stable_benchmark").unwrap();
    assert!(
        regression.is_none(),
        "Should not detect regression for stable performance"
    );

    println!(
        "✓ Stable performance analysis passed (slope: {:.6})",
        trend.slope
    );
}

/// Test comprehensive performance reporting
#[test]
fn test_comprehensive_performance_reporting() {
    println!("Testing Comprehensive Performance Reporting...");

    let config = BenchmarkConfig {
        iterations: 2,
        enable_memory: true,
        enable_scaling: false,
        detailed_reporting: true,
        ..Default::default()
    };

    let suite = BenchmarkSuite::with_config(config);
    let results = suite.run_full_suite();

    assert!(
        results.is_ok(),
        "Benchmark suite should execute successfully"
    );

    let results = results.unwrap();
    assert!(!results.is_empty(), "Should have benchmark results");

    // Create analyzer and generate reports
    let mut analyzer = PerformanceAnalyzer::with_default_config();

    // Add current results to analyzer
    for result in &results {
        analyzer.add_result(&result.name, result.performance.clone());
    }

    let reports = analyzer.generate_report(&results).unwrap();

    assert_eq!(
        reports.len(),
        results.len(),
        "Should have report for each result"
    );

    for report in reports {
        // Each report should have current metrics
        assert!(
            report.current_metrics.mean > 0.0,
            "Should have valid current metrics"
        );

        // Check recommendations are generated
        assert!(
            !report.recommendations.is_empty(),
            "Should have recommendations for {}",
            report.benchmark_name
        );

        println!(
            "  📊 {}: {:.3}ms ± {:.3}ms ({} recommendations)",
            report.benchmark_name,
            report.current_metrics.mean * 1000.0,
            report.current_metrics.std_dev * 1000.0,
            report.recommendations.len()
        );
    }

    println!("✓ Comprehensive performance reporting passed");
}

/// Test statistical analysis robustness
#[test]
fn test_statistical_analysis_robustness() {
    println!("Testing Statistical Analysis Robustness...");

    let mut analyzer = PerformanceAnalyzer::new(RegressionConfig {
        degradation_threshold: 10.0, // High threshold to avoid false positives
        min_samples: 3,
        confidence_level: 0.95,
        lookback_window: 20,
    });

    // Test with minimal data
    analyzer.add_result(
        "minimal_test",
        PerformanceMetrics {
            mean: 1.0,
            std_dev: 0.1,
            min: 0.9,
            max: 1.1,
            median: 1.0,
            samples: 10,
        },
    );

    let trend_result = analyzer.analyze_trend("minimal_test");
    assert!(trend_result.is_err(), "Should reject insufficient data");

    // Add more data
    analyzer.add_result(
        "minimal_test",
        PerformanceMetrics {
            mean: 1.05,
            std_dev: 0.1,
            min: 0.95,
            max: 1.15,
            median: 1.05,
            samples: 10,
        },
    );
    analyzer.add_result(
        "minimal_test",
        PerformanceMetrics {
            mean: 1.10,
            std_dev: 0.1,
            min: 1.0,
            max: 1.20,
            median: 1.10,
            samples: 10,
        },
    );

    let trend = analyzer.analyze_trend("minimal_test").unwrap();
    assert!(
        trend.r_squared >= 0.0 && trend.r_squared <= 1.0,
        "R-squared should be valid"
    );
    assert!(
        trend.p_value >= 0.0 && trend.p_value <= 1.0,
        "P-value should be valid"
    );

    println!("✓ Statistical analysis robustness passed");
}

/// Test performance metrics calculations
#[test]
fn test_performance_metrics_calculations() {
    println!("Testing Performance Metrics Calculations...");

    let metrics = PerformanceMetrics {
        mean: 1.0,
        std_dev: 0.2,
        min: 0.8,
        max: 1.2,
        median: 1.0,
        samples: 10,
    };

    // Test coefficient of variation
    let cv = metrics.coefficient_of_variation();
    assert!(
        (cv - 0.2).abs() < 1e-10,
        "Coefficient of variation should be 0.2: {}",
        cv
    );

    // Test stability
    assert!(
        metrics.is_stable(0.25),
        "Should be stable with CV < threshold"
    );
    assert!(
        !metrics.is_stable(0.15),
        "Should not be stable with CV > threshold"
    );

    // Test edge cases
    let zero_std = PerformanceMetrics {
        mean: 1.0,
        std_dev: 0.0,
        min: 1.0,
        max: 1.0,
        median: 1.0,
        samples: 1,
    };

    assert_eq!(
        zero_std.coefficient_of_variation(),
        0.0,
        "Zero std dev should give CV = 0"
    );

    let zero_mean = PerformanceMetrics {
        mean: 0.0,
        std_dev: 1.0,
        min: -1.0,
        max: 1.0,
        median: 0.0,
        samples: 10,
    };

    // CV is undefined for zero mean, should return 0
    assert_eq!(
        zero_mean.coefficient_of_variation(),
        0.0,
        "Zero mean should return CV = 0"
    );

    println!("✓ Performance metrics calculations passed");
}

/// Integration test for performance analysis pipeline
#[test]
fn test_performance_analysis_pipeline() {
    println!("🧪 Testing Performance Analysis Pipeline Integration...");

    // Create analyzer with custom config
    let config = RegressionConfig {
        degradation_threshold: 5.0,
        min_samples: 4,
        confidence_level: 0.90,
        lookback_window: 15,
    };

    let mut analyzer = PerformanceAnalyzer::new(config);

    // Simulate a realistic benchmark history with gradual improvement
    let history = vec![
        // Initial performance (baseline)
        PerformanceMetrics {
            mean: 2.0,
            std_dev: 0.1,
            min: 1.9,
            max: 2.1,
            median: 2.0,
            samples: 10,
        },
        // Some degradation (possible regression)
        PerformanceMetrics {
            mean: 2.1,
            std_dev: 0.1,
            min: 2.0,
            max: 2.2,
            median: 2.1,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 2.15,
            std_dev: 0.1,
            min: 2.05,
            max: 2.25,
            median: 2.15,
            samples: 10,
        },
        // Recovery and improvement
        PerformanceMetrics {
            mean: 2.0,
            std_dev: 0.1,
            min: 1.9,
            max: 2.1,
            median: 2.0,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.9,
            std_dev: 0.1,
            min: 1.8,
            max: 2.0,
            median: 1.9,
            samples: 10,
        },
        PerformanceMetrics {
            mean: 1.8,
            std_dev: 0.1,
            min: 1.7,
            max: 1.9,
            median: 1.8,
            samples: 10,
        },
        // Stable improved performance
        PerformanceMetrics {
            mean: 1.75,
            std_dev: 0.05,
            min: 1.7,
            max: 1.8,
            median: 1.75,
            samples: 10,
        },
    ];

    for (i, metric) in history.iter().enumerate() {
        analyzer.add_result("cfd_solver_benchmark", *metric);
        println!(
            "  Added data point {}: {:.3}ms ± {:.3}ms",
            i + 1,
            metric.mean * 1000.0,
            metric.std_dev * 1000.0
        );
    }

    // Analyze trend
    let trend = analyzer.analyze_trend("cfd_solver_benchmark").unwrap();
    println!(
        "  📈 Trend Analysis: slope={:.6}, R²={:.3}, p={:.4}, type={:?}",
        trend.slope, trend.r_squared, trend.p_value, trend.trend_type
    );

    // Check for regressions
    let regression = analyzer.detect_regression("cfd_solver_benchmark").unwrap();
    if let Some(alert) = regression {
        println!(
            "  🚨 Regression Alert: {:.2}% degradation (confidence: {:.1}%)",
            alert.degradation_rate,
            alert.confidence * 100.0
        );
    } else {
        println!("  ✅ No performance regression detected");
    }

    // Run benchmark suite and generate reports
    let benchmark_config = BenchmarkConfig {
        iterations: 1, // Quick test
        enable_memory: false,
        enable_scaling: false,
        detailed_reporting: false,
        ..Default::default()
    };

    let suite = BenchmarkSuite::new(benchmark_config);
    if let Ok(results) = suite.run_full_suite() {
        let reports = analyzer.generate_report(&results).unwrap();

        println!("  📊 Generated {} performance reports:", reports.len());
        for report in &reports {
            let perf_ms = report.current_metrics.mean * 1000.0;
            let std_ms = report.current_metrics.std_dev * 1000.0;
            println!(
                "    • {}: {:.2}ms ± {:.2}ms",
                report.benchmark_name, perf_ms, std_ms
            );

            if !report.recommendations.is_empty() {
                println!("      Recommendations: {}", report.recommendations.len());
            }
        }
    }

    println!("✅ Performance analysis pipeline integration passed");
}

/// Property-based tests for performance analysis
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Test trend analysis with various performance data
        #[test]
        fn test_trend_analysis_properties(
            means in proptest::collection::vec(0.1f64..10.0, 5..20),
            std_devs in proptest::collection::vec(0.01f64..1.0, 5..20)
        ) {
            prop_assume!(means.len() == std_devs.len());
            prop_assume!(means.len() >= 5);

            let mut analyzer = PerformanceAnalyzer::with_default_config();

            // Create metrics with varying means and constant std dev
            for (i, (&mean, &std_dev)) in means.iter().zip(std_devs.iter()).enumerate() {
                let metric = PerformanceMetrics {
                    mean,
                    std_dev,
                    min: mean - std_dev,
                    max: mean + std_dev,
                    median: mean,
                    samples: 10,
                };
                analyzer.add_result("property_test", metric);
            }

            let trend_result = analyzer.analyze_trend("property_test");

            // Should either succeed or fail with insufficient data
            if trend_result.is_ok() {
                let trend = trend_result.unwrap();
                prop_assert!(trend.r_squared >= 0.0 && trend.r_squared <= 1.0);
                prop_assert!(trend.p_value >= 0.0 && trend.p_value <= 1.0);
                prop_assert!(matches!(trend.trend_type,
                    TrendType::Improving | TrendType::Degrading | TrendType::Stable | TrendType::Volatile));
            }
        }

        /// Test regression detection sensitivity
        #[test]
        fn test_regression_detection_properties(
            threshold in 1.0f64..20.0,
            confidence in 0.8f64..0.99
        ) {
            let config = RegressionConfig {
                degradation_threshold: threshold,
                min_samples: 5,
                confidence_level: confidence,
                lookback_window: 10,
            };

            let mut analyzer = PerformanceAnalyzer::new(config);

            // Add stable data
            for i in 0..6 {
                let metric = PerformanceMetrics {
                    mean: 1.0,
                    std_dev: 0.1,
                    min: 0.9,
                    max: 1.1,
                    median: 1.0,
                    samples: 10,
                };
                analyzer.add_result("stable_test", metric);
            }

            let regression = analyzer.detect_regression("stable_test").unwrap();
            prop_assert!(regression.is_none(), "Stable data should not trigger regression");
        }
    }
}
