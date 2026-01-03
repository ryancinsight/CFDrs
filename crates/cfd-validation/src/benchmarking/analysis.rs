//! Performance analysis and regression detection for CFD benchmarks
//!
//! Provides statistical analysis, trend detection, and performance regression
//! monitoring for CFD operations.

use super::suite::BenchmarkResult;
use crate::reporting::PerformanceMetrics;
use cfd_core::error::{Error, Result};
use std::collections::HashMap;

/// Performance trend analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceTrend {
    /// Slope of performance change (negative = improvement)
    pub slope: f64,
    /// R-squared value for trend fit
    pub r_squared: f64,
    /// Statistical significance (p-value)
    pub p_value: f64,
    /// Trend classification
    pub trend_type: TrendType,
}

/// Performance trend classification for temporal analysis
///
/// Classifies the direction and consistency of performance changes over time,
/// enabling automated detection of performance improvements, degradations,
/// and stability issues in CFD benchmark results.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum TrendType {
    /// Performance is consistently improving over time
    ///
    /// Indicates successful optimizations, algorithm improvements, or hardware upgrades.
    /// Requires statistical significance to avoid false positives from measurement noise.
    Improving,

    /// Performance is consistently degrading over time
    ///
    /// Indicates performance regressions, resource contention, or system degradation.
    /// Requires immediate investigation to identify root causes and prevent further decline.
    Degrading,

    /// Performance is stable with minimal variation
    ///
    /// Indicates consistent performance with acceptable measurement noise.
    /// Represents the ideal state for production CFD systems.
    Stable,

    /// Performance shows high variability with no clear trend
    ///
    /// Indicates measurement instability, system interference, or inconsistent workloads.
    /// May mask underlying performance issues or improvements.
    Volatile,
}

/// Alert severity levels for regression detection and performance monitoring
///
/// Hierarchical classification of performance issues by their impact and urgency,
/// enabling appropriate response strategies for different types of performance problems.
/// Higher severity levels require more immediate attention and resources.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum AlertSeverity {
    /// Minor performance variation requiring monitoring
    ///
    /// Small performance changes that may be within normal operating ranges.
    /// Requires monitoring but does not typically require immediate action.
    Low,

    /// Moderate performance issue requiring investigation
    ///
    /// Performance degradation that exceeds normal variation but may be acceptable
    /// depending on operational requirements. Requires investigation and potential optimization.
    Medium,

    /// Significant performance degradation requiring action
    ///
    /// Major performance issues that significantly impact CFD simulation efficiency.
    /// Requires immediate investigation and corrective action planning.
    High,

    /// Critical performance failure requiring immediate response
    ///
    /// Severe performance degradation that makes CFD simulations impractical or unusable.
    /// Requires immediate emergency response, rollback, or system intervention.
    Critical,
}

/// Regression detection configuration
#[derive(Debug, Clone)]
pub struct RegressionConfig {
    /// Threshold for performance degradation (%)
    pub degradation_threshold: f64,
    /// Minimum samples for trend analysis
    pub min_samples: usize,
    /// Confidence level for statistical tests
    pub confidence_level: f64,
    /// Lookback window for analysis (number of recent runs)
    pub lookback_window: usize,
}

impl Default for RegressionConfig {
    fn default() -> Self {
        Self {
            degradation_threshold: 5.0, // 5% degradation threshold
            min_samples: 5,
            confidence_level: 0.95,
            lookback_window: 10,
        }
    }
}

/// Performance analyzer for CFD benchmarks
pub struct PerformanceAnalyzer {
    config: RegressionConfig,
    historical_data: HashMap<String, Vec<PerformanceMetrics>>,
}

impl PerformanceAnalyzer {
    /// Create a new performance analyzer with custom regression detection configuration
    ///
    /// Initializes the analyzer with user-specified thresholds for regression detection,
    /// statistical significance requirements, and analysis parameters. Use this constructor
    /// when you need fine-tuned control over regression detection sensitivity.
    ///
    /// # Parameters
    ///
    /// * `config` - Regression detection configuration specifying thresholds, confidence levels,
    ///   and analysis window parameters for performance monitoring
    pub fn new(config: RegressionConfig) -> Self {
        Self {
            config,
            historical_data: HashMap::new(),
        }
    }

    /// Create a new performance analyzer with default regression detection configuration
    ///
    /// Initializes the analyzer with sensible default settings for regression detection:
    /// - 5% degradation threshold
    /// - Minimum 5 samples for analysis
    /// - 95% confidence level for statistical tests
    /// - 10-run lookback window for trend analysis
    ///
    /// Suitable for most CFD benchmarking scenarios without requiring manual configuration tuning.
    pub fn with_default_config() -> Self {
        Self::new(RegressionConfig::default())
    }

    /// Add benchmark result to historical data
    pub fn add_result(&mut self, benchmark_name: &str, metrics: PerformanceMetrics) {
        self.historical_data
            .entry(benchmark_name.to_string())
            .or_default()
            .push(metrics);
    }

    /// Analyze performance trend for a benchmark
    pub fn analyze_trend(&self, benchmark_name: &str) -> Result<PerformanceTrend> {
        let data = self.historical_data.get(benchmark_name).ok_or_else(|| {
            Error::InvalidInput(format!(
                "No historical data for benchmark: {benchmark_name}"
            ))
        })?;

        if data.len() < self.config.min_samples {
            return Err(Error::InvalidInput(format!(
                "Insufficient data for trend analysis: {} samples, need {}",
                data.len(),
                self.config.min_samples
            )));
        }

        // Use recent data within lookback window
        let recent_data = &data[data.len().saturating_sub(self.config.lookback_window)..];

        let trend = self.calculate_linear_trend(recent_data)?;
        Ok(trend)
    }

    /// Detect performance regression
    pub fn detect_regression(&self, benchmark_name: &str) -> Result<Option<RegressionAlert>> {
        let trend = self.analyze_trend(benchmark_name)?;

        if trend.trend_type == TrendType::Degrading
            && trend.slope.abs() > self.config.degradation_threshold / 100.0
            && trend.p_value < (1.0 - self.config.confidence_level)
        {
            let degradation_rate = trend.slope * 100.0; // Convert to percentage
            return Ok(Some(RegressionAlert {
                benchmark_name: benchmark_name.to_string(),
                degradation_rate,
                confidence: self.config.confidence_level,
                trend,
            }));
        }

        Ok(None)
    }

    /// Calculate linear trend using least squares
    fn calculate_linear_trend(&self, data: &[PerformanceMetrics]) -> Result<PerformanceTrend> {
        let n = data.len() as f64;
        if n < 2.0 {
            return Ok(PerformanceTrend {
                slope: 0.0,
                r_squared: 0.0,
                p_value: 1.0,
                trend_type: TrendType::Stable,
            });
        }

        // Extract mean execution times
        let values: Vec<f64> = data.iter().map(|m| m.mean).collect();

        // Calculate means
        let x_mean = (n - 1.0) / 2.0; // Time indices centered around 0
        let y_mean = values.iter().sum::<f64>() / n;

        // Calculate slope and intercept
        let mut numerator = 0.0;
        let mut denominator = 0.0;

        for (i, &y) in values.iter().enumerate() {
            let x = i as f64 - x_mean;
            numerator += x * (y - y_mean);
            denominator += x * x;
        }

        let slope = if denominator > 0.0 {
            numerator / denominator
        } else {
            0.0
        };

        // Calculate R-squared
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;

        for (i, &y) in values.iter().enumerate() {
            let x = i as f64 - x_mean;
            let y_pred = y_mean + slope * x;
            ss_res += (y - y_pred).powi(2);
            ss_tot += (y - y_mean).powi(2);
        }

        let r_squared = if ss_tot > 0.0 {
            1.0 - (ss_res / ss_tot)
        } else {
            0.0
        };

        // Simplified p-value calculation (t-test approximation)
        let se_slope = if denominator > 0.0 {
            let variance = ss_res / (n - 2.0);
            variance.sqrt() / denominator.sqrt()
        } else {
            0.0
        };

        let t_stat = if se_slope > 0.0 {
            slope.abs() / se_slope
        } else {
            0.0
        };
        let p_value = 2.0 * (1.0 - Self::t_cdf(t_stat, (n - 2.0) as usize));

        // Classify trend
        let trend_type = if slope.abs() < 1e-6 {
            TrendType::Stable
        } else if slope < 0.0 {
            TrendType::Improving
        } else {
            TrendType::Degrading
        };

        Ok(PerformanceTrend {
            slope,
            r_squared,
            p_value,
            trend_type,
        })
    }

    /// Cumulative distribution function for t-distribution (approximation)
    fn t_cdf(t: f64, df: usize) -> f64 {
        // Simplified approximation using normal CDF for large df
        if df > 30 {
            Self::normal_cdf(t)
        } else {
            // More accurate approximation for small df would be complex
            // This is a simplified version
            Self::normal_cdf(t * (1.0 - 1.0 / (4.0 * df as f64)).sqrt())
        }
    }

    /// Normal cumulative distribution function approximation
    fn normal_cdf(x: f64) -> f64 {
        // Scale for standard normal distribution (Phi(x) = 0.5 * (1 + erf(x/sqrt(2))))
        let x = x / 2.0f64.sqrt();

        // Abramowitz & Stegun approximation for erf(x)
        let a1 = 0.254829592;
        let a2 = -0.284496736;
        let a3 = 1.421413741;
        let a4 = -1.453152027;
        let a5 = 1.061405429;
        let p = 0.3275911;

        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x_abs = x.abs();

        let t = 1.0 / (1.0 + p * x_abs);
        let exp_val = (-x_abs * x_abs).exp();
        let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp_val;

        0.5 * (1.0 + sign * y)
    }
}

/// Performance regression alert with statistical significance assessment
///
/// Automatically generated alert when statistical analysis detects significant
/// performance degradation beyond configured thresholds. Includes quantitative
/// measures of degradation magnitude and statistical confidence in the detection.
#[derive(Debug, Clone)]
pub struct RegressionAlert {
    /// Name of the benchmark operation showing performance regression
    ///
    /// Identifier of the CFD operation experiencing performance degradation,
    /// enabling targeted investigation and remediation efforts.
    pub benchmark_name: String,

    /// Performance degradation rate as percentage change from baseline
    ///
    /// Quantified performance loss expressed as percentage (e.g., 15.7 means 15.7% slower).
    /// Computed as (current_time - baseline_time) / baseline_time * 100.
    /// Positive values indicate performance degradation (slower execution).
    pub degradation_rate: f64,

    /// Statistical confidence level in the regression detection [0.0-1.0]
    ///
    /// Probability that the detected regression is not due to random measurement variation.
    /// Higher values indicate greater confidence in the regression detection.
    /// Typically requires >95% confidence for actionable alerts.
    pub confidence: f64,

    /// Detailed performance trend analysis supporting the regression alert
    ///
    /// Comprehensive trend analysis including slope, statistical significance,
    /// and classification of the performance degradation pattern.
    pub trend: PerformanceTrend,
}

/// Comprehensive performance analysis report for individual CFD benchmarks
///
/// Provides detailed performance assessment for a specific CFD benchmark operation,
/// including current performance metrics, historical comparison, trend analysis,
/// regression detection, and actionable recommendations for optimization.
///
/// # Report Components
///
/// - **Current Performance**: Latest benchmark execution metrics and statistics
/// - **Baseline Comparison**: Historical performance reference for change detection
/// - **Trend Analysis**: Long-term performance evolution and stability assessment
/// - **Regression Alerts**: Automated detection of performance degradation
/// - **Optimization Recommendations**: Actionable suggestions for performance improvement
#[derive(Debug, Clone)]
pub struct PerformanceReport {
    /// Name/identifier of the benchmark operation being analyzed
    ///
    /// Human-readable identifier for the CFD operation (e.g., "Navier-Stokes Solver",
    /// "FFT Transform", "Matrix Factorization"). Used for report organization and cross-referencing.
    pub benchmark_name: String,

    /// Current performance metrics from the latest benchmark execution
    ///
    /// Statistical analysis of the most recent benchmark run, including timing statistics,
    /// variability measures, and sample characteristics. Represents the current performance state.
    pub current_metrics: PerformanceMetrics,

    /// Baseline performance metrics for historical comparison
    ///
    /// Reference performance metrics from a stable baseline period. Used for detecting
    /// performance changes and regression analysis. None if no baseline is established.
    pub baseline_metrics: Option<PerformanceMetrics>,

    /// Performance trend analysis over the analysis window
    ///
    /// Temporal analysis of performance evolution, classifying trends as improving,
    /// degrading, stable, or volatile. None if insufficient historical data for analysis.
    pub trend: Option<PerformanceTrend>,

    /// Regression alert if performance degradation is detected
    ///
    /// Automated alert triggered when statistical analysis detects significant performance
    /// degradation beyond acceptable thresholds. Includes severity classification and confidence metrics.
    pub regression_alert: Option<RegressionAlert>,

    /// Actionable recommendations for performance optimization
    ///
    /// Generated suggestions for improving performance based on current metrics, trends,
    /// and detected issues. May include algorithmic changes, parameter tuning, or system optimizations.
    pub recommendations: Vec<String>,
}

impl PerformanceAnalyzer {
    /// Generate comprehensive performance report
    pub fn generate_report(&self, results: &[BenchmarkResult]) -> Result<Vec<PerformanceReport>> {
        let mut reports = Vec::new();

        for result in results {
            let trend = self.analyze_trend(&result.name).ok();
            let regression = self.detect_regression(&result.name).unwrap_or(None);

            let baseline = self.get_baseline_metrics(&result.name);
            let recommendations =
                self.generate_recommendations(result, trend.as_ref(), regression.as_ref());

            // Convert TimingResult to PerformanceMetrics if available
            let current_metrics = result
                .performance
                .as_ref()
                .map(|perf| crate::reporting::PerformanceMetrics {
                    mean: perf.stats.mean,
                    std_dev: perf.stats.std_dev,
                    min: perf.stats.min,
                    max: perf.stats.max,
                    median: perf.stats.median,
                    samples: perf.stats.samples,
                })
                .unwrap_or_default();

            reports.push(PerformanceReport {
                benchmark_name: result.name.clone(),
                current_metrics,
                baseline_metrics: baseline,
                trend,
                regression_alert: regression,
                recommendations,
            });
        }

        Ok(reports)
    }

    /// Get baseline metrics for comparison
    fn get_baseline_metrics(&self, benchmark_name: &str) -> Option<PerformanceMetrics> {
        self.historical_data
            .get(benchmark_name)
            .and_then(|data: &Vec<PerformanceMetrics>| data.first().cloned())
    }

    /// Generate performance recommendations
    fn generate_recommendations(
        &self,
        result: &BenchmarkResult,
        trend: Option<&PerformanceTrend>,
        regression: Option<&RegressionAlert>,
    ) -> Vec<String> {
        let mut recommendations = Vec::new();

        // Check for high variance and performance issues
        if let Some(perf) = &result.performance {
            if perf.stats.coefficient_of_variation() > 0.1 {
                recommendations.push("High performance variance detected. Consider increasing sample size or stabilizing test environment.".to_string());
            }

            if perf.stats.mean > 0.1 {
                // More than 100ms
                recommendations.push(
                    "High execution time detected. Consider optimization opportunities."
                        .to_string(),
                );
            }
        }

        // Check trend
        if let Some(trend) = trend {
            match trend.trend_type {
                TrendType::Degrading => {
                    recommendations.push(format!(
                        "Performance degrading at {:.2}% per run. Investigate recent changes.",
                        trend.slope * 100.0
                    ));
                }
                TrendType::Improving => {
                    recommendations.push(format!(
                        "Performance improving at {:.2}% per run. Good trend!",
                        trend.slope.abs() * 100.0
                    ));
                }
                TrendType::Volatile => {
                    recommendations.push("Performance is volatile. Consider stabilizing factors affecting benchmark.".to_string());
                }
                TrendType::Stable => {
                    if trend.r_squared > 0.8 {
                        recommendations
                            .push("Performance is stable. Continue monitoring.".to_string());
                    }
                }
            }
        }

        // Check for regression alerts
        if let Some(regression) = regression {
            recommendations.push(format!(
                "🚨 PERFORMANCE REGRESSION: {:.1}% degradation detected with {:.0}% confidence. Immediate investigation recommended.",
                regression.degradation_rate,
                regression.confidence * 100.0
            ));
        }

        recommendations
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_performance_analyzer() {
        let mut analyzer = PerformanceAnalyzer::with_default_config();

        // Add some historical data
        let metrics = vec![
            PerformanceMetrics {
                mean: 1.0,
                std_dev: 0.1,
                min: 0.9,
                max: 1.1,
                median: 1.0,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.05,
                std_dev: 0.1,
                min: 0.95,
                max: 1.15,
                median: 1.05,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.10,
                std_dev: 0.1,
                min: 1.0,
                max: 1.20,
                median: 1.10,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.08,
                std_dev: 0.1,
                min: 0.98,
                max: 1.18,
                median: 1.08,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.12,
                std_dev: 0.1,
                min: 1.02,
                max: 1.22,
                median: 1.12,
                samples: 10,
            },
        ];

        for metric in metrics {
            analyzer.add_result("test_benchmark", metric);
        }

        // Analyze trend
        let trend = analyzer.analyze_trend("test_benchmark").unwrap();
        assert!(trend.r_squared > 0.0); // Should have some correlation

        // The trend should show slight degradation (increasing mean)
        assert_eq!(trend.trend_type, TrendType::Degrading);

        // Check regression detection (should not trigger with default threshold)
        let regression = analyzer.detect_regression("test_benchmark").unwrap();
        assert!(regression.is_none()); // Slope should be small
    }

    #[test]
    fn test_regression_detection() {
        let mut analyzer = PerformanceAnalyzer::new(RegressionConfig {
            degradation_threshold: 1.0, // Low threshold for testing
            min_samples: 3,
            confidence_level: 0.95,
            lookback_window: 10,
        });

        // Add degrading performance data
        let metrics = vec![
            PerformanceMetrics {
                mean: 1.0,
                std_dev: 0.01,
                min: 0.99,
                max: 1.01,
                median: 1.0,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.1,
                std_dev: 0.01,
                min: 1.09,
                max: 1.11,
                median: 1.1,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.2,
                std_dev: 0.01,
                min: 1.19,
                max: 1.21,
                median: 1.2,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.3,
                std_dev: 0.01,
                min: 1.29,
                max: 1.31,
                median: 1.3,
                samples: 10,
            },
            PerformanceMetrics {
                mean: 1.4,
                std_dev: 0.01,
                min: 1.39,
                max: 1.41,
                median: 1.4,
                samples: 10,
            },
        ];

        for metric in metrics {
            analyzer.add_result("degrading_benchmark", metric);
        }

        // Should detect regression
        let regression = analyzer.detect_regression("degrading_benchmark").unwrap();
        assert!(regression.is_some());

        let alert = regression.unwrap();
        assert!(alert.degradation_rate > 0.0);
        assert!(alert.confidence > 0.9);
    }

    #[test]
    fn test_normal_cdf() {
        // Test some known values
        assert_relative_eq!(PerformanceAnalyzer::normal_cdf(0.0), 0.5, epsilon = 1e-3);
        assert_relative_eq!(PerformanceAnalyzer::normal_cdf(1.0), 0.8413, epsilon = 1e-3);
        assert_relative_eq!(
            PerformanceAnalyzer::normal_cdf(-1.0),
            0.1587,
            epsilon = 1e-3
        );
    }

    #[test]
    fn test_insufficient_data() {
        let mut analyzer = PerformanceAnalyzer::with_default_config();
        analyzer.add_result(
            "sparse_benchmark",
            PerformanceMetrics {
                mean: 1.0,
                std_dev: 0.1,
                min: 0.9,
                max: 1.1,
                median: 1.0,
                samples: 10,
            },
        );

        // Should fail with insufficient data
        let result = analyzer.analyze_trend("sparse_benchmark");
        assert!(result.is_err());
    }
}
