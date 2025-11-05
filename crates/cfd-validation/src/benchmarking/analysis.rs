//! Performance analysis and regression detection for CFD benchmarks
//!
//! Provides statistical analysis, trend detection, and performance regression
//! monitoring for CFD operations.

use crate::reporting::PerformanceMetrics;
use super::suite::BenchmarkResult;
use super::memory::MemoryStats;
use super::performance::TimingResult;
use super::scaling::ScalingResult;
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

/// Performance trend classification
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum TrendType {
    Improving,
    Degrading,
    Stable,
    Volatile,
}

/// Alert severity levels
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum AlertSeverity {
    Low,
    Medium,
    High,
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
    pub fn new(config: RegressionConfig) -> Self {
        Self {
            config,
            historical_data: HashMap::new(),
        }
    }

    pub fn with_default_config() -> Self {
        Self::new(RegressionConfig::default())
    }

    /// Add benchmark result to historical data
    pub fn add_result(&mut self, benchmark_name: &str, metrics: PerformanceMetrics) {
        self.historical_data
            .entry(benchmark_name.to_string())
            .or_insert_with(Vec::new)
            .push(metrics);
    }

    /// Analyze performance trend for a benchmark
    pub fn analyze_trend(&self, benchmark_name: &str) -> Result<PerformanceTrend> {
        let data = self.historical_data.get(benchmark_name)
            .ok_or_else(|| Error::InvalidInput(format!("No historical data for benchmark: {}", benchmark_name)))?;

        if data.len() < self.config.min_samples {
            return Err(Error::InvalidInput(format!(
                "Insufficient data for trend analysis: {} samples, need {}",
                data.len(), self.config.min_samples
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

        if trend.trend_type == TrendType::Degrading &&
           trend.slope.abs() > self.config.degradation_threshold / 100.0 &&
           trend.p_value < (1.0 - self.config.confidence_level) {

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

        let slope = if denominator > 0.0 { numerator / denominator } else { 0.0 };

        // Calculate R-squared
        let mut ss_res = 0.0;
        let mut ss_tot = 0.0;

        for (i, &y) in values.iter().enumerate() {
            let x = i as f64 - x_mean;
            let y_pred = y_mean + slope * x;
            ss_res += (y - y_pred).powi(2);
            ss_tot += (y - y_mean).powi(2);
        }

        let r_squared = if ss_tot > 0.0 { 1.0 - (ss_res / ss_tot) } else { 0.0 };

        // Simplified p-value calculation (t-test approximation)
        let se_slope = if denominator > 0.0 {
            let variance = ss_res / (n - 2.0);
            variance.sqrt() / denominator.sqrt()
        } else {
            0.0
        };

        let t_stat = if se_slope > 0.0 { slope.abs() / se_slope } else { 0.0 };
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
        // Abramowitz & Stegun approximation
        let a1 =  0.254829592;
        let a2 = -0.284496736;
        let a3 =  1.421413741;
        let a4 = -1.453152027;
        let a5 =  1.061405429;
        let p  =  0.3275911;

        let sign = if x < 0.0 { -1.0 } else { 1.0 };
        let x = x.abs();

        let t = 1.0 / (1.0 + p * x);
        let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

        0.5 * (1.0 + sign * y)
    }
}

/// Performance regression alert
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RegressionAlert {
    pub benchmark_name: String,
    pub degradation_rate: f64, // percentage
    pub confidence: f64,
    pub trend: PerformanceTrend,
}

/// Comprehensive performance report
#[derive(Debug, Clone)]
pub struct PerformanceReport {
    pub benchmark_name: String,
    pub current_metrics: PerformanceMetrics,
    pub baseline_metrics: Option<PerformanceMetrics>,
    pub trend: Option<PerformanceTrend>,
    pub regression_alert: Option<RegressionAlert>,
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
            let recommendations = self.generate_recommendations(&result, trend.as_ref(), regression.as_ref());

            // Convert TimingResult to PerformanceMetrics if available
            let current_metrics = result.performance.as_ref().map(|perf| {
                crate::reporting::PerformanceMetrics {
                    mean: perf.stats.mean,
                    std_dev: perf.stats.std_dev,
                    min: perf.stats.min,
                    max: perf.stats.max,
                    median: perf.stats.median,
                    samples: perf.stats.samples,
                }
            }).unwrap_or_default();

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
        self.historical_data.get(benchmark_name)
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

            if perf.stats.mean > 0.1 { // More than 100ms
                recommendations.push("High execution time detected. Consider optimization opportunities.".to_string());
            }
        }

        // Check trend
        if let Some(trend) = trend {
            match trend.trend_type {
                TrendType::Degrading => {
                    recommendations.push(format!("Performance degrading at {:.2}% per run. Investigate recent changes.", trend.slope * 100.0));
                }
                TrendType::Improving => {
                    recommendations.push(format!("Performance improving at {:.2}% per run. Good trend!", trend.slope.abs() * 100.0));
                }
                TrendType::Volatile => {
                    recommendations.push("Performance is volatile. Consider stabilizing factors affecting benchmark.".to_string());
                }
                TrendType::Stable => {
                    if trend.r_squared > 0.8 {
                        recommendations.push("Performance is stable. Continue monitoring.".to_string());
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
            PerformanceMetrics { mean: 1.0, std_dev: 0.1, min: 0.9, max: 1.1, median: 1.0, samples: 10 },
            PerformanceMetrics { mean: 1.05, std_dev: 0.1, min: 0.95, max: 1.15, median: 1.05, samples: 10 },
            PerformanceMetrics { mean: 1.10, std_dev: 0.1, min: 1.0, max: 1.20, median: 1.10, samples: 10 },
            PerformanceMetrics { mean: 1.08, std_dev: 0.1, min: 0.98, max: 1.18, median: 1.08, samples: 10 },
            PerformanceMetrics { mean: 1.12, std_dev: 0.1, min: 1.02, max: 1.22, median: 1.12, samples: 10 },
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
            PerformanceMetrics { mean: 1.0, std_dev: 0.01, min: 0.99, max: 1.01, median: 1.0, samples: 10 },
            PerformanceMetrics { mean: 1.02, std_dev: 0.01, min: 1.01, max: 1.03, median: 1.02, samples: 10 },
            PerformanceMetrics { mean: 1.04, std_dev: 0.01, min: 1.03, max: 1.05, median: 1.04, samples: 10 },
            PerformanceMetrics { mean: 1.06, std_dev: 0.01, min: 1.05, max: 1.07, median: 1.06, samples: 10 },
            PerformanceMetrics { mean: 1.08, std_dev: 0.01, min: 1.07, max: 1.09, median: 1.08, samples: 10 },
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
        assert_relative_eq!(PerformanceAnalyzer::normal_cdf(-1.0), 0.1587, epsilon = 1e-3);
    }

    #[test]
    fn test_insufficient_data() {
        let mut analyzer = PerformanceAnalyzer::with_default_config();
        analyzer.add_result("sparse_benchmark",
            PerformanceMetrics { mean: 1.0, std_dev: 0.1, min: 0.9, max: 1.1, median: 1.0, samples: 10 });

        // Should fail with insufficient data
        let result = analyzer.analyze_trend("sparse_benchmark");
        assert!(result.is_err());
    }
}
