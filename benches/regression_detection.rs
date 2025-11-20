//! Regression detection and performance alerting for CFD operations
//!
//! This module provides automated regression detection including:
//! - Performance baseline tracking and comparison
//! - Statistical regression analysis
//! - Automated alerting for performance regressions
//! - Historical performance trend analysis
//! - Performance validation against requirements

use super::{BenchmarkConfig, PerformanceMetrics};
use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

/// Regression detection configuration
#[derive(Debug, Clone)]
pub struct RegressionConfig {
    /// Baseline file path for comparison
    pub baseline_path: Option<PathBuf>,
    /// Performance regression threshold (%)
    pub regression_threshold_pct: f64,
    /// Minimum statistical significance (p-value threshold)
    pub significance_threshold: f64,
    /// Historical data retention period (days)
    pub retention_days: u64,
    /// Alert configuration
    pub alerting: AlertConfig,
}

#[derive(Debug, Clone)]
pub struct AlertConfig {
    /// Enable console alerts
    pub console_alerts: bool,
    /// Enable file-based alerts
    pub file_alerts: bool,
    /// Alert output directory
    pub alert_dir: PathBuf,
    /// Critical regression threshold (%)
    pub critical_threshold_pct: f64,
}

impl Default for RegressionConfig {
    fn default() -> Self {
        Self {
            baseline_path: Some(PathBuf::from("performance_baselines.json")),
            regression_threshold_pct: 10.0,
            significance_threshold: 0.05,
            retention_days: 90,
            alerting: AlertConfig::default(),
        }
    }
}

impl Default for AlertConfig {
    fn default() -> Self {
        Self {
            console_alerts: true,
            file_alerts: true,
            alert_dir: PathBuf::from("performance_alerts"),
            critical_threshold_pct: 25.0,
        }
    }
}

/// Performance baseline data
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceBaseline {
    pub timestamp: u64,
    pub benchmarks: HashMap<String, BenchmarkBaseline>,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BenchmarkBaseline {
    pub name: String,
    pub mean_time_ns: u64,
    pub std_dev_ns: u64,
    pub throughput: f64,
    pub samples: usize,
    pub confidence_interval: (u64, u64), // (lower, upper) bounds in ns
}

/// Regression analysis result
#[derive(Debug, Clone)]
pub struct RegressionAnalysis {
    pub benchmark_name: String,
    pub current_performance: BenchmarkResult,
    pub baseline_performance: Option<BenchmarkBaseline>,
    pub regression_detected: bool,
    pub regression_magnitude: f64,     // percentage change
    pub statistical_significance: f64, // p-value
    pub trend: PerformanceTrend,
    pub alert_level: AlertLevel,
}

/// Performance trend analysis
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PerformanceTrend {
    Improving,
    Degrading,
    Stable,
    Volatile,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlertLevel {
    None,
    Warning,
    Critical,
}

/// Current benchmark result
#[derive(Debug, Clone)]
pub struct BenchmarkResult {
    pub name: String,
    pub mean_time: Duration,
    pub std_dev: Duration,
    pub throughput: f64,
    pub samples: usize,
    pub timestamp: u64,
}

/// Regression detector
pub struct RegressionDetector {
    config: RegressionConfig,
    historical_data: Vec<PerformanceBaseline>,
}

impl RegressionDetector {
    pub fn new(config: RegressionConfig) -> Self {
        let historical_data = Self::load_historical_data(&config);
        Self {
            config,
            historical_data,
        }
    }

    /// Load historical performance data
    fn load_historical_data(config: &RegressionConfig) -> Vec<PerformanceBaseline> {
        let mut data = Vec::new();

        if let Some(baseline_path) = &config.baseline_path {
            if baseline_path.exists() {
                match fs::read_to_string(baseline_path) {
                    Ok(content) => {
                        if let Ok(baseline) = serde_json::from_str::<PerformanceBaseline>(&content)
                        {
                            data.push(baseline);
                        }
                    }
                    Err(e) => eprintln!("Failed to load baseline data: {}", e),
                }
            }
        }

        // Load historical data directory if it exists
        let historical_dir = Path::new("performance_history");
        if historical_dir.exists() {
            if let Ok(entries) = fs::read_dir(historical_dir) {
                for entry in entries.flatten() {
                    if let Ok(content) = fs::read_to_string(entry.path()) {
                        if let Ok(baseline) = serde_json::from_str::<PerformanceBaseline>(&content)
                        {
                            // Check if data is within retention period
                            let current_time = SystemTime::now()
                                .duration_since(UNIX_EPOCH)
                                .unwrap_or_default()
                                .as_secs();

                            if current_time.saturating_sub(baseline.timestamp)
                                <= config.retention_days * 24 * 3600
                            {
                                data.push(baseline);
                            }
                        }
                    }
                }
            }
        }

        data.sort_by_key(|b| b.timestamp);
        data
    }

    /// Analyze current results against baseline
    pub fn analyze_regression(
        &self,
        current_results: &[BenchmarkResult],
    ) -> Vec<RegressionAnalysis> {
        let mut analyses = Vec::new();

        for current in current_results {
            let baseline = self.find_baseline(&current.name);
            let analysis = self.analyze_single_benchmark(current, baseline.as_ref());

            if analysis.regression_detected {
                self.handle_alert(&analysis);
            }

            analyses.push(analysis);
        }

        analyses
    }

    /// Find baseline for a specific benchmark
    fn find_baseline(&self, benchmark_name: &str) -> Option<BenchmarkBaseline> {
        // Use most recent baseline that contains this benchmark
        for baseline in self.historical_data.iter().rev() {
            if let Some(benchmark) = baseline.benchmarks.get(benchmark_name) {
                return Some(benchmark.clone());
            }
        }
        None
    }

    /// Analyze regression for a single benchmark
    fn analyze_single_benchmark(
        &self,
        current: &BenchmarkResult,
        baseline: Option<&BenchmarkBaseline>,
    ) -> RegressionAnalysis {
        let regression_detected;
        let regression_magnitude;
        let statistical_significance;
        let trend;

        if let Some(baseline) = baseline {
            // Calculate regression metrics
            let current_time_ns = current.mean_time.as_nanos() as f64;
            let baseline_time_ns = baseline.mean_time_ns as f64;

            regression_magnitude =
                ((current_time_ns - baseline_time_ns) / baseline_time_ns) * 100.0;

            // Simple statistical significance test (t-test approximation)
            let current_std = current.std_dev.as_nanos() as f64;
            let baseline_std = baseline.std_dev_ns as f64;

            let se = (current_std.powi(2) / current.samples as f64
                + baseline_std.powi(2) / baseline.samples as f64)
                .sqrt();

            let t_stat = (current_time_ns - baseline_time_ns) / se;
            statistical_significance = if t_stat.abs() > 2.0 {
                0.01
            } else if t_stat.abs() > 1.96 {
                0.05
            } else {
                0.1
            };

            regression_detected = regression_magnitude.abs() > self.config.regression_threshold_pct
                && statistical_significance < self.config.significance_threshold;

            // Determine trend
            trend = if regression_magnitude < -5.0 {
                PerformanceTrend::Improving
            } else if regression_magnitude > 5.0 {
                PerformanceTrend::Degrading
            } else {
                PerformanceTrend::Stable
            };
        } else {
            // No baseline available
            regression_detected = false;
            regression_magnitude = 0.0;
            statistical_significance = 1.0;
            trend = PerformanceTrend::Stable;
        }

        let alert_level = if regression_detected {
            if regression_magnitude.abs() > self.config.alerting.critical_threshold_pct {
                AlertLevel::Critical
            } else {
                AlertLevel::Warning
            }
        } else {
            AlertLevel::None
        };

        RegressionAnalysis {
            benchmark_name: current.name.clone(),
            current_performance: current.clone(),
            baseline_performance: baseline.cloned(),
            regression_detected,
            regression_magnitude,
            statistical_significance,
            trend,
            alert_level,
        }
    }

    /// Handle performance alerts
    fn handle_alert(&self, analysis: &RegressionAnalysis) {
        let alert_message = format!(
            "PERFORMANCE {}: {} - {:.1}% change (p={:.3})",
            match analysis.alert_level {
                AlertLevel::Critical => "REGRESSION",
                AlertLevel::Warning => "WARNING",
                AlertLevel::None => "INFO",
            },
            analysis.benchmark_name,
            analysis.regression_magnitude,
            analysis.statistical_significance
        );

        if self.config.alerting.console_alerts {
            match analysis.alert_level {
                AlertLevel::Critical => eprintln!("🚨 {}", alert_message),
                AlertLevel::Warning => println!("⚠️  {}", alert_message),
                AlertLevel::None => println!("ℹ️  {}", alert_message),
            }
        }

        if self.config.alerting.file_alerts {
            self.save_alert_to_file(analysis, &alert_message);
        }
    }

    /// Save alert to file
    fn save_alert_to_file(&self, analysis: &RegressionAnalysis, message: &str) {
        if let Err(e) = fs::create_dir_all(&self.config.alerting.alert_dir) {
            eprintln!("Failed to create alert directory: {}", e);
            return;
        }

        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();

        let filename = format!("alert_{}_{}.txt", analysis.benchmark_name, timestamp);
        let filepath = self.config.alerting.alert_dir.join(filename);

        let alert_content = format!(
            "Timestamp: {}\nBenchmark: {}\nAlert Level: {:?}\nRegression: {:.2}%\nStatistical Significance: {:.3}\nTrend: {:?}\n\nMessage: {}\n",
            timestamp,
            analysis.benchmark_name,
            analysis.alert_level,
            analysis.regression_magnitude,
            analysis.statistical_significance,
            analysis.trend,
            message
        );

        if let Err(e) = fs::write(&filepath, alert_content) {
            eprintln!("Failed to write alert file {}: {}", filepath.display(), e);
        }
    }

    /// Update baseline with current results
    pub fn update_baseline(&mut self, current_results: &[BenchmarkResult]) {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs();

        let mut benchmarks = HashMap::new();

        for result in current_results {
            let baseline = BenchmarkBaseline {
                name: result.name.clone(),
                mean_time_ns: result.mean_time.as_nanos() as u64,
                std_dev_ns: result.std_dev.as_nanos() as u64,
                throughput: result.throughput,
                samples: result.samples,
                confidence_interval: self.calculate_confidence_interval(result),
            };
            benchmarks.insert(result.name.clone(), baseline);
        }

        let new_baseline = PerformanceBaseline {
            timestamp,
            benchmarks,
        };

        // Save to baseline file
        if let Some(baseline_path) = &self.config.baseline_path {
            if let Ok(json) = serde_json::to_string_pretty(&new_baseline) {
                if let Err(e) = fs::write(baseline_path, json) {
                    eprintln!(
                        "Failed to save baseline to {}: {}",
                        baseline_path.display(),
                        e
                    );
                } else {
                    println!("Baseline updated: {}", baseline_path.display());
                }
            }
        }

        // Save to historical data
        let historical_dir = Path::new("performance_history");
        if let Err(e) = fs::create_dir_all(historical_dir) {
            eprintln!("Failed to create historical data directory: {}", e);
        } else {
            let historical_file = historical_dir.join(format!("baseline_{}.json", timestamp));
            if let Ok(json) = serde_json::to_string_pretty(&new_baseline) {
                let _ = fs::write(historical_file, json); // Ignore errors for historical data
            }
        }

        self.historical_data.push(new_baseline);
    }

    /// Calculate confidence interval for benchmark result
    fn calculate_confidence_interval(&self, result: &BenchmarkResult) -> (u64, u64) {
        let mean_ns = result.mean_time.as_nanos() as f64;
        let std_dev_ns = result.std_dev.as_nanos() as f64;
        let samples = result.samples as f64;

        // 95% confidence interval using t-distribution approximation
        let t_value = 1.96; // Approximately 95% CI for large samples
        let margin = t_value * std_dev_ns / samples.sqrt();

        let lower = (mean_ns - margin).max(0.0) as u64;
        let upper = (mean_ns + margin) as u64;

        (lower, upper)
    }
}

/// Run comprehensive regression detection benchmark
pub fn benchmark_regression_detection(c: &mut Criterion, config: &BenchmarkConfig) {
    let regression_config = RegressionConfig::default();
    let mut detector = RegressionDetector::new(regression_config);

    let mut group = c.benchmark_group("regression_detection");

    // Run current benchmarks
    let current_results = run_current_benchmarks(config);

    // Analyze regressions
    let analyses = detector.analyze_regression(&current_results);

    // Report results
    println!("\n=== Regression Detection Results ===");
    let mut regressions_found = 0;
    for analysis in &analyses {
        if analysis.regression_detected {
            regressions_found += 1;
            println!(
                "🚨 REGRESSION: {} - {:.1}% (p={:.3})",
                analysis.benchmark_name,
                analysis.regression_magnitude,
                analysis.statistical_significance
            );
        } else {
            println!("✅ STABLE: {}", analysis.benchmark_name);
        }
    }

    if regressions_found == 0 {
        println!("🎉 No performance regressions detected!");
    } else {
        println!(
            "⚠️  {} performance regressions detected!",
            regressions_found
        );
    }

    // Update baseline with current results
    detector.update_baseline(&current_results);

    group.finish();
}

/// Run current benchmark suite and collect results
fn run_current_benchmarks(config: &BenchmarkConfig) -> Vec<BenchmarkResult> {
    let mut results = Vec::new();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    // Simulate running various benchmarks
    for &size in &config.problem_sizes {
        // Memory allocation benchmark
        let result = run_memory_allocation_benchmark(size);
        results.push(BenchmarkResult {
            name: format!("memory_allocation_{}", size),
            mean_time: result.0,
            std_dev: result.1,
            throughput: calculate_throughput(result.0, size),
            samples: 100,
            timestamp,
        });

        // CFD computation benchmark
        let result = run_cfd_computation_benchmark(size);
        results.push(BenchmarkResult {
            name: format!("cfd_computation_{}", size),
            mean_time: result.0,
            std_dev: result.1,
            throughput: calculate_throughput(result.0, size),
            samples: 100,
            timestamp,
        });
    }

    results
}

/// Helper functions for benchmark execution
fn run_memory_allocation_benchmark(size: usize) -> (Duration, Duration) {
    let mut times = Vec::new();

    for _ in 0..10 {
        let start = std::time::Instant::now();
        let data = black_box(vec![0.0f64; size * size]);
        let _sum = data.iter().sum::<f64>();
        drop(data);
        times.push(start.elapsed());
    }

    let mean = times.iter().sum::<Duration>() / times.len() as u32;
    let variance = times
        .iter()
        .map(|&t| {
            let diff = if t > mean { t - mean } else { mean - t };
            diff.as_secs_f64().powi(2)
        })
        .sum::<f64>()
        / times.len() as f64;

    let std_dev = Duration::from_secs_f64(variance.sqrt());
    (mean, std_dev)
}

fn run_cfd_computation_benchmark(size: usize) -> (Duration, Duration) {
    let mut times = Vec::new();

    for _ in 0..10 {
        let start = std::time::Instant::now();
        use nalgebra::DMatrix;

        let mut field = DMatrix::<f64>::zeros(size, size);
        for i in 0..size.min(20) {
            for j in 0..size.min(20) {
                field[(i, j)] = (i as f64 * j as f64).sin();
            }
        }
        let _result = black_box(field);
        times.push(start.elapsed());
    }

    let mean = times.iter().sum::<Duration>() / times.len() as u32;
    let variance = times
        .iter()
        .map(|&t| {
            let diff = if t > mean { t - mean } else { mean - t };
            diff.as_secs_f64().powi(2)
        })
        .sum::<f64>()
        / times.len() as f64;

    let std_dev = Duration::from_secs_f64(variance.sqrt());
    (mean, std_dev)
}

fn calculate_throughput(duration: Duration, problem_size: usize) -> f64 {
    let operations = (problem_size * problem_size) as f64;
    operations / duration.as_secs_f64()
}

/// Generate regression detection recommendations
pub fn generate_regression_recommendations(analyses: &[RegressionAnalysis]) -> Vec<String> {
    let mut recommendations = Vec::new();

    let regressions: Vec<_> = analyses.iter().filter(|a| a.regression_detected).collect();

    if regressions.is_empty() {
        recommendations
            .push("✅ No performance regressions detected. Performance is stable.".to_string());
        return recommendations;
    }

    recommendations.push(format!(
        "⚠️ {} performance regressions detected. Immediate investigation required.",
        regressions.len()
    ));

    let critical_regressions: Vec<_> = regressions
        .iter()
        .filter(|a| matches!(a.alert_level, AlertLevel::Critical))
        .collect();

    if !critical_regressions.is_empty() {
        recommendations.push(format!(
            "🚨 {} critical regressions found. Performance has degraded significantly.",
            critical_regressions.len()
        ));
    }

    // Analyze trends
    let improving: Vec<_> = analyses
        .iter()
        .filter(|a| matches!(a.trend, PerformanceTrend::Improving))
        .collect();

    if !improving.is_empty() {
        recommendations.push(format!(
            "📈 {} benchmarks show performance improvements.",
            improving.len()
        ));
    }

    let degrading: Vec<_> = analyses
        .iter()
        .filter(|a| matches!(a.trend, PerformanceTrend::Degrading))
        .collect();

    if !degrading.is_empty() {
        recommendations.push(format!(
            "📉 {} benchmarks show performance degradation.",
            degrading.len()
        ));
    }

    recommendations
        .push("💡 Consider reviewing recent code changes and running profiling tools.".to_string());

    recommendations
}
