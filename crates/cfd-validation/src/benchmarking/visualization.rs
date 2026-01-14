//! Performance visualization and reporting tools
//!
//! This module provides comprehensive visualization capabilities for performance
//! benchmarking results, including charts, graphs, and HTML reports.
//!
//! Architecture follows domain-driven design with clear separation:
//! - Chart generation for performance metrics
//! - HTML report generation with embedded charts
//! - Interactive visualization using Chart.js
//! - Regression alert visualization and analysis

use super::analysis::RegressionAlert;
use super::scaling::ScalingResult;
use super::{BenchmarkResult, BenchmarkStatus};
use std::collections::HashMap;
use std::fmt;
use std::fmt::Write as _;

/// Chart types for CFD performance data visualization
///
/// Defines the available chart types for visualizing benchmark results,
/// performance trends, and scaling analysis in CFD validation reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChartType {
    /// Line chart for temporal performance trends and time series data
    ///
    /// Ideal for showing performance evolution over time, regression analysis,
    /// and historical benchmark comparisons. Best for continuous data visualization.
    Line,

    /// Bar chart for categorical performance comparisons
    ///
    /// Suitable for comparing performance across different CFD algorithms,
    /// problem sizes, or hardware configurations. Effective for discrete comparisons.
    Bar,

    /// Scatter plot for correlation analysis and outlier detection
    ///
    /// Used for analyzing relationships between performance variables,
    /// identifying correlations, and detecting measurement outliers in CFD benchmarks.
    Scatter,

    /// Histogram for performance distribution analysis
    ///
    /// Visualizes the distribution of performance measurements, helping identify
    /// measurement stability, variability, and potential performance clusters
    /// in CFD benchmark results.
    Histogram,
}

/// Configuration parameters for CFD performance visualization
///
/// Defines the visual properties and layout parameters for generating
/// performance charts, graphs, and data visualizations in CFD validation reports.
/// Controls both aesthetic and functional aspects of benchmark result presentation.
#[derive(Debug, Clone)]
pub struct VisualizationConfig {
    /// Chart width in pixels for output image or display
    ///
    /// Determines the horizontal resolution of generated visualizations.
    /// Higher values provide more detail but increase file size and rendering time.
    pub width: usize,

    /// Chart height in pixels for output image or display
    ///
    /// Determines the vertical resolution of generated visualizations.
    /// Should be proportional to width for optimal aspect ratio in CFD reports.
    pub height: usize,

    /// Chart title text displayed at the top of the visualization
    ///
    /// Descriptive title explaining the benchmark data being visualized.
    /// Should clearly indicate the CFD operation, algorithm, or analysis type.
    pub title: String,

    /// Label for the X-axis describing the independent variable
    ///
    /// Typically represents time, problem size, processor count, or other
    /// benchmark parameters in CFD performance analysis.
    pub x_label: String,

    /// Label for the Y-axis describing the dependent variable
    ///
    /// Usually represents performance metrics like execution time, throughput,
    /// memory usage, or efficiency measures in CFD benchmarks.
    pub y_label: String,

    /// Whether to display grid lines on the chart background
    ///
    /// Grid lines aid in reading precise values from performance charts.
    /// Recommended for detailed performance analysis but can be disabled
    /// for cleaner presentation in executive summaries.
    pub show_grid: bool,

    /// Color palette for chart elements (lines, bars, points)
    ///
    /// Array of color specifications (hex codes, names) for different data series.
    /// Should provide sufficient contrast for accessibility and clear differentiation
    /// of multiple CFD benchmark results in the same visualization.
    pub colors: Vec<String>,
}

impl Default for VisualizationConfig {
    fn default() -> Self {
        Self {
            width: 800,
            height: 600,
            title: "Performance Benchmark Results".to_string(),
            x_label: "Problem Size".to_string(),
            y_label: "Execution Time (ms)".to_string(),
            show_grid: true,
            colors: vec![
                "#1f77b4".to_string(), // Blue
                "#ff7f0e".to_string(), // Orange
                "#2ca02c".to_string(), // Green
                "#d62728".to_string(), // Red
                "#9467bd".to_string(), // Purple
            ],
        }
    }
}

/// Complete chart data structure for CFD performance visualization
///
/// Contains all data series and metadata needed to render a complete
/// performance chart for CFD benchmark analysis and validation reporting.
pub struct ChartData {
    /// X-axis labels for categorical or discrete data points
    ///
    /// Labels corresponding to each data point on the independent axis.
    /// For time-series data, these might be timestamps or iteration numbers.
    /// For scaling analysis, these could be processor counts or problem sizes.
    pub labels: Vec<String>,

    /// Data series to be plotted on the chart
    ///
    /// Collection of datasets representing different CFD benchmarks, algorithms,
    /// or performance metrics. Each dataset can have its own color and styling.
    pub datasets: Vec<Dataset>,
}

/// Individual data series for CFD performance visualization
///
/// Represents a single performance metric or benchmark result series
/// that can be plotted on a chart with consistent styling and labeling.
/// Multiple datasets can be combined in a single chart for comparative analysis.
#[derive(Debug, Clone)]
pub struct Dataset {
    /// Human-readable label for this data series
    ///
    /// Descriptive name that appears in legends and identifies the CFD operation,
    /// algorithm, or performance metric being visualized (e.g., "Navier-Stokes CPU",
    /// "Memory Usage", "GPU Throughput").
    pub label: String,

    /// Performance measurement values for this series [units vary by metric]
    ///
    /// Array of numerical values representing performance measurements.
    /// Could be execution times, throughput rates, memory usage, or other
    /// CFD performance metrics. Must correspond 1:1 with chart labels.
    pub data: Vec<f64>,

    /// Color specification for rendering this data series
    ///
    /// Color used to draw lines, bars, or points for this dataset.
    /// Should be specified as hex color code (e.g., "#FF0000") or named color.
    /// Must provide sufficient contrast for accessibility and clear differentiation.
    pub color: String,
}

/// HTML report generator
pub struct HtmlReportGenerator {
    config: VisualizationConfig,
}

impl HtmlReportGenerator {
    /// Create a new HTML report generator with custom visualization configuration
    ///
    /// Initializes the report generator with user-specified chart dimensions,
    /// styling preferences, and output parameters for CFD performance visualization.
    /// Use this constructor when you need fine control over report appearance
    /// and chart formatting for specific CFD validation requirements.
    ///
    /// # Parameters
    ///
    /// * `config` - Visualization configuration specifying chart dimensions,
    ///   colors, labels, and other display parameters for performance reports
    pub fn new(config: VisualizationConfig) -> Self {
        Self { config }
    }

    /// Create a new HTML report generator with default visualization configuration
    ///
    /// Initializes the report generator with sensible defaults optimized for
    /// CFD performance analysis:
    /// - 800x600 pixel charts for good readability
    /// - Grid lines enabled for precise value reading
    /// - Standard color palette for accessibility
    /// - Appropriate labels for CFD performance metrics
    ///
    /// Suitable for most CFD benchmarking scenarios without requiring manual styling configuration.
    pub fn with_default_config() -> Self {
        Self::new(VisualizationConfig::default())
    }

    /// Generate comprehensive HTML performance report
    pub fn generate_report(
        &self,
        results: &[BenchmarkResult],
        _scaling_results: Option<&[ScalingResult]>,
        alerts: Option<&[RegressionAlert]>,
    ) -> String {
        let mut html = String::new();

        // HTML header
        html.push_str(&self.generate_html_header());

        // Title and summary
        html.push_str(&Self::generate_summary_section(
            results,
            alerts.unwrap_or(&[]),
        ));

        // Performance charts with Chart.js
        html.push_str(&Self::generate_performance_chart(results));

        // Scaling analysis
        // Use results which contain scaling info directly
        html.push_str(&Self::generate_scaling_chart(results));

        // Regression alerts (if available)
        if let Some(alerts) = alerts {
            html.push_str(&Self::generate_alerts_section(alerts));
        }

        // Detailed results table
        html.push_str(&Self::generate_results_table(results));

        // Recommendations
        html.push_str(&Self::generate_recommendations(
            results,
            alerts.unwrap_or(&[]),
        ));

        // HTML footer
        html.push_str(&Self::generate_html_footer());

        html
    }

    fn generate_html_header(&self) -> String {
        format!(
            r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{} - Performance Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
    <script src="https://cdn.jsdelivr.net/npm/chartjs-adapter-date-fns@3.0.0"></script>
    <script src="https://cdn.jsdelivr.net/npm/date-fns@3.0.0"></script>
    <style>
        :root {{
            --primary-gradient: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            --success-color: #28a745;
            --warning-color: #ffc107;
            --danger-color: #dc3545;
            --info-color: #17a2b8;
            --light-bg: #f8f9fa;
            --shadow: 0 4px 6px rgba(0,0,0,0.1);
            --border-radius: 12px;
        }}

        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: var(--primary-gradient);
            min-height: 100vh;
            line-height: 1.6;
            color: #333;
        }}

        .container {{
            max-width: 1400px;
            margin: 0 auto 20px auto;
            background: white;
            padding: 40px;
            border-radius: var(--border-radius);
            box-shadow: var(--shadow);
        }}

        .header {{
            text-align: center;
            margin-bottom: 50px;
            position: relative;
        }}

        .header::after {{
            content: '';
            position: absolute;
            bottom: -10px;
            left: 50%;
            transform: translateX(-50%);
            width: 100px;
            height: 4px;
            background: var(--primary-gradient);
            border-radius: 2px;
        }}

        h1 {{
            color: #2c3e50;
            font-size: 3em;
            font-weight: 300;
            margin-bottom: 10px;
            text-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}

        h2 {{
            color: #34495e;
            font-size: 2em;
            font-weight: 400;
            margin: 50px 0 20px 0;
            padding-bottom: 15px;
            border-bottom: 3px solid #3498db;
            position: relative;
        }}

        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 25px;
            margin: 40px 0;
        }}

        .summary-card {{
            background: linear-gradient(135deg, #667eea, #764ba2);
            color: white;
            padding: 30px;
            border-radius: var(--border-radius);
            text-align: center;
            box-shadow: var(--shadow);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
            border: 2px solid transparent;
        }}

        .summary-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }}

        .summary-card .value {{
            font-size: 3em;
            font-weight: bold;
            margin-bottom: 10px;
            text-shadow: 0 1px 2px rgba(0,0,0,0.3);
        }}

        .summary-card .label {{
            font-size: 1.1em;
            opacity: 0.9;
            font-weight: 500;
        }}

        .chart-container {{
            margin: 40px 0;
            padding: 30px;
            background: var(--light-bg);
            border-radius: var(--border-radius);
            border: 1px solid #e9ecef;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}

        .chart-wrapper {{
            position: relative;
            height: 500px;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            padding: 20px;
            box-shadow: inset 0 1px 3px rgba(0,0,0,0.1);
        }}

        .alert {{
            padding: 20px 25px;
            margin: 15px 0;
            border-radius: 10px;
            border: none;
            box-shadow: var(--shadow);
            display: flex;
            align-items: center;
            font-weight: 500;
            position: relative;
            overflow: hidden;
        }}

        .alert::before {{
            content: '';
            position: absolute;
            left: 0;
            top: 0;
            bottom: 0;
            width: 4px;
        }}

        .alert-critical {{
            background: linear-gradient(135deg, #ff6b6b, #ee5a52);
            color: white;
        }}
        .alert-critical::before {{ background: #dc3545; }}

        .alert-high {{
            background: linear-gradient(135deg, #ff9f43, #ee8c3a);
            color: white;
        }}
        .alert-high::before {{ background: #fd7e14; }}

        .alert-medium {{
            background: linear-gradient(135deg, #ffd93d, #ffc93c);
            color: #856404;
        }}
        .alert-medium::before {{ background: #ffc107; }}

        .alert-low {{
            background: linear-gradient(135deg, #6bcf7f, #4ecdc4);
            color: white;
        }}
        .alert-low::before {{ background: #20c997; }}

        .alert-content {{
            flex: 1;
        }}

        .alert strong {{
            display: block;
            margin-bottom: 5px;
            font-size: 1.1em;
        }}

        .alert small {{
            opacity: 0.8;
            font-size: 0.9em;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 30px 0;
            font-size: 0.95em;
            box-shadow: var(--shadow);
            border-radius: var(--border-radius);
            overflow: hidden;
            background: white;
        }}

        th, td {{
            padding: 18px 15px;
            text-align: left;
            border-bottom: 1px solid #e9ecef;
        }}

        th {{
            background: var(--primary-gradient);
            color: white;
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 0.5px;
            position: sticky;
            top: 0;
        }}

        tr:nth-child(even) {{
            background: #f8f9fa;
        }}

        tr:hover {{
            background: #e8f4fd;
            transition: background-color 0.3s ease;
        }}

        .status-indicator {{
            display: inline-block;
            width: 12px;
            height: 12px;
            border-radius: 50%;
            margin-right: 10px;
            vertical-align: middle;
        }}

        .status-passed {{ background: var(--success-color); }}
        .status-failed {{ background: var(--danger-color); }}
        .status-regression {{ background: var(--warning-color); }}
        .status-skipped {{ background: #6c757d; }}

        .recommendations {{
            background: linear-gradient(135deg, #a8e6cf, #dcedc8);
            padding: 30px;
            border-radius: var(--border-radius);
            margin: 40px 0;
            border: 2px solid #c8e6c9;
        }}

        .recommendations h2 {{
            color: #2e7d32;
            margin-top: 0;
            font-size: 1.8em;
        }}

        .recommendation-item {{
            margin: 15px 0;
            padding: 15px 20px;
            background: rgba(255, 255, 255, 0.9);
            border-radius: 8px;
            border-left: 4px solid #4caf50;
            font-weight: 500;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            transition: transform 0.2s ease;
        }}

        .recommendation-item:hover {{
            transform: translateX(5px);
        }}

        .performance-chart {{
            margin: 20px 0;
            padding: 20px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}

        .chart-legend {{
            display: flex;
            flex-wrap: wrap;
            gap: 20px;
            margin: 15px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 6px;
        }}

        .legend-item {{
            display: flex;
            align-items: center;
            font-size: 0.9em;
        }}

        .legend-color {{
            width: 16px;
            height: 16px;
            border-radius: 3px;
            margin-right: 8px;
            border: 1px solid rgba(0,0,0,0.1);
        }}

        @media (max-width: 1024px) {{
            .container {{
                margin: 10px;
                padding: 25px;
            }}

            .summary-grid {{
                grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
                gap: 20px;
            }}

            h1 {{
                font-size: 2.5em;
            }}

            h2 {{
                font-size: 1.6em;
            }}

            .chart-wrapper {{
                height: 400px;
            }}
        }}

        @media (max-width: 768px) {{
            body {{
                padding: 10px;
            }}

            .container {{
                padding: 20px;
            }}

            .summary-grid {{
                grid-template-columns: 1fr;
            }}

            table {{
                font-size: 0.8em;
            }}

            th, td {{
                padding: 12px 8px;
            }}

            .alert {{
                padding: 15px 20px;
                flex-direction: column;
                align-items: flex-start;
            }}

            .chart-wrapper {{
                height: 300px;
            }}
        }}

        @keyframes fadeIn {{
            from {{ opacity: 0; transform: translateY(20px); }}
            to {{ opacity: 1; transform: translateY(0); }}
        }}

        .chart-container, .summary-card, .alert, .recommendation-item {{
            animation: fadeIn 0.6s ease-out;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{}</h1>
            <p>Generated on {}</p>
        </div>"#,
            self.config.title,
            self.config.title,
            chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
        )
    }

    fn generate_summary_section(results: &[BenchmarkResult], alerts: &[RegressionAlert]) -> String {
        let total_operations = results.len();
        let success_rate = Self::calculate_success_rate(results);
        let total_alerts = alerts.len();
        let critical_alerts = alerts.iter().filter(|a| a.degradation_rate > 25.0).count();

        format!(
            r#"<div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Operations</h3>
                    <p style="font-size: 2em; margin: 0;">{}</p>
                </div>
                <div class="summary-card">
                    <h3>Success Rate</h3>
                    <p style="font-size: 2em; margin: 0; color: {};">{:.1}%</p>
                </div>
                <div class="summary-card">
                    <h3>Critical Alerts</h3>
                    <p style="font-size: 2em; margin: 0; color: {};">{}</p>
                </div>
                <div class="summary-card">
                    <h3>Total Alerts</h3>
                    <p style="font-size: 2em; margin: 0;">{}</p>
                </div>
            </div>"#,
            total_operations,
            if success_rate >= 95.0 {
                "#28a745"
            } else if success_rate >= 80.0 {
                "#ffc107"
            } else {
                "#dc3545"
            },
            success_rate,
            if critical_alerts > 0 {
                "#dc3545"
            } else {
                "#28a745"
            },
            critical_alerts,
            total_alerts
        )
    }

    fn generate_alerts_section(alerts: &[RegressionAlert]) -> String {
        if alerts.is_empty() {
            return r#"<div class="chart-container">
                <h2>Performance Alerts</h2>
                <p style="color: #28a745;">✅ No performance regressions detected</p>
            </div>"#
                .to_string();
        }

        let mut alerts_html = String::from(
            r#"<div class="chart-container">
            <h2>Performance Alerts</h2>"#,
        );

        for alert in alerts {
            // Determine severity based on degradation rate (architectural mapping)
            let (alert_class, severity_text) = if alert.degradation_rate > 25.0 {
                ("alert-critical", "Critical")
            } else if alert.degradation_rate > 15.0 {
                ("alert-high", "High")
            } else if alert.degradation_rate > 5.0 {
                ("alert-medium", "Medium")
            } else {
                ("alert-low", "Low")
            };

            write!(
                &mut alerts_html,
                r#"<div class="alert {}">
                    <strong>{}:</strong> {}<br>
                    <small>Degradation: {:.1}%, Confidence: {:.1}%</small>
                </div>"#,
                alert_class,
                severity_text,
                alert.benchmark_name,
                alert.degradation_rate,
                alert.confidence * 100.0
            )
            .expect("writing to String cannot fail");
        }

        alerts_html.push_str("</div>");
        alerts_html
    }

    fn generate_results_table(results: &[BenchmarkResult]) -> String {
        let mut table_html = r#"<div class="chart-container">
            <h2>Detailed Results</h2>
            <table>
                <thead>
                    <tr>
                        <th>Operation</th>
                        <th>Duration</th>
                        <th>Performance</th>
                        <th>Memory</th>
                        <th>Scaling</th>
                        <th>Status</th>
                    </tr>
                </thead>
                <tbody>"#
            .to_string();

        for result in results {
            // Extract performance metrics
            let perf_time = result.performance.as_ref().map_or_else(
                || "N/A".to_string(),
                |p| {
                    format!(
                        "{:.3}ms ± {:.3}ms",
                        p.stats.mean * 1000.0,
                        p.stats.std_dev * 1000.0
                    )
                },
            );

            // Extract memory metrics (use peak memory usage for footprint analysis)
            let memory_usage = result.memory.as_ref().map_or_else(
                || "N/A".to_string(),
                |m| format!("{:.1}MB", m.peak_allocated as f64 / 1_048_576.0),
            );

            // Extract scaling info
            let scaling_info = result
                .scaling
                .as_ref()
                .and_then(|s| {
                    // Get average speedup across all measurements
                    if s.speedup_factors.is_empty() {
                        None
                    } else {
                        let total: f64 = s.speedup_factors.values().sum();
                        let avg = total / s.speedup_factors.len() as f64;
                        Some(format!("{avg:.2}x avg speedup"))
                    }
                })
                .unwrap_or_else(|| "N/A".to_string());

            // Determine status color and text
            let (status_color, status_text) = match (&result.status, &result.regression_detected) {
                (BenchmarkStatus::Passed, None) => ("#28a745", "✅ Passed"),
                (BenchmarkStatus::Passed, Some(reg)) if *reg > 25.0 => {
                    ("#dc3545", "🚨 Major Regression")
                }
                (BenchmarkStatus::Passed, Some(reg)) if *reg > 10.0 => {
                    ("#ffc107", "⚠️ Minor Regression")
                }
                (BenchmarkStatus::Passed, Some(_)) => ("#17a2b8", "⚠️ Regression"),
                (BenchmarkStatus::Failed, _) => ("#dc3545", "❌ Failed"),
                (BenchmarkStatus::Regression, _) => ("#dc3545", "📉 Regression"),
                (BenchmarkStatus::Skipped, _) => ("#6c757d", "⏭️ Skipped"),
                (BenchmarkStatus::Running, _) => ("#007bff", "⏳ Running"),
            };

            write!(
                &mut table_html,
                r#"<tr>
                    <td>{}</td>
                    <td>{:.3}ms</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td>{}</td>
                    <td style="color: {};">{}</td>
                </tr>"#,
                result.name,
                result.duration.as_secs_f64() * 1000.0,
                perf_time,
                memory_usage,
                scaling_info,
                status_color,
                status_text
            )
            .expect("writing to String cannot fail");
        }

        table_html.push_str(r"</tbody></table></div>");
        table_html
    }

    fn generate_recommendations(results: &[BenchmarkResult], alerts: &[RegressionAlert]) -> String {
        let mut recommendations = Vec::new();

        // Analyze results for recommendations with architectural purity
        let critical_alerts = alerts.iter().filter(|a| a.degradation_rate > 25.0).count();

        if critical_alerts > 0 {
            recommendations.push(format!(
                "🚨 Address {critical_alerts} critical performance regressions immediately"
            ));
        }

        let failed_tests = results
            .iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Failed))
            .count();

        if failed_tests > 0 {
            recommendations.push(format!("❌ Fix {failed_tests} failed benchmark operations"));
        }

        // Check for general regressions (any regression > 5% but < 25%)
        let minor_regressions = alerts
            .iter()
            .filter(|a| a.degradation_rate > 5.0 && a.degradation_rate <= 25.0)
            .count();

        if minor_regressions > 0 {
            recommendations.push(format!(
                "⚠️ Review {minor_regressions} minor performance regressions"
            ));
        }

        // Memory efficiency analysis
        let high_memory_usage = results
            .iter()
            .filter_map(|r| r.memory.as_ref())
            .filter(|m| m.total_allocated > 100_000_000) // 100MB
            .count();

        if high_memory_usage > 0 {
            recommendations.push(format!(
                "🧠 Optimize memory usage in {high_memory_usage} operations (>100MB)"
            ));
        }

        // Scaling efficiency analysis
        let poor_scaling = results
            .iter()
            .filter_map(|r| r.scaling.as_ref())
            .filter(|s| {
                // Calculate average efficiency
                if s.parallel_efficiency.is_empty() {
                    false
                } else {
                    let total: f64 = s.parallel_efficiency.values().sum();
                    let avg = total / s.parallel_efficiency.len() as f64;
                    avg < 0.7
                }
            })
            .count();

        if poor_scaling > 0 {
            recommendations.push(format!(
                "⚡ Improve parallel scaling efficiency in {poor_scaling} operations"
            ));
        }

        if recommendations.is_empty() {
            recommendations.push(
                "✅ Performance is stable - continue monitoring for architectural consistency"
                    .to_string(),
            );
        }

        let mut rec_html = r#"<div class="recommendations">
            <h2>Architectural Recommendations</h2>"#
            .to_string();

        for rec in recommendations {
            write!(
                &mut rec_html,
                r#"<div class="recommendation-item">• {rec}</div>"#
            )
            .expect("writing to String cannot fail");
        }

        rec_html.push_str("</div>");
        rec_html
    }

    /// Generate Chart.js performance chart
    fn generate_performance_chart(results: &[BenchmarkResult]) -> String {
        let mut labels = Vec::new();
        let mut execution_times = Vec::new();
        let mut memory_usage = Vec::new();

        for result in results {
            labels.push(format!("\"{}\"", result.name));

            // Extract execution time (convert to milliseconds)
            execution_times.push(format!("{:.3}", result.duration.as_secs_f64() * 1000.0));

            // Extract memory usage (convert to MB) - use peak memory usage
            let mem_mb = result
                .memory
                .as_ref()
                .map_or(0.0, |m| m.peak_allocated as f64 / 1_048_576.0);
            memory_usage.push(format!("{mem_mb:.1}"));
        }

        format!(
            r#"
        <div class="chart-container">
            <h2>Performance Metrics Overview</h2>
            <div class="chart-legend">
                <div class="legend-item">
                    <div class="legend-color" style="background: #3498db;"></div>
                    <span>Execution Time (ms)</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background: #e74c3c;"></div>
                    <span>Memory Usage (MB)</span>
                </div>
            </div>
            <div class="chart-wrapper">
                <canvas id="performanceChart"></canvas>
            </div>
        </div>
        <script>
            const ctx = document.getElementById('performanceChart').getContext('2d');
            new Chart(ctx, {{
                type: 'bar',
                data: {{
                    labels: [{}],
                    datasets: [{{
                        label: 'Execution Time (ms)',
                        data: [{}],
                        backgroundColor: 'rgba(52, 152, 219, 0.8)',
                        borderColor: 'rgba(52, 152, 219, 1)',
                        borderWidth: 1,
                        yAxisID: 'y'
                    }}, {{
                        label: 'Memory Usage (MB)',
                        data: [{}],
                        backgroundColor: 'rgba(231, 76, 60, 0.8)',
                        borderColor: 'rgba(231, 76, 60, 1)',
                        borderWidth: 1,
                        yAxisID: 'y1'
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    interaction: {{
                        mode: 'index',
                        intersect: false,
                    }},
                    plugins: {{
                        title: {{
                            display: true,
                            text: 'CFD Performance Analysis',
                            font: {{
                                size: 16,
                                weight: 'bold'
                            }}
                        }},
                        legend: {{
                            display: true,
                            position: 'top'
                        }},
                        tooltip: {{
                            callbacks: {{
                                label: function(context) {{
                                    let label = context.dataset.label || '';
                                    if (label) {{
                                        label += ': ';
                                    }}
                                    if (context.parsed.y !== null) {{
                                        label += context.parsed.y.toLocaleString();
                                        if (context.datasetIndex === 0) label += ' ms';
                                        else if (context.datasetIndex === 1) label += ' MB';
                                        else if (context.datasetIndex === 2) label += ' ops/sec';
                                    }}
                                    return label;
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            display: true,
                            title: {{
                                display: true,
                                text: 'Benchmark Operations'
                            }}
                        }},
                        y: {{
                            type: 'linear',
                            display: true,
                            position: 'left',
                            title: {{
                                display: true,
                                text: 'Execution Time (ms)'
                            }},
                            grid: {{
                                drawOnChartArea: false,
                            }},
                        }},
                        y1: {{
                            type: 'linear',
                            display: true,
                            position: 'right',
                            title: {{
                                display: true,
                                text: 'Memory Usage (MB)'
                            }},
                            grid: {{
                                drawOnChartArea: false,
                            }},
                        }}
                    }}
                }}
            }});
        </script>"#,
            labels.join(","),
            execution_times.join(","),
            memory_usage.join(",")
        )
    }

    /// Generate scaling analysis chart (Speedup vs Processors)
    fn generate_scaling_chart(results: &[BenchmarkResult]) -> String {
        let scaling_benchmarks: Vec<&BenchmarkResult> =
            results.iter().filter(|r| r.scaling.is_some()).collect();

        if scaling_benchmarks.is_empty() {
            return String::new();
        }

        // Collect all unique processor counts to define the X-axis
        let mut all_processor_counts = Vec::new();
        for result in &scaling_benchmarks {
            if let Some(scaling) = &result.scaling {
                for &count in &scaling.processor_counts {
                    if !all_processor_counts.contains(&count) {
                        all_processor_counts.push(count);
                    }
                }
            }
        }
        all_processor_counts.sort_unstable();

        if all_processor_counts.is_empty() {
            return String::new();
        }

        // Prepare datasets
        let mut datasets = Vec::new();
        let colors = [
            "#4e79a7", "#f28e2c", "#e15759", "#76b7b2", "#59a14f", "#edc949", "#af7aa1", "#ff9da7",
            "#9c755f", "#bab0ac",
        ];

        for (i, result) in scaling_benchmarks.iter().enumerate() {
            if let Some(scaling) = &result.scaling {
                let mut data = Vec::new();
                for &proc_count in &all_processor_counts {
                    // Look up speedup for this processor count and the benchmark's problem size
                    let key = (result.problem_size, proc_count);
                    if let Some(&speedup) = scaling.speedup_factors.get(&key) {
                        data.push(Some(speedup));
                    } else {
                        data.push(None);
                    }
                }

                datasets.push(serde_json::json!({
                    "label": format!("{} (N={})", result.name, result.problem_size),
                    "data": data,
                    "borderColor": colors[i % colors.len()],
                    "backgroundColor": colors[i % colors.len()],
                    "fill": false,
                    "tension": 0.4
                }));
            }
        }

        let chart_data = serde_json::json!({
            "labels": all_processor_counts,
            "datasets": datasets
        });

        let chart_data_json = serde_json::to_string(&chart_data).unwrap_or_else(|e| {
            serde_json::json!({
                "labels": [],
                "datasets": [],
                "error": e.to_string()
            })
            .to_string()
        });

        format!(
            r#"<div class="chart-container">
                <div class="chart-header">
                    <h2>Parallel Scaling Performance</h2>
                    <p class="chart-subtitle">Speedup vs Processor Count</p>
                </div>
                <div class="chart-wrapper">
                    <canvas id="scalingChart"></canvas>
                </div>
            </div>
            <script>
                const scalingCtx = document.getElementById('scalingChart').getContext('2d');
                new Chart(scalingCtx, {{
                    type: 'line',
                    data: {chart_data_json},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        interaction: {{
                            mode: 'index',
                            intersect: false,
                        }},
                        plugins: {{
                            title: {{
                                display: false
                            }},
                            legend: {{
                                display: true,
                                position: 'top'
                            }},
                            tooltip: {{
                                callbacks: {{
                                    label: function(context) {{
                                        let label = context.dataset.label || '';
                                        if (label) {{
                                            label += ': ';
                                        }}
                                        if (context.parsed.y !== null) {{
                                            label += context.parsed.y.toFixed(2) + 'x';
                                        }}
                                        return label;
                                    }}
                                }}
                            }}
                        }},
                        scales: {{
                            x: {{
                                title: {{
                                    display: true,
                                    text: 'Processor Count'
                                }}
                            }},
                            y: {{
                                title: {{
                                    display: true,
                                    text: 'Speedup Factor'
                                }},
                                min: 0
                            }}
                        }}
                    }}
                }});
            </script>"#,
        )
    }

    /// Calculate success rate across all results
    fn calculate_success_rate(results: &[BenchmarkResult]) -> f64 {
        if results.is_empty() {
            return 0.0;
        }

        let successful = results
            .iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Passed))
            .count();

        (successful as f64 / results.len() as f64) * 100.0
    }

    fn generate_html_footer() -> String {
        r"</div></body></html>".to_string()
    }
}

/// Performance dashboard generator
pub struct PerformanceDashboard {
    generator: HtmlReportGenerator,
}

impl PerformanceDashboard {
    /// Create a new performance dashboard generator with default configuration
    pub fn new() -> Self {
        Self {
            generator: HtmlReportGenerator::with_default_config(),
        }
    }

    /// Generate complete performance dashboard
    pub fn generate_dashboard(
        &self,
        results: &[BenchmarkResult],
        scaling_results: Option<&[ScalingResult]>,
        alerts: Option<&[RegressionAlert]>,
        output_path: &str,
    ) -> std::io::Result<()> {
        let html_content = self
            .generator
            .generate_report(results, scaling_results, alerts);

        std::fs::write(output_path, html_content)?;

        println!("Performance dashboard generated: {output_path}");
        println!("Open in browser to view interactive charts and detailed analysis");

        Ok(())
    }

    /// Generate summary report (text format)
    pub fn generate_summary_report(&self, results: &[BenchmarkResult]) -> String {
        let mut report = String::new();

        report.push_str("CFD Performance Benchmark Summary\n");
        report.push_str("==================================\n\n");

        writeln!(
            &mut report,
            "Total Operations Benchmarked: {}",
            results.len()
        )
        .expect("writing to String cannot fail");
        writeln!(
            &mut report,
            "Report Generated: {}\n",
            chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
        )
        .expect("writing to String cannot fail");

        // Summary statistics
        let total_time: f64 = results.iter().map(|r| r.duration.as_secs_f64()).sum();

        let avg_time = total_time / results.len() as f64;
        let max_time = results
            .iter()
            .map(|r| r.duration.as_secs_f64())
            .fold(0.0, f64::max);

        report.push_str("Performance Statistics:\n");
        writeln!(
            &mut report,
            "  Average Execution Time: {:.3} ms",
            avg_time * 1000.0
        )
        .expect("writing to String cannot fail");
        writeln!(
            &mut report,
            "  Maximum Execution Time: {:.3} ms",
            max_time * 1000.0
        )
        .expect("writing to String cannot fail");
        writeln!(&mut report, "  Total Execution Time: {total_time:.3} s\n")
            .expect("writing to String cannot fail");

        // Operation breakdown
        report.push_str("Operation Breakdown:\n");
        let mut operations: HashMap<String, Vec<&BenchmarkResult>> = HashMap::new();

        for result in results {
            operations
                .entry(result.name.clone())
                .or_default()
                .push(result);
        }

        for (operation, op_results) in operations {
            let op_avg_time = op_results
                .iter()
                .map(|r| r.duration.as_secs_f64())
                .sum::<f64>()
                / op_results.len() as f64;

            writeln!(
                &mut report,
                "  {}: {:.3} ms ({} runs)",
                operation,
                op_avg_time * 1000.0,
                op_results.len()
            )
            .expect("writing to String cannot fail");
        }

        report
    }
}

impl Default for PerformanceDashboard {
    fn default() -> Self {
        Self::new()
    }
}

/// Export formats for CFD benchmark results and performance reports
///
/// Defines the available output formats for exporting benchmark data,
/// analysis results, and performance visualizations from CFD validation runs.
/// Each format is optimized for different use cases in CFD development workflows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    /// Interactive HTML report with embedded charts and visualizations
    ///
    /// Generates comprehensive web-based reports with interactive performance charts,
    /// scaling analysis graphs, and detailed CFD benchmark results. Ideal for
    /// sharing results with stakeholders and detailed performance analysis.
    Html,

    /// Machine-readable JSON format for programmatic analysis
    ///
    /// Structured data format suitable for automated processing, CI/CD integration,
    /// and further analysis with external tools. Contains all raw benchmark data,
    /// statistical analysis, and metadata in a parseable format.
    Json,

    /// Comma-separated values for spreadsheet analysis
    ///
    /// Tabular format optimized for import into Excel, Google Sheets, or other
    /// spreadsheet applications. Contains summarized benchmark results and key
    /// performance metrics in a human-readable table format.
    Csv,

    /// Plain text format for command-line and log analysis
    ///
    /// Simple text-based format suitable for terminal output, log files, and
    /// basic reporting. Contains formatted benchmark summaries and key metrics
    /// without complex formatting or visualizations.
    Text,
}

/// Benchmark result exporter
pub struct BenchmarkExporter;

impl BenchmarkExporter {
    /// Export benchmark results to specified format
    pub fn export_results(
        results: &[BenchmarkResult],
        format: ExportFormat,
        output_path: &str,
    ) -> std::io::Result<()> {
        match format {
            ExportFormat::Html => {
                let dashboard = PerformanceDashboard::new();
                dashboard.generate_dashboard(results, None, None, output_path)
            }
            ExportFormat::Json => {
                let json = serde_json::to_string_pretty(results)?;
                std::fs::write(output_path, json)
            }
            ExportFormat::Csv => {
                let csv_content = Self::results_to_csv(results);
                std::fs::write(output_path, csv_content)
            }
            ExportFormat::Text => {
                let dashboard = PerformanceDashboard::new();
                let text_report = dashboard.generate_summary_report(results);
                std::fs::write(output_path, text_report)
            }
        }
    }

    /// Convert benchmark results to CSV format
    fn results_to_csv(results: &[BenchmarkResult]) -> String {
        let mut csv = String::from("benchmark_name,problem_size,duration_ms,memory_usage_mb,throughput_ops_per_sec,status,regression_detected\n");

        for result in results {
            let throughput = if result.duration.as_secs_f64() > 0.0 {
                result.problem_size as f64 / result.duration.as_secs_f64()
            } else {
                0.0
            };

            let memory_mb = match result.memory.as_ref() {
                Some(memory) => {
                    let mut memory_mb = String::new();
                    let _ = write!(
                        &mut memory_mb,
                        "{:.1}",
                        memory.total_allocated as f64 / 1024.0 / 1024.0
                    );
                    memory_mb
                }
                None => "N/A".to_string(),
            };

            let _ = writeln!(
                &mut csv,
                "{},{},{:.3},{},{:.2},{:?},{}",
                result.name,
                result.problem_size,
                result.duration.as_secs_f64() * 1000.0,
                memory_mb,
                throughput,
                result.status,
                result.regression_detected.is_some()
            );
        }

        csv
    }
}

/// Performance comparison utilities
pub struct PerformanceComparator;

impl PerformanceComparator {
    /// Compare current results against baseline
    pub fn compare_against_baseline(
        current: &[BenchmarkResult],
        baseline: &[BenchmarkResult],
    ) -> Vec<PerformanceComparison> {
        let mut comparisons = Vec::new();

        // Create lookup map for baseline results
        let baseline_map: HashMap<(String, usize), &BenchmarkResult> = baseline
            .iter()
            .map(|r| ((r.name.clone(), r.problem_size), r))
            .collect();

        for current_result in current {
            let key = (current_result.name.clone(), current_result.problem_size);

            if let Some(baseline_result) = baseline_map.get(&key) {
                let current_time = current_result.duration.as_secs_f64();
                let baseline_time = baseline_result.duration.as_secs_f64();

                let time_change = ((current_time - baseline_time) / baseline_time) * 100.0;
                let is_improvement = time_change < 0.0;

                comparisons.push(PerformanceComparison {
                    benchmark_name: current_result.name.clone(),
                    problem_size: current_result.problem_size,
                    baseline_time,
                    current_time,
                    time_change_percent: time_change,
                    is_improvement,
                });
            }
        }

        comparisons
    }
}

/// Performance comparison result
#[derive(Debug, Clone)]
pub struct PerformanceComparison {
    /// Name of the benchmark
    pub benchmark_name: String,
    /// Size of the problem (e.g., number of grid points)
    pub problem_size: usize,
    /// Execution time in baseline (seconds)
    pub baseline_time: f64,
    /// Execution time in current run (seconds)
    pub current_time: f64,
    /// Percentage change in execution time
    pub time_change_percent: f64,
    /// Whether the change represents an improvement (lower time)
    pub is_improvement: bool,
}

impl fmt::Display for PerformanceComparison {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let change_symbol = if self.is_improvement { "↓" } else { "↑" };
        write!(
            f,
            "{} (size: {}): {:.3}ms → {:.3}ms ({}{:.1}%)",
            self.benchmark_name,
            self.problem_size,
            self.baseline_time * 1000.0,
            self.current_time * 1000.0,
            change_symbol,
            self.time_change_percent.abs()
        )
    }
}
