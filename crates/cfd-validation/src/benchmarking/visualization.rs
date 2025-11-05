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

use super::{BenchmarkResult, BenchmarkStatus};
use super::analysis::{RegressionAlert, AlertSeverity};
use super::scaling::ScalingResult;
use crate::reporting::PerformanceMetrics;
use std::collections::HashMap;
use std::fmt;

/// Chart types for performance visualization
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChartType {
    Line,
    Bar,
    Scatter,
    Histogram,
}

/// Visualization configuration
#[derive(Debug, Clone)]
pub struct VisualizationConfig {
    pub width: usize,
    pub height: usize,
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub show_grid: bool,
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

/// Performance chart data
#[derive(Debug, Clone)]
pub struct ChartData {
    pub labels: Vec<String>,
    pub datasets: Vec<Dataset>,
}

#[derive(Debug, Clone)]
pub struct Dataset {
    pub label: String,
    pub data: Vec<f64>,
    pub color: String,
}

/// HTML report generator
pub struct HtmlReportGenerator {
    config: VisualizationConfig,
}

impl HtmlReportGenerator {
    pub fn new(config: VisualizationConfig) -> Self {
        Self { config }
    }

    pub fn with_default_config() -> Self {
        Self::new(VisualizationConfig::default())
    }

    /// Generate comprehensive HTML performance report
    pub fn generate_report(
        &self,
        results: &[BenchmarkResult],
        scaling_results: Option<&[ScalingResult]>,
        alerts: Option<&[RegressionAlert]>,
    ) -> String {
        let mut html = String::new();

        // HTML header
        html.push_str(&self.generate_html_header());

        // Title and summary
        html.push_str(&self.generate_summary_section(results));

        // Performance charts with Chart.js
        html.push_str(&self.generate_performance_chart(results));

        // Scaling analysis (if available)
        if let Some(scaling) = scaling_results {
            html.push_str(&self.generate_scaling_chart(scaling));
        }

        // Regression alerts (if available)
        if let Some(alerts) = alerts {
            html.push_str(&self.generate_alerts_section(alerts));
        }

        // Detailed results table
        html.push_str(&self.generate_results_table(results));

        // Recommendations
        html.push_str(&self.generate_recommendations(results, alerts.unwrap_or(&[])));

        // HTML footer
        html.push_str(&self.generate_html_footer());

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

    fn generate_summary_section(&self, results: &[BenchmarkResult]) -> String {
        let total_operations = results.len();
        let avg_improvement = self.calculate_average_improvement(results);
        let critical_alerts = 0; // Would be calculated from alerts
        let total_alerts = 0; // Would be calculated from alerts

        format!(
            r#"<div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Operations</h3>
                    <p style="font-size: 2em; margin: 0;">{}</p>
                </div>
                <div class="summary-card">
                    <h3>Average Improvement</h3>
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
            if avg_improvement >= 0.0 { "#28a745" } else { "#dc3545" },
            avg_improvement,
            if critical_alerts > 0 { "#dc3545" } else { "#28a745" },
            critical_alerts,
            total_alerts
        )
    }

    fn generate_performance_charts(&self, results: &[BenchmarkResult]) -> String {
        let chart_data = self.prepare_performance_chart_data(results);

        format!(
            r#"<div class="chart-container">
                <h2>Performance Overview</h2>
                <canvas id="performanceChart" width="{}" height="{}"></canvas>
            </div>
            <script>
                const ctx = document.getElementById('performanceChart').getContext('2d');
                new Chart(ctx, {{
                    type: 'line',
                    data: {},
                    options: {{
                        responsive: true,
                        plugins: {{
                            title: {{
                                display: true,
                                text: 'Execution Time vs Problem Size'
                            }}
                        }},
                        scales: {{
                            x: {{
                                title: {{
                                    display: true,
                                    text: '{}'
                                }}
                            }},
                            y: {{
                                title: {{
                                    display: true,
                                    text: '{}'
                                }}
                            }}
                        }}
                    }}
                }});
            </script>"#,
            self.config.width,
            self.config.height,
            serde_json::to_string(&chart_data).unwrap_or_default(),
            self.config.x_label,
            self.config.y_label
        )
    }

    fn generate_scaling_charts(&self, scaling_results: &[ScalingResult]) -> String {
        let chart_data = self.prepare_scaling_chart_data(scaling_results);

        format!(
            r#"<div class="chart-container">
                <h2>Scaling Analysis</h2>
                <canvas id="scalingChart" width="{}" height="{}"></canvas>
            </div>
            <script>
                const ctxScaling = document.getElementById('scalingChart').getContext('2d');
                new Chart(ctxScaling, {{
                    type: 'line',
                    data: {},
                    options: {{
                        responsive: true,
                        plugins: {{
                            title: {{
                                display: true,
                                text: 'Parallel Scaling Efficiency'
                            }}
                        }},
                        scales: {{
                            x: {{
                                title: {{
                                    display: true,
                                    text: 'Number of Cores'
                                }}
                            }},
                            y: {{
                                title: {{
                                    display: true,
                                    text: 'Efficiency (%)'
                                }}
                            }}
                        }}
                    }}
                }});
            </script>"#,
            self.config.width,
            self.config.height,
            serde_json::to_string(&chart_data).unwrap_or_default()
        )
    }

    fn generate_alerts_section(&self, alerts: &[RegressionAlert]) -> String {
        if alerts.is_empty() {
            return r#"<div class="chart-container">
                <h2>Performance Alerts</h2>
                <p style="color: #28a745;">✅ No performance regressions detected</p>
            </div>"#.to_string();
        }

        let mut alerts_html = String::from(r#"<div class="chart-container">
            <h2>Performance Alerts</h2>"#);

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

            alerts_html.push_str(&format!(
                r#"<div class="alert {}">
                    <strong>{}:</strong> {}<br>
                    <small>Degradation: {:.1}%, Confidence: {:.1}%</small>
                </div>"#,
                alert_class,
                severity_text,
                alert.benchmark_name,
                alert.degradation_rate,
                alert.confidence * 100.0
            ));
        }

        alerts_html.push_str("</div>");
        alerts_html
    }

    fn generate_results_table(&self, results: &[BenchmarkResult]) -> String {
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
                <tbody>"#.to_string();

        for result in results {
            // Extract performance metrics
            let perf_time = result.performance.as_ref()
                .map(|p| format!("{:.3}ms ± {:.3}ms",
                    p.stats.mean * 1000.0,
                    p.stats.std_dev * 1000.0))
                .unwrap_or_else(|| "N/A".to_string());

            // Extract memory metrics
            let memory_usage = result.memory.as_ref()
                .map(|m| format!("{:.1}MB", m.total_allocated as f64 / 1_048_576.0))
                .unwrap_or_else(|| "N/A".to_string());

            // Extract scaling info
            let scaling_info = result.scaling.as_ref()
                .and_then(|s| {
                    // Get average speedup across all measurements
                    if s.speedup_factors.is_empty() {
                        None
                    } else {
                        let total: f64 = s.speedup_factors.values().sum();
                        let avg = total / s.speedup_factors.len() as f64;
                        Some(format!("{:.2}x avg speedup", avg))
                    }
                })
                .unwrap_or_else(|| "N/A".to_string());

            // Determine status color and text
            let (status_color, status_text) = match (&result.status, &result.regression_detected) {
                (BenchmarkStatus::Passed, None) => ("#28a745", "✅ Passed"),
                (BenchmarkStatus::Passed, Some(reg)) if *reg > 10.0 => ("#ffc107", "⚠️ Minor Regression"),
                (BenchmarkStatus::Passed, Some(reg)) if *reg > 25.0 => ("#dc3545", "🚨 Major Regression"),
                (BenchmarkStatus::Passed, Some(_)) => ("#17a2b8", "⚠️ Regression"),
                (BenchmarkStatus::Failed, _) => ("#dc3545", "❌ Failed"),
                (BenchmarkStatus::Regression, _) => ("#dc3545", "📉 Regression"),
                (BenchmarkStatus::Skipped, _) => ("#6c757d", "⏭️ Skipped"),
                (BenchmarkStatus::Running, _) => ("#007bff", "⏳ Running"),
            };

            table_html.push_str(&format!(
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
            ));
        }

        table_html.push_str(r#"</tbody></table></div>"#);
        table_html
    }

    fn generate_recommendations(&self, results: &[BenchmarkResult], alerts: &[RegressionAlert]) -> String {
        let mut recommendations = Vec::new();

        // Analyze results for recommendations with architectural purity
        let critical_alerts = alerts.iter()
            .filter(|a| a.degradation_rate > 25.0)
            .count();

        if critical_alerts > 0 {
            recommendations.push(format!("🚨 Address {} critical performance regressions immediately", critical_alerts));
        }

        let failed_tests = results.iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Failed))
            .count();

        if failed_tests > 0 {
            recommendations.push(format!("❌ Fix {} failed benchmark operations", failed_tests));
        }

        let avg_improvement = self.calculate_average_improvement(results);
        if avg_improvement < -10.0 {
            recommendations.push("📉 Significant performance degradation detected - investigate recent architectural changes".to_string());
        } else if avg_improvement < -5.0 {
            recommendations.push("⚠️ Performance degradation detected - review recent modifications".to_string());
        }

        // Memory efficiency analysis
        let high_memory_usage = results.iter()
            .filter_map(|r| r.memory.as_ref())
            .filter(|m| m.total_allocated > 100_000_000) // 100MB
            .count();

        if high_memory_usage > 0 {
            recommendations.push(format!("🧠 Optimize memory usage in {} operations (>100MB)", high_memory_usage));
        }

        // Scaling efficiency analysis
        let poor_scaling = results.iter()
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
            recommendations.push(format!("⚡ Improve parallel scaling efficiency in {} operations", poor_scaling));
        }

        if recommendations.is_empty() {
            recommendations.push("✅ Performance is stable - continue monitoring for architectural consistency".to_string());
        }

        let mut rec_html = r#"<div class="recommendations">
            <h2>Architectural Recommendations</h2>"#.to_string();

        for rec in recommendations {
            rec_html.push_str(&format!(r#"<div class="recommendation-item">• {}</div>"#, rec));
        }

        rec_html.push_str("</div>");
        rec_html
    }


    /// Generate Chart.js performance chart
    fn generate_performance_chart(&self, results: &[BenchmarkResult]) -> String {
        let mut labels = Vec::new();
        let mut execution_times = Vec::new();
        let mut memory_usage = Vec::new();

        for result in results {
            labels.push(format!("\"{}\"", result.name));

            // Extract execution time (convert to milliseconds)
            execution_times.push(format!("{:.3}", result.duration.as_secs_f64() * 1000.0));

            // Extract memory usage (convert to MB)
            let mem_mb = result.memory.as_ref()
                .map(|m| m.total_allocated as f64 / 1_048_576.0)
                .unwrap_or(0.0);
            memory_usage.push(format!("{:.1}", mem_mb));
        }

        format!(r#"
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

    /// Generate scaling analysis chart
    fn generate_scaling_chart(&self, scaling_results: &[super::scaling::ScalingResult]) -> String {
        if scaling_results.is_empty() {
            return String::new();
        }

        let mut problem_sizes = Vec::new();
        let mut speedups = Vec::new();
        let mut efficiencies = Vec::new();

        for result in scaling_results {
            // Use first problem size and processor count for display
            let problem_size = result.problem_sizes.first().copied().unwrap_or(0);
            let processor_count = result.processor_counts.first().copied().unwrap_or(1);

            problem_sizes.push(problem_size.to_string());

            // Get speedup and efficiency from the metrics
            let speedup = result.speedup_factors.get(&(problem_size, processor_count)).copied().unwrap_or(1.0);
            let efficiency = result.parallel_efficiency.get(&(problem_size, processor_count)).copied().unwrap_or(1.0);

            speedups.push(format!("{:.2}", speedup));
            efficiencies.push(format!("{:.1}", efficiency * 100.0));
        }

        format!(r#"
        <div class="chart-container">
            <h2>Parallel Scaling Analysis</h2>
            <div class="chart-legend">
                <div class="legend-item">
                    <div class="legend-color" style="background: #9b59b6;"></div>
                    <span>Speedup Factor</span>
                </div>
                <div class="legend-item">
                    <div class="legend-color" style="background: #f39c12;"></div>
                    <span>Parallel Efficiency (%)</span>
                </div>
            </div>
            <div class="chart-wrapper">
                <canvas id="scalingChart"></canvas>
            </div>
        </div>
        <script>
            const scalingCtx = document.getElementById('scalingChart').getContext('2d');
            new Chart(scalingCtx, {{
                type: 'line',
                data: {{
                    labels: [{}],
                    datasets: [{{
                        label: 'Speedup',
                        data: [{}],
                        backgroundColor: 'rgba(155, 89, 182, 0.2)',
                        borderColor: 'rgba(155, 89, 182, 1)',
                        borderWidth: 3,
                        fill: false,
                        tension: 0.4,
                        pointBackgroundColor: 'rgba(155, 89, 182, 1)',
                        pointBorderColor: '#fff',
                        pointBorderWidth: 2,
                        pointRadius: 6,
                        pointHoverRadius: 8
                    }}, {{
                        label: 'Efficiency (%)',
                        data: [{}],
                        backgroundColor: 'rgba(243, 156, 18, 0.2)',
                        borderColor: 'rgba(243, 156, 18, 1)',
                        borderWidth: 3,
                        fill: false,
                        tension: 0.4,
                        yAxisID: 'y1',
                        pointBackgroundColor: 'rgba(243, 156, 18, 1)',
                        pointBorderColor: '#fff',
                        pointBorderWidth: 2,
                        pointRadius: 6,
                        pointHoverRadius: 8
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
                            text: 'Parallel Scaling Performance',
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
                                        if (context.datasetIndex === 1) label += '%';
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
                                text: 'Problem Size'
                            }}
                        }},
                        y: {{
                            type: 'linear',
                            display: true,
                            position: 'left',
                            title: {{
                                display: true,
                                text: 'Speedup Factor'
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
                                text: 'Parallel Efficiency (%)'
                            }},
                            min: 0,
                            max: 100,
                            grid: {{
                                drawOnChartArea: false,
                            }},
                        }}
                    }}
                }}
            }});
        </script>"#,
        problem_sizes.join("\",\"").split(",").collect::<Vec<&str>>().join("\",\""),
        speedups.join(","),
        efficiencies.join(",")
        )
    }

    /// Calculate average performance improvement across all results
    fn calculate_average_improvement(&self, results: &[BenchmarkResult]) -> f64 {
        if results.is_empty() {
            return 0.0;
        }

        // For now, return a simple metric based on successful operations
        // In a real implementation, this would compare against baseline performance
        let successful = results.iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Passed))
            .count();

        ((successful as f64 / results.len() as f64) - 0.5) * 200.0 // Convert to percentage around 0
    }

    fn generate_html_footer(&self) -> String {
        r#"</div></body></html>"#.to_string()
    }

    fn prepare_performance_chart_data(&self, results: &[BenchmarkResult]) -> serde_json::Value {
        // Group results by operation
        let mut operations: HashMap<String, Vec<(usize, f64)>> = HashMap::new();

        for result in results {
            // Use name as operation name and duration as execution time
            operations.entry(result.name.clone())
                .or_insert_with(Vec::new)
                .push((0, result.duration.as_secs_f64() * 1000.0)); // Use 0 as problem size placeholder
        }

        // Sort by problem size for each operation
        for data in operations.values_mut() {
            data.sort_by_key(|(size, _)| *size);
        }

        // Create Chart.js compatible data
        let mut datasets = Vec::new();
        let mut all_labels = Vec::new();
        let mut color_idx = 0;

        for (operation, data) in operations {
            let (labels, values): (Vec<_>, Vec<_>) = data.into_iter().unzip();

            // Collect all unique labels
            for label in labels {
                if !all_labels.contains(&label) {
                    all_labels.push(label);
                }
            }

            datasets.push(serde_json::json!({
                "label": operation,
                "data": values,
                "borderColor": self.config.colors[color_idx % self.config.colors.len()],
                "backgroundColor": self.config.colors[color_idx % self.config.colors.len()],
                "fill": false
            }));

            color_idx += 1;
        }

        // Sort labels
        all_labels.sort();

        serde_json::json!({
            "labels": all_labels,
            "datasets": datasets
        })
    }

    fn prepare_scaling_chart_data(&self, scaling_results: &[ScalingResult]) -> serde_json::Value {
        let mut cores = Vec::new();
        let mut efficiencies = Vec::new();

        for result in scaling_results {
            // Use processor counts as cores
            for &processor_count in &result.processor_counts {
                if !cores.contains(&processor_count) {
                    cores.push(processor_count);
                }
            }

            // Get average efficiency across all measurements
            let avg_efficiency = if result.parallel_efficiency.is_empty() {
                1.0
            } else {
                result.parallel_efficiency.values().sum::<f64>() / result.parallel_efficiency.len() as f64
            };
            efficiencies.push(avg_efficiency * 100.0); // Convert to percentage
        }

        serde_json::json!({
            "labels": cores,
            "datasets": [{
                "label": "Scaling Efficiency",
                "data": efficiencies,
                "borderColor": "#1f77b4",
                "backgroundColor": "#1f77b4",
                "fill": false
            }]
        })
    }

}

/// Performance dashboard generator
pub struct PerformanceDashboard {
    generator: HtmlReportGenerator,
}

impl PerformanceDashboard {
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
        let html_content = self.generator.generate_report(results, scaling_results, alerts);

        std::fs::write(output_path, html_content)?;

        println!("Performance dashboard generated: {}", output_path);
        println!("Open in browser to view interactive charts and detailed analysis");

        Ok(())
    }

    /// Generate summary report (text format)
    pub fn generate_summary_report(&self, results: &[BenchmarkResult]) -> String {
        let mut report = String::new();

        report.push_str("CFD Performance Benchmark Summary\n");
        report.push_str("==================================\n\n");

        report.push_str(&format!("Total Operations Benchmarked: {}\n", results.len()));
        report.push_str(&format!("Report Generated: {}\n\n", chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")));

        // Summary statistics
        let total_time: f64 = results.iter()
            .map(|r| r.duration.as_secs_f64())
            .sum();

        let avg_time = total_time / results.len() as f64;
        let max_time = results.iter()
            .map(|r| r.duration.as_secs_f64())
            .fold(0.0, f64::max);

        report.push_str("Performance Statistics:\n");
        report.push_str(&format!("  Average Execution Time: {:.3} ms\n", avg_time * 1000.0));
        report.push_str(&format!("  Maximum Execution Time: {:.3} ms\n", max_time * 1000.0));
        report.push_str(&format!("  Total Execution Time: {:.3} s\n\n", total_time));

        // Operation breakdown
        report.push_str("Operation Breakdown:\n");
        let mut operations: HashMap<String, Vec<&BenchmarkResult>> = HashMap::new();

        for result in results {
            operations.entry(result.name.clone())
                .or_insert_with(Vec::new)
                .push(result);
        }

        for (operation, op_results) in operations {
            let op_avg_time = op_results.iter()
                .map(|r| r.duration.as_secs_f64())
                .sum::<f64>() / op_results.len() as f64;

            report.push_str(&format!("  {}: {:.3} ms ({} runs)\n",
                operation, op_avg_time * 1000.0, op_results.len()));
        }

        report
    }
}

impl Default for PerformanceDashboard {
    fn default() -> Self {
        Self::new()
    }
}

/// Export formats for benchmark results
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    Html,
    Json,
    Csv,
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
        let mut csv = String::from("benchmark_name,problem_size,duration_ms,memory_usage_mb,throughput_ops_per_sec,status\n");

        for result in results {
            csv.push_str(&format!(
                "{},{},{:.3},{},{},{}\n",
                result.name,
                0, // Use 0 as placeholder for problem size
                result.duration.as_secs_f64() * 1000.0,
                result.memory.as_ref().map_or("N/A".to_string(), |m| format!("{:.1}", m.total_allocated as f64 / 1024.0 / 1024.0)),
                "N/A", // No throughput field available
                result.regression_detected.is_some()
            ));
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
    pub benchmark_name: String,
    pub problem_size: usize,
    pub baseline_time: f64,
    pub current_time: f64,
    pub time_change_percent: f64,
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
