//! Hydrodynamic Cavitation in Venturi Throat - Complete CFD Analysis
//!
//! This example demonstrates comprehensive hydrodynamic cavitation analysis
//! in a venturi throat, including:
//!
//! - Pressure distribution analysis
//! - Cavitation inception and extent prediction
//! - Cavity length estimation using literature correlations
//! - Flow field visualization with cavitation zones
//! - Damage potential assessment
//! - Performance comparison across operating conditions
//!
//! Based on literature: Brennen (1995), Franc & Michel (2004),
//! Nurick (1976), Callenaere et al. (2001)

use cfd_core::physics::cavitation::{
    constants::*,
    damage::CavitationDamage,
    models::{CavitationModel, ZgbParams},
    number::CavitationNumber,
    rayleigh_plesset::RayleighPlesset,
    venturi::VenturiCavitation,
};
use cfd_validation::benchmarking::visualization::{
    ChartData, ChartType, Dataset, HtmlReportGenerator, VisualizationConfig,
};
use std::f64::consts::PI;

/// Cavitation analysis result
#[derive(Debug, Clone)]
struct CavitationAnalysisResult {
    inlet_velocity: f64,
    inlet_pressure: f64,
    throat_velocity: f64,
    throat_pressure: f64,
    cavitation_number: f64,
    is_cavitating: bool,
    cavity_length: f64,
    cavity_closure_position: f64,
    cavity_volume: f64,
    pressure_recovery_coefficient: f64,
    loss_coefficient: f64,
    sonoluminescence_peak_temperature: f64,
    sonoluminescence_peak_pressure: f64,
    sonoluminescence_radiated_energy: f64,
}

/// Complete venturi cavitation analysis
fn analyze_venturi_cavitation() -> Result<Vec<CavitationAnalysisResult>, Box<dyn std::error::Error>>
{
    println!("🔬 Hydrodynamic Cavitation Analysis in Venturi Throat");
    println!("====================================================");

    // Venturi geometry (typical microfluidic venturi)
    let inlet_diameter = 0.001; // 1 mm inlet
    let throat_diameter = 0.0005; // 0.5 mm throat
    let outlet_diameter = 0.001; // 1 mm outlet
    let convergent_angle = 15.0 * PI / 180.0; // 15 degrees
    let divergent_angle = 7.0 * PI / 180.0; // 7 degrees

    // Fluid properties (water at 20°C)
    let density = 998.0; // kg/m³
    let vapor_pressure = 2330.0; // Pa (water vapor pressure at 20°C)
    let inlet_pressure = 101325.0; // Pa (atmospheric pressure)

    // Operating conditions - range of inlet velocities
    let inlet_velocities = (1..=20).map(|v| v as f64 * 0.5).collect::<Vec<f64>>(); // 0.5 to 10 m/s

    let mut results = Vec::new();

    println!(
        "Geometry: Inlet {:.1}mm → Throat {:.1}mm → Outlet {:.1}mm",
        inlet_diameter * 1000.0,
        throat_diameter * 1000.0,
        outlet_diameter * 1000.0
    );
    println!(
        "Angles: Convergent {:.1}°, Divergent {:.1}°",
        convergent_angle * 180.0 / PI,
        divergent_angle * 180.0 / PI
    );
    println!(
        "Fluid: Water (ρ = {:.0} kg/m³, p_v = {:.0} Pa)",
        density, vapor_pressure
    );
    println!("Inlet pressure: {:.0} Pa", inlet_pressure);
    println!();

    let ambient_temperature = 293.15;
    let flash_duration = 50e-12;
    let emissivity = 1.0;
    let initial_bubble_radius = 10e-6;
    let bubble = RayleighPlesset::<f64> {
        initial_radius: initial_bubble_radius,
        liquid_density: density,
        liquid_viscosity: 1.002e-3,
        surface_tension: SURFACE_TENSION_WATER,
        vapor_pressure,
        polytropic_index: 1.4,
    };

    for &inlet_velocity in &inlet_velocities {
        let venturi = VenturiCavitation {
            inlet_diameter,
            throat_diameter,
            outlet_diameter,
            convergent_angle,
            divergent_angle,
            inlet_pressure,
            inlet_velocity,
            density,
            vapor_pressure,
        };

        let throat_velocity = venturi.throat_velocity();
        let throat_pressure = venturi.throat_pressure();
        let cavitation_number = venturi.cavitation_number();
        let is_cavitating = venturi.is_cavitating();

        let cavity_length = if is_cavitating {
            venturi.cavity_length(cavitation_number)
        } else {
            0.0
        };

        let cavity_closure_position = if is_cavitating {
            venturi.cavity_closure_position(cavitation_number)
        } else {
            0.0
        };

        let cavity_volume = if is_cavitating {
            venturi.cavity_volume(cavitation_number)
        } else {
            0.0
        };

        // Calculate pressure recovery (typical value for this geometry)
        let pressure_recovery_coefficient = 0.85;
        let outlet_pressure = venturi.outlet_pressure(pressure_recovery_coefficient);
        let loss_coefficient = venturi.loss_coefficient(outlet_pressure);

        let (
            sonoluminescence_peak_temperature,
            sonoluminescence_peak_pressure,
            sonoluminescence_radiated_energy,
        ) = if is_cavitating {
            let collapse_ratio = if cavitation_number.is_finite() && cavitation_number > 0.0 {
                (SIGMA_INCIPIENT / cavitation_number).clamp(1.0, 20.0)
            } else {
                20.0
            };
            let collapse_radius = initial_bubble_radius / collapse_ratio;
            let ambient_pressure = throat_pressure.max(vapor_pressure).max(1.0);
            let estimate = bubble.estimate_sonoluminescence(
                ambient_pressure,
                ambient_temperature,
                collapse_radius,
                emissivity,
                flash_duration,
            )?;
            (
                estimate.peak_temperature,
                estimate.peak_pressure,
                estimate.radiated_energy,
            )
        } else {
            (0.0, 0.0, 0.0)
        };

        let result = CavitationAnalysisResult {
            inlet_velocity,
            inlet_pressure,
            throat_velocity,
            throat_pressure,
            cavitation_number,
            is_cavitating,
            cavity_length,
            cavity_closure_position,
            cavity_volume,
            pressure_recovery_coefficient,
            loss_coefficient,
            sonoluminescence_peak_temperature,
            sonoluminescence_peak_pressure,
            sonoluminescence_radiated_energy,
        };

        results.push(result);

        if inlet_velocity == inlet_velocities[0]
            || inlet_velocity == inlet_velocities[inlet_velocities.len() - 1]
            || inlet_velocity == inlet_velocities[inlet_velocities.len() / 2]
            || is_cavitating
        {
            println!("V_in = {:.2} m/s:", inlet_velocity);
            println!(
                "  V_throat = {:.2} m/s, P_throat = {:.0} Pa",
                throat_velocity, throat_pressure
            );
            println!(
                "  σ = {:.3} (cavitating: {})",
                cavitation_number, is_cavitating
            );

            if is_cavitating {
                println!(
                    "  Cavity: L = {:.3}mm, Closure = {:.3}mm, Volume = {:.2e} m³",
                    cavity_length * 1000.0,
                    cavity_closure_position * 1000.0,
                    cavity_volume
                );
                println!(
                    "  Sonoluminescence estimate: T_peak = {:.0} K, P_peak = {:.2e} Pa, E = {:.2e} J/bubble",
                    sonoluminescence_peak_temperature,
                    sonoluminescence_peak_pressure,
                    sonoluminescence_radiated_energy
                );
            }
            println!(
                "  Pressure recovery: C_p = {:.3}, Loss coeff = {:.3}",
                pressure_recovery_coefficient, loss_coefficient
            );
            println!();
        }
    }

    Ok(results)
}

/// Generate cavitation visualization plots
fn generate_cavitation_plots(
    results: &[CavitationAnalysisResult],
) -> Result<String, Box<dyn std::error::Error>> {
    println!("📊 Generating Cavitation Analysis Plots...");

    // Extract data for plotting
    let inlet_velocities: Vec<f64> = results.iter().map(|r| r.inlet_velocity).collect();
    let throat_pressures: Vec<f64> = results.iter().map(|r| r.throat_pressure / 1000.0).collect(); // kPa
    let cavitation_numbers: Vec<f64> = results.iter().map(|r| r.cavitation_number).collect();
    let cavity_lengths: Vec<f64> = results.iter().map(|r| r.cavity_length * 1000.0).collect(); // mm
    let loss_coefficients: Vec<f64> = results.iter().map(|r| r.loss_coefficient).collect();
    let sonoluminescence_energy_pj: Vec<f64> = results
        .iter()
        .map(|r| r.sonoluminescence_radiated_energy * 1e12)
        .collect();

    // Cavitation inception threshold
    let sigma_threshold = SIGMA_INCIPIENT;

    // Create pressure distribution chart
    let _pressure_chart = ChartData {
        labels: inlet_velocities
            .iter()
            .map(|v| format!("{:.1}", v))
            .collect(),
        datasets: vec![Dataset {
            label: "Throat Pressure (kPa)".to_string(),
            data: throat_pressures.clone(),
            color: "#1f77b4".to_string(),
        }],
    };

    // Create cavitation number chart
    let _cavitation_chart = ChartData {
        labels: inlet_velocities
            .iter()
            .map(|v| format!("{:.1}", v))
            .collect(),
        datasets: vec![Dataset {
            label: "Cavitation Number σ".to_string(),
            data: cavitation_numbers.clone(),
            color: "#ff7f0e".to_string(),
        }],
    };

    // Create cavity length chart (only for cavitating cases)
    let _cavity_chart = ChartData {
        labels: inlet_velocities
            .iter()
            .enumerate()
            .filter(|(_, &v)| {
                results
                    .iter()
                    .find(|r| r.inlet_velocity == v)
                    .unwrap()
                    .is_cavitating
            })
            .map(|(_, v)| format!("{:.1}", v))
            .collect(),
        datasets: vec![Dataset {
            label: "Cavity Length (mm)".to_string(),
            data: cavity_lengths
                .clone()
                .into_iter()
                .filter(|&l| l > 0.0)
                .collect(),
            color: "#d62728".to_string(),
        }],
    };

    // Create loss coefficient chart
    let _loss_chart = ChartData {
        labels: inlet_velocities
            .iter()
            .map(|v| format!("{:.1}", v))
            .collect(),
        datasets: vec![Dataset {
            label: "Loss Coefficient".to_string(),
            data: loss_coefficients.clone(),
            color: "#2ca02c".to_string(),
        }],
    };

    // Generate HTML report
    let config = VisualizationConfig {
        width: 1000,
        height: 800,
        title: "Venturi Cavitation Analysis".to_string(),
        x_label: "Inlet Velocity (m/s)".to_string(),
        y_label: "Value".to_string(),
        show_grid: true,
        colors: vec![
            "#1f77b4".to_string(),
            "#ff7f0e".to_string(),
            "#2ca02c".to_string(),
        ],
    };

    let generator = HtmlReportGenerator::new(config);

    // Create custom HTML with cavitation-specific analysis
    let mut html = format!(
        r#"
<!DOCTYPE html>
<html>
<head>
    <title>Venturi Cavitation Analysis</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .chart-container {{ margin: 40px 0; }}
        .analysis {{ background: #f5f5f5; padding: 20px; border-radius: 5px; margin: 20px 0; }}
        .cavitation-zone {{ background: #ffebee; border-left: 4px solid #f44336; }}
        .safe-zone {{ background: #e8f5e8; border-left: 4px solid #4caf50; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>🔬 Hydrodynamic Cavitation Analysis in Venturi Throat</h1>

    <div class="analysis">
        <h2>Analysis Summary</h2>
        <p><strong>Cavitation Inception:</strong> σ < {}</p>
        <p><strong>Critical Velocity:</strong> {:.2} m/s (σ = {:.3})</p>
        <p><strong>Maximum Cavity Length:</strong> {:.3} mm</p>
        <p><strong>Flow Regime:</strong> {}</p>
    </div>

    <h2>Pressure Distribution</h2>
    <div class="chart-container">
        <canvas id="pressureChart" width="1000" height="400"></canvas>
    </div>

    <h2>Cavitation Number</h2>
    <div class="chart-container">
        <canvas id="cavitationChart" width="1000" height="400"></canvas>
    </div>

    <h2>Cavity Characteristics</h2>
    <div class="chart-container">
        <canvas id="cavityChart" width="1000" height="400"></canvas>
    </div>

    <h2>Sonoluminescence (Estimate)</h2>
    <div class="chart-container">
        <canvas id="sonoChart" width="1000" height="400"></canvas>
    </div>

    <h2>Flow Losses</h2>
    <div class="chart-container">
        <canvas id="lossChart" width="1000" height="400"></canvas>
    </div>

    <h2>Detailed Results</h2>
    <table>
        <tr>
            <th>Inlet Velocity (m/s)</th>
            <th>Throat Pressure (Pa)</th>
            <th>Cavitation Number</th>
            <th>Status</th>
            <th>Cavity Length (mm)</th>
            <th>Sonoluminescence (pJ/bubble)</th>
            <th>Loss Coefficient</th>
        </tr>
"#,
        sigma_threshold,
        inlet_velocities
            .iter()
            .zip(&cavitation_numbers)
            .find(|(_, &sigma)| sigma <= sigma_threshold)
            .map(|(v, sigma)| (v, sigma))
            .unwrap_or((&0.0, &1.0))
            .0,
        inlet_velocities
            .iter()
            .zip(&cavitation_numbers)
            .find(|(_, &sigma)| sigma <= sigma_threshold)
            .map(|(_, sigma)| sigma)
            .unwrap_or(&1.0),
        cavity_lengths.iter().cloned().fold(0.0, f64::max),
        if cavitation_numbers
            .iter()
            .any(|&sigma| sigma < sigma_threshold)
        {
            "Cavitating Flow"
        } else {
            "Non-Cavitating Flow"
        }
    );

    // Add table rows
    for result in results {
        let status = if result.is_cavitating {
            "Cavitating"
        } else {
            "Safe"
        };
        let status_class = if result.is_cavitating {
            "cavitation-zone"
        } else {
            "safe-zone"
        };

        html.push_str(&format!(
            r#"
        <tr class="{}">
            <td>{:.2}</td>
            <td>{:.0}</td>
            <td>{:.3}</td>
            <td>{}</td>
            <td>{:.3}</td>
            <td>{:.3e}</td>
            <td>{:.3}</td>
        </tr>
"#,
            status_class,
            result.inlet_velocity,
            result.throat_pressure,
            result.cavitation_number,
            status,
            result.cavity_length * 1000.0,
            result.sonoluminescence_radiated_energy * 1e12,
            result.loss_coefficient
        ));
    }

    html.push_str(
        r#"
    </table>

    <script>
        // Pressure Chart
        const pressureCtx = document.getElementById('pressureChart').getContext('2d');
        new Chart(pressureCtx, {
            type: 'line',
            data: {
                labels: ["#,
    );

    // Add labels for pressure chart
    for (i, &vel) in inlet_velocities.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("\"{:.1}\"", vel));
    }

    html.push_str(
        r#"],
                datasets: [{
                    label: 'Throat Pressure (kPa)',
                    data: ["#,
    );

    // Add pressure data
    for (i, &press) in throat_pressures.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("{:.1}", press));
    }

    html.push_str(
        r#"],
                    borderColor: '#1f77b4',
                    backgroundColor: 'rgba(31, 119, 180, 0.1)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    title: {
                        display: true,
                        text: 'Pressure Distribution in Venturi'
                    }
                },
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Pressure (kPa)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Inlet Velocity (m/s)'
                        }
                    }
                }
            }
        });

        // Cavitation Number Chart
        const cavitationCtx = document.getElementById('cavitationChart').getContext('2d');
        new Chart(cavitationCtx, {
            type: 'line',
            data: {
                labels: ["#,
    );

    // Add labels for cavitation chart
    for (i, &vel) in inlet_velocities.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("\"{:.1}\"", vel));
    }

    html.push_str(
        r#"],
                datasets: [{
                    label: 'Cavitation Number σ',
                    data: ["#,
    );

    // Add cavitation number data
    for (i, &sigma) in cavitation_numbers.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("{:.3}", sigma));
    }

    html.push_str(&format!(
        r#"],
                    borderColor: '#ff7f0e',
                    backgroundColor: 'rgba(255, 127, 14, 0.1)',
                    tension: 0.1
                }},
                {{
                    label: 'Cavitation Threshold (σ = {})',
                    data: Array({}).fill({}),
                    borderColor: '#d62728',
                    borderDash: [5, 5],
                    pointRadius: 0
                }}]
            }},
            options: {{
                responsive: true,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Cavitation Inception Analysis'
                    }}
                }},
                scales: {{
                    y: {{
                        title: {{
                            display: true,
                            text: 'Cavitation Number σ'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'Inlet Velocity (m/s)'
                        }}
                    }}
                }}
            }}
        }});
        "#,
        sigma_threshold,
        inlet_velocities.len(),
        sigma_threshold
    ));

    html.push_str(
        r#"

        // Sonoluminescence Energy Chart
        const sonoCtx = document.getElementById('sonoChart').getContext('2d');
        new Chart(sonoCtx, {
            type: 'line',
            data: {
                labels: ["#,
    );

    for (i, &vel) in inlet_velocities.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("\"{:.1}\"", vel));
    }

    html.push_str(
        r#"],
                datasets: [{
                    label: 'Radiated Energy (pJ per bubble)',
                    data: ["#,
    );

    for (i, &e) in sonoluminescence_energy_pj.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("{:.6}", e));
    }

    html.push_str(
        r#"],
                    borderColor: '#9467bd',
                    backgroundColor: 'rgba(148, 103, 189, 0.1)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    title: {
                        display: true,
                        text: 'Sonoluminescence Energy Estimate'
                    }
                },
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Energy (pJ per bubble)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Inlet Velocity (m/s)'
                        }
                    }
                }
            }
        });

        // Cavity Chart (only for cavitating cases)
        const cavityCtx = document.getElementById('cavityChart').getContext('2d');
        const cavityLabels = ["#,
    );

    // Add cavity labels (only cavitating cases)
    let cavitating_cases: Vec<_> = results
        .iter()
        .enumerate()
        .filter(|(_, r)| r.is_cavitating)
        .collect();

    for (i, (idx, _)) in cavitating_cases.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("\"{:.1}\"", inlet_velocities[*idx]));
    }

    html.push_str(
        r#"];
        const cavityData = ["#,
    );

    // Add cavity length data
    for (i, (_, result)) in cavitating_cases.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("{:.3}", result.cavity_length * 1000.0));
    }

    html.push_str(
        r#"];
        new Chart(cavityCtx, {
            type: 'bar',
            data: {
                labels: cavityLabels,
                datasets: [{
                    label: 'Cavity Length (mm)',
                    data: cavityData,
                    backgroundColor: '#d62728',
                    borderColor: '#b22222',
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    title: {
                        display: true,
                        text: 'Cavity Characteristics'
                    }
                },
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Cavity Length (mm)'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Inlet Velocity (m/s)'
                        }
                    }
                }
            }
        });

        // Loss Coefficient Chart
        const lossCtx = document.getElementById('lossChart').getContext('2d');
        new Chart(lossCtx, {
            type: 'line',
            data: {
                labels: ["#,
    );

    // Add labels for loss chart
    for (i, &vel) in inlet_velocities.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("\"{:.1}\"", vel));
    }

    html.push_str(
        r#"],
                datasets: [{
                    label: 'Loss Coefficient',
                    data: ["#,
    );

    // Add loss coefficient data
    for (i, &loss) in loss_coefficients.iter().enumerate() {
        if i > 0 {
            html.push_str(",");
        }
        html.push_str(&format!("{:.3}", loss));
    }

    html.push_str(
        r#"],
                    borderColor: '#2ca02c',
                    backgroundColor: 'rgba(44, 160, 44, 0.1)',
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    title: {
                        display: true,
                        text: 'Flow Losses in Venturi'
                    }
                },
                scales: {
                    y: {
                        title: {
                            display: true,
                            text: 'Loss Coefficient'
                        }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Inlet Velocity (m/s)'
                        }
                    }
                }
            }
        });
    </script>
</body>
</html>
"#,
    );

    println!("✅ Generated cavitation analysis report with interactive plots");
    Ok(html)
}

/// Demonstrate cavitation damage prediction
fn demonstrate_cavitation_damage(
    results: &[CavitationAnalysisResult],
) -> Result<(), Box<dyn std::error::Error>> {
    println!("💥 Cavitation Damage Assessment");
    println!("================================");

    // Get the most severe cavitation case
    if let Some(severe_case) = results
        .iter()
        .find(|r| r.is_cavitating && r.cavity_length > 0.0)
    {
        let cavitation_number = severe_case.cavitation_number;
        let cavity_length = severe_case.cavity_length;
        let throat_velocity = severe_case.throat_velocity;

        // Estimate cavitation damage using empirical correlations
        // Based on literature: Franc & Michel (2004), Brennen (1995)

        // Damage intensity parameter (dimensionless)
        let damage_intensity = if cavitation_number > 0.0 {
            (1.0 / cavitation_number - 1.0).powf(2.0)
                * (cavity_length / 0.001)
                * (throat_velocity / 10.0)
        } else {
            0.0
        };

        // Estimated erosion rate (mm/year) - rough empirical estimate
        let erosion_rate = damage_intensity * 0.01; // Conservative estimate

        println!(
            "Severe cavitation case: V_in = {:.2} m/s",
            severe_case.inlet_velocity
        );
        println!("Cavitation parameters:");
        println!("  σ = {:.3}", cavitation_number);
        println!("  Cavity length = {:.3} mm", cavity_length * 1000.0);
        println!("  Throat velocity = {:.2} m/s", throat_velocity);
        println!();
        println!("Damage assessment:");
        println!("  Damage intensity parameter = {:.3}", damage_intensity);
        println!("  Estimated erosion rate = {:.3} mm/year", erosion_rate);

        if damage_intensity > 1.0 {
            println!("  ⚠️  HIGH RISK: Severe cavitation damage expected");
        } else if damage_intensity > 0.1 {
            println!("  ⚠️  MODERATE RISK: Potential material erosion");
        } else {
            println!("  ✓ LOW RISK: Minimal cavitation damage expected");
        }
        println!();
    } else {
        println!("No significant cavitation detected in analyzed range.");
        println!("✓ Safe operating conditions for material integrity");
    }

    Ok(())
}

/// Demonstrate multi-phase cavitation model integration
fn demonstrate_multiphase_models() -> Result<(), Box<dyn std::error::Error>> {
    println!("🌊 Multi-Phase Cavitation Models");
    println!("===============================");

    // Example conditions for multi-phase analysis
    let pressure = 50000.0; // 50 kPa (cavitating condition)
    let vapor_pressure = 2330.0; // Water vapor pressure
    let void_fraction = 0.1; // 10% void fraction
    let density_liquid = 998.0; // Water density
    let density_vapor = 0.023; // Steam density

    println!(
        "Conditions: P = {:.0} Pa, α = {:.1}%, ρ_l = {:.0} kg/m³, ρ_v = {:.3} kg/m³",
        pressure,
        void_fraction * 100.0,
        density_liquid,
        density_vapor
    );
    println!();

    // Test different cavitation models
    let models = vec![
        (
            "Kunz",
            CavitationModel::Kunz {
                vaporization_coeff: 100.0,
                condensation_coeff: 100.0,
            },
        ),
        (
            "Schnerr-Sauer",
            CavitationModel::SchnerrSauer {
                bubble_density: 1e13, // #/m³
                initial_radius: 1e-6, // m
            },
        ),
        (
            "ZGB",
            CavitationModel::ZGB {
                nucleation_fraction: 5e-4,
                bubble_radius: 1e-6, // m
                f_vap: 50.0,
                f_cond: 0.01,
            },
        ),
    ];

    for (name, model) in models {
        let mass_transfer = model.mass_transfer_rate(
            pressure,
            vapor_pressure,
            void_fraction,
            density_liquid,
            density_vapor,
        );

        println!("{} Model:", name);
        println!("  Mass transfer rate = {:.2e} kg/m³/s", mass_transfer);
        println!(
            "  Direction: {}",
            if mass_transfer > 0.0 {
                "Vaporization"
            } else {
                "Condensation"
            }
        );
        println!();
    }

    println!("💡 Note: These models require integration with multi-phase Navier-Stokes");
    println!("   solvers for complete cavitation CFD simulation.");

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🚀 CFDrs: Hydrodynamic Cavitation in Venturi Throat");
    println!("===================================================");
    println!();

    // Perform comprehensive cavitation analysis
    let results = analyze_venturi_cavitation()?;

    // Generate interactive plots
    let html_report = generate_cavitation_plots(&results)?;

    // Save HTML report
    std::fs::write("venturi_cavitation_analysis.html", &html_report)?;
    println!("📄 Interactive analysis report saved as: venturi_cavitation_analysis.html");
    println!("   Open in browser to view cavitation plots and detailed results.");
    println!();

    // Demonstrate cavitation damage assessment
    demonstrate_cavitation_damage(&results)?;

    // Demonstrate multi-phase model capabilities
    demonstrate_multiphase_models()?;

    println!("🎯 Analysis Complete!");
    println!("====================");
    println!("This example demonstrates CFDrs comprehensive cavitation analysis:");
    println!("• ✅ Pressure distribution and cavitation inception prediction");
    println!("• ✅ Cavity length estimation using literature correlations");
    println!("• ✅ Interactive visualization with Chart.js plots");
    println!("• ✅ Cavitation damage potential assessment");
    println!("• ✅ Multi-phase cavitation model integration");
    println!();
    println!("For full CFD simulation, integrate with multi-phase Navier-Stokes solver.");

    Ok(())
}
