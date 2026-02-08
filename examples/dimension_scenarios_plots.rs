use cfd_3d::vof::{
    AdvectionMethod, CavitationVofConfig, CavitationVofSolver, InterfaceReconstruction, VofConfig,
};
use cfd_core::physics::cavitation::{models::CavitationModel, venturi::VenturiCavitation};
use nalgebra::{DMatrix, Vector3};

struct DatasetSpec {
    label: String,
    data: Vec<f64>,
    color: String,
}

struct ChartSpec {
    canvas_id: String,
    heading: String,
    y_label: String,
    labels: Vec<String>,
    datasets: Vec<DatasetSpec>,
}

fn write_html(path: &std::path::Path, html: &str) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    std::fs::write(path, html)
}

fn chart_js_html(title: &str, charts: &[ChartSpec]) -> String {
    let mut html = String::new();

    html.push_str("<!DOCTYPE html><html><head><meta charset=\"utf-8\">");
    html.push_str(&format!("<title>{}</title>", title));
    html.push_str("<script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>");
    html.push_str("<style>");
    html.push_str("body{font-family:Arial,sans-serif;margin:32px;}");
    html.push_str(".chart{margin:24px 0;}");
    html.push_str("canvas{max-width:1200px;}");
    html.push_str("</style></head><body>");
    html.push_str(&format!("<h1>{}</h1>", title));

    for chart in charts {
        html.push_str("<div class=\"chart\">");
        html.push_str(&format!("<h2>{}</h2>", chart.heading));
        html.push_str(&format!(
            "<canvas id=\"{}\" width=\"1200\" height=\"450\"></canvas>",
            chart.canvas_id
        ));
        html.push_str(&format!("<div><strong>Y:</strong> {}</div>", chart.y_label));
        html.push_str("</div>");
    }

    html.push_str("<script>");
    for chart in charts {
        html.push_str(&format!(
            "const ctx_{id}=document.getElementById('{id}').getContext('2d');",
            id = chart.canvas_id
        ));

        html.push_str(&format!(
            "new Chart(ctx_{id},{{type:'line',data:{{labels:[",
            id = chart.canvas_id
        ));
        for (i, label) in chart.labels.iter().enumerate() {
            if i > 0 {
                html.push(',');
            }
            html.push('"');
            html.push_str(&label.replace('"', "\\\""));
            html.push('"');
        }
        html.push_str("],datasets:[");

        for (i, dataset) in chart.datasets.iter().enumerate() {
            if i > 0 {
                html.push(',');
            }
            html.push('{');
            html.push_str(&format!("label:'{}',", dataset.label.replace('\'', "\\'")));
            html.push_str("data:[");
            for (j, v) in dataset.data.iter().enumerate() {
                if j > 0 {
                    html.push(',');
                }
                html.push_str(&format!("{:.12}", v));
            }
            html.push_str("],");
            html.push_str(&format!("borderColor:'{}',", dataset.color));
            html.push_str(&format!("backgroundColor:'{}22',", dataset.color));
            html.push_str("tension:0.15,pointRadius:0,fill:false}");
        }

        html.push_str("]} ,options:{responsive:true,plugins:{legend:{display:true,position:'top'}},scales:{x:{title:{display:true,text:'Index'}},y:{title:{display:true,text:'");
        html.push_str(&chart.y_label.replace('\'', "\\'"));
        html.push_str("'}}}}});");
    }
    html.push_str("</script>");

    html.push_str("</body></html>");
    html
}

fn venturi_1d_sweep() -> (Vec<String>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let inlet_diameter = 0.001;
    let throat_diameter = 0.0005;
    let outlet_diameter = 0.001;
    let convergent_angle = 15.0 * std::f64::consts::PI / 180.0;
    let divergent_angle = 7.0 * std::f64::consts::PI / 180.0;
    let density = 998.0;
    let vapor_pressure = 2330.0;
    let inlet_pressure = 101325.0;

    let base = VenturiCavitation::<f64> {
        inlet_diameter,
        throat_diameter,
        outlet_diameter,
        convergent_angle,
        divergent_angle,
        inlet_pressure,
        inlet_velocity: 1.0,
        density,
        vapor_pressure,
    };

    let inlet_velocities = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 12.5, 15.0];
    let mut labels = Vec::with_capacity(inlet_velocities.len());
    let mut throat_pressures_kpa = Vec::with_capacity(inlet_velocities.len());
    let mut sigmas = Vec::with_capacity(inlet_velocities.len());
    let mut cavity_lengths_mm = Vec::with_capacity(inlet_velocities.len());

    for &vin in &inlet_velocities {
        let mut v = base.clone();
        v.inlet_velocity = vin;
        let p_throat = v.throat_pressure();
        let sigma = v.cavitation_number();
        let cav_len = v.cavity_length(sigma);

        labels.push(format!("{:.1}", vin));
        throat_pressures_kpa.push(p_throat / 1000.0);
        sigmas.push(sigma);
        cavity_lengths_mm.push(cav_len * 1000.0);
    }

    (labels, throat_pressures_kpa, sigmas, cavity_lengths_mm)
}

fn make_vof_config() -> CavitationVofConfig {
    CavitationVofConfig {
        vof_config: VofConfig {
            surface_tension_coefficient: 0.072,
            interface_compression: 0.1,
            reconstruction_method: InterfaceReconstruction::PLIC,
            advection_method: AdvectionMethod::Geometric,
            max_iterations: 10,
            tolerance: 1e-6,
            cfl_number: 0.3,
            enable_compression: false,
        },
        cavitation_model: CavitationModel::Kunz {
            vaporization_coeff: 100.0,
            condensation_coeff: 100.0,
        },
        damage_model: None,
        bubble_dynamics: None,
        inception_threshold: 0.5,
        max_void_fraction: 0.8,
        relaxation_time: 1e-6,
        vapor_pressure: 2330.0,
        liquid_density: 998.0,
        liquid_viscosity: 0.001, // Water at ~20°C [Pa·s]
        vapor_density: 0.023,
        sound_speed: 1500.0,
    }
}

struct VenturiLikeFieldParams {
    dx: f64,
    inlet_velocity: f64,
    throat_velocity: f64,
    inlet_pressure: f64,
    density: f64,
}

fn build_venturi_like_fields(
    nx: usize,
    ny: usize,
    nz: usize,
    params: VenturiLikeFieldParams,
) -> (Vec<Vector3<f64>>, DMatrix<f64>, DMatrix<f64>) {
    let grid_size = nx * ny * nz;
    let mut velocity = vec![Vector3::<f64>::zeros(); grid_size];
    let mut pressure = DMatrix::<f64>::zeros(nx, ny * nz);
    let density_field = DMatrix::<f64>::from_element(nx, ny * nz, params.density);

    let x0 = 0.5 * (nx as f64 - 1.0) * params.dx;
    let w = 0.12 * (nx as f64 - 1.0) * params.dx;
    let bump = (params.throat_velocity - params.inlet_velocity).max(0.0);

    for i in 0..nx {
        let x = i as f64 * params.dx;
        let g = if w > 0.0 {
            (-(x - x0).powi(2) / (2.0 * w * w)).exp()
        } else {
            0.0
        };
        let ux = params.inlet_velocity + bump * g;
        let p = params.inlet_pressure
            - 0.5 * params.density * (ux * ux - params.inlet_velocity * params.inlet_velocity);

        for j in 0..ny {
            for k in 0..nz {
                let idx = (k * ny + j) * nx + i;
                let col = j + k * ny;
                velocity[idx] = Vector3::new(ux, 0.0, 0.0);
                pressure[(i, col)] = p;
            }
        }
    }

    (velocity, pressure, density_field)
}

fn run_vof_scenario(
    nx: usize,
    ny: usize,
    nz: usize,
    steps: usize,
    dt: f64,
) -> (Vec<String>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let config = make_vof_config();
    let mut solver = CavitationVofSolver::new(nx, ny, nz, config).expect("solver init");

    let field_params = VenturiLikeFieldParams {
        dx: 0.01,
        inlet_velocity: 2.0,
        throat_velocity: 15.0,
        inlet_pressure: 101_325.0,
        density: 998.0,
    };
    let (velocity, pressure, density) = build_venturi_like_fields(nx, ny, nz, field_params);

    let mut labels = Vec::with_capacity(steps + 1);
    let mut cav_frac = Vec::with_capacity(steps + 1);
    let mut total_void = Vec::with_capacity(steps + 1);
    let mut max_void = Vec::with_capacity(steps + 1);

    for step in 0..=steps {
        let stats = solver.cavitation_statistics();
        labels.push(format!("{}", step));
        cav_frac.push(stats.cavitation_fraction);
        total_void.push(stats.total_void_fraction);
        max_void.push(stats.max_void_fraction);

        if step < steps {
            solver
                .step(dt, &velocity, &pressure, &density)
                .expect("step failed");
        }
    }

    (labels, cav_frac, total_void, max_void)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = std::path::Path::new("target").join("plots");

    let (labels, throat_p_kpa, sigma, cavity_mm) = venturi_1d_sweep();
    let charts_1d = vec![
        ChartSpec {
            canvas_id: "p_throat".to_string(),
            heading: "Throat Pressure vs Inlet Velocity".to_string(),
            y_label: "Pressure (kPa)".to_string(),
            labels: labels.clone(),
            datasets: vec![DatasetSpec {
                label: "Throat Pressure".to_string(),
                data: throat_p_kpa,
                color: "#1f77b4".to_string(),
            }],
        },
        ChartSpec {
            canvas_id: "sigma".to_string(),
            heading: "Cavitation Number vs Inlet Velocity".to_string(),
            y_label: "Cavitation Number σ".to_string(),
            labels: labels.clone(),
            datasets: vec![DatasetSpec {
                label: "σ".to_string(),
                data: sigma,
                color: "#ff7f0e".to_string(),
            }],
        },
        ChartSpec {
            canvas_id: "cavity".to_string(),
            heading: "Nurick Cavity Length vs Inlet Velocity".to_string(),
            y_label: "Cavity Length (mm)".to_string(),
            labels: labels.clone(),
            datasets: vec![DatasetSpec {
                label: "Cavity Length".to_string(),
                data: cavity_mm,
                color: "#d62728".to_string(),
            }],
        },
    ];
    let html_1d = chart_js_html("1D Venturi Cavitation Plots", &charts_1d);
    write_html(&out_dir.join("venturi_1d.html"), &html_1d)?;

    let (labels_2d, cav_2d, total_void_2d, max_void_2d) = run_vof_scenario(60, 40, 1, 80, 1e-5);
    let charts_2d = vec![
        ChartSpec {
            canvas_id: "cav_frac_2d".to_string(),
            heading: "Cavitation Fraction vs Step".to_string(),
            y_label: "Cavitating Fraction".to_string(),
            labels: labels_2d.clone(),
            datasets: vec![DatasetSpec {
                label: "Cavitation Fraction".to_string(),
                data: cav_2d,
                color: "#1f77b4".to_string(),
            }],
        },
        ChartSpec {
            canvas_id: "void_2d".to_string(),
            heading: "Void Fraction Metrics vs Step".to_string(),
            y_label: "Void Fraction (sum/max)".to_string(),
            labels: labels_2d.clone(),
            datasets: vec![
                DatasetSpec {
                    label: "Total Void Fraction".to_string(),
                    data: total_void_2d,
                    color: "#ff7f0e".to_string(),
                },
                DatasetSpec {
                    label: "Max Void Fraction".to_string(),
                    data: max_void_2d,
                    color: "#2ca02c".to_string(),
                },
            ],
        },
    ];
    let html_2d = chart_js_html("2D Cavitation-VOF Plots (nz=1)", &charts_2d);
    write_html(&out_dir.join("cavitation_2d.html"), &html_2d)?;

    let (labels_3d, cav_3d, total_void_3d, max_void_3d) = run_vof_scenario(30, 20, 10, 80, 1e-5);
    let charts_3d = vec![
        ChartSpec {
            canvas_id: "cav_frac_3d".to_string(),
            heading: "Cavitation Fraction vs Step".to_string(),
            y_label: "Cavitating Fraction".to_string(),
            labels: labels_3d.clone(),
            datasets: vec![DatasetSpec {
                label: "Cavitation Fraction".to_string(),
                data: cav_3d,
                color: "#1f77b4".to_string(),
            }],
        },
        ChartSpec {
            canvas_id: "void_3d".to_string(),
            heading: "Void Fraction Metrics vs Step".to_string(),
            y_label: "Void Fraction (sum/max)".to_string(),
            labels: labels_3d.clone(),
            datasets: vec![
                DatasetSpec {
                    label: "Total Void Fraction".to_string(),
                    data: total_void_3d,
                    color: "#ff7f0e".to_string(),
                },
                DatasetSpec {
                    label: "Max Void Fraction".to_string(),
                    data: max_void_3d,
                    color: "#2ca02c".to_string(),
                },
            ],
        },
    ];
    let html_3d = chart_js_html("3D Cavitation-VOF Plots", &charts_3d);
    write_html(&out_dir.join("cavitation_3d.html"), &html_3d)?;

    println!("Wrote plots to: {}", out_dir.display());
    println!(" - {}", out_dir.join("venturi_1d.html").display());
    println!(" - {}", out_dir.join("cavitation_2d.html").display());
    println!(" - {}", out_dir.join("cavitation_3d.html").display());

    Ok(())
}
