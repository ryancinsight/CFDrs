//! Medical-Grade Millifluidic CFD Screening — 2D Bifurcation Solver
//!
//! 2D counterpart of `cfd-1d::medical_millifluidic_screening`. Replaces the
//! lumped 1D network with a full 2D Navier-Stokes bifurcation solver, giving
//! spatially resolved velocity and pressure fields for a Murray's-law-sized
//! blood bifurcation within a 96-well plate footprint.
//!
//! Comprehensive analysis:
//! 1. Murray's law bifurcation schematic (cfd-schematics)
//! 2. 2D CFD with Casson non-Newtonian blood (BifurcationSolver2D)
//! 3. Hemolysis index (Giersiepen–Wurzinger)
//! 4. Wall shear stress estimation
//! 5. FDA shear limit screening
//! 6. Cavitation risk assessment
//! 7. Flow rate distribution
//! 8. Pressure distribution
//! 9. Render 6 overlays + 1 plain schematic + JSON export
//!
//! Run with:
//! `cargo run -p cfd-2d --example medical_screening_2d`

use cfd_2d::solvers::bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
use cfd_2d::solvers::ns_fvm_2d::{BloodModel, SIMPLEConfig};
use cfd_core::physics::cavitation::CavitationNumber;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::hemolysis::HemolysisModel;
use cfd_schematics::geometry::types::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

/// FDA conservative wall shear stress limit for whole-blood devices [Pa].
const FDA_WSS_LIMIT: f64 = 150.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🏥 Medical-Grade 2D CFD Screening — Bifurcation");
    println!("================================================\n");

    // ── 0. Output directory ──────────────────────────────────────────────────
    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("medical_screening");
    fs::create_dir_all(&out)?;

    // ── Phase 1: Design (cfd-schematics) ─────────────────────────────────────
    println!("📐 Phase 1: Schematic Design (Murray's Law Bifurcation)\n");

    // Murray's law: D_parent³ = D_d1³ + D_d2³  →  r_d = r_p / 2^(1/3)
    let r_parent = 1.0_f64;                          // mm
    let r_daughter = r_parent / 2.0_f64.powf(1.0 / 3.0); // ≈ 0.794 mm
    let w_parent = r_parent * 2.0;
    let w_daughter = r_daughter * 2.0;

    let angle_deg = 30.0_f64;
    let angle_rad = angle_deg.to_radians();
    let l_parent = 6.0_f64;   // mm — parent channel length
    let l_daughter = 5.0_f64;  // mm — daughter channel length

    let jx = l_parent;
    let jy = 5.0_f64;

    let d1_end_x = jx + l_daughter * angle_rad.cos();
    let d1_end_y = jy + l_daughter * angle_rad.sin();
    let d2_end_x = jx + l_daughter * angle_rad.cos();
    let d2_end_y = jy - l_daughter * angle_rad.sin();

    let nodes: Vec<Node> = vec![
        Node { id: 0, point: (0.0, jy),             metadata: None },
        Node { id: 1, point: (jx, jy),              metadata: None },
        Node { id: 2, point: (d1_end_x, d1_end_y),  metadata: None },
        Node { id: 3, point: (d2_end_x, d2_end_y),  metadata: None },
    ];

    let channels: Vec<Channel> = vec![
        Channel { id: 0, from_node: 0, to_node: 1, width: w_parent,   height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 1, from_node: 1, to_node: 2, width: w_daughter, height: 1.0, channel_type: ChannelType::Straight, metadata: None },
        Channel { id: 2, from_node: 1, to_node: 3, width: w_daughter, height: 1.0, channel_type: ChannelType::Straight, metadata: None },
    ];

    let system = ChannelSystem {
        box_dims: (d1_end_x + 2.0, jy * 2.0),
        nodes,
        channels,
        box_outline: vec![],
    };

    let renderer = create_plotters_renderer();
    let schematic_cfg = RenderConfig {
        width: 1000,
        height: 700,
        title: format!(
            "96-Well Bifurcation — Murray's Law (r_p={:.1}mm, r_d={:.2}mm, θ={:.0}°)",
            r_parent, r_daughter, angle_deg
        ),
        ..Default::default()
    };

    let schematic_path = out.join("schematic.png");
    renderer.render_system(&system, schematic_path.to_str().unwrap(), &schematic_cfg)?;
    println!("  ✅ Schematic → {}", schematic_path.display());
    println!("  Murray's law: r_p³ = {:.4}, r_d1³+r_d2³ = {:.4}",
        r_parent.powi(3), r_daughter.powi(3) * 2.0);

    // ── Phase 2: 2D CFD Simulation ───────────────────────────────────────────
    println!("\n⚙️  Phase 2: 2D Navier-Stokes Simulation (Casson Blood)\n");

    let geom = BifurcationGeometry::new_symmetric(
        w_parent * 1e-3,    // parent width [m]
        l_parent * 1e-3,    // parent length [m]
        w_daughter * 1e-3,  // daughter width [m]
        l_daughter * 1e-3,  // daughter length [m]
        angle_rad,
    );

    let casson = CassonBlood::<f64>::normal_blood();
    let density = 1060.0_f64;
    let blood = BloodModel::Casson(casson);
    let u_inlet = 0.1_f64; // 100 mm/s

    let mut config = SIMPLEConfig::default();
    config.max_iterations = 3000;
    config.tolerance = 1e-4;
    config.alpha_u = 0.5;
    config.alpha_p = 0.2;

    println!("  Parent:   w={:.2}mm, L={:.1}mm", w_parent, l_parent);
    println!("  Daughter: w={:.3}mm, L={:.1}mm (Murray's law)", w_daughter, l_daughter);
    println!("  Inlet:    u={:.0} mm/s", u_inlet * 1e3);
    println!("  Fluid:    Casson blood, ρ={:.0} kg/m³", density);
    println!("  Device:   96-well plate footprint (ANSI/SLAS 1-2004)");

    let mut solver = BifurcationSolver2D::new(geom, blood, density, 80, 50, config);
    let sol = solver.solve(u_inlet)?;

    println!("\n  📊 2D Simulation Results:");
    println!("    Q_parent    = {:.4e} m²/s", sol.q_parent);
    println!("    Q_daughter1 = {:.4e} m²/s", sol.q_daughter1);
    println!("    Q_daughter2 = {:.4e} m²/s", sol.q_daughter2);
    println!("    Mass balance error = {:.2e}", sol.mass_balance_error);

    let symmetry_err = (sol.q_daughter1 - sol.q_daughter2).abs() / sol.q_parent.abs().max(1e-20);
    println!("    Symmetry error = {:.2e}", symmetry_err);

    // ── Phase 3: Medical Analysis ────────────────────────────────────────────
    println!("\n🔬 Phase 3: Medical Analysis\n");

    let hemolysis_model = HemolysisModel::giersiepen_standard();
    let vapor_pressure = 6300.0; // Pa — blood at 37°C

    // Channel analysis data: (channel_id, width_m, length_m, flow_rate)
    let channel_info: Vec<(usize, f64, f64, f64)> = vec![
        (0, w_parent * 1e-3,   l_parent * 1e-3,   sol.q_parent),
        (1, w_daughter * 1e-3, l_daughter * 1e-3, sol.q_daughter1),
        (2, w_daughter * 1e-3, l_daughter * 1e-3, sol.q_daughter2),
    ];

    let height_m: f64 = 1.0e-3; // 1 mm channel height

    let mut edge_hemolysis = HashMap::<usize, f64>::new();
    let mut edge_shear = HashMap::<usize, f64>::new();
    let mut edge_flow = HashMap::<usize, f64>::new();
    let mut edge_velocity = HashMap::<usize, f64>::new();
    let mut edge_pressure = HashMap::<usize, f64>::new();
    let mut edge_cavitation = HashMap::<usize, f64>::new();
    let mut node_pressure = HashMap::<usize, f64>::new();
    let mut fda_violations: Vec<(usize, f64)> = Vec::new();

    let mut max_hi = 0.0_f64;
    let mut max_wss = 0.0_f64;

    // Estimate total pressure drop from 2D flow rate using Hele-Shaw relation
    // ΔP = (12 · µ_eff · L · Q) / (w · h³)
    let mu_eff = 3.5e-3; // approximate Casson viscosity at moderate shear [Pa·s]

    // Node pressure estimates (inlet → junction → outlets)
    let dp_parent = 12.0 * mu_eff * (l_parent * 1e-3) * sol.q_parent.abs()
        / ((w_parent * 1e-3) * height_m.powi(3));
    let dp_daughter = 12.0 * mu_eff * (l_daughter * 1e-3) * sol.q_daughter1.abs()
        / ((w_daughter * 1e-3) * height_m.powi(3));
    let p_inlet = dp_parent + dp_daughter;
    let p_junction = dp_daughter;
    let p_outlet = 0.0;

    node_pressure.insert(0, p_inlet);
    node_pressure.insert(1, p_junction);
    node_pressure.insert(2, p_outlet);
    node_pressure.insert(3, p_outlet);

    for &(ch_id, w, l, q) in &channel_info {
        let area = w * height_m;
        let velocity = q.abs() / area;
        let d_h = 2.0 * w * height_m / (w + height_m);
        let shear_rate = 8.0 * velocity / d_h;

        // Wall shear stress (Casson approximation at computed shear rate)
        let wall_shear = mu_eff * shear_rate;

        // Exposure time
        let exposure_time = if velocity > 1e-12 { l / velocity } else { 0.0 };

        // Hemolysis index
        let hi = hemolysis_model.damage_index(wall_shear, exposure_time).unwrap_or(0.0);

        // Cavitation
        let local_p = if ch_id == 0 {
            (p_inlet + p_junction) / 2.0
        } else {
            (p_junction + p_outlet) / 2.0
        };
        let cav = CavitationNumber {
            reference_pressure: local_p,
            vapor_pressure,
            density,
            velocity,
        };
        let sigma = cav.calculate();

        // FDA screening
        if wall_shear > FDA_WSS_LIMIT {
            fda_violations.push((ch_id, wall_shear));
        }

        edge_hemolysis.insert(ch_id, hi);
        edge_shear.insert(ch_id, wall_shear);
        edge_flow.insert(ch_id, q.abs());
        edge_velocity.insert(ch_id, velocity);
        edge_pressure.insert(ch_id, local_p);
        edge_cavitation.insert(ch_id, 1.0 / sigma.max(1e-12));

        max_hi = max_hi.max(hi);
        max_wss = max_wss.max(wall_shear);
    }

    println!("  Max hemolysis index:    {:.4e}", max_hi);
    println!("  Max wall shear stress:  {:.2} Pa", max_wss);
    println!("  FDA violations:         {} of {} channels", fda_violations.len(), channel_info.len());
    println!("  FDA limit (τ_w):        {:.0} Pa", FDA_WSS_LIMIT);

    if !fda_violations.is_empty() {
        println!("\n  ⚠ FDA Shear Limit Violations:");
        for &(id, wss) in &fda_violations {
            println!("    Channel {id}: τ_w = {:.2} Pa ({:.1}× limit)", wss, wss / FDA_WSS_LIMIT);
        }
    }

    // ── Phase 4: Overlay Visualization ───────────────────────────────────────
    println!("\n🎨 Phase 4: Overlay Visualization\n");

    let overlays: Vec<(&str, AnalysisOverlay, String)> = vec![
        (
            "hemolysis_index.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Hemolysis Index (Giersiepen)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_hemolysis.clone())
            .with_node_data(node_pressure.clone()),
            format!("Hemolysis Index — Giersiepen (max HI = {:.2e})", max_hi),
        ),
        (
            "wall_shear_stress.png",
            AnalysisOverlay::new(AnalysisField::WallShearStress, ColormapKind::BlueRed)
                .with_edge_data(edge_shear.clone())
                .with_node_data(node_pressure.clone()),
            format!("Wall Shear Stress (max τ_w = {:.2} Pa, FDA limit: {:.0} Pa)", max_wss, FDA_WSS_LIMIT),
        ),
        (
            "fda_screening.png",
            {
                let fda_ratio: HashMap<usize, f64> = edge_shear
                    .iter()
                    .map(|(&k, &v)| (k, (v / FDA_WSS_LIMIT).min(3.0)))
                    .collect();
                AnalysisOverlay::new(
                    AnalysisField::Custom("FDA Shear Exceedance Ratio".into()),
                    ColormapKind::BlueRed,
                )
                .with_edge_data(fda_ratio)
                .with_node_data(node_pressure.clone())
            },
            format!("FDA Screening (limit: {:.0} Pa, {} violations)", FDA_WSS_LIMIT, fda_violations.len()),
        ),
        (
            "cavitation_risk.png",
            AnalysisOverlay::new(
                AnalysisField::Custom("Cavitation Risk (1/σ)".into()),
                ColormapKind::BlueRed,
            )
            .with_edge_data(edge_cavitation.clone())
            .with_node_data(node_pressure.clone()),
            "Cavitation Risk (1/σ — red = high risk)".to_string(),
        ),
        (
            "flow_rate.png",
            AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
                .with_edge_data(edge_flow.clone())
                .with_node_data(node_pressure.clone()),
            "Flow Rate Distribution [m²/s]".to_string(),
        ),
        (
            "pressure.png",
            AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::Viridis)
                .with_edge_data(edge_pressure.clone())
                .with_node_data(node_pressure.clone()),
            "Pressure Distribution [Pa]".to_string(),
        ),
    ];

    for (filename, overlay, title) in &overlays {
        let path = out.join(filename);
        let cfg = RenderConfig {
            width: 1000,
            height: 700,
            title: title.clone(),
            show_axes: true,
            show_grid: false,
            ..Default::default()
        };
        renderer.render_analysis(&system, path.to_str().unwrap(), &cfg, overlay)?;
        println!("  ✓ {filename}");
    }

    // ── Phase 5: JSON Export ─────────────────────────────────────────────────
    println!("\n📄 Phase 5: JSON Export\n");

    let channel_report: Vec<serde_json::Value> = channel_info
        .iter()
        .map(|&(ch_id, _w, _l, _q)| {
            serde_json::json!({
                "channel_id": ch_id,
                "hemolysis_index": edge_hemolysis.get(&ch_id),
                "wall_shear_stress_pa": edge_shear.get(&ch_id),
                "flow_rate_m2_s": edge_flow.get(&ch_id),
                "velocity_m_s": edge_velocity.get(&ch_id),
                "local_pressure_pa": edge_pressure.get(&ch_id),
                "cavitation_risk_inv_sigma": edge_cavitation.get(&ch_id),
                "fda_violation": fda_violations.iter().any(|&(id, _)| id == ch_id)
            })
        })
        .collect();

    let results = serde_json::json!({
        "analysis": "medical_screening_2d",
        "solver": "BifurcationSolver2D (Navier-Stokes FVM, Casson)",
        "device": {
            "application": "Sonodynamic Therapy (SDT)",
            "standard": "ANSI/SLAS 1-2004 (96-well)",
            "topology": "Murray's law symmetric bifurcation"
        },
        "geometry": {
            "parent_width_mm": w_parent,
            "parent_length_mm": l_parent,
            "daughter_width_mm": w_daughter,
            "daughter_length_mm": l_daughter,
            "branch_angle_deg": angle_deg,
            "murrays_law_verification": {
                "r_parent_cubed": r_parent.powi(3),
                "sum_r_daughter_cubed": r_daughter.powi(3) * 2.0
            }
        },
        "fluid": {
            "model": "Casson blood",
            "density_kg_m3": density,
            "effective_viscosity_pa_s": mu_eff,
            "vapor_pressure_pa": vapor_pressure
        },
        "boundary_conditions": {
            "inlet_velocity_m_s": u_inlet
        },
        "simulation_results": {
            "q_parent_m2_s": sol.q_parent,
            "q_daughter1_m2_s": sol.q_daughter1,
            "q_daughter2_m2_s": sol.q_daughter2,
            "mass_balance_error": sol.mass_balance_error,
            "symmetry_error": symmetry_err
        },
        "medical_analysis": {
            "max_hemolysis_index": max_hi,
            "max_wall_shear_stress_pa": max_wss,
            "fda_wss_limit_pa": FDA_WSS_LIMIT,
            "fda_violations": fda_violations.len(),
            "hemolysis_model": "Giersiepen-Wurzinger"
        },
        "channel_data": channel_report
    });

    let json_path = out.join("results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("  ✅ JSON → {}", json_path.display());

    println!("\n✅ Medical screening 2D example complete.");
    Ok(())
}
