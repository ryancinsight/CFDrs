//! Cross-Fidelity Branching Validation Example
//!
//! Uses `cfd-schematics` to compare solutions between `cfd-1d` (Hagen-Poiseuille)
//! and `cfd-2d` (Depth-averaged LBM) for standard bifurcation and trifurcation
//! millifluidic designs.
//!
//! # Architecture
//! - Generates `NetworkBlueprint` using presets (`bifurcation_rect`, `trifurcation_rect`)
//! - Calls `process_blueprint_with_reference_trace` with `run_2d_reference: true`
//! - Outputs a comparative table of mass conservation and pressure drops
//! - Renders an analysis overlay of flow rates for visual verification
//!
//! # Execution
//! `cargo run --example cross_fidelity_branching --release`

use cfd_3d::blueprint_integration::{
    process_blueprint_with_reference_trace, Blueprint3dProcessingConfig,
};
use cfd_schematics::interface::presets::{
    bifurcation_rect, pentafurcation_rect, quadfurcation_rect, trifurcation_rect,
};
use cfd_schematics::visualizations::{
    create_plotters_renderer, AnalysisField, AnalysisOverlay, ColormapKind, RenderConfig,
    SchematicRenderer,
};
use cfd_1d::physics::cell_separation::{CellProperties, CellSeparationModel};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("🔬 Cross-Fidelity Branching Validation");
    println!("=====================================");

    let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("cross_fidelity_branching");
    fs::create_dir_all(&out_dir)?;

    // Shared Fluid & Flow Parameters
    // Blood-like properties for physiological relevance
    let density_kg_m3 = 1060.0;
    let viscosity_pa_s = 3.5e-3;
    let total_flow_rate_m3_s = 1.0e-7; // 100 µL/s or 0.1 mL/s

    // We will evaluate two cases:
    // 1. Symmetric Bifurcation
    // 2. Symmetric Trifurcation
    let bp_bif = bifurcation_rect(
        "Symmetric Bifurcation",
        60e-3,  // parent length
        10e-3,  // daughter length
        4.0e-3, // parent width
        4.0e-3, // daughter width
        4.0e-3, // height
    );

    let bp_tri = trifurcation_rect(
        "Symmetric Trifurcation",
        60e-3,  // parent length
        10e-3,  // daughter length
        4.0e-3, // parent width
        4.0e-3, // daughter width
        4.0e-3, // height
    );

    let bp_quad = quadfurcation_rect(
        "Symmetric Quadfurcation",
        60e-3, 10e-3, 4.0e-3, 4.0e-3, 4.0e-3,
    );

    let bp_penta = pentafurcation_rect(
        "Symmetric Pentafurcation",
        60e-3, 10e-3, 4.0e-3, 4.0e-3, 4.0e-3,
    );

    let cell_cancer = CellProperties::mcf7_breast_cancer();
    let cell_rbc = CellProperties::red_blood_cell();

    for (bp, filename_prefix) in [
        (bp_bif, "bifurcation_rect"),
        (bp_tri, "trifurcation_rect"),
        (bp_quad, "quadfurcation_rect"),
    ] {
        println!("\n▶ Evaluating: {}", bp.name);

        // ── 1. Process Blueprint (1D + 2D) ──────────────────────────────────
        let config = Blueprint3dProcessingConfig {
            density_kg_m3,
            viscosity_pa_s,
            total_flow_rate_m3_s,
            // 2D grid resolution — kept moderate for example execution time
            two_d_grid_nx: 64,
            two_d_grid_ny: 16,
            two_d_tolerance: 1e-9,
            run_2d_reference: true, // Crucial: explicitly solve in 2D
            ..Default::default()
        };

        let result = process_blueprint_with_reference_trace(&bp, &config)
            .expect("Integration pipeline failed");

        // ── 2. Comparative Analysis Table ───────────────────────────────────
        println!("  • 1D vs 2D Results Table\n");
        println!(
            "  {:<15} | {:<16} | {:<16} | {:<16} | {:<10} | {:<14} | {:<14} | {:<10}",
            "Channel", "1D Flow [m³/s]", "2D Conv Flow", "2D Err [%]", "1D dP [Pa]", "1D Sep Eff [%]", "2D Sep Eff [%]", "2D Shear [Pa]"
        );
        println!("  {:-<15}-+-{:-<16}-+-{:-<16}-+-{:-<16}-+-{:-<10}-+-{:-<14}-+-{:-<14}-+-{:-<10}", "", "", "", "", "", "", "", "");

        let mut edge_flow_map = HashMap::new();
        let mut node_flow_map = HashMap::new();

        for trace in &result.channel_traces {
            let flow_1d = trace.reference_flow_rate_m3_s;
            let dp_1d = trace.reference_pressure_drop_pa;

            // Optional 2D results
            let flow_2d_str = trace
                .two_d_outlet_flow_error_pct
                .map(|e| {
                    let flow_2d = flow_1d * (1.0 + e / 100.0);
                    format!("{:.2e}", flow_2d)
                })
                .unwrap_or_else(|| "N/A".to_string());

            let err_2d_str = trace
                .two_d_outlet_flow_error_pct
                .map(|e| format!("{:>6.2}%", e))
                .unwrap_or_else(|| "N/A".to_string());

            let shear_2d_str = trace
                .two_d_field_wall_shear_mean_pa
                .map(|s| format!("{:>6.2}", s))
                .unwrap_or_else(|| "N/A".to_string());

            let mut sep_eff_str = "N/A".to_string();
            // Look up channel geometry to evaluate cell separation efficiency
            if let Some(geom) = bp.channels.iter().find(|ch| ch.id.as_str() == trace.channel_id.as_str()) {
                if let cfd_schematics::domain::model::CrossSectionSpec::Rectangular { width_m, height_m } = geom.cross_section {
                    if geom.length_m > 0.0 && width_m > 0.0 && height_m > 0.0 {
                        let sep_model = CellSeparationModel::new(width_m, height_m, geom.length_m, None);
                        if let Some(analysis) = sep_model.analyze(
                            &cell_cancer,
                            &cell_rbc,
                            density_kg_m3,
                            viscosity_pa_s,
                            trace.reference_mean_velocity_m_s,
                        ) {
                            sep_eff_str = format!("{:>6.2}%", analysis.separation_efficiency * 100.0);
                        }
                    }
                }
            }

            let sep_2d_str = trace
                .two_d_field_separation_efficiency_pct
                .map(|s| format!("{:>6.2}%", s))
                .unwrap_or_else(|| "N/A".to_string());

            println!(
                "  {:<15} | {:<16.2e} | {:<16} | {:<16} | {:<10.1} | {:<14} | {:<14} | {:<10}",
                trace.channel_id, flow_1d, flow_2d_str, err_2d_str, dp_1d, sep_eff_str, sep_2d_str, shear_2d_str
            );

            // Mapping for Visualization (assuming channel IDs end in an index for basic visualization,
            // or we use the topology order. In process_blueprint, channel ordering matches blueprint.channels).
        }

        // Map flow rates backward to edge indices for `cfd-schematics` renderer
        for (i, ch) in bp.channels.iter().enumerate() {
            if let Some(trace) = result
                .channel_traces
                .iter()
                .find(|t| t.channel_id == ch.id.as_str())
            {
                edge_flow_map.insert(i, trace.reference_flow_rate_m3_s);
                
                // Hacky node mapping: from/to node IDs to indices
                let mut from_idx = 0;
                let mut to_idx = 0;
                for (j, n) in bp.nodes.iter().enumerate() {
                    if n.id == ch.from {
                        from_idx = j;
                    }
                    if n.id == ch.to {
                        to_idx = j;
                    }
                }
                node_flow_map.insert(from_idx, trace.reference_flow_rate_m3_s);
                node_flow_map.insert(to_idx, trace.reference_flow_rate_m3_s);
            }
        }

        // ── 3. Overlay Visualization ────────────────────────────────────────
        let renderer = create_plotters_renderer();
        let overlay = AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis)
            .with_edge_data(edge_flow_map)
            .with_node_data(node_flow_map);

        let overlay_cfg = RenderConfig {
            width: 1000,
            height: 600,
            title: format!("{} — 1D Flow Rate Distribution [m³/s]", bp.name),
            ..Default::default()
        };

        let overlay_path = out_dir.join(format!("{}_flow_overlay.svg", filename_prefix));
        renderer.render_analysis(
            &bp,
            overlay_path.to_str().unwrap(),
            &overlay_cfg,
            &overlay,
        )?;
        
        let overlay_path_png = out_dir.join(format!("{}_flow_overlay.png", filename_prefix));
        renderer.render_analysis(
            &bp,
            overlay_path_png.to_str().unwrap(),
            &overlay_cfg,
            &overlay,
        )?;

        println!("  ✅ Rendered overlay → {}", overlay_path.display());
    }

    println!("\n✅ Cross-fidelity branching validation complete.");
    Ok(())
}
