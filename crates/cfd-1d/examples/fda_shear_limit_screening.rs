//! End-to-end FDA-oriented blood shear limit screening example.
//!
//! This example:
//! 1) Builds and solves a simple 1D channel network defined via `ChannelSpec`,
//! 2) Runs flow analysis,
//! 3) Flags components that exceed configured blood shear limits.
//!
//! Network topology and physical geometry are defined entirely via
//! `cfd-schematics` `ChannelSpec` — no `cfd_1d::channel` imports required.
//!
//! Run with:
//! `cargo run -p cfd-1d --example fda_shear_limit_screening`

use cfd_1d::{
    BloodShearLimits, EdgeProperties, Network, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
    NetworkBuilder,
};
use cfd_schematics::domain::model::ChannelSpec;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Approximate whole-blood Newtonianized properties for this screening demo.
    let blood_like = cfd_core::physics::fluid::ConstantPropertyFluid::<f64>::new(
        "Blood-like (screening)".to_string(),
        1060.0,
        0.0035,
        3770.0,
        0.52,
        1570.0,
    );

    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let edge = builder.connect_with_pipe(inlet, outlet, "treatment_channel".to_string());
    let graph = builder.build()?;

    let mut network = Network::new(graph, blood_like);

    // Intentionally narrow channel to demonstrate a potential shear-limit flag.
    // Hagen-Poiseuille: R = 128·μ·L / (π·D⁴)
    let diameter = 120e-6_f64;
    let length = 2.0e-3_f64;
    let mu = 0.0035_f64;
    let resistance = 128.0 * mu * length / (std::f64::consts::PI * diameter.powi(4));

    let channel_spec = ChannelSpec::new_pipe(
        "treatment_channel",
        "inlet",
        "outlet",
        length,
        diameter,
        resistance,
        0.0,
    );

    // EdgeProperties::from(&ChannelSpec) populates ChannelGeometry from CrossSectionSpec
    network.add_edge_properties(edge, EdgeProperties::from(&channel_spec));
    network.update_resistances()?;

    // Pressure-driven setup.
    network.set_pressure(inlet, 12_000.0);
    network.set_pressure(outlet, 0.0);

    let mut analyzer = NetworkAnalyzerOrchestrator::<f64>::new();
    let result: NetworkAnalysisResult<f64> = analyzer.analyze(&mut network)?;

    let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood();
    let violations = result
        .flow_analysis
        .flag_fda_shear_limit_violations(&limits);

    println!("FDA-oriented shear screening summary");
    println!(
        "Configured max wall shear stress: {:.1} Pa",
        limits.max_wall_shear_stress_pa
    );
    println!("Components analyzed: {}", result.flow_analysis.component_flows.len());

    if violations.is_empty() {
        println!("No shear-limit violations flagged.");
    } else {
        println!("Violations flagged: {}", violations.len());
        for v in violations {
            println!(
                "- {}: tau_w={:.2} Pa (limit {:.2} Pa, x{:.2})",
                v.component_id,
                v.wall_shear_stress_pa,
                v.stress_limit_pa,
                v.stress_exceedance_ratio
            );
        }
    }

    Ok(())
}
