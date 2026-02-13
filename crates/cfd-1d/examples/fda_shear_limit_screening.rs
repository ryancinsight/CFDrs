//! End-to-end FDA-oriented blood shear limit screening example.
//!
//! This example:
//! 1) Builds and solves a simple 1D channel network,
//! 2) Runs flow analysis,
//! 3) Flags components that exceed configured blood shear limits.
//!
//! Run with:
//! `cargo run -p cfd-1d --example fda_shear_limit_screening`

use cfd_1d::{
    BloodShearLimits, EdgeProperties, Network, NetworkAnalysisResult, NetworkAnalyzerOrchestrator,
    NetworkBuilder,
};
use cfd_1d::channel::{ChannelGeometry, ChannelType, CrossSection, SurfaceProperties, Wettability};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Approximate whole-blood Newtonianized properties for this screening demo.
    // (Use your validated blood model/rheology calibration for production studies.)
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
    let diameter = 120e-6;
    let length = 2.0e-3;
    let area = std::f64::consts::PI * (diameter * 0.5_f64).powi(2);
    let geometry = ChannelGeometry {
        channel_type: ChannelType::Straight,
        length,
        cross_section: CrossSection::Circular { diameter },
        surface: SurfaceProperties {
            roughness: 0.0,
            contact_angle: None,
            surface_energy: None,
            wettability: Wettability::Hydrophilic,
        },
        variations: Vec::new(),
    };

    network.add_edge_properties(
        edge,
        EdgeProperties {
            id: "treatment_channel".to_string(),
            component_type: cfd_1d::network::ComponentType::Pipe,
            length,
            area,
            hydraulic_diameter: Some(diameter),
            resistance: 0.0,
            geometry: Some(geometry),
            properties: std::collections::HashMap::new(),
        },
    );

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
