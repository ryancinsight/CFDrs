//! 1D Blood Flow through a TPMS-Inspired Channel Network
//!
//! Demonstrates simulating blood flow through a millifluidic network whose
//! topology derives from a TPMS (Triply Periodic Minimal Surface) Gyroid
//! lattice embedded in a 96-well plate footprint.
//!
//! The TPMS lattice partitions the internal cavity into interconnected
//! channels. This example models a representative sub-network as a
//! 1D lumped-parameter system using the Hagen-Poiseuille resistance model
//! with Carreau-Yasuda non-Newtonian blood rheology.
//!
//! Network topology (Gyroid-inspired multi-path bifurcation):
//!
//! ```text
//!                 ┌─ ch_upper_1 ─ J2 ─ ch_upper_2 ─┐
//!   Inlet ─ ch_in ─ J1                               J3 ─ ch_out ─ Outlet
//!                 └─ ch_lower_1 ─ J4 ─ ch_lower_2 ─┘
//! ```
//!
//! This mirrors how a Gyroid lattice naturally bifurcates flow into parallel
//! paths through its bicontinuous pore network, then recombines at the outlet.
//!
//! # Run
//!
//! ```sh
//! cargo run -p cfd-1d --example tpms_blood_1d
//! ```

use cfd_1d::domain::network::{EdgeProperties, Network, NetworkBuilder};
use cfd_1d::solver::core::{NetworkProblem, NetworkSolver, SolverConfig};
use cfd_core::compute::solver::Solver;
use cfd_core::physics::fluid::non_newtonian::CarreauYasuda;
use cfd_core::physics::fluid::FluidTrait;
use cfd_schematics::domain::model::{ChannelSpec, NodeKind, NodeSpec};

/// Hagen-Poiseuille resistance for a circular tube: R = 128·μ·L / (π·D⁴)
fn hagen_poiseuille_resistance(mu: f64, length_m: f64, diameter_m: f64) -> f64 {
    128.0 * mu * length_m / (std::f64::consts::PI * diameter_m.powi(4))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("╔══════════════════════════════════════════════════════════╗");
    println!("║  1D Blood Flow in TPMS-Inspired Channel Network        ║");
    println!("║  Gyroid lattice → lumped Hagen-Poiseuille solver        ║");
    println!("╚══════════════════════════════════════════════════════════╝");
    println!();

    // ── 1. Fluid: Carreau-Yasuda Blood ───────────────────────────────────────
    //
    // μ(γ̇) = μ_∞ + (μ₀ − μ_∞) · [1 + (λ·γ̇)^a]^((n−1)/a)
    //
    // Parameters (validated reference):
    //   ρ = 1060 kg/m³, μ₀ = 0.056 Pa·s, μ_∞ = 0.0035 Pa·s
    //   λ = 3.313 s, n = 0.3568, a = 2.0
    let blood = CarreauYasuda::<f64>::blood();
    println!("Fluid Model: {}", blood.name());
    println!("  Density       : {} kg/m³", blood.density);
    println!("  μ₀ (zero-shear) : {:.4} Pa·s", blood.viscosity_zero);
    println!("  μ_∞ (high-shear): {:.4} Pa·s", blood.viscosity_inf);
    println!();

    // ── 2. Network Topology ──────────────────────────────────────────────────
    //
    // TPMS Gyroid lattice with λ = 4 mm creates pore sizes of ~2 mm.
    // We model a representative cross-section as a multi-path network
    // with two parallel branches (upper and lower Gyroid channels).
    //
    // Physical dimensions (mm → m):
    //   Inlet channel   : D = 4.0 mm, L = 10.0 mm (port connection)
    //   Upper branch ch1: D = 2.0 mm, L = 50.0 mm (Gyroid pore path)
    //   Upper branch ch2: D = 2.0 mm, L = 50.0 mm
    //   Lower branch ch1: D = 2.5 mm, L = 45.0 mm (wider Gyroid pore)
    //   Lower branch ch2: D = 2.5 mm, L = 45.0 mm
    //   Outlet channel  : D = 4.0 mm, L = 10.0 mm (port connection)
    let nodes = vec![
        NodeSpec::new("Inlet", NodeKind::Inlet),
        NodeSpec::new("J1", NodeKind::Junction), // bifurcation
        NodeSpec::new("J2", NodeKind::Junction), // upper mid-point
        NodeSpec::new("J3", NodeKind::Junction), // recombination
        NodeSpec::new("J4", NodeKind::Junction), // lower mid-point
        NodeSpec::new("Outlet", NodeKind::Outlet),
    ];

    // Use mid-shear apparent viscosity for initial resistance estimates.
    // The Picard solver iterates to self-consistent viscosity.
    let mu_ref = 3.5e-3; // Pa·s (blood at ~100 s⁻¹)

    let channels = vec![
        // Inlet port → first bifurcation
        ChannelSpec::new_pipe(
            "ch_in",
            "Inlet",
            "J1",
            0.010,
            0.004,
            hagen_poiseuille_resistance(mu_ref, 0.010, 0.004),
            0.0,
        ),
        // Upper branch: J1 → J2 → J3
        ChannelSpec::new_pipe(
            "ch_upper_1",
            "J1",
            "J2",
            0.050,
            0.002,
            hagen_poiseuille_resistance(mu_ref, 0.050, 0.002),
            0.0,
        ),
        ChannelSpec::new_pipe(
            "ch_upper_2",
            "J2",
            "J3",
            0.050,
            0.002,
            hagen_poiseuille_resistance(mu_ref, 0.050, 0.002),
            0.0,
        ),
        // Lower branch: J1 → J4 → J3
        ChannelSpec::new_pipe(
            "ch_lower_1",
            "J1",
            "J4",
            0.045,
            0.0025,
            hagen_poiseuille_resistance(mu_ref, 0.045, 0.0025),
            0.0,
        ),
        ChannelSpec::new_pipe(
            "ch_lower_2",
            "J4",
            "J3",
            0.045,
            0.0025,
            hagen_poiseuille_resistance(mu_ref, 0.045, 0.0025),
            0.0,
        ),
        // Recombination → Outlet port
        ChannelSpec::new_pipe(
            "ch_out",
            "J3",
            "Outlet",
            0.010,
            0.004,
            hagen_poiseuille_resistance(mu_ref, 0.010, 0.004),
            0.0,
        ),
    ];

    println!("Network Topology (Gyroid-inspired):");
    for ch in &channels {
        println!(
            "  {} : {} → {}, D_h = {:.1} mm, L = {:.1} mm, R = {:.2e} Pa·s/m³",
            ch.id.as_str(),
            ch.from.as_str(),
            ch.to.as_str(),
            ch.cross_section.hydraulic_diameter() * 1e3,
            ch.length_m * 1e3,
            ch.resistance,
        );
    }
    println!();

    // ── 3. Build Network ─────────────────────────────────────────────────────
    let mut builder = NetworkBuilder::<f64>::new();

    let mut node_indices = std::collections::HashMap::new();
    for node_spec in &nodes {
        let idx = match node_spec.kind {
            NodeKind::Inlet => builder.add_inlet(node_spec.id.as_str().to_string()),
            NodeKind::Outlet => builder.add_outlet(node_spec.id.as_str().to_string()),
            NodeKind::Junction | NodeKind::Reservoir => {
                builder.add_junction(node_spec.id.as_str().to_string())
            }
        };
        node_indices.insert(node_spec.id.as_str().to_string(), idx);
    }

    let mut edge_indices = std::collections::HashMap::new();
    for ch in &channels {
        let from = node_indices[ch.from.as_str()];
        let to = node_indices[ch.to.as_str()];
        let eidx = builder.connect_with_pipe(from, to, ch.id.as_str().to_string());
        edge_indices.insert(ch.id.as_str().to_string(), eidx);
    }

    let graph = builder.build()?;
    let mut network = Network::new(graph, blood.clone());

    // Assign edge properties from channel specs
    for ch in &channels {
        let props = EdgeProperties::from(ch);
        network.add_edge_properties(edge_indices[ch.id.as_str()], props);
    }

    // ── 4. Boundary Conditions ───────────────────────────────────────────────
    //
    // Inlet: volumetric flow rate Q = 2 mL/min = 3.33e-8 m³/s
    //   (typical millifluidic blood processing rate)
    // Outlet: atmospheric pressure P = 101325 Pa
    let inlet_idx = node_indices["Inlet"];
    let outlet_idx = node_indices["Outlet"];
    let q_inlet = 2.0 / 60.0 * 1e-6; // 2 mL/min → m³/s

    network.set_neumann_flow(inlet_idx, q_inlet);
    network.set_pressure(outlet_idx, 101_325.0);

    // Initial flow rate guesses (split evenly between branches)
    network.set_flow_rate(edge_indices["ch_in"], q_inlet);
    network.set_flow_rate(edge_indices["ch_upper_1"], q_inlet / 2.0);
    network.set_flow_rate(edge_indices["ch_upper_2"], q_inlet / 2.0);
    network.set_flow_rate(edge_indices["ch_lower_1"], q_inlet / 2.0);
    network.set_flow_rate(edge_indices["ch_lower_2"], q_inlet / 2.0);
    network.set_flow_rate(edge_indices["ch_out"], q_inlet);
    network.update_resistances()?;

    println!("Boundary Conditions:");
    println!("  Inlet flow  : {:.2} mL/min ({:.2e} m³/s)", 2.0, q_inlet);
    println!("  Outlet press: 101325 Pa (atmospheric)");
    println!();

    // ── 5. Solve ─────────────────────────────────────────────────────────────
    let config = SolverConfig {
        tolerance: 1e-6,
        max_iterations: 100,
    };
    let solver = NetworkSolver::<f64, CarreauYasuda<f64>>::with_config(config);
    println!("Solving (Anderson-accelerated Picard iteration) ...");
    let solution = solver.solve(&NetworkProblem::new(network))?;

    // ── 6. Results ───────────────────────────────────────────────────────────
    println!();
    println!("═══════════════════════════════════════════════════════");
    println!("  NODE PRESSURES");
    println!("═══════════════════════════════════════════════════════");
    for idx in solution.graph.node_indices() {
        let n = solution.graph.node_weight(idx).unwrap();
        let p = *solution.pressures.get(idx.index()).unwrap_or(&0.0);
        println!(
            "  {:10} : {:>10.2} Pa  ({:>7.2} mmHg)",
            n.id,
            p,
            p / 133.322
        );
    }

    println!();
    println!("═══════════════════════════════════════════════════════");
    println!("  CHANNEL FLOW RATES & HEMODYNAMICS");
    println!("═══════════════════════════════════════════════════════");
    for idx in solution.graph.edge_indices() {
        let e = solution.graph.edge_weight(idx).unwrap();
        let q = *solution.flow_rates.get(idx.index()).unwrap_or(&0.0);
        let props = solution.properties.get(&idx);

        let (velocity, shear_rate) = if let Some(p) = props {
            let vel = q.abs() / p.area;
            let sr = p.hydraulic_diameter.map(|d| 8.0 * vel / d).unwrap_or(0.0);
            (vel, sr)
        } else {
            (0.0, 0.0)
        };

        let viscosity = blood
            .viscosity_at_shear(shear_rate, 310.15, 101325.0)
            .unwrap_or(blood.viscosity_inf);
        let wall_shear = viscosity * shear_rate;

        println!("  {}", e.id);
        println!("    Flow rate         : {:.4} mL/min", q.abs() * 1e6 * 60.0);
        println!("    Mean velocity     : {:.2} cm/s", velocity * 100.0);
        println!("    Wall shear rate   : {:.1} s⁻¹", shear_rate);
        println!("    Apparent viscosity: {:.3} mPa·s", viscosity * 1000.0);
        println!("    Wall shear stress : {:.3} Pa", wall_shear);
    }

    // ── 7. Flow Distribution Analysis ────────────────────────────────────────
    let q_upper = solution
        .flow_rates
        .get(edge_indices["ch_upper_1"].index())
        .unwrap_or(&0.0)
        .abs();
    let q_lower = solution
        .flow_rates
        .get(edge_indices["ch_lower_1"].index())
        .unwrap_or(&0.0)
        .abs();
    let q_total = q_upper + q_lower;

    println!();
    println!("═══════════════════════════════════════════════════════");
    println!("  TPMS LATTICE FLOW DISTRIBUTION");
    println!("═══════════════════════════════════════════════════════");
    println!(
        "  Upper branch (D=2.0mm): {:.1}%",
        q_upper / q_total * 100.0
    );
    println!(
        "  Lower branch (D=2.5mm): {:.1}%",
        q_lower / q_total * 100.0
    );
    println!();
    println!("  The lower branch carries more flow due to its larger");
    println!("  diameter (lower Hagen-Poiseuille resistance, R ∝ D⁻⁴).");
    println!("  This is analogous to how different-sized Gyroid pores");
    println!("  naturally partition flow by their hydraulic resistance.");

    // ── 8. JSON Export ───────────────────────────────────────────────────────
    let output_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("tpms_blood_1d");
    std::fs::create_dir_all(&output_dir)?;

    let results = serde_json::json!({
        "simulation": "1D TPMS Blood Flow",
        "fluid": "Carreau-Yasuda Blood",
        "inlet_flow_mL_per_min": 2.0,
        "tpms_surface": "Gyroid",
        "tpms_period_mm": 4.0,
        "upper_branch_flow_fraction": q_upper / q_total,
        "lower_branch_flow_fraction": q_lower / q_total,
    });

    let json_path = output_dir.join("results.json");
    std::fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;
    println!("  Results exported to: {}", json_path.display());

    Ok(())
}
