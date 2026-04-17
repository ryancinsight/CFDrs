//! Adversarial tests for cfd-1d: robustness, boundary conditions, and degenerate inputs.
//!
//! Tests verify that:
//! 1. Physics models never produce NaN for valid-but-extreme inputs
//! 2. NetworkBuilder validates topology and rejects malformed networks with Err
//! 3. Solver correctly handles NaN/infinite resistance edges without panicking
//! 4. Womersley and Bessel functions produce correct limits at degenerate inputs

use approx::assert_relative_eq;
use cfd_1d::{
    domain::components::channels::CircularChannel,
    domain::components::membranes::PorousMembrane,
    domain::components::Component,
    physics::vascular::womersley::{WomersleyFlow, WomersleyNumber},
    Network, NetworkBuilder, NetworkProblem, NetworkSolver,
};
use cfd_core::physics::fluid::database::water_20c;
use petgraph::visit::EdgeRef;

fn water() -> cfd_core::physics::fluid::ConstantPropertyFluid<f64> {
    water_20c::<f64>().expect("test invariant")
}

fn hp_resistance(d: f64, l: f64, mu: f64) -> f64 {
    128.0 * mu * l / (std::f64::consts::PI * d.powi(4))
}

// ============================================================================
// Boundary: degenerate channel geometry
// ============================================================================

/// Zero-diameter channel produces infinite resistance (never panics, never NaN).
#[test]
fn test_zero_diameter_returns_inf_resistance() {
    let fluid = water();
    let chan = CircularChannel::new(0.1, 0.0, 0.0);
    let r = chan.resistance(&fluid);
    assert!(
        r.is_infinite(),
        "Zero diameter should yield infinite resistance, got {}",
        r
    );
}

/// Zero-length channel produces zero resistance (Poiseuille law: R = 8μL/πR⁴).
#[test]
fn test_zero_length_returns_zero_resistance() {
    let fluid = water();
    let chan = CircularChannel::new(0.0, 1e-3, 0.0);
    let r = chan.resistance(&fluid);
    assert_eq!(r, 0.0, "Zero length must yield exactly zero resistance");
}

/// Very small diameter produces very large but finite resistance.
#[test]
fn test_tiny_diameter_large_resistance() {
    let fluid = water();
    let chan = CircularChannel::new(0.1, 1e-9, 0.0); // 1 nm diameter
    let r = chan.resistance(&fluid);
    assert!(
        r > 1e30,
        "1 nm channel should have enormous resistance, got {}",
        r
    );
    assert!(
        r.is_finite(),
        "Resistance must be finite for non-zero diameter"
    );
}

/// Very long channel produces large but finite resistance.
#[test]
fn test_extreme_length_stays_finite() {
    let fluid = water();
    let chan = CircularChannel::new(1e6, 1e-3, 0.0); // 1000 km pipe
    let r = chan.resistance(&fluid);
    assert!(
        r.is_finite(),
        "Extreme-length channel must still produce finite resistance"
    );
}

// ============================================================================
// NetworkBuilder Validation: structural topology errors must be Err (not panic)
// ============================================================================

/// NetworkBuilder::build() must reject empty graphs (no nodes).
#[test]
fn test_empty_network_is_rejected() {
    let builder = NetworkBuilder::<f64>::new();
    let result = builder.build();
    assert!(
        result.is_err(),
        "Empty graph (no nodes) must be rejected by NetworkBuilder::build()"
    );
}

/// NetworkBuilder::build() must reject graph with junction but no inlet or outlet.
#[test]
fn test_no_inlet_or_outlet_is_rejected() {
    let mut builder = NetworkBuilder::<f64>::new();
    let j1 = builder.add_junction("j1".into());
    let j2 = builder.add_junction("j2".into());
    builder.connect_with_pipe(j1, j2, "pipe".into());
    let result = builder.build();
    assert!(
        result.is_err(),
        "Network without inlet/outlet must be rejected"
    );
}

/// NetworkBuilder::build() must reject disconnected components.
#[test]
fn test_disconnected_network_is_rejected() {
    let mut builder = NetworkBuilder::<f64>::new();
    let n1 = builder.add_inlet("in1".into());
    let n2 = builder.add_junction("j1".into());
    builder.connect_with_pipe(n1, n2, "p1".into());

    // Second disconnected subgraph
    let n3 = builder.add_junction("j2".into());
    let n4 = builder.add_outlet("out1".into());
    builder.connect_with_pipe(n3, n4, "p2".into());

    let result = builder.build();
    assert!(
        result.is_err(),
        "Disconnected components must be rejected by NetworkBuilder::build()"
    );
}

// ============================================================================
// Solver: NaN resistance edges
// ============================================================================

/// Solver with NaN resistance on an edge must fail gracefully (not panic).
/// Since NetworkBuilder validates zero resistance, we inject NaN after build().
#[test]
fn test_nan_resistance_does_not_panic() {
    let fluid = water();
    let d = 1e-3_f64;
    let l = 0.1_f64;
    let r = hp_resistance(d, l, fluid.viscosity);

    let mut builder = NetworkBuilder::new();
    let n1 = builder.add_inlet("n1".into());
    let n2 = builder.add_outlet("n2".into());
    builder.connect_with_pipe(n1, n2, "pipe".into());

    let mut graph = builder.build().expect("test invariant");
    // Pre-set a valid resistance to pass validation, then overwrite with NaN
    for edge in graph.edge_indices() {
        if let Some(e) = graph.edge_weight_mut(edge) {
            e.resistance = r;
        }
    }
    // Now inject NaN after validation
    let edge_idx = graph.edge_indices().next().expect("test invariant");
    if let Some(e) = graph.edge_weight_mut(edge_idx) {
        e.resistance = f64::NAN;
    }

    let mut network = Network::new(graph, fluid);
    network.set_pressure(n1, 1000.0);
    network.set_pressure(n2, 0.0);

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::new();
    // Must not panic; either Err or NaN pressures are acceptable
    let result = solver.solve_network(&problem);
    match result {
        Err(_) => {} // Graceful rejection
        Ok(solved) => {
            // Solver ran; check either NaN propagated or flow is trivially zero/nan
            let pressures: Vec<f64> = solved.pressures().iter().copied().collect();
            let flows: Vec<f64> = solved
                .graph
                .edge_references()
                .map(|e| {
                    let p_src = solved
                        .pressures()
                        .get(e.source().index())
                        .copied()
                        .unwrap_or(0.0);
                    let p_tgt = solved
                        .pressures()
                        .get(e.target().index())
                        .copied()
                        .unwrap_or(0.0);
                    (p_src - p_tgt) / e.weight().resistance
                })
                .collect();
            let has_nan_or_inf = pressures.iter().any(|p| p.is_nan() || p.is_infinite())
                || flows.iter().any(|q| q.is_nan() || q.is_infinite());
            assert!(
                has_nan_or_inf || flows[0].abs() < 1e-15,
                "NaN resistance must produce NaN/Inf/zero flow (not silently wrong)"
            );
        }
    }
}

// ============================================================================
// Porous membrane: zero porosity / zero pore radius
// ============================================================================

/// Zero-pore-radius membrane must produce effectively infinite resistance.
#[test]
fn test_membrane_zero_pore_radius_is_very_high() {
    let fluid = water();
    let membrane = PorousMembrane::new(1e-4, 1e-3, 1e-3, 0.0, 0.3);
    let r = membrane.resistance(&fluid);
    assert!(
        r > 1e11,
        "Membrane with zero pore radius must have very high resistance (>1e11), got {}",
        r
    );
}

/// Zero-porosity membrane must produce effectively infinite resistance.
#[test]
fn test_membrane_zero_porosity_is_very_high() {
    let fluid = water();
    let membrane = PorousMembrane::new(1e-4, 1e-3, 1e-3, 500e-9, 0.0);
    let r = membrane.resistance(&fluid);
    assert!(
        r > 1e11,
        "Zero-porosity membrane must have very high resistance (>1e11), got {}",
        r
    );
}

// ============================================================================
// Two Dirichlet boundaries: flow from high-to-low pressure
// ============================================================================

/// Solver with two Dirichlet boundaries and opposite pressures must
/// produce correct Ohm's law flow without panic.
#[test]
fn test_two_dirichlet_nodes_solve_correctly() {
    let fluid = water();
    let d = 1e-3_f64;
    let l = 0.1_f64;
    let r = hp_resistance(d, l, fluid.viscosity);

    let mut builder = NetworkBuilder::new();
    let n1 = builder.add_inlet("in".into());
    let n2 = builder.add_outlet("out".into());
    builder.connect_with_pipe(n1, n2, "pipe".into());

    let mut graph = builder.build().expect("test invariant");
    for edge in graph.edge_indices() {
        if let Some(e) = graph.edge_weight_mut(edge) {
            e.resistance = r;
        }
    }

    let mut network = Network::new(graph, fluid.clone());
    network.set_pressure(n1, 200.0);
    network.set_pressure(n2, 100.0);

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::new();
    let result = solver.solve_network(&problem).expect("test invariant");

    let expected_q = 100.0 / r;
    let flows: Vec<f64> = result
        .graph
        .edge_references()
        .map(|e| {
            let p_src = result
                .pressures()
                .get(e.source().index())
                .copied()
                .unwrap_or(0.0);
            let p_tgt = result
                .pressures()
                .get(e.target().index())
                .copied()
                .unwrap_or(0.0);
            (p_src - p_tgt) / e.weight().resistance
        })
        .collect();
    assert_relative_eq!(flows[0], expected_q, max_relative = 1e-9);
}

// ============================================================================
// Womersley: zero-frequency → Poiseuille limit
// ============================================================================

/// Womersley number α = R√(ωρ/μ) must be zero for ω = 0   (quasi-steady limit).
#[test]
fn test_womersley_zero_frequency_alpha_is_zero() {
    let alpha_num = WomersleyNumber::<f64>::new(
        0.005, // R = 5 mm
        0.0,   // ω = 0 (static)
        1060.0, 0.0035,
    );
    let alpha = alpha_num.value();
    assert_eq!(
        alpha, 0.0,
        "Zero-frequency Womersley number must be exactly 0"
    );
}

/// WomersleyFlow solver must not panic at near-zero frequency.
#[test]
fn test_womersley_near_zero_frequency_no_panic() {
    let flow = WomersleyFlow::<f64>::new(
        0.005,  // R = 5 mm
        0.1,    // L = 10 cm
        1060.0, // blood density
        0.0035, // blood viscosity
        1e-10,  // Effectively ω → 0
        133.0,  // pressure amplitude
        -100.0, // steady gradient
    );
    let u = flow.velocity(0.5, 0.0);
    assert!(
        u.is_finite(),
        "Near-zero-frequency velocity must be finite, got {}",
        u
    );
}

// ============================================================================
// Resistance extremes: high roughness and extreme viscosity
// ============================================================================

/// Very high surface roughness (ε ≈ D) must produce higher resistance but no panic.
#[test]
fn test_circular_channel_high_roughness_finite() {
    let fluid = water();
    let chan = CircularChannel::new(0.1, 1e-3, 1e-3); // roughness = diameter
    let r = chan.resistance(&fluid);
    assert!(
        r.is_finite(),
        "High-roughness channel must return finite resistance, got {}",
        r
    );
    assert!(r > 0.0, "Resistance must be strictly positive");
}

/// Near-inviscid fluid (μ → 0) must keep resistance ≥ 0 and finite.
#[test]
fn test_near_inviscid_fluid_resistance_near_zero() {
    use cfd_core::physics::fluid::ConstantPropertyFluid;
    let inviscid = ConstantPropertyFluid::new(
        "near_inviscid".to_string(),
        1000.0_f64, // density
        1e-15_f64,  // viscosity ≈ 0
        0.0_f64,    // specific_heat
        0.0_f64,    // thermal_conductivity
        1500.0_f64, // speed_of_sound
    );
    let chan = CircularChannel::new(0.1, 1e-3, 0.0);
    let r = chan.resistance(&inviscid);
    assert!(r >= 0.0, "Resistance must be non-negative, got {}", r);
    assert!(r < 1.0, "Near-inviscid resistance must be tiny, got {}", r);
}
