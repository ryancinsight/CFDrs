//! Cross-fidelity validation: public cfd-1d branching API vs 2D junction solvers.
//!
//! These tests close a remaining phase-1 gap: existing branching validations
//! compare 2D junction solvers against analytical split laws, but they do not
//! explicitly exercise the real `cfd-1d` public solve path (`solve_reference_trace`)
//! on the same bifurcation contracts.

use cfd_2d::network::solve_reference_trace;
use cfd_2d::solvers::{BifurcationGeometry, BifurcationSolver2D};
use cfd_core::physics::fluid::BloodModel;
use cfd_math::pressure_velocity::SIMPLEConfig;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId,
    NodeKind, NodeSpec,
};
use cfd_schematics::geometry::metadata::GeometryAuthoringProvenance;

const MU: f64 = 3.5e-3;
const RHO: f64 = 1060.0;
const DEPTH_M: f64 = 0.5e-3;

fn rectangular_bifurcation_blueprint(
    name: &str,
    parent_width_m: f64,
    parent_length_m: f64,
    daughter1_width_m: f64,
    daughter2_width_m: f64,
    daughter_length_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new_with_explicit_positions(name);
    bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (-1.0, 0.0)));
    bp.add_node(NodeSpec::new_at("junction", NodeKind::Junction, (0.0, 0.0)));
    bp.add_node(NodeSpec::new_at("daughter1_outlet", NodeKind::Outlet, (1.0, 1.0)));
    bp.add_node(NodeSpec::new_at("daughter2_outlet", NodeKind::Outlet, (1.0, -1.0)));

    let rectangular_section = |width_m: f64| CrossSectionSpec::Rectangular {
        width_m,
        height_m: DEPTH_M,
    };

    bp.add_channel(ChannelSpec {
        id: EdgeId::new("parent"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("inlet"),
        to: NodeId::new("junction"),
        length_m: parent_length_m,
        cross_section: rectangular_section(parent_width_m),
        channel_shape: ChannelShape::Straight,
        resistance: 0.0,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        visual_role: None,
        path: vec![(-1.0, 0.0), (0.0, 0.0)],
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    });
    bp.add_channel(ChannelSpec {
        id: EdgeId::new("daughter1"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("junction"),
        to: NodeId::new("daughter1_outlet"),
        length_m: daughter_length_m,
        cross_section: rectangular_section(daughter1_width_m),
        channel_shape: ChannelShape::Straight,
        resistance: 0.0,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        visual_role: None,
        path: vec![(0.0, 0.0), (1.0, 1.0)],
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    });
    bp.add_channel(ChannelSpec {
        id: EdgeId::new("daughter2"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("junction"),
        to: NodeId::new("daughter2_outlet"),
        length_m: daughter_length_m,
        cross_section: rectangular_section(daughter2_width_m),
        channel_shape: ChannelShape::Straight,
        resistance: 0.0,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        visual_role: None,
        path: vec![(0.0, 0.0), (1.0, -1.0)],
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    });
    bp.metadata
        .get_or_insert_with(Default::default)
        .insert(GeometryAuthoringProvenance::selective_wrapper());
    bp
}

fn channel_flow(trace: &cfd_2d::network::NetworkReferenceTrace<f64>, channel_id: &str) -> f64 {
    trace
        .channel_traces
        .iter()
        .find(|channel| channel.channel_id == channel_id)
        .map(|channel| channel.flow_rate_m3_s.abs())
        .unwrap_or_else(|| panic!("Missing 1D channel trace for '{channel_id}'"))
}

/// ## Theorem (Public 1D/2D Symmetric Bifurcation Agreement)
///
/// For a symmetric rectangular bifurcation with identical daughter branches,
/// the public `cfd-1d` solve path and the 2D junction solver must predict the
/// same equal-flow partition up to the junction-loss tolerance already accepted
/// by the 2D suite.
#[test]
fn cross_fidelity_public_api_symmetric_bifurcation() {
    let parent_width_m = 2.0e-3_f64;
    let parent_length_m = 6.0e-3_f64;
    let daughter_width_m = 1.5e-3_f64;
    let daughter_length_m = 6.0e-3_f64;
    let angle_rad = 0.3_f64;
    let inlet_velocity_m_s = 0.005_f64;
    let target_total_flow_m3_s = inlet_velocity_m_s * parent_width_m * DEPTH_M;

    let blueprint = rectangular_bifurcation_blueprint(
        "public-symmetric-bifurcation",
        parent_width_m,
        parent_length_m,
        daughter_width_m,
        daughter_width_m,
        daughter_length_m,
    );
    let trace_1d = solve_reference_trace::<f64>(&blueprint, RHO, MU, target_total_flow_m3_s)
        .expect("1D reference trace");

    let trace_mass_error =
        ((trace_1d.total_inlet_flow_m3_s - trace_1d.total_outlet_flow_m3_s) / trace_1d.total_inlet_flow_m3_s)
            .abs();
    assert!(trace_mass_error < 1e-4, "1D mass error {trace_mass_error} exceeds 0.01%");

    let q1_1d = channel_flow(&trace_1d, "daughter1");
    let q2_1d = channel_flow(&trace_1d, "daughter2");
    let share1_1d = q1_1d / (q1_1d + q2_1d).max(1e-30);

    let geometry = BifurcationGeometry {
        parent_width: parent_width_m,
        parent_length: parent_length_m,
        daughter1_width: daughter_width_m,
        daughter1_length: daughter_length_m,
        daughter1_angle: angle_rad,
        daughter2_width: daughter_width_m,
        daughter2_length: daughter_length_m,
        daughter2_angle: -angle_rad,
    };
    let config = SIMPLEConfig {
        max_iterations: 500,
        ..SIMPLEConfig::default()
    };
    let mut solver = BifurcationSolver2D::new(
        geometry,
        BloodModel::Newtonian(MU),
        RHO,
        60,
        40,
        config,
    );
    let result = solver.solve(inlet_velocity_m_s).expect("2D bifurcation solve");

    let q1_2d = result.q_daughter1.abs();
    let q2_2d = result.q_daughter2.abs();
    let share1_2d = q1_2d / (q1_2d + q2_2d).max(1e-30);

    assert!(result.mass_balance_error < 0.10, "2D mass error {} exceeds 10%", result.mass_balance_error);
    assert!(
        (share1_2d - share1_1d).abs() < 0.10,
        "Symmetric bifurcation daughter-1 share mismatch: 1D={share1_1d:.3}, 2D={share1_2d:.3}"
    );
}

/// ## Theorem (Public 1D/2D Asymmetric Bifurcation Ordering)
///
/// For an asymmetric rectangular bifurcation with identical daughter lengths,
/// both the public `cfd-1d` solve path and the 2D junction solver must route a
/// larger fraction of the total flow into the wider daughter, and their wider-
/// daughter shares must remain within a junction-scale tolerance band.
#[test]
fn cross_fidelity_public_api_asymmetric_bifurcation() {
    let parent_width_m = 2.0e-3_f64;
    let parent_length_m = 10.0e-3_f64;
    let daughter1_width_m = 1.5e-3_f64;
    let daughter2_width_m = 0.75e-3_f64;
    let daughter_length_m = 20.0e-3_f64;
    let angle_rad = 0.2_f64;
    let inlet_velocity_m_s = 0.005_f64;
    let target_total_flow_m3_s = inlet_velocity_m_s * parent_width_m * DEPTH_M;

    let blueprint = rectangular_bifurcation_blueprint(
        "public-asymmetric-bifurcation",
        parent_width_m,
        parent_length_m,
        daughter1_width_m,
        daughter2_width_m,
        daughter_length_m,
    );
    let trace_1d = solve_reference_trace::<f64>(&blueprint, RHO, MU, target_total_flow_m3_s)
        .expect("1D reference trace");

    let q1_1d = channel_flow(&trace_1d, "daughter1");
    let q2_1d = channel_flow(&trace_1d, "daughter2");
    let share1_1d = q1_1d / (q1_1d + q2_1d).max(1e-30);

    let geometry = BifurcationGeometry {
        parent_width: parent_width_m,
        parent_length: parent_length_m,
        daughter1_width: daughter1_width_m,
        daughter1_length: daughter_length_m,
        daughter1_angle: angle_rad,
        daughter2_width: daughter2_width_m,
        daughter2_length: daughter_length_m,
        daughter2_angle: -angle_rad,
    };
    let config = SIMPLEConfig {
        max_iterations: 1000,
        ..SIMPLEConfig::default()
    };
    let mut solver = BifurcationSolver2D::new(
        geometry,
        BloodModel::Newtonian(MU),
        RHO,
        100,
        60,
        config,
    );
    let result = solver.solve(inlet_velocity_m_s).expect("2D bifurcation solve");

    let q1_2d = result.q_daughter1.abs();
    let q2_2d = result.q_daughter2.abs();
    let share1_2d = q1_2d / (q1_2d + q2_2d).max(1e-30);

    assert!(share1_1d > 0.5, "1D wider daughter must receive more than half the flow");
    assert!(share1_2d > 0.5, "2D wider daughter must receive more than half the flow");
    assert!(result.mass_balance_error < 0.10, "2D mass error {} exceeds 10%", result.mass_balance_error);
    assert!(
        (share1_2d - share1_1d).abs() < 0.20,
        "Asymmetric bifurcation wider-daughter share mismatch: 1D={share1_1d:.3}, 2D={share1_2d:.3}"
    );
}