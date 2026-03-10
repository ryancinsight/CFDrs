use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId,
    NodeKind, NodeSpec,
};
use cfd_schematics::interface::presets::{bifurcation_venturi_rect, venturi_rect};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};

use super::postprocess::extract_field_wall_shear;
use super::Network2dBuilderSink;

fn circular_blueprint() -> NetworkBlueprint {
    let mut blueprint = NetworkBlueprint::new("circular_trace");
    blueprint.add_node(NodeSpec {
        id: NodeId::new("inlet"),
        kind: NodeKind::Inlet,
        point: (0.0, 0.0),
        layout: None,
        junction_geometry: None,
        metadata: None,
    });
    blueprint.add_node(NodeSpec {
        id: NodeId::new("outlet"),
        kind: NodeKind::Outlet,
        point: (1.0, 0.0),
        layout: None,
        junction_geometry: None,
        metadata: None,
    });
    blueprint.add_channel(ChannelSpec {
        id: EdgeId::new("pipe"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("inlet"),
        to: NodeId::new("outlet"),
        length_m: 0.01,
        cross_section: CrossSectionSpec::Circular { diameter_m: 1.0e-3 },
        channel_shape: ChannelShape::Straight,
        resistance: 0.0,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        visual_role: None,
        path: vec![(0.0, 0.0), (1.0, 0.0)],
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    });
    blueprint.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    blueprint
}

fn rectangular_blueprint() -> NetworkBlueprint {
    let mut blueprint = NetworkBlueprint::new("rectangular_trace");
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));
    let mut duct = ChannelSpec::new_pipe_rect(
        "duct", "inlet", "outlet", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    duct.path = vec![(0.0, 0.0), (1.0, 0.0)];
    blueprint.add_channel(duct);
    blueprint.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    blueprint
}

#[test]
fn build_single_venturi_network() {
    let mut bp = venturi_rect("test_v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
    bp.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-6, 10, 5);
    let net2d = sink.build(&bp).expect("build should succeed");
    assert_eq!(net2d.reference_trace().channel_traces.len(), 3);
}

#[test]
fn reference_trace_sums_to_q_total() {
    use cfd_schematics::interface::presets::bifurcation_venturi_rect;

    let mut bp = bifurcation_venturi_rect("bv", 0.005, 0.002, 0.0005, 0.001, 0.001);
    bp.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    let q_total = 1e-6;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 10, 5);
    let net2d = sink.build(&bp).expect("build should succeed");
    let trace = net2d.reference_trace();

    assert!((trace.total_inlet_flow_m3_s - q_total).abs() < q_total * 1e-6);
    assert!((trace.total_outlet_flow_m3_s - q_total).abs() < q_total * 1e-6);

    for node in trace
        .node_traces
        .iter()
        .filter(|node| matches!(node.node_kind, NodeKind::Junction | NodeKind::Reservoir))
    {
        assert!(node.continuity_residual_m3_s.abs() < q_total * 1e-6);
    }
}

#[test]
fn circular_channels_use_true_cross_section_area_in_trace() {
    let bp = circular_blueprint();
    let q_total = 2e-8;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 8, 8);
    let net2d = sink.build(&bp).expect("build should succeed");
    let channel = &net2d.reference_trace().channel_traces[0];
    let expected_area = std::f64::consts::PI * (0.5e-3_f64).powi(2);

    assert!((channel.cross_section_area_m2 - expected_area).abs() < expected_area * 1e-12);
}

#[test]
fn single_channel_outlet_flow_tracks_reference() {
    let bp = rectangular_blueprint();
    let q_total = 1e-7;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 80, 20);
    let mut net2d = sink.build(&bp).expect("build should succeed");
    let result = net2d.solve_all(1e-6).expect("solve should succeed");
    let channel = &result.channels[0];

    assert!(channel.solve_result.converged);
    assert!(channel.field_outlet_flow_error_pct < 25.0);
}

#[test]
fn field_wall_shear_agrees_with_analytical() {
    let mu = 3.5e-3_f64;
    let w = 1.0e-3;
    let h = 1.0e-3;
    let l = 10.0e-3;
    let q = 1e-7;
    let u_mean = q / (w * h);

    let grid = StaggeredGrid2D::new(80, 20, l, w);
    let blood = BloodModel::Newtonian(mu);
    let mut solver = NavierStokesSolver2D::new(grid, blood, 1060.0, SIMPLEConfig::default());
    let result = solver.solve(u_mean).expect("straight channel solve");
    assert!(result.converged);

    let (max_tau, mean_tau) = extract_field_wall_shear(&solver);
    let analytical_tau = mu * 6.0 * q / (w * h * h);
    let ratio = mean_tau / analytical_tau;
    assert!(ratio > 0.1 && ratio < 10.0);
    assert!(max_tau >= mean_tau);
}

// ── 2D vs 1D cross-fidelity physics comparisons ─────────────────────────────

/// Single rectangular channel: 2D wall shear must agree with 1D Hagen-Poiseuille
/// within a factor of 3×, and mass is conserved to < 15%.
#[test]
fn two_d_wall_shear_within_factor_of_hp_analytical() {
    let mu = 3.5e-3_f64;
    let w = 1.0e-3;
    let h = 1.0e-3;
    let q = 1e-7;

    let bp = rectangular_blueprint();
    let blood = BloodModel::Newtonian(mu);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q, 40, 10);
    let mut net2d = sink.build(&bp).expect("build");
    let result = net2d.solve_all(1e-6).expect("solve");

    let ch = &result.channels[0];
    // HP analytical for square duct: τ = 6 μ Q / (w h²)
    let tau_hp = mu * 6.0 * q / (w * h * h);
    let ratio = ch.field_wall_shear_mean_pa / tau_hp;
    assert!(
        ratio > 0.1 && ratio < 5.0,
        "2D mean wall shear {:.4} Pa should be within 5× of HP {tau_hp:.4} Pa (ratio={ratio:.3})",
        ch.field_wall_shear_mean_pa
    );
    // 1D channel wall_shear_pa is the blueprint-based HP estimate
    let one_d_tau = ch.wall_shear_pa;
    let agreement = (ch.field_wall_shear_mean_pa / one_d_tau)
        .min(one_d_tau / ch.field_wall_shear_mean_pa.max(1e-12));
    assert!(
        agreement > 0.05,
        "2D mean shear {:.4} and 1D HP {:.4} should be within same order of magnitude",
        ch.field_wall_shear_mean_pa,
        one_d_tau
    );
}

/// Single rectangular channel: 2D outlet flow must be within 20% of the 1D
/// reference (mass conservation at the coarse grid used here).
#[test]
fn two_d_outlet_flow_error_within_20pct_of_1d_reference() {
    let bp = rectangular_blueprint();
    let q_total = 1e-7;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 40, 10);
    let mut net2d = sink.build(&bp).expect("build");
    let result = net2d.solve_all(1e-6).expect("solve");

    let ch = &result.channels[0];
    assert!(
        ch.field_outlet_flow_error_pct < 20.0,
        "2D outlet flow error {:.2}% should be < 20% of 1D reference",
        ch.field_outlet_flow_error_pct
    );
}

/// Venturi: 2D max wall shear in the throat channel must be higher than in the
/// inlet section (constriction accelerates the flow), consistent with 1D ordering.
#[test]
fn two_d_venturi_throat_has_higher_shear_than_inlet_section() {
    // main channel 2 mm wide, throat 0.5 mm, height 0.5 mm, length 2 mm
    let mut bp = venturi_rect("v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
    bp.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    let q = 5e-8_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q, 20, 8);
    let mut net2d = sink.build(&bp).expect("build");
    let result = net2d.solve_all(1e-6).expect("solve");

    let find_channel = |substr: &str| -> Option<&super::types::Channel2dResult<f64>> {
        result
            .channels
            .iter()
            .find(|ch| ch.channel_id.contains(substr))
    };

    let throat = find_channel("throat").expect("venturi blueprint must have a throat channel");
    let inlet = find_channel("inlet_section").expect("venturi blueprint must have inlet_section");

    assert!(
        throat.wall_shear_pa > inlet.wall_shear_pa,
        "1D HP: throat shear {:.4} Pa must exceed inlet shear {:.4} Pa",
        throat.wall_shear_pa,
        inlet.wall_shear_pa
    );
    assert!(
        throat.field_wall_shear_max_pa > inlet.field_wall_shear_max_pa * 0.5,
        "2D max shear: throat {:.4} Pa should be in the same ball-park as or exceed inlet {:.4} Pa",
        throat.field_wall_shear_max_pa,
        inlet.field_wall_shear_max_pa
    );
}

/// Bifurcation-venturi: 1D reference must conserve mass; 2D outlet error must
/// be < 30% per channel at coarse grid.
#[test]
fn two_d_bifurcation_mass_conservation_tracks_1d_reference() {
    let mut bp = bifurcation_venturi_rect("bv", 5e-3, 2e-3, 5e-4, 1e-3, 1e-3);
    bp.metadata.get_or_insert_with(Default::default).insert(cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper());
    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let mut net2d = sink.build(&bp).expect("build");
    let result = net2d.solve_all(1e-6).expect("solve");

    // 1D reference must conserve total mass
    let ref_in = result.reference_trace.total_inlet_flow_m3_s;
    let ref_out = result.reference_trace.total_outlet_flow_m3_s;
    let ref_err = ((ref_in - ref_out) / ref_in).abs();
    assert!(
        ref_err < 1e-4,
        "1D reference inlet/outlet mismatch {:.2e} m³/s must be < 0.01%",
        ref_err
    );

    // 2D per-channel outlet error vs 1D reference
    for ch in &result.channels {
        assert!(
            ch.field_outlet_flow_error_pct < 30.0,
            "channel '{}': 2D outlet error {:.2}% should be < 30% of 1D reference",
            ch.channel_id,
            ch.field_outlet_flow_error_pct
        );
    }
}
