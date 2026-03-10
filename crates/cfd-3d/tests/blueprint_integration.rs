use cfd_3d::blueprint_integration::{
    process_blueprint_with_reference_trace, Blueprint3dProcessingConfig,
};
use cfd_schematics::interface::presets::{
    bifurcation_serpentine_rect, double_bifurcation_serpentine_rect, serpentine_chain,
    symmetric_bifurcation, trifurcation_serpentine_rect, venturi_chain,
};

fn reference_flow<'a>(
    trace: &'a cfd_3d::blueprint_integration::Blueprint3dTrace,
    channel_id: &str,
) -> f64 {
    trace
        .channel_traces
        .iter()
        .find(|channel| channel.channel_id == channel_id)
        .unwrap_or_else(|| panic!("missing channel trace for '{channel_id}'"))
        .reference_flow_rate_m3_s
}

fn continuity_residual<'a>(
    trace: &'a cfd_3d::blueprint_integration::Blueprint3dTrace,
    node_id: &str,
) -> f64 {
    trace
        .node_traces
        .iter()
        .find(|node| node.node_id == node_id)
        .unwrap_or_else(|| panic!("missing node trace for '{node_id}'"))
        .continuity_residual_m3_s
}

fn no_chip_trace_config() -> Blueprint3dProcessingConfig {
    Blueprint3dProcessingConfig {
        mesh: cfd_mesh::application::pipeline::PipelineConfig {
            include_chip_body: false,
            ..Default::default()
        },
        total_flow_rate_m3_s: 1.0e-9,
        run_2d_reference: false,
        ..Default::default()
    }
}

#[test]
fn blueprint_trace_captures_mesh_and_reference_contract() {
    let blueprint = symmetric_bifurcation("trace_bif", 0.010, 0.010, 0.004, 0.003);
    let config = no_chip_trace_config();

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should build for symmetric bifurcation");

    assert_eq!(trace.channel_traces.len(), blueprint.channels.len());
    assert_eq!(trace.node_traces.len(), blueprint.nodes.len());
    assert!(trace.volume_trace.fluid_mesh_volume_mm3 > 0.0);
    assert!(trace
        .channel_traces
        .iter()
        .all(|channel| channel.reference_flow_rate_m3_s.abs() > 0.0));
}

#[test]
fn blueprint_trace_tracks_serpentine_connectors_and_node_continuity() {
    let blueprint = serpentine_chain("trace_s", 4, 0.010, 0.004);
    let flow_rate = 1.0e-9;
    let config = Blueprint3dProcessingConfig {
        total_flow_rate_m3_s: flow_rate,
        ..no_chip_trace_config()
    };

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should preserve serpentine connector diagnostics");

    assert_eq!(trace.channel_traces.len(), blueprint.channels.len());
    assert!(trace
        .layout_segments
        .iter()
        .any(|segment| segment.is_synthetic_connector));
    assert!(trace.volume_trace.synthetic_connector_volume_mm3 > 0.0);
    assert!(trace
        .layout_segments
        .iter()
        .all(|segment| { !segment.is_synthetic_connector || segment.source_channel_id.is_none() }));
    assert!(trace
        .node_traces
        .iter()
        .all(|node| node.continuity_residual_m3_s.abs() < flow_rate * 1.0e-6));
}

#[test]
fn blueprint_trace_balances_bifurcation_serpentine_arms() {
    let blueprint = bifurcation_serpentine_rect("trace_bs", 0.010, 3, 0.008, 0.004, 0.004);
    let flow_rate = 1.0e-9;
    let config = Blueprint3dProcessingConfig {
        total_flow_rate_m3_s: flow_rate,
        ..no_chip_trace_config()
    };

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should preserve symmetric bifurcation serpentine balance");

    assert_eq!(trace.channel_traces.len(), blueprint.channels.len());
    assert_eq!(trace.volume_trace.synthetic_connector_volume_mm3, 0.0);
    assert!(trace
        .layout_segments
        .iter()
        .all(|segment| !segment.is_synthetic_connector && segment.source_channel_id.is_some()));
    assert!((reference_flow(&trace, "inlet_section") - flow_rate).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "trunk_out") - flow_rate).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "segment_1") - flow_rate * 0.5).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm2_seg_1") - flow_rate * 0.5).abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "split_jn").abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "merge_jn").abs() < flow_rate * 1.0e-6);
}

#[test]
fn blueprint_trace_balances_double_bifurcation_serpentine_arms() {
    let blueprint = double_bifurcation_serpentine_rect("trace_dbs", 0.012, 3, 0.008, 0.004, 0.004);
    let flow_rate = 1.0e-9;
    let config = Blueprint3dProcessingConfig {
        total_flow_rate_m3_s: flow_rate,
        ..no_chip_trace_config()
    };

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should preserve double bifurcation serpentine balance");

    assert_eq!(trace.channel_traces.len(), blueprint.channels.len());
    assert_eq!(trace.volume_trace.synthetic_connector_volume_mm3, 0.0);
    assert!(trace
        .layout_segments
        .iter()
        .all(|segment| !segment.is_synthetic_connector && segment.source_channel_id.is_some()));
    assert!((reference_flow(&trace, "inlet_section") - flow_rate).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "branch_a_in") - flow_rate * 0.5).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "branch_b_in") - flow_rate * 0.5).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm1_seg_1") - flow_rate * 0.25).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm2_seg_1") - flow_rate * 0.25).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm3_seg_1") - flow_rate * 0.25).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm4_seg_1") - flow_rate * 0.25).abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "split_1").abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "merge_1").abs() < flow_rate * 1.0e-6);
}

#[test]
fn blueprint_trace_balances_trifurcation_serpentine_arms() {
    let blueprint = trifurcation_serpentine_rect("trace_ts", 0.012, 3, 0.008, 0.004, 0.004);
    let flow_rate = 1.0e-9;
    let config = Blueprint3dProcessingConfig {
        total_flow_rate_m3_s: flow_rate,
        ..no_chip_trace_config()
    };

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should preserve trifurcation serpentine balance");

    assert_eq!(trace.channel_traces.len(), blueprint.channels.len());
    assert_eq!(trace.volume_trace.synthetic_connector_volume_mm3, 0.0);
    assert!(trace
        .layout_segments
        .iter()
        .all(|segment| !segment.is_synthetic_connector && segment.source_channel_id.is_some()));
    assert!((reference_flow(&trace, "inlet_section") - flow_rate).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "trunk_out") - flow_rate).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "segment_1") - flow_rate / 3.0).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm2_seg_1") - flow_rate / 3.0).abs() < flow_rate * 1.0e-6);
    assert!((reference_flow(&trace, "arm3_seg_1") - flow_rate / 3.0).abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "split_jn").abs() < flow_rate * 1.0e-6);
    assert!(continuity_residual(&trace, "merge_jn").abs() < flow_rate * 1.0e-6);
}

#[test]
fn blueprint_trace_can_attach_two_d_comparison() {
    let blueprint = venturi_chain("trace_v", 0.030, 0.004, 0.002);
    let config = Blueprint3dProcessingConfig {
        mesh: cfd_mesh::application::pipeline::PipelineConfig {
            include_chip_body: false,
            ..Default::default()
        },
        two_d_grid_nx: 12,
        two_d_grid_ny: 6,
        two_d_tolerance: 1.0e-4,
        total_flow_rate_m3_s: 1.0e-9,
        ..Default::default()
    };

    let trace = process_blueprint_with_reference_trace(&blueprint, &config)
        .expect("3D blueprint trace should attach cfd-2d comparison for venturi chain");

    assert_eq!(trace.two_d_converged_count, Some(blueprint.channels.len()));
    assert!(trace.channel_traces.iter().all(|channel| {
        channel.two_d_outlet_flow_error_pct.is_some() && channel.two_d_converged == Some(true)
    }));
}
