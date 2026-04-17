use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId,
    NodeKind, NodeSpec,
};
use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
use cfd_schematics::geometry::metadata::{
    BranchBoundaryMetadata, JunctionFamily, JunctionGeometryMetadata,
};
use cfd_schematics::interface::presets::{
    bifurcation_rect, bifurcation_venturi_rect, primitive_selective_split_tree_rect,
    serpentine_rect, venturi_rect,
};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};

use super::postprocess::extract_field_wall_shear;
use super::Network2dBuilderSink;

fn circular_blueprint() -> NetworkBlueprint {
    let mut blueprint = NetworkBlueprint::new_with_explicit_positions("circular_trace");
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
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(
            cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
        );
    blueprint
}

fn rectangular_blueprint() -> NetworkBlueprint {
    let mut blueprint = NetworkBlueprint::new_with_explicit_positions("rectangular_trace");
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));
    let mut duct =
        ChannelSpec::new_pipe_rect("duct", "inlet", "outlet", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    duct.path = vec![(0.0, 0.0), (1.0, 0.0)];
    blueprint.add_channel(duct);
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(
            cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
        );
    blueprint
}

fn annotate_with_provenance(mut blueprint: NetworkBlueprint) -> NetworkBlueprint {
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(
            cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
        );
    blueprint
}

fn projected_network(
    blueprint: &NetworkBlueprint,
    grid_nx: usize,
    grid_ny: usize,
) -> super::types::Network2DSolver<f64> {
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-6, grid_nx, grid_ny);
    sink.build(blueprint).expect("build should succeed")
}

fn assert_projection_summary(network: &super::types::Network2DSolver<f64>) {
    let projection = network.projection_summary_ref();
    assert_eq!(projection.channel_count(), network.channels.len());
    assert!(projection.total_fluid_cell_count > 0);
    assert!(projection.mean_fluid_fraction > 0.0);
    assert!(projection.mean_fluid_fraction <= 1.0);

    for (entry, summary) in network
        .channels
        .iter()
        .zip(projection.channel_summaries.iter())
    {
        assert_eq!(entry.projection.channel_id, summary.channel_id);
        assert_eq!(entry.projection.fluid_cell_count, summary.fluid_cell_count);
        assert!(
            summary.fluid_cell_count > 0,
            "{} should project fluid cells",
            summary.channel_id
        );
        assert!(summary.fluid_fraction > 0.0);
        assert!(summary.fluid_fraction <= 1.0);
        assert!(summary.grid_length_m > 0.0);
        assert!(summary.grid_width_m > 0.0);
        assert!(summary.path_length_m > 0.0);
    }
}

#[test]
fn straight_projection_is_rasterized_into_fluid_cells() {
    let blueprint = annotate_with_provenance(rectangular_blueprint());
    let network = projected_network(&blueprint, 32, 12);
    assert_projection_summary(&network);

    let channel = &network.channels[0];
    assert_eq!(channel.projection.path_span_y_m, 0.0);
    assert!(channel.projection.path_span_x_m > 0.0);
}

#[test]
fn serpentine_projection_preserves_path_bend_span() {
    let blueprint = annotate_with_provenance(serpentine_rect("serp", 4, 8.0e-3, 2.0e-3, 1.0e-3));
    let network = projected_network(&blueprint, 32, 12);
    assert_projection_summary(&network);

    assert!(
        network
            .channels
            .iter()
            .any(|entry| entry.projection.path_span_y_m > 0.0),
        "serpentine projection should preserve a non-zero y-span"
    );
}

#[test]
fn venturi_projection_marks_throat_channels() {
    let blueprint = annotate_with_provenance(venturi_rect("vent", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3));
    let network = projected_network(&blueprint, 32, 12);
    assert_projection_summary(&network);

    assert!(
        network.channels.iter().any(|entry| entry.is_venturi_throat),
        "venturi preset should mark at least one throat channel"
    );
}

#[test]
fn bifurcation_projection_covers_multiple_channels() {
    let blueprint = annotate_with_provenance(bifurcation_rect(
        "bif", 4.0e-3, 2.0e-3, 1.0e-3, 0.7e-3, 1.0e-3,
    ));
    let network = projected_network(&blueprint, 32, 12);
    assert_projection_summary(&network);

    assert!(
        network.channels.len() > 1,
        "bifurcation blueprint should produce multiple 2D channel solves"
    );
}

#[test]
fn selective_tree_projection_is_available_for_primitives() {
    let blueprint = annotate_with_provenance(primitive_selective_split_tree_rect(
        "pst",
        (127.76, 85.47),
        &[PrimitiveSelectiveSplitKind::Tri],
        2.0e-3,
        0.45,
        0.45,
        0.68,
        0.5e-3,
        1.0e-3,
        1.0e-3,
        false,
        0,
        None,
    ));
    let network = projected_network(&blueprint, 24, 10);
    assert_projection_summary(&network);

    assert!(
        network.channels.len() >= 3,
        "selective-tree blueprint should produce multiple projected channels"
    );
}

#[test]
fn projected_solve_returns_projection_metadata() {
    let blueprint = annotate_with_provenance(rectangular_blueprint());
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-7, 24, 10);
    let mut network = sink.build(&blueprint).expect("build should succeed");

    let projected = network
        .solve_projected(1e-6)
        .expect("projected solve should succeed");

    assert_eq!(
        projected.result.channels.len(),
        projected.projection.channel_count()
    );
    assert!(projected.projection.total_fluid_cell_count > 0);
    assert!(projected.result.channels[0].solve_result.converged);
    assert_eq!(
        projected.result.channels[0].projection.channel_id,
        projected.projection.channel_summaries[0].channel_id
    );
    assert_eq!(
        projected.result.channels[0].projection.fluid_cell_count,
        projected.projection.channel_summaries[0].fluid_cell_count
    );
    assert!(projected.result.channels[0].field_effective_resistance_pa_s_per_m3 > 0.0);
    assert!(projected.result.channels[0].field_pressure_drop_pa > 0.0);
}

#[test]
fn coupled_multi_channel_solve_converges_on_junction_network() {
    let mut blueprint = bifurcation_venturi_rect("coupled_bv", 5e-3, 2e-3, 5e-4, 1e-3, 1e-3);
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(
            cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
        );

    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let mut network = sink.build(&blueprint).expect("build should succeed");

    let coupled = network
        .solve_coupled(1e-4)
        .expect("coupled solve should succeed");

    assert!(coupled.coupling_iterations > 0);
    assert_eq!(
        coupled.result.channels.len(),
        coupled.projection.channel_count()
    );
    assert!(coupled
        .result
        .channels
        .iter()
        .all(|channel| channel.solve_result.converged));
    assert!(coupled
        .result
        .channels
        .iter()
        .all(|channel| channel.field_effective_resistance_pa_s_per_m3 > 0.0));
    assert!(
        (coupled.result.reference_trace.total_inlet_flow_m3_s - q_total).abs() < q_total * 1e-6
    );
}

#[test]
fn coupled_solve_honors_branch_boundary_metadata() {
    let mut blueprint = annotate_with_provenance(NetworkBlueprint::new_with_explicit_positions(
        "branch_boundary_metadata",
    ));
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(NodeSpec::new_at("split", NodeKind::Junction, (0.5, 0.0)));
    blueprint.add_node(
        NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0))
            .with_metadata(BranchBoundaryMetadata::pressure(0.25)),
    );

    let mut upstream =
        ChannelSpec::new_pipe_rect("upstream", "inlet", "split", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "split",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    blueprint.add_channel(downstream);

    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let mut network = sink.build(&blueprint).expect("build should succeed");

    let coupled = network
        .solve_coupled(1e-4)
        .expect("coupled solve should succeed");

    assert!(coupled.coupling_iterations > 0);

    let scale = coupled.result.reference_trace.pressure_scale_factor;
    let inlet = coupled
        .result
        .reference_trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "inlet")
        .expect("inlet node trace must exist");
    let split = coupled
        .result
        .reference_trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "split")
        .expect("split node trace must exist");
    let outlet = coupled
        .result
        .reference_trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "outlet")
        .expect("outlet node trace must exist");

    assert!((inlet.pressure_pa - scale).abs() < scale * 1e-6);
    assert!((outlet.pressure_pa - scale * 0.25).abs() < scale * 1e-6);
    assert!(split.continuity_residual_m3_s.abs() < q_total * 1e-6);
}

#[test]
fn coupled_solve_honors_negative_branch_flow_sink_metadata() {
    let mut blueprint = annotate_with_provenance(NetworkBlueprint::new_with_explicit_positions(
        "branch_flow_sink_metadata_coupled",
    ));
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(
        NodeSpec::new_at("split", NodeKind::Junction, (0.5, 0.0))
            .with_metadata(BranchBoundaryMetadata::flow_rate(-2.0e-7)),
    );
    blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut upstream =
        ChannelSpec::new_pipe_rect("upstream", "inlet", "split", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "split",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    blueprint.add_channel(downstream);

    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let mut network = sink.build(&blueprint).expect("build should succeed");

    let coupled = network
        .solve_coupled(1e-4)
        .expect("coupled solve should succeed");

    assert!(coupled.coupling_iterations > 0);
    assert_eq!(coupled.result.channels.len(), 2);

    let split = coupled
        .result
        .reference_trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "split")
        .expect("split node trace must exist");

    assert!(split.prescribed_boundary_flow_m3_s < 0.0);
    assert!(split.continuity_residual_m3_s.abs() < q_total * 1e-6);
}

#[test]
fn reference_trace_records_branch_flow_metadata() {
    let mut blueprint = annotate_with_provenance(NetworkBlueprint::new_with_explicit_positions(
        "branch_flow_metadata",
    ));
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(
        NodeSpec::new_at("split", NodeKind::Junction, (0.5, 0.0))
            .with_metadata(BranchBoundaryMetadata::flow_rate(2.0e-7)),
    );
    blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut upstream =
        ChannelSpec::new_pipe_rect("upstream", "inlet", "split", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "split",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    blueprint.add_channel(downstream);

    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let network = sink.build(&blueprint).expect("build should succeed");
    let trace = network.reference_trace();
    let scale = trace.pressure_scale_factor;

    let split = trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "split")
        .expect("split node trace must exist");

    assert!((split.prescribed_boundary_flow_m3_s - scale * 2.0e-7).abs() < scale * 1e-6);
    assert!(split.continuity_residual_m3_s.abs() < q_total * 1e-6);
}

#[test]
fn reference_trace_records_negative_branch_flow_metadata() {
    let mut blueprint = annotate_with_provenance(NetworkBlueprint::new_with_explicit_positions(
        "branch_flow_sink_metadata",
    ));
    blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    blueprint.add_node(
        NodeSpec::new_at("split", NodeKind::Junction, (0.5, 0.0))
            .with_metadata(BranchBoundaryMetadata::flow_rate(-2.0e-7)),
    );
    blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut upstream =
        ChannelSpec::new_pipe_rect("upstream", "inlet", "split", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "split",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    blueprint.add_channel(downstream);

    let q_total = 1e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 20, 8);
    let network = sink.build(&blueprint).expect("build should succeed");
    let trace = network.reference_trace();
    let scale = trace.pressure_scale_factor;

    let split = trace
        .node_traces
        .iter()
        .find(|node| node.node_id == "split")
        .expect("split node trace must exist");

    assert!((split.prescribed_boundary_flow_m3_s + scale * 2.0e-7).abs() < scale * 1e-6);
    assert!(split.prescribed_boundary_flow_m3_s < 0.0);
    assert!(split.continuity_residual_m3_s.abs() < q_total * 1e-6);
}

#[test]
fn branch_metadata_biases_split_and_merge_coupling_weights() {
    let mut split_blueprint = NetworkBlueprint::new_with_explicit_positions("split_weighting");
    split_blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    split_blueprint.add_node(
        NodeSpec::new_at("junction", NodeKind::Junction, (0.5, 0.0)).with_metadata(
            JunctionGeometryMetadata {
                junction_family: JunctionFamily::Bifurcation,
                branch_angles_deg: vec![45.0, 60.0],
                merge_angles_deg: vec![10.0],
            },
        ),
    );
    split_blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut upstream = ChannelSpec::new_pipe_rect(
        "upstream", "inlet", "junction", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    split_blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "junction",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    split_blueprint.add_channel(downstream);

    let split_weights = super::coupled::build_channel_coupling_weights::<f64>(&split_blueprint);
    assert_eq!(split_weights.len(), 2);
    assert!(split_weights[1] > split_weights[0]);
    assert!(split_weights[1] > 1.0);
    assert!(split_weights[0] < 1.0);

    let mut split_boundary_blueprint =
        NetworkBlueprint::new_with_explicit_positions("split_weighting_boundary");
    split_boundary_blueprint.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    split_boundary_blueprint.add_node(
        NodeSpec::new_at("junction", NodeKind::Junction, (0.5, 0.0))
            .with_metadata(JunctionGeometryMetadata {
                junction_family: JunctionFamily::Bifurcation,
                branch_angles_deg: vec![45.0, 60.0],
                merge_angles_deg: vec![10.0],
            })
            .with_metadata(BranchBoundaryMetadata::flow_rate(2.0e-7)),
    );
    split_boundary_blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut upstream = ChannelSpec::new_pipe_rect(
        "upstream", "inlet", "junction", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    upstream.path = vec![(0.0, 0.0), (0.5, 0.0)];
    split_boundary_blueprint.add_channel(upstream);

    let mut downstream = ChannelSpec::new_pipe_rect(
        "downstream",
        "junction",
        "outlet",
        0.01,
        1.0e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    downstream.path = vec![(0.5, 0.0), (1.0, 0.0)];
    split_boundary_blueprint.add_channel(downstream);

    let split_boundary_weights =
        super::coupled::build_channel_coupling_weights::<f64>(&split_boundary_blueprint);
    assert!(split_boundary_weights[1] > split_weights[1]);
    assert!(split_boundary_weights[0] < split_weights[0]);

    let mut merge_blueprint = NetworkBlueprint::new_with_explicit_positions("merge_weighting");
    merge_blueprint.add_node(NodeSpec::new_at("inlet_a", NodeKind::Inlet, (0.0, 0.2)));
    merge_blueprint.add_node(NodeSpec::new_at("inlet_b", NodeKind::Inlet, (0.0, -0.2)));
    merge_blueprint.add_node(
        NodeSpec::new_at("merge", NodeKind::Junction, (0.5, 0.0)).with_metadata(
            JunctionGeometryMetadata {
                junction_family: JunctionFamily::Merge,
                branch_angles_deg: vec![25.0],
                merge_angles_deg: vec![70.0],
            },
        ),
    );
    merge_blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut inlet_a = ChannelSpec::new_pipe_rect(
        "inlet_a", "inlet_a", "merge", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    inlet_a.path = vec![(0.0, 0.2), (0.5, 0.0)];
    merge_blueprint.add_channel(inlet_a);

    let mut inlet_b = ChannelSpec::new_pipe_rect(
        "inlet_b", "inlet_b", "merge", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    inlet_b.path = vec![(0.0, -0.2), (0.5, 0.0)];
    merge_blueprint.add_channel(inlet_b);

    let mut outlet =
        ChannelSpec::new_pipe_rect("outlet", "merge", "outlet", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    outlet.path = vec![(0.5, 0.0), (1.0, 0.0)];
    merge_blueprint.add_channel(outlet);

    let merge_weights = super::coupled::build_channel_coupling_weights::<f64>(&merge_blueprint);
    assert_eq!(merge_weights.len(), 3);
    assert!(merge_weights[0] > merge_weights[2]);
    assert!(merge_weights[1] > merge_weights[2]);
    assert!(merge_weights[2] < 1.0);

    let mut merge_boundary_blueprint =
        NetworkBlueprint::new_with_explicit_positions("merge_weighting_boundary");
    merge_boundary_blueprint.add_node(NodeSpec::new_at("inlet_a", NodeKind::Inlet, (0.0, 0.2)));
    merge_boundary_blueprint.add_node(NodeSpec::new_at("inlet_b", NodeKind::Inlet, (0.0, -0.2)));
    merge_boundary_blueprint.add_node(
        NodeSpec::new_at("merge", NodeKind::Junction, (0.5, 0.0))
            .with_metadata(JunctionGeometryMetadata {
                junction_family: JunctionFamily::Merge,
                branch_angles_deg: vec![25.0],
                merge_angles_deg: vec![70.0],
            })
            .with_metadata(BranchBoundaryMetadata::flow_rate(-2.0e-7)),
    );
    merge_boundary_blueprint.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));

    let mut inlet_a = ChannelSpec::new_pipe_rect(
        "inlet_a", "inlet_a", "merge", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    inlet_a.path = vec![(0.0, 0.2), (0.5, 0.0)];
    merge_boundary_blueprint.add_channel(inlet_a);

    let mut inlet_b = ChannelSpec::new_pipe_rect(
        "inlet_b", "inlet_b", "merge", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0,
    );
    inlet_b.path = vec![(0.0, -0.2), (0.5, 0.0)];
    merge_boundary_blueprint.add_channel(inlet_b);

    let mut outlet =
        ChannelSpec::new_pipe_rect("outlet", "merge", "outlet", 0.01, 1.0e-3, 1.0e-3, 0.0, 0.0);
    outlet.path = vec![(0.5, 0.0), (1.0, 0.0)];
    merge_boundary_blueprint.add_channel(outlet);

    let merge_boundary_weights =
        super::coupled::build_channel_coupling_weights::<f64>(&merge_boundary_blueprint);
    assert!(merge_boundary_weights[0] > merge_weights[0]);
    assert!(merge_boundary_weights[1] > merge_weights[1]);
    assert!(merge_boundary_weights[2] < merge_weights[2]);
}

#[test]
fn separation_tracking_is_disabled_by_default() {
    let blueprint = annotate_with_provenance(rectangular_blueprint());
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-7, 24, 10);
    let mut network = sink.build(&blueprint).expect("build should succeed");
    let result = network.solve_all(1e-6).expect("solve should succeed");

    assert!(result
        .channels
        .iter()
        .all(|channel| channel.field_separation_efficiency_pct.is_none()));
}

#[test]
fn build_single_venturi_network() {
    let mut bp = venturi_rect("test_v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
    bp.metadata.get_or_insert_with(Default::default).insert(
        cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
    );
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let sink = Network2dBuilderSink::new(blood, 1060.0, 1e-6, 10, 5);
    let net2d = sink.build(&bp).expect("build should succeed");
    assert_eq!(net2d.reference_trace().channel_traces.len(), 3);
}

#[test]
fn reference_trace_sums_to_q_total() {
    let mut bp = bifurcation_venturi_rect("bv", 0.005, 0.002, 0.0005, 0.001, 0.001);
    bp.metadata.get_or_insert_with(Default::default).insert(
        cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
    );
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
    bp.metadata.get_or_insert_with(Default::default).insert(
        cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
    );
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
    bp.metadata.get_or_insert_with(Default::default).insert(
        cfd_schematics::geometry::metadata::GeometryAuthoringProvenance::selective_wrapper(),
    );
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
