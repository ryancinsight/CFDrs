use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::physics::resistance::models::{
    FlowConditions, SerpentineCrossSection, SerpentineModel,
};
use cfd_core::physics::fluid::database::water_20c;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, NetworkBlueprint, NodeKind, NodeSpec,
};
use cfd_schematics::geometry::generator::{
    create_selective_tree_geometry, CenterSerpentinePathSpec, SelectiveTreeRequest,
    SelectiveTreeTopology,
};
use cfd_schematics::geometry::metadata::{
    JunctionFamily, JunctionGeometryMetadata, VenturiGeometryMetadata,
};

fn throat_edge_coefficients(bp: &NetworkBlueprint) -> (f64, f64) {
    let fluid = water_20c::<f64>().expect("water database entry must exist");
    let network = network_from_blueprint(bp, fluid).expect("network build must succeed");
    let throat = network
        .graph
        .edge_references()
        .find(|edge| edge.weight().id == "throat_section")
        .expect("throat edge must exist");
    (throat.weight().resistance, throat.weight().quad_coeff)
}

fn branch_quad_coeff(bp: &NetworkBlueprint, edge_id: &str) -> f64 {
    let fluid = water_20c::<f64>().expect("water database entry must exist");
    let network = network_from_blueprint(bp, fluid).expect("network build must succeed");
    network
        .graph
        .edge_references()
        .find(|edge| edge.weight().id == edge_id)
        .expect("branch edge must exist")
        .weight()
        .quad_coeff
}

fn edge_coefficients(bp: &NetworkBlueprint, edge_id: &str) -> (f64, f64) {
    let fluid = water_20c::<f64>().expect("water database entry must exist");
    let network = network_from_blueprint(bp, fluid).expect("network build must succeed");
    let edge = network
        .graph
        .edge_references()
        .find(|candidate| candidate.weight().id == edge_id)
        .expect("edge must exist");
    (edge.weight().resistance, edge.weight().quad_coeff)
}

fn blueprint_channel<'a>(bp: &'a NetworkBlueprint, edge_id: &str) -> &'a ChannelSpec {
    bp.channels
        .iter()
        .find(|channel| channel.id.as_str() == edge_id)
        .expect("channel must exist")
}

fn selective_serpentine_blueprint(spec: CenterSerpentinePathSpec) -> NetworkBlueprint {
    let request = SelectiveTreeRequest {
        name: format!("selective-serp-{}", spec.segments),
        box_dims_mm: (127.76, 85.47),
        trunk_length_m: 12.0e-3,
        branch_length_m: 10.0e-3,
        hybrid_branch_length_m: 8.0e-3,
        main_width_m: 1.2e-3,
        throat_width_m: 0.4e-3,
        throat_length_m: 3.0e-3,
        channel_height_m: 0.5e-3,
        topology: SelectiveTreeTopology::CascadeCenterTrifurcation {
            n_levels: 3,
            center_frac: 0.45,
            venturi_treatment_enabled: false,
            center_serpentine: Some(spec),
        },
    };
    create_selective_tree_geometry(&request)
}

fn serpentine_analysis(channel: &ChannelSpec) -> (usize, f64, f64, f64) {
    let (segments, bend_radius_m) = match channel.channel_shape {
        ChannelShape::Serpentine {
            segments,
            bend_radius_m,
        } => (segments, bend_radius_m),
        _ => panic!("channel must carry inferred serpentine metadata"),
    };

    let cross_section = match channel.cross_section {
        CrossSectionSpec::Circular { diameter_m } => SerpentineCrossSection::Circular {
            diameter: diameter_m,
        },
        CrossSectionSpec::Rectangular { width_m, height_m } => {
            SerpentineCrossSection::Rectangular {
                width: width_m,
                height: height_m,
            }
        }
    };

    let fluid = water_20c::<f64>().expect("water database entry must exist");
    let model = SerpentineModel::new(channel.length_m, segments, cross_section, bend_radius_m);
    let analysis = model
        .analyze(&fluid, &FlowConditions::new(0.25))
        .expect("serpentine analysis must succeed");
    (
        segments,
        bend_radius_m,
        analysis.dean_number,
        analysis.dp_bends,
    )
}

fn venturi_blueprint(
    convergent_half_angle_deg: f64,
    divergent_half_angle_deg: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new("venturi-angle-test");
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "inlet_section",
        "inlet",
        "contraction",
        5.0e-3,
        1.0e-3,
        0.5e-3,
        1.0,
        0.0,
    ));
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "throat_section",
            "contraction",
            "throat",
            8.0e-4,
            2.5e-4,
            0.5e-3,
            1.0,
            0.0,
        )
        .with_venturi_geometry(VenturiGeometryMetadata {
            throat_width_m: 2.5e-4,
            throat_height_m: 0.5e-3,
            throat_length_m: 8.0e-4,
            inlet_width_m: 1.0e-3,
            outlet_width_m: 1.0e-3,
            convergent_half_angle_deg,
            divergent_half_angle_deg,
            throat_position: 0.5,
        })
        .with_metadata(VenturiGeometryMetadata {
            throat_width_m: 2.5e-4,
            throat_height_m: 0.5e-3,
            throat_length_m: 8.0e-4,
            inlet_width_m: 1.0e-3,
            outlet_width_m: 1.0e-3,
            convergent_half_angle_deg,
            divergent_half_angle_deg,
            throat_position: 0.5,
        }),
    );
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "diffuser_section",
        "throat",
        "outlet",
        5.0e-3,
        1.0e-3,
        0.5e-3,
        1.0,
        0.0,
    ));
    bp
}

fn split_blueprint(branch_angles_deg: Vec<f64>) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new("junction-angle-test");
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(
        NodeSpec::new("split", NodeKind::Junction)
            .with_junction_geometry(JunctionGeometryMetadata {
                junction_family: JunctionFamily::Trifurcation,
                branch_angles_deg: branch_angles_deg.clone(),
                merge_angles_deg: Vec::new(),
            })
            .with_metadata(JunctionGeometryMetadata {
                junction_family: JunctionFamily::Trifurcation,
                branch_angles_deg,
                merge_angles_deg: Vec::new(),
            }),
    );
    bp.add_node(NodeSpec::new("run_out", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("branch_up", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("branch_dn", NodeKind::Outlet));
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "inlet_section",
        "inlet",
        "split",
        6.0e-3,
        1.2e-3,
        0.6e-3,
        1.0,
        0.0,
    ));
    for (edge_id, node_id) in [
        ("run_edge", "run_out"),
        ("branch_up_edge", "branch_up"),
        ("branch_dn_edge", "branch_dn"),
    ] {
        bp.add_channel(ChannelSpec::new_pipe_rect(
            edge_id, "split", node_id, 6.0e-3, 0.8e-3, 0.6e-3, 1.0, 0.0,
        ));
    }
    bp
}

fn bifurcation_flow_blueprint(branch_width_m: f64) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new("branch-width-test");
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(
        NodeSpec::new("split", NodeKind::Junction).with_junction_geometry(
            JunctionGeometryMetadata {
                junction_family: JunctionFamily::Bifurcation,
                branch_angles_deg: vec![-35.0, 35.0],
                merge_angles_deg: Vec::new(),
            },
        ),
    );
    bp.add_node(NodeSpec::new("wide_out", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("narrow_out", NodeKind::Outlet));
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "inlet_section",
        "inlet",
        "split",
        8.0e-3,
        1.2e-3,
        0.6e-3,
        1.0,
        0.0,
    ));
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "treatment_edge",
        "split",
        "wide_out",
        9.0e-3,
        branch_width_m,
        0.6e-3,
        1.0,
        0.0,
    ));
    bp.add_channel(ChannelSpec::new_pipe_rect(
        "bypass_edge",
        "split",
        "narrow_out",
        9.0e-3,
        0.45e-3,
        0.6e-3,
        1.0,
        0.0,
    ));
    bp
}

fn predicted_branch_split(bp: &NetworkBlueprint, inlet_flow_m3_s: f64) -> (f64, f64, f64, f64) {
    let fluid = water_20c::<f64>().expect("water database entry must exist");
    let network = network_from_blueprint(bp, fluid).expect("network build must succeed");
    let split = network
        .graph
        .node_indices()
        .find(|idx| {
            network
                .graph
                .node_weight(*idx)
                .is_some_and(|node| node.id == "split")
        })
        .expect("split node must exist");
    let mut treatment = None::<(f64, f64)>;
    let mut bypass = None::<(f64, f64)>;
    for edge in network.graph.edges(split) {
        match edge.weight().id.as_str() {
            "treatment_edge" => {
                let resistance = edge.weight().resistance.max(1.0e-12);
                treatment = Some((resistance, 1.0 / resistance));
            }
            "bypass_edge" => {
                let resistance = edge.weight().resistance.max(1.0e-12);
                bypass = Some((resistance, 1.0 / resistance));
            }
            _ => {}
        }
    }
    let (r_treatment, g_treatment) = treatment.expect("treatment edge must exist");
    let (r_bypass, g_bypass) = bypass.expect("bypass edge must exist");
    let conductance_sum = g_treatment + g_bypass;
    let q_treatment = inlet_flow_m3_s * g_treatment / conductance_sum;
    let q_bypass = inlet_flow_m3_s * g_bypass / conductance_sum;
    (r_treatment, q_treatment, r_bypass, q_bypass)
}

#[test]
fn venturi_angles_change_resistance_coefficients() {
    let shallow = venturi_blueprint(6.0, 7.0);
    let steep = venturi_blueprint(22.0, 24.0);

    let (r_shallow, k_shallow) = throat_edge_coefficients(&shallow);
    let (r_steep, k_steep) = throat_edge_coefficients(&steep);

    assert!(
        (r_shallow - r_steep).abs() > 1e-9 || (k_shallow - k_steep).abs() > 1e-9,
        "venturi angle metadata must perturb 1D coefficients"
    );
}

#[test]
fn junction_angles_change_branch_minor_loss() {
    let narrow = split_blueprint(vec![-15.0, 0.0, 15.0]);
    let wide = split_blueprint(vec![-60.0, 0.0, 60.0]);

    let k_narrow = branch_quad_coeff(&narrow, "branch_up_edge");
    let k_wide = branch_quad_coeff(&wide, "branch_up_edge");

    assert!(
        k_wide > k_narrow,
        "wider branch angles should increase branch minor-loss coefficient"
    );
}

#[test]
fn branch_width_changes_resistance_and_flow_split() {
    let narrow = bifurcation_flow_blueprint(0.55e-3);
    let wide = bifurcation_flow_blueprint(0.95e-3);

    let (r_treat_narrow, q_treat_narrow, _r_bypass_narrow, q_bypass_narrow) =
        predicted_branch_split(&narrow, 1.0e-9);
    let (r_treat_wide, q_treat_wide, _r_bypass_wide, q_bypass_wide) =
        predicted_branch_split(&wide, 1.0e-9);

    assert!(
        r_treat_wide < r_treat_narrow,
        "wider treatment branch should reduce hydraulic resistance"
    );
    assert!(
        q_treat_wide > q_treat_narrow,
        "wider treatment branch should attract more predicted branch flow"
    );
    assert!(
        q_bypass_wide < q_bypass_narrow,
        "widening one branch should reduce the predicted bypass flow split"
    );
}

#[test]
fn venturi_throat_width_and_length_change_coefficients() {
    let short_wide = venturi_blueprint(8.0, 8.0);
    let mut long_narrow = venturi_blueprint(8.0, 8.0);
    if let Some(throat) = long_narrow
        .channels
        .iter_mut()
        .find(|channel| channel.id.as_str() == "throat_section")
    {
        throat.length_m = 1.6e-3;
        throat.cross_section = cfd_schematics::domain::model::CrossSectionSpec::Rectangular {
            width_m: 1.8e-4,
            height_m: 0.5e-3,
        };
        throat.venturi_geometry = Some(VenturiGeometryMetadata {
            throat_width_m: 1.8e-4,
            throat_height_m: 0.5e-3,
            throat_length_m: 1.6e-3,
            inlet_width_m: 1.0e-3,
            outlet_width_m: 1.0e-3,
            convergent_half_angle_deg: 8.0,
            divergent_half_angle_deg: 8.0,
            throat_position: 0.5,
        });
    }

    let (r_short_wide, k_short_wide) = throat_edge_coefficients(&short_wide);
    let (r_long_narrow, k_long_narrow) = throat_edge_coefficients(&long_narrow);

    assert!(
        r_long_narrow > r_short_wide,
        "longer, narrower throats must increase linear throat resistance"
    );
    assert!(
        (r_long_narrow - r_short_wide).abs() > 1.0e-9
            || (k_long_narrow - k_short_wide).abs() > 1.0e-9,
        "changing throat geometry must perturb at least one throat-loss coefficient"
    );
}

#[test]
fn inferred_selective_serpentine_metadata_changes_1d_losses() {
    let tight = selective_serpentine_blueprint(CenterSerpentinePathSpec {
        segments: 9,
        bend_radius_m: 0.45e-3,
    });
    let gentle = selective_serpentine_blueprint(CenterSerpentinePathSpec {
        segments: 5,
        bend_radius_m: 2.2e-3,
    });

    let tight_channel = blueprint_channel(&tight, "center_lv0");
    let gentle_channel = blueprint_channel(&gentle, "center_lv0");
    let (tight_segments, tight_radius_m, tight_dean, tight_dp_bends) =
        serpentine_analysis(tight_channel);
    let (gentle_segments, gentle_radius_m, gentle_dean, gentle_dp_bends) =
        serpentine_analysis(gentle_channel);
    let (tight_r, tight_k) = edge_coefficients(&tight, "center_lv0");
    let (gentle_r, gentle_k) = edge_coefficients(&gentle, "center_lv0");

    assert!(
        tight_segments > gentle_segments,
        "tighter selective routing should infer more serpentine segments"
    );
    assert!(
        tight_radius_m < gentle_radius_m,
        "tighter selective routing should infer a smaller bend radius"
    );
    assert!(
        tight_dean > gentle_dean,
        "smaller inferred bend radii should raise Dean number on the representative selective lane"
    );
    assert!(
        tight_dp_bends > gentle_dp_bends,
        "more turns with tighter curvature should increase bend-loss pressure drop"
    );
    assert!(
        (tight_r - gentle_r).abs() > 1.0e-12 || (tight_k - gentle_k).abs() > 1.0e-12,
        "network_from_blueprint must remain sensitive to inferred selective serpentine metadata"
    );
}
