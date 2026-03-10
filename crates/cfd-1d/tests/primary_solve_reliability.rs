use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::{NetworkProblem, NetworkSolver, SolveFailureReason, SolverConfig};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
use cfd_schematics::geometry::metadata::VenturiGeometryMetadata;
use cfd_schematics::interface::presets::primitive_selective_split_tree_rect;

#[derive(Debug, Clone)]
struct PrimitiveSolveCase {
    sequence: Vec<PrimitiveSelectiveSplitKind>,
    first_trifurcation_center_frac: f64,
    later_trifurcation_center_frac: f64,
    bifurcation_treatment_frac: f64,
    throat_width_m: f64,
    venturi_treatment_enabled: bool,
    treatment_branch_throat_count: u8,
    inlet_flow_m3_s: f64,
}

fn primitive_tree(
    sequence: &[PrimitiveSelectiveSplitKind],
    first_trifurcation_center_frac: f64,
    later_trifurcation_center_frac: f64,
    bifurcation_treatment_frac: f64,
    throat_width_m: f64,
    venturi_treatment_enabled: bool,
    treatment_branch_throat_count: u8,
) -> cfd_schematics::NetworkBlueprint {
    primitive_selective_split_tree_rect(
        "primary-solve-reliability",
        (127.76, 85.47),
        sequence,
        4.0e-3,
        first_trifurcation_center_frac,
        later_trifurcation_center_frac,
        bifurcation_treatment_frac,
        throat_width_m,
        120.0e-6,
        1.0e-3,
        venturi_treatment_enabled,
        treatment_branch_throat_count,
        None,
    )
}

fn solve_blueprint(
    blueprint: &cfd_schematics::NetworkBlueprint,
    inlet_flow_m3_s: f64,
    max_iterations: usize,
) -> Result<
    (
        cfd_1d::Network<f64, CassonBlood<f64>>,
        cfd_1d::PrimarySolveDiagnostics,
    ),
    cfd_1d::PrimarySolveError,
> {
    let blood = CassonBlood::<f64>::normal_blood();
    let mut network = network_from_blueprint(blueprint, blood).expect("network build must succeed");
    let inlet_nodes: Vec<_> = network
        .graph
        .node_indices()
        .filter(|idx| {
            network
                .graph
                .node_weight(*idx)
                .is_some_and(|node| node.node_type == cfd_1d::NodeType::Inlet)
        })
        .collect();
    let outlet_nodes: Vec<_> = network
        .graph
        .node_indices()
        .filter(|idx| {
            network
                .graph
                .node_weight(*idx)
                .is_some_and(|node| node.node_type == cfd_1d::NodeType::Outlet)
        })
        .collect();
    for inlet in &inlet_nodes {
        network.set_neumann_flow(*inlet, inlet_flow_m3_s / inlet_nodes.len() as f64);
    }
    for outlet in &outlet_nodes {
        network.set_pressure(*outlet, 0.0);
    }

    let solver = NetworkSolver::<f64, CassonBlood<f64>>::with_config(SolverConfig {
        tolerance: 1.0e-8,
        max_iterations,
    });
    solver.solve_network_with_diagnostics(&NetworkProblem::new(network))
}

fn representative_case_search_with_failure_detail(
    sequence: &[PrimitiveSelectiveSplitKind],
) -> (
    Option<(
        PrimitiveSolveCase,
        cfd_1d::Network<f64, CassonBlood<f64>>,
        cfd_1d::PrimarySolveDiagnostics,
    )>,
    String,
) {
    let first_tri_fracs = [0.33, 0.40, 0.45];
    let later_tri_fracs = [0.25, 0.33, 0.50];
    let bi_treat_fracs = [0.55, 0.68];
    let throat_widths_m = [35.0e-6, 45.0e-6, 120.0e-6];
    let inlet_flows_m3_s = [80.0 / 6.0e7, 150.0 / 6.0e7, 300.0 / 6.0e7];
    let mut last_failure = String::from("no attempts executed");

    for &venturi_treatment_enabled in &[false, true] {
        for &throat_width_m in &throat_widths_m {
            for &first_trifurcation_center_frac in &first_tri_fracs {
                for &later_trifurcation_center_frac in &later_tri_fracs {
                    for &bifurcation_treatment_frac in &bi_treat_fracs {
                        for &inlet_flow_m3_s in &inlet_flows_m3_s {
                            let case = PrimitiveSolveCase {
                                sequence: sequence.to_vec(),
                                first_trifurcation_center_frac,
                                later_trifurcation_center_frac,
                                bifurcation_treatment_frac,
                                throat_width_m,
                                venturi_treatment_enabled,
                                treatment_branch_throat_count: 1,
                                inlet_flow_m3_s,
                            };
                            let blueprint = primitive_tree(
                                &case.sequence,
                                case.first_trifurcation_center_frac,
                                case.later_trifurcation_center_frac,
                                case.bifurcation_treatment_frac,
                                case.throat_width_m,
                                case.venturi_treatment_enabled,
                                case.treatment_branch_throat_count,
                            );
                            match solve_blueprint(&blueprint, case.inlet_flow_m3_s, 1000) {
                                Ok((network, diagnostics)) => {
                                    return (Some((case, network, diagnostics)), last_failure);
                                }
                                Err(err) => {
                                    last_failure = format!(
                                        "reason={:?}, detail={}, residual={:?}, delta={:?}, iterations={}",
                                        err.reason,
                                        err.diagnostics
                                            .failure_detail
                                            .as_deref()
                                            .unwrap_or("n/a"),
                                        err.diagnostics.last_residual_norm,
                                        err.diagnostics.last_solution_change_norm,
                                        err.diagnostics.picard_iterations,
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    (None, last_failure)
}

#[test]
fn representative_primitive_selective_trees_primary_converge() {
    let sequences = [
        vec![
            PrimitiveSelectiveSplitKind::Bi,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Bi,
        ],
        vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
    ];

    for sequence in sequences {
        let (found, failure_detail) = representative_case_search_with_failure_detail(&sequence);
        let Some((case, network, diagnostics)) = found else {
            panic!(
                "failed to find a primary-converged representative case for primitive sequence {:?}; last failure: {}",
                sequence,
                failure_detail,
            );
        };

        assert!(
            diagnostics.picard_iterations >= 1,
            "primary solve must record iteration count for {:?}",
            case.sequence
        );
        assert!(
            diagnostics
                .last_residual_norm
                .is_some_and(|residual| residual.is_finite() && residual >= 0.0),
            "primary solve must report a finite residual for {:?}",
            case.sequence
        );
        assert!(
            !diagnostics.degraded_geometry_for_recovery,
            "trusted primary solves must not degrade geometry for {:?}",
            case.sequence
        );

        for edge in network.graph.edge_weights() {
            assert!(
                edge.resistance.is_finite() && edge.resistance > 0.0,
                "effective linear resistance must stay positive and finite for edge '{}' in {:?}",
                edge.id,
                case.sequence
            );
            assert!(
                edge.quad_coeff.is_finite() && edge.quad_coeff >= 0.0,
                "effective quadratic coefficient must stay finite and non-negative for edge '{}' in {:?}",
                edge.id,
                case.sequence
            );
        }
    }
}

#[test]
fn aggressive_primitive_selective_case_primary_converges_without_recovery() {
    let blueprint = primitive_tree(
        &[
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        0.45,
        0.55,
        0.68,
        30.0e-6,
        true,
        4,
    );
    let (network, diagnostics) = solve_blueprint(&blueprint, 500.0 / 6.0e7, 1000).expect(
        "aggressive but valid primitive selective trees should stay on the trusted primary path",
    );

    assert!(diagnostics.picard_iterations >= 1);
    assert!(
        diagnostics.linear_solver_method.is_some(),
        "primary-converged aggressive cases must record which linear solver was chosen"
    );
    assert!(!diagnostics.degraded_geometry_for_recovery);
    for edge in network.graph.edge_weights() {
        assert!(edge.resistance.is_finite() && edge.resistance > 0.0);
        assert!(edge.quad_coeff.is_finite() && edge.quad_coeff >= 0.0);
    }
}

#[test]
fn constrained_iteration_budget_classifies_max_iterations_exceeded() {
    let blueprint = primitive_tree(
        &[
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        0.45,
        0.55,
        0.68,
        35.0e-6,
        true,
        4,
    );
    let err = solve_blueprint(&blueprint, 260.0 / 6.0e7, 1)
        .expect_err("single-iteration budget must force a classified primary failure");

    assert_eq!(err.reason, SolveFailureReason::MaxIterationsExceeded);
    assert_eq!(err.diagnostics.picard_iterations, 1);
    assert!(
        err.diagnostics.linear_solver_method.is_some(),
        "classified failures after matrix assembly should record the chosen linear solver"
    );
}

#[test]
fn missing_branching_junction_metadata_is_rejected_before_solve() {
    let mut blueprint = primitive_tree(
        &[
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Bi,
        ],
        0.45,
        0.55,
        0.68,
        35.0e-6,
        true,
        2,
    );
    let branching_node = blueprint
        .nodes
        .iter()
        .position(|node| {
            node.junction_geometry.is_some()
                && blueprint
                    .channels
                    .iter()
                    .filter(|channel| channel.from == node.id || channel.to == node.id)
                    .count()
                    >= 3
        })
        .expect("primitive selective tree must carry explicit branching metadata");
    blueprint.nodes[branching_node].junction_geometry = None;

    let err = network_from_blueprint(&blueprint, CassonBlood::<f64>::normal_blood())
        .expect_err("branching nodes missing explicit geometry must be rejected");

    assert!(
        err.to_string()
            .contains("missing explicit junction metadata"),
        "unexpected error when junction metadata is missing: {err}"
    );
}

#[test]
fn missing_venturi_geometry_is_rejected_before_solve() {
    let mut blueprint = primitive_tree(
        &[
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        0.45,
        0.55,
        0.68,
        35.0e-6,
        true,
        2,
    );
    let throat_channel = blueprint
        .channels
        .iter_mut()
        .find(|channel| channel.venturi_geometry.is_some())
        .expect("primitive selective venturi tree must contain at least one throat channel");
    throat_channel.venturi_geometry = None;
    if let Some(metadata) = throat_channel.metadata.as_mut() {
        metadata.remove::<VenturiGeometryMetadata>();
    }

    let err = network_from_blueprint(&blueprint, CassonBlood::<f64>::normal_blood())
        .expect_err("channels declaring venturi throats must require explicit venturi geometry");

    assert!(
        err.to_string()
            .contains("declares venturi throats but has no explicit venturi geometry"),
        "unexpected error when venturi geometry is missing: {err}"
    );
}
