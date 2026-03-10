use cfd_1d::domain::network::network_from_blueprint;
use cfd_1d::{NetworkProblem, NetworkSolver};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
use cfd_schematics::interface::presets::primitive_selective_split_tree_rect;
use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};

#[derive(Debug, Clone)]
struct PrimitiveSolveCase {
    name: &'static str,
    sequence: Vec<PrimitiveSelectiveSplitKind>,
    first_trifurcation_center_frac: f64,
    later_trifurcation_center_frac: f64,
    bifurcation_treatment_frac: f64,
    throat_width_m: f64,
    venturi_treatment_enabled: bool,
    throat_count: u8,
    inlet_flow_m3_s: f64,
}

fn primitive_tree(case: &PrimitiveSolveCase) -> cfd_schematics::NetworkBlueprint {
    primitive_selective_split_tree_rect(
        case.name,
        (127.76, 85.47),
        &case.sequence,
        4.0e-3,
        case.first_trifurcation_center_frac,
        case.later_trifurcation_center_frac,
        case.bifurcation_treatment_frac,
        case.throat_width_m,
        120.0e-6,
        1.0e-3,
        case.venturi_treatment_enabled,
        case.throat_count,
        None,
    )
}

fn representative_primary_cases() -> Vec<PrimitiveSolveCase> {
    let sequences = [
        (
            "Bi->Tri",
            vec![
                PrimitiveSelectiveSplitKind::Bi,
                PrimitiveSelectiveSplitKind::Tri,
            ],
        ),
        (
            "Tri->Bi",
            vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Bi,
            ],
        ),
        (
            "Tri->Tri",
            vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
        ),
        (
            "Tri->Tri->Tri",
            vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
        ),
    ];
    let first_tri_fracs = [0.33, 0.40, 0.45];
    let later_tri_fracs = [0.25, 0.33, 0.50];
    let bi_treat_fracs = [0.55, 0.68];
    let throat_widths_m = [35.0e-6, 45.0e-6, 120.0e-6];
    let inlet_flows_m3_s = [80.0 / 6.0e7, 150.0 / 6.0e7, 300.0 / 6.0e7];
    let solver = NetworkSolver::<f64, CassonBlood<f64>>::new();

    sequences
        .into_iter()
        .map(|(label, sequence)| {
            for &venturi_treatment_enabled in &[false, true] {
                for &throat_width_m in &throat_widths_m {
                    for &first_trifurcation_center_frac in &first_tri_fracs {
                        for &later_trifurcation_center_frac in &later_tri_fracs {
                            for &bifurcation_treatment_frac in &bi_treat_fracs {
                                for &inlet_flow_m3_s in &inlet_flows_m3_s {
                                    let case = PrimitiveSolveCase {
                                        name: label,
                                        sequence: sequence.clone(),
                                        first_trifurcation_center_frac,
                                        later_trifurcation_center_frac,
                                        bifurcation_treatment_frac,
                                        throat_width_m,
                                        venturi_treatment_enabled,
                                        throat_count: 1,
                                        inlet_flow_m3_s,
                                    };
                                    let blueprint = primitive_tree(&case);
                                    let network = prepared_network(&blueprint, inlet_flow_m3_s);
                                    if solver
                                        .solve_network_with_diagnostics(&NetworkProblem::new(
                                            network,
                                        ))
                                        .is_ok()
                                    {
                                        return case;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            panic!("no representative primary-converged benchmark case found for {label}");
        })
        .collect()
}

fn representative_blueprints() -> Vec<(&'static str, cfd_schematics::NetworkBlueprint)> {
    representative_primary_cases()
        .into_iter()
        .map(|case| (case.name, primitive_tree(&case)))
        .collect()
}

fn prepared_network(
    blueprint: &cfd_schematics::NetworkBlueprint,
    inlet_flow_m3_s: f64,
) -> cfd_1d::Network<f64, CassonBlood<f64>> {
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
    network
}

fn bench_network_from_blueprint(c: &mut Criterion) {
    let mut group = c.benchmark_group("primitive_selective_network_from_blueprint");
    let blueprints = representative_blueprints();
    let blood = CassonBlood::<f64>::normal_blood();

    for (label, blueprint) in &blueprints {
        group.bench_with_input(BenchmarkId::from_parameter(label), blueprint, |b, bp| {
            b.iter(|| black_box(network_from_blueprint(black_box(bp), black_box(blood))).unwrap());
        });
    }

    group.finish();
}

fn bench_solve_network(c: &mut Criterion) {
    let mut group = c.benchmark_group("primitive_selective_solve_network");
    let solver = NetworkSolver::<f64, CassonBlood<f64>>::new();
    let representative = representative_primary_cases();

    for case in &representative {
        let blueprint = primitive_tree(case);
        let network = prepared_network(&blueprint, case.inlet_flow_m3_s);
        group.bench_with_input(
            BenchmarkId::new("primary", case.name),
            &network,
            |b, template| {
                b.iter_batched(
                    || template.clone(),
                    |network| {
                        black_box(
                            solver
                                .solve_network_with_diagnostics(&NetworkProblem::new(network))
                                .expect("representative primitive selective tree must converge"),
                        )
                    },
                    BatchSize::SmallInput,
                );
            },
        );
    }

    let difficult = primitive_tree(&PrimitiveSolveCase {
        name: "bench-difficult-tri-tri-tri",
        sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        first_trifurcation_center_frac: 0.45,
        later_trifurcation_center_frac: 0.55,
        bifurcation_treatment_frac: 0.68,
        throat_width_m: 45.0e-6,
        venturi_treatment_enabled: true,
        throat_count: 4,
        inlet_flow_m3_s: 260.0 / 6.0e7,
    });
    let difficult_network = prepared_network(&difficult, 260.0 / 6.0e7);
    group.bench_function("difficult/Tri->Tri->Tri", |b| {
        b.iter_batched(
            || difficult_network.clone(),
            |network| {
                black_box(solver.solve_network_with_diagnostics(&NetworkProblem::new(network)))
            },
            BatchSize::SmallInput,
        );
    });

    let outliers = [
        (
            "outlier/Bi->Tri",
            primitive_tree(&PrimitiveSolveCase {
                name: "bench-outlier-bi-tri",
                sequence: vec![
                    PrimitiveSelectiveSplitKind::Bi,
                    PrimitiveSelectiveSplitKind::Tri,
                ],
                first_trifurcation_center_frac: 0.45,
                later_trifurcation_center_frac: 0.40,
                bifurcation_treatment_frac: 0.55,
                throat_width_m: 120.0e-6,
                venturi_treatment_enabled: false,
                throat_count: 1,
                inlet_flow_m3_s: 80.0 / 6.0e7,
            }),
        ),
        (
            "outlier/Tri->Bi",
            primitive_tree(&PrimitiveSolveCase {
                name: "bench-outlier-tri-bi",
                sequence: vec![
                    PrimitiveSelectiveSplitKind::Tri,
                    PrimitiveSelectiveSplitKind::Bi,
                ],
                first_trifurcation_center_frac: 0.40,
                later_trifurcation_center_frac: 0.50,
                bifurcation_treatment_frac: 0.60,
                throat_width_m: 120.0e-6,
                venturi_treatment_enabled: false,
                throat_count: 1,
                inlet_flow_m3_s: 80.0 / 6.0e7,
            }),
        ),
        (
            "outlier/Tri->Tri",
            primitive_tree(&PrimitiveSolveCase {
                name: "bench-outlier-tri-tri",
                sequence: vec![
                    PrimitiveSelectiveSplitKind::Tri,
                    PrimitiveSelectiveSplitKind::Tri,
                ],
                first_trifurcation_center_frac: 0.45,
                later_trifurcation_center_frac: 0.55,
                bifurcation_treatment_frac: 0.68,
                throat_width_m: 120.0e-6,
                venturi_treatment_enabled: false,
                throat_count: 1,
                inlet_flow_m3_s: 80.0 / 6.0e7,
            }),
        ),
    ];
    for (label, blueprint) in outliers {
        let network = prepared_network(&blueprint, 80.0 / 6.0e7);
        group.bench_function(label, |b| {
            b.iter_batched(
                || network.clone(),
                |network| {
                    black_box(solver.solve_network_with_diagnostics(&NetworkProblem::new(network)))
                },
                BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, bench_network_from_blueprint, bench_solve_network);
criterion_main!(benches);
