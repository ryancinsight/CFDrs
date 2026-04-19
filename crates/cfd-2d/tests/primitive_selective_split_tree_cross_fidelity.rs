use cfd_2d::network::{solve_reference_trace, Network2dBuilderSink};
use cfd_core::physics::fluid::BloodModel;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{NetworkBlueprint, NodeKind};
use cfd_schematics::geometry::generator::PrimitiveSelectiveSplitKind;
use cfd_schematics::interface::presets::primitive_selective_split_tree_rect;

const PLATE_DIMS_MM: (f64, f64) = (127.76, 85.47);

#[derive(Clone, Copy)]
struct ComparisonCase {
    label: &'static str,
    split_kind: PrimitiveSelectiveSplitKind,
    max_flow_error_pct: f64,
    mean_flow_error_pct: f64,
}

fn count_inlets_outlets(blueprint: &NetworkBlueprint) -> (usize, usize) {
    let inlets = blueprint
        .nodes
        .iter()
        .filter(|node| node.kind == NodeKind::Inlet)
        .count();
    let outlets = blueprint
        .nodes
        .iter()
        .filter(|node| node.kind == NodeKind::Outlet)
        .count();
    (inlets, outlets)
}

fn build_selective_tree_case(case: ComparisonCase) -> NetworkBlueprint {
    primitive_selective_split_tree_rect(
        case.label,
        PLATE_DIMS_MM,
        &[case.split_kind],
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
    )
}

fn compare_case(case: ComparisonCase) {
    let q_total = 1.0e-6_f64;
    let blood = BloodModel::Newtonian(3.5e-3_f64);
    let blueprint = build_selective_tree_case(case);

    let (inlets, outlets) = count_inlets_outlets(&blueprint);
    assert_eq!(inlets, 1, "{} must retain a single inlet", case.label);
    assert_eq!(outlets, 1, "{} must retain a single outlet", case.label);
    assert!(
        (blueprint.box_dims.0 - PLATE_DIMS_MM.0).abs() < 1e-9
            && (blueprint.box_dims.1 - PLATE_DIMS_MM.1).abs() < 1e-9,
        "{} must stay within the 96-well plate footprint",
        case.label
    );
    solve_reference_trace::<f64>(&blueprint, 1060.0, 3.5e-3_f64, q_total)
        .expect("reference trace should build");

    let sink = Network2dBuilderSink::new(blood, 1060.0, q_total, 24, 10);
    let mut network = sink.build(&blueprint).expect("build should succeed");
    let projection = network.projection_summary_ref();
    assert!(
        projection
            .channel_summaries
            .iter()
            .all(|summary| summary.fluid_cell_count > 0),
        "{} all projected channels must retain fluid occupancy",
        case.label
    );
    let result = network.solve_all(1e-6).expect("solve should succeed");

    let reference = result.reference_trace.clone();
    let inlet_outlet_error =
        (reference.total_inlet_flow_m3_s - reference.total_outlet_flow_m3_s).abs() / q_total;

    assert!(
        inlet_outlet_error < 1e-8,
        "{} reference trace must conserve mass",
        case.label
    );
    assert_eq!(
        result.channels.len(),
        reference.channel_traces.len(),
        "{} must produce a 1D reference entry for every 2D channel",
        case.label
    );
    assert!(
        result.converged_count > 0,
        "{} 2D channel solves must converge at least once",
        case.label
    );
    assert!(
        result.max_field_outlet_flow_error_pct <= case.max_flow_error_pct,
        "{} max 2D/1D outlet-flow error {:.2}% exceeds {:.2}%",
        case.label,
        result.max_field_outlet_flow_error_pct,
        case.max_flow_error_pct
    );
    assert!(
        result.mean_field_outlet_flow_error_pct <= case.mean_flow_error_pct,
        "{} mean 2D/1D outlet-flow error {:.2}% exceeds {:.2}%",
        case.label,
        result.mean_field_outlet_flow_error_pct,
        case.mean_flow_error_pct
    );
}

#[test]
fn primitive_selective_split_tree_1d_2d_cross_fidelity() {
    let cases = [
        ComparisonCase {
            label: "Bi",
            split_kind: PrimitiveSelectiveSplitKind::Bi,
            max_flow_error_pct: 50.0,
            mean_flow_error_pct: 35.0,
        },
        ComparisonCase {
            label: "Tri",
            split_kind: PrimitiveSelectiveSplitKind::Tri,
            max_flow_error_pct: 50.0,
            mean_flow_error_pct: 35.0,
        },
        ComparisonCase {
            label: "Quad",
            split_kind: PrimitiveSelectiveSplitKind::Quad,
            max_flow_error_pct: 50.0,
            mean_flow_error_pct: 35.0,
        },
        ComparisonCase {
            label: "Penta",
            split_kind: PrimitiveSelectiveSplitKind::Penta,
            max_flow_error_pct: 50.0,
            mean_flow_error_pct: 35.0,
        },
    ];

    for case in cases {
        compare_case(case);
    }
}
