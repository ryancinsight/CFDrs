use super::*;
use crate::domain::network::{Network, NetworkBuilder};
use aequitas::systems::si::quantities::{HydraulicResistance, Pressure};
use cfd_core::physics::fluid::database::water_20c;

/// Empty-history slice is not stalled — fallback should not start.
#[test]
fn empty_history_is_not_stalling() {
    let budget = FallbackBudget::default();
    assert!(!budget.is_stalling(&[]));
}

/// Brief history (shorter than `stall_window`) is not yet pattern-classified.
#[test]
fn short_history_is_not_stalling() {
    let budget = FallbackBudget::default();
    let shorter: Vec<f64> = (0..budget.stall_window - 1).map(|_| 1.0).collect();
    assert!(!budget.is_stalling(&shorter));
}

/// The derived recovery work stays within the caller's finite iteration
/// budget, including the bounded Krylov work per Newton step.
#[test]
fn derived_budget_is_bounded() {
    for max_iterations in [1, 64, 200, 1_000, 10_000] {
        let budget = FallbackBudget::for_max_iterations(max_iterations);
        assert!(budget.warmup_iterations < max_iterations || max_iterations == 0);
        assert!(
            budget.warmup_iterations + budget.max_newton_iterations * budget.max_krylov_iterations
                <= max_iterations.max(1)
        );
    }
}

/// Bounded fixed-amplitude cycle is the canonical stalled trajectory.
#[test]
fn bounded_amplitude_cycle_is_stalling() {
    let budget = FallbackBudget::default();
    // max/min = 1.0 → ratio 1 ≤ stall_ratio (10) so this is stalled.
    let cycle: Vec<f64> = (0..budget.stall_window)
        .map(|i| if i % 2 == 0 { 1.0e-3 } else { 1.5e-3 })
        .collect();
    assert!(budget.is_stalling(&cycle));
}

/// Decaying envelope (typical converging Picard) is NOT stalled.
#[test]
fn decaying_envelope_is_not_stalling() {
    let budget = FallbackBudget::default();
    // max/min = 1e0 / 1e-12 = 1e12 ≫ stall_ratio (10) → not stalled.
    let decaying: Vec<f64> = (0..budget.stall_window)
        .map(|i| 10f64.powi(-(i as i32)))
        .collect();
    assert!(!budget.is_stalling(&decaying));
}

/// Oscillation at a slightly larger envelope (factor 5x) is still stalled
/// within Walker-Ni bounds.
#[test]
fn narrow_envelope_oscillation_is_stalling() {
    let budget = FallbackBudget::default();
    let narrow: Vec<f64> = (0..budget.stall_window)
        .map(|i| if i % 2 == 0 { 1.0e-3 } else { 5.0e-3 })
        .collect();
    // ratio 5 ≤ 10 → stalled
    assert!(budget.is_stalling(&narrow));
}

/// Wide-amplitude oscillation (factor 20x) is treated as productive.
#[test]
fn wide_envelope_oscillation_is_not_stalling() {
    let budget = FallbackBudget::default();
    let wide: Vec<f64> = (0..budget.stall_window)
        .map(|i| if i % 2 == 0 { 1.0e-3 } else { 2.0e-2 })
        .collect();
    // ratio 20 > 10 → not stalled
    assert!(!budget.is_stalling(&wide));
}

/// Transient non-finite samples are skipped, not interpreted as stagnation.
#[test]
fn nan_and_zero_samples_are_skipped() {
    let budget = FallbackBudget::default();
    let mut cycle = vec![0.0_f64; budget.stall_window];
    cycle[0] = f64::NAN;
    cycle[1] = 0.0;
    // All NaN/0 sample → no usable history, not stalled.
    assert!(!budget.is_stalling(&cycle));
}

/// A real two-node linear network exercises the full residual callback,
/// matrix assembly, linear solve, and JFNK convergence contract.
#[test]
fn linear_network_fallback_recovers_pressure_solution() {
    let fluid = water_20c::<f64>().expect("water properties are valid");
    let mut builder = NetworkBuilder::new();
    let inlet = builder.add_inlet("inlet".into());
    let outlet = builder.add_outlet("outlet".into());
    builder.connect_with_pipe(inlet, outlet, "pipe".into());
    let mut graph = builder.build().expect("network graph is valid");
    let edge = graph
        .edge_indices()
        .next()
        .expect("network contains the pipe edge");
    graph[edge].resistance = HydraulicResistance::from_base(2.0);

    let mut network = Network::new(graph, fluid);
    network.set_pressure(inlet, Pressure::from_base(10.0));
    network.set_pressure(outlet, Pressure::from_base(0.0));

    let n = network.node_count();
    let mut workspace = SolverWorkspace::new(n);
    MatrixAssembler::classify_boundary_conditions_into(
        &network,
        &mut workspace.dirichlet_values,
        &mut workspace.neumann_sources,
    )
    .expect("boundary conditions are valid");
    let assembler = MatrixAssembler::new();
    let matrix = assembler
        .assemble_into(&network, &mut workspace)
        .expect("network matrix is valid");
    let mut warm_solution = Array1::from_elem([n], 0.0);
    let initial_solution = LinearSystemSolver::new()
        .with_method(LinearSolverMethod::ConjugateGradient)
        .with_tolerance(1.0e-10)
        .with_max_iterations(20)
        .solve_with_initial_guess(&matrix, &workspace.rhs, &mut warm_solution)
        .expect("linear warm start is valid");
    let budget = FallbackBudget {
        warmup_iterations: 1,
        stall_window: 1,
        stall_ratio: 10.0,
        max_newton_iterations: 10,
        max_krylov_iterations: 10,
    };

    let solution = jfnk_fallback(
        &assembler,
        &network,
        workspace.dirichlet_values,
        workspace.neumann_sources,
        LinearSolverMethod::ConjugateGradient,
        1.0e-10,
        &budget,
        &initial_solution,
    )
    .expect("JFNK recovers the linear network root");

    assert!((solution[inlet.index()] - 10.0).abs() < 1.0e-8);
    assert!(solution[outlet.index()].abs() < 1.0e-8);
}
