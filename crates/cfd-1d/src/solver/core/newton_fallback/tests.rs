use super::*;
use crate::domain::network::{Network, NetworkBuilder};
use aequitas::systems::si::quantities::{HydraulicResistance, Pressure};
use cfd_core::physics::fluid::database::water_20c;
use petgraph::graph::NodeIndex;

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

/// Inlet held at 10 Pa and outlet at 0 Pa, joined by one pipe of the given
/// linear resistance.
fn pipe_network(
    resistance: f64,
) -> (
    Network<f64, impl FluidTrait<f64> + Clone>,
    NodeIndex,
    NodeIndex,
) {
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
    graph[edge].resistance = HydraulicResistance::from_base(resistance);

    let mut network = Network::new(graph, fluid);
    network.set_pressure(inlet, Pressure::from_base(10.0));
    network.set_pressure(outlet, Pressure::from_base(0.0));
    (network, inlet, outlet)
}

/// Workspace with the network's boundary conditions classified.
fn classified_workspace<F: FluidTrait<f64>>(network: &Network<f64, F>) -> SolverWorkspace<f64> {
    let mut workspace = SolverWorkspace::new(network.node_count());
    MatrixAssembler::classify_boundary_conditions_into(
        network,
        &mut workspace.dirichlet_values,
        &mut workspace.neumann_sources,
    )
    .expect("boundary conditions are valid");
    workspace
}

/// Budget with an explicit Newton limit and a fixed Krylov limit.
fn budget(max_newton_iterations: usize) -> FallbackBudget {
    FallbackBudget {
        warmup_iterations: 1,
        stall_window: 1,
        stall_ratio: 10.0,
        max_newton_iterations,
        max_krylov_iterations: 10,
    }
}

/// A real two-node linear network exercises the full residual callback,
/// matrix assembly, linear solve, and JFNK convergence contract.
#[test]
fn linear_network_fallback_recovers_pressure_solution() {
    let (network, inlet, outlet) = pipe_network(2.0);
    let n = network.node_count();
    let mut workspace = classified_workspace(&network);
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

    let solution = jfnk_fallback(
        &assembler,
        &network,
        workspace.dirichlet_values,
        workspace.neumann_sources,
        LinearSolverMethod::ConjugateGradient,
        1.0e-10,
        &budget(10),
        &initial_solution,
    )
    .expect("JFNK recovers the linear network root");

    assert!((solution[inlet.index()] - 10.0).abs() < 1.0e-8);
    assert!(solution[outlet.index()].abs() < 1.0e-8);
}

/// A warm start whose length differs from the node count is rejected as a
/// typed configuration error before any residual evaluation.
#[test]
fn warm_start_length_mismatch_is_invalid_configuration() {
    let (network, _, _) = pipe_network(2.0);
    let workspace = classified_workspace(&network);
    let warm_solution = Array1::from_elem([network.node_count() + 1], 0.0);

    let err = jfnk_fallback(
        &MatrixAssembler::new(),
        &network,
        workspace.dirichlet_values,
        workspace.neumann_sources,
        LinearSolverMethod::ConjugateGradient,
        1.0e-10,
        &budget(10),
        &warm_solution,
    )
    .expect_err("a mismatched warm start has no defined residual");

    assert!(
        matches!(
            &err,
            Error::InvalidConfiguration(message)
                if message.contains("warm-start length 3 does not match network node count 2")
        ),
        "unexpected error: {err:?}"
    );
}

/// A physically invalid edge coefficient raised inside the residual
/// callback reaches the caller as the original typed error.
#[test]
fn invalid_edge_resistance_propagates_as_invalid_configuration() {
    let (network, _, _) = pipe_network(-2.0);
    let workspace = classified_workspace(&network);
    let warm_solution = Array1::from_elem([network.node_count()], 0.0);

    let err = jfnk_fallback(
        &MatrixAssembler::new(),
        &network,
        workspace.dirichlet_values,
        workspace.neumann_sources,
        LinearSolverMethod::ConjugateGradient,
        1.0e-10,
        &budget(10),
        &warm_solution,
    )
    .expect_err("negative resistance has no admissible linearization");

    assert!(
        matches!(
            &err,
            Error::InvalidConfiguration(message) if message.contains("negative resistance")
        ),
        "unexpected error: {err:?}"
    );
}

/// A zero-iteration Newton budget from an unconverged warm start reports
/// the typed non-convergence variant carrying that budget.
#[test]
fn exhausted_newton_budget_is_typed_non_convergence() {
    let (network, _, _) = pipe_network(2.0);
    let workspace = classified_workspace(&network);
    let warm_solution = Array1::from_elem([network.node_count()], 0.0);

    let err = jfnk_fallback(
        &MatrixAssembler::new(),
        &network,
        workspace.dirichlet_values,
        workspace.neumann_sources,
        LinearSolverMethod::ConjugateGradient,
        1.0e-10,
        &budget(0),
        &warm_solution,
    )
    .expect_err("the zero-iteration budget cannot reach the root from zero pressure");

    assert!(
        matches!(
            err,
            Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { max: 0 })
        ),
        "unexpected error: {err:?}"
    );
}
