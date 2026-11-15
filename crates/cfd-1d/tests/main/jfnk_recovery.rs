//! Solver-level witnesses for the bounded JFNK recovery.
//!
//! Network: inlet 100 Pa, edge A (`R = 0.2`, `k = 1`), junction, edge B
//! (`R = 0.1`, `k = 0`), outlet 0 Pa.
//!
//! # Derived: the lagged Picard map alternates slowly
//!
//! Each edge is linearized as `r = R + k |Q_prev|`. One assembly gives the
//! divider `p = 100 r_B / (r_A + r_B)` and both edges the flow
//! `Q = 100 / (r_A + r_B)`, so `Q_{m+1} = 100 / (0.3 + Q_m)` with fixed point
//! `Q* = 9.851`, `p* = 0.985` and multiplier `-Q* / (0.3 + Q*) = -0.9704`: the
//! iterates alternate and shrink by only 3 % per step.
//!
//! # Observed, not derived
//!
//! The loop also mixes iterates (Anderson, then damping), so its decay is not
//! the undamped multiplier. At `tolerance = 1e-8`, `max_iterations = 300` the
//! pressure change norm is 3.34e-2 at iteration 136 and 7.63e-3 at iteration
//! 199: a window ratio of 4.4, below `stall_ratio = 10`, so the detector fires
//! at iteration 199. No analytic window bound is claimed.
//!
//! # Derived: what one JFNK trajectory leaves behind
//!
//! JFNK stops at `x_s` with `||F(x_s)|| < atol = tolerance`, `F(x) = G(x) - x`
//! on the network snapshot. The solver installs `x_s` with the flows `F`
//! computed, so the next pressure change is at most `||F(x_s)||`: a change
//! below `tolerance` right after the stall iteration proves a converged
//! trajectory ran. The flows are not in the JFNK state, so the cycle resumes
//! (`ATLAS-CFDRS-JFNK-RECOVERY-110`).
//!
//! # Derived and observed: budget and refire
//!
//! `FallbackBudget::for_max_iterations(n)` gives `warmup = min(200, n - 1)`,
//! `remaining = n - warmup`, `krylov = min(50, remaining)` and
//! `newton = max(1, remaining / krylov)`: one Newton step for `n` of 201, 250
//! and 267. The check at iteration `i` reads iterations `i - 63 ..= i`. The
//! jump leaves 1e-11 at iteration 200 and the cycle resumes with 3.9e-2 at 201
//! (observed); windows ending at 200 to 264 contain one of them (derived).
//! Iteration 265 reads 202 to 265, whose observed ratio is below 10: the first
//! iteration at which a detector without the once-only latch fires again.

use aequitas::systems::si::quantities::{
    HydraulicResistance, Pressure, QuadraticHydraulicResistance,
};
use cfd_1d::domain::network::{Edge, EdgeType, Network, NetworkBuilder};
use cfd_1d::{NetworkProblem, NetworkSolver, PrimarySolveError, SolveFailureReason, SolverConfig};
use cfd_core::error::{ConvergenceErrorKind, Error};
use cfd_core::physics::fluid::newtonian::ConstantPropertyFluid;

/// Tolerance of the stalled solve; also JFNK's `atol`.
const TOLERANCE: f64 = 1.0e-8;

/// Picard iterations run when the detector fires at iteration 199.
const STALLED_ITERATIONS: usize = 200;

fn stalling_chain() -> Network<f64, ConstantPropertyFluid<f64>> {
    let mut builder = NetworkBuilder::<f64>::new();
    let inlet = builder.add_inlet("inlet".to_string());
    let junction = builder.add_junction("junction".to_string());
    let outlet = builder.add_outlet("outlet".to_string());
    let first = builder.add_edge(
        inlet,
        junction,
        Edge::new("first".to_string(), EdgeType::Pipe),
    );
    let second = builder.add_edge(
        junction,
        outlet,
        Edge::new("second".to_string(), EdgeType::Pipe),
    );
    let graph = builder.build().expect("series graph is valid");
    let fluid = ConstantPropertyFluid::<f64>::water_20c().expect("water properties are valid");
    let mut network = Network::new(graph, fluid);
    for (edge, resistance, quadratic) in [(first, 0.2, 1.0), (second, 0.1, 0.0)] {
        let weight = network
            .graph
            .edge_weight_mut(edge)
            .expect("series edge exists");
        weight.resistance = HydraulicResistance::from_base(resistance);
        weight.quad_coeff = QuadraticHydraulicResistance::from_base(quadratic);
    }
    network.set_pressure(inlet, Pressure::from_base(100.0));
    network.set_pressure(outlet, Pressure::from_base(0.0));
    network
}

/// Solve the stalling chain, which converges in no budget used here.
fn failed_solve(max_iterations: usize, tolerance: f64) -> PrimarySolveError {
    let solver = NetworkSolver::<f64, ConstantPropertyFluid<f64>>::with_config(SolverConfig {
        tolerance,
        max_iterations,
        require_flow_convergence: true,
    });
    match solver.solve_network_with_diagnostics(&NetworkProblem::new(stalling_chain())) {
        Ok((_, diagnostics)) => panic!("the chain must not converge: {diagnostics:?}"),
        Err(error) => error,
    }
}

fn source_error(error: &PrimarySolveError) -> &Error {
    std::error::Error::source(error)
        .and_then(|source| source.downcast_ref::<Error>())
        .expect("primary solve error carries a cfd-core error source")
}

fn final_change(error: &PrimarySolveError) -> f64 {
    error
        .diagnostics
        .last_solution_change_norm
        .expect("the final iteration records its change norm")
}

/// With one iteration left after the stall, the final pressure change is below
/// `tolerance` only if a converged JFNK trajectory was installed.
#[test]
fn stalled_picard_cycle_runs_a_converged_jfnk_trajectory() {
    let error = failed_solve(STALLED_ITERATIONS + 1, TOLERANCE);

    assert_eq!(error.reason, SolveFailureReason::MaxIterationsExceeded);
    assert!(error.diagnostics.fallback_used);
    assert_eq!(error.diagnostics.picard_iterations, STALLED_ITERATIONS + 1);
    let change = final_change(&error);
    assert!(
        change < TOLERANCE,
        "pressure change after the JFNK jump must be below atol: {change:e}"
    );
}

/// The trajectory runs once: an unlatched detector would fire again at
/// iteration 265 and end a 267-iteration solve below `tolerance`.
#[test]
fn jfnk_trajectory_is_not_repeated_when_the_window_refills() {
    let error = failed_solve(267, TOLERANCE);

    assert_eq!(error.reason, SolveFailureReason::MaxIterationsExceeded);
    assert_eq!(error.diagnostics.picard_iterations, 267);
    let change = final_change(&error);
    assert!(
        change > TOLERANCE,
        "the Picard cycle, not a second JFNK jump, ends the solve: {change:e}"
    );
}

/// An exhausted one-step Newton budget fails the solve at the stall iteration
/// with that budget in the typed source: `||F||` stays near 1e-11 (observed).
#[test]
fn exhausted_jfnk_budget_fails_the_solve_at_the_stall_iteration() {
    let error = failed_solve(250, 1.0e-13);

    assert_eq!(error.reason, SolveFailureReason::MaxIterationsExceeded);
    assert!(error.diagnostics.fallback_used);
    assert_eq!(error.diagnostics.picard_iterations, STALLED_ITERATIONS);
    assert!(matches!(
        source_error(&error),
        Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { max: 1 })
    ));
}
