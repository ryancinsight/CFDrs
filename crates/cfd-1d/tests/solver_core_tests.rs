//! Tests for solver core modules: network construction, Ohm's law, KCL, error paths,
//! and ConvergenceChecker invariants.
//!
//! All tests use the correct `NetworkBuilder → Network → NetworkProblem → NetworkSolver` chain.

use approx::assert_relative_eq;
use cfd_1d::solver::core::ConvergenceChecker;
use cfd_1d::{Network, NetworkBuilder, NetworkProblem, NetworkSolver};
use cfd_core::physics::fluid::database::water_20c;
use nalgebra::DVector;
use petgraph::visit::EdgeRef;

fn water() -> cfd_core::physics::fluid::ConstantPropertyFluid<f64> {
    water_20c::<f64>().expect("test invariant")
}

fn hp_resistance(d: f64, l: f64, mu: f64) -> f64 {
    128.0 * mu * l / (std::f64::consts::PI * d.powi(4))
}

// ============================================================
// Ohm's Law: Q = ΔP / R for a single edge
// ============================================================

/// For a 2-node network with a single HP resistor, Q = ΔP / R exactly.
#[test]
fn test_solver_ohms_law_single_edge() {
    let fluid = water();
    let d = 1e-3_f64;
    let l = 0.1_f64;
    let r = hp_resistance(d, l, fluid.viscosity);
    let p_in = 1000.0_f64;

    let mut builder = NetworkBuilder::new();
    let inlet = builder.add_inlet("inlet".into());
    let outlet = builder.add_outlet("outlet".into());
    builder.connect_with_pipe(inlet, outlet, "pipe".into());

    let mut graph = builder.build().expect("test invariant");
    for edge in graph.edge_indices() {
        if let Some(e) = graph.edge_weight_mut(edge) {
            e.resistance = r;
        }
    }

    let mut network = Network::new(graph, fluid.clone());
    network.set_pressure(inlet, p_in);
    network.set_pressure(outlet, 0.0);

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::new();
    let solved = solver.solve_network(&problem).expect("test invariant");

    let flows: Vec<f64> = solved
        .graph
        .edge_references()
        .map(|e| {
            let p_from = solved.pressures().get(e.source().index()).copied().unwrap_or(0.0);
            let p_to = solved.pressures().get(e.target().index()).copied().unwrap_or(0.0);
            (p_from - p_to) / e.weight().resistance
        })
        .collect();

    assert_eq!(flows.len(), 1);
    assert_relative_eq!(flows[0], p_in / r, max_relative = 1e-9);
}

// ============================================================
// KCL: mass conservation at Y-junction
// ============================================================

/// KCL: Q_in = Q_out1 + Q_out2 at a symmetric T-junction.
#[test]
fn test_solver_y_junction_kcl() {
    let fluid = water();
    let d = 1e-3_f64;
    let l = 0.05_f64;
    let r = hp_resistance(d, l, fluid.viscosity);

    let mut builder = NetworkBuilder::new();
    let n_in = builder.add_inlet("in".into());
    let n_mid = builder.add_junction("mid".into());
    let n_out1 = builder.add_outlet("out1".into());
    let n_out2 = builder.add_outlet("out2".into());
    builder.connect_with_pipe(n_in, n_mid, "edge0".into());
    builder.connect_with_pipe(n_mid, n_out1, "edge1".into());
    builder.connect_with_pipe(n_mid, n_out2, "edge2".into());

    let mut graph = builder.build().expect("test invariant");
    for edge in graph.edge_indices() {
        if let Some(e) = graph.edge_weight_mut(edge) {
            e.resistance = r;
        }
    }

    let mut network = Network::new(graph, fluid.clone());
    network.set_pressure(n_in, 1000.0);
    network.set_pressure(n_out1, 0.0);
    network.set_pressure(n_out2, 0.0);

    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::new();
    let solved = solver.solve_network(&problem).expect("test invariant");

    let mut q_into_mid = 0.0_f64;
    let mut q_out_of_mid = 0.0_f64;
    for e in solved.graph.edge_references() {
        let p_from = solved.pressures().get(e.source().index()).copied().unwrap_or(0.0);
        let p_to = solved.pressures().get(e.target().index()).copied().unwrap_or(0.0);
        let q = (p_from - p_to) / e.weight().resistance;
        if e.target() == n_mid {
            q_into_mid += q;
        }
        if e.source() == n_mid {
            q_out_of_mid += q;
        }
    }
    assert_relative_eq!(q_into_mid, q_out_of_mid, max_relative = 1e-9);
}

// ============================================================
// Error paths: singular system (no Dirichlet BCs)
// ============================================================

/// Network missing Dirichlet BCs prevents a unique solution.
///
/// The solver must NOT silently return a plausible solution for an unconstrained
/// (rank-deficient) system. It must either return an `Err` or panic.
/// We accept both — the critical invariant is "no silent wrong answer".
#[test]
fn test_solver_no_dirichlet_bc_singular_system() {
    let fluid = water();
    let d = 1e-3_f64;
    let l = 0.1_f64;
    let r = hp_resistance(d, l, fluid.viscosity);

    // Build a valid topology (inlet+outlet required by builder) but do NOT set pressures,
    // so the assembled Laplacian has no Dirichlet rows → singular system.
    let mut builder = NetworkBuilder::new();
    let inlet = builder.add_inlet("n0".into());
    let outlet = builder.add_outlet("n1".into());
    builder.connect_with_pipe(inlet, outlet, "edge".into());

    let mut graph = builder.build().expect("test invariant");
    for edge in graph.edge_indices() {
        if let Some(e) = graph.edge_weight_mut(edge) {
            e.resistance = r;
        }
    }

    // Intentionally skip set_pressure() → no Dirichlet BC → null space of Laplacian is non-trivial
    let network = Network::new(graph, fluid);
    let problem = NetworkProblem::new(network);
    let solver = NetworkSolver::new();

    let outcome = std::panic::catch_unwind(|| solver.solve_network(&problem));
    match outcome {
        Ok(Ok(_)) => panic!("Solver must not succeed on a singular (unconstrained) system"),
        Ok(Err(_)) => { /* Graceful error — ideal */ }
        Err(_) => { /* Panic — acceptable for an unconstrained Laplacian */ }
    }
}

// ============================================================
// ConvergenceChecker invariants
// ============================================================

#[test]
fn test_convergence_nan_returns_diverged() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    let sol = DVector::from_vec(vec![1.0, f64::NAN]);
    assert!(checker.check(&sol).is_err());
}

#[test]
fn test_convergence_inf_returns_diverged() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    let sol = DVector::from_vec(vec![f64::INFINITY, 1.0]);
    assert!(checker.check(&sol).is_err());
}

#[test]
fn test_convergence_max_iterations() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    assert!(checker.max_iterations_reached(1000));
    assert!(!checker.max_iterations_reached(999));
}

#[test]
fn test_convergence_dual_identical_converged() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    let x = DVector::from_vec(vec![1.0_f64, 2.0, 3.0]);
    let converged = checker
        .has_converged_dual(&x, &x, 0.0, 1.0)
        .expect("test invariant");
    assert!(converged);
}

#[test]
fn test_convergence_dual_large_change_not_converged() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    let x_old = DVector::from_vec(vec![0.0_f64, 0.0]);
    let x_new = DVector::from_vec(vec![100.0_f64, 100.0]);
    let converged = checker
        .has_converged_dual(&x_new, &x_old, 50.0, 1.0)
        .expect("test invariant");
    assert!(!converged);
}
