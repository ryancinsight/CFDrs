//! Modularized network solver for 1D CFD analysis
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies.
//!
//! ## Mathematical Foundation
//!
//! ### Theorem: Network Conservation Principle (Mass Conservation at Nodes)
//!
//! **Statement**: For incompressible fluid flow in a network of channels, the conservation
//! of mass at each interior node requires that the sum of all inflow rates equals the
//! sum of all outflow rates.
//!
//! **Mathematical Formulation**:
//! ```text
//! ∑_{edges connected to node i} Q_{edge} = 0    for all interior nodes i
//! ```
//!
//! where Q_{edge} is positive for flow into the node and negative for flow out of the node.
//!
//! **Kirchhoff's Current Law Analogy**: This principle is analogous to Kirchhoff's current
//! law in electrical circuits, where the sum of currents entering a node equals the sum
//! leaving the node.
//!
//! **Proof**: Conservation of mass requires that mass cannot accumulate at nodes in
//! steady-state incompressible flow. For each node i:
//!
//! ∂ρ/∂t + ∇·(ρu) = 0
//!
//! For steady-state (∂ρ/∂t = 0) and incompressible (ρ = constant) flow:
//!
//! ∇·u = 0
//!
//! Integrating over a control volume around node i and applying the divergence theorem:
//!
//! ∮_S u·dS = 0
//!
//! Where S is the surface enclosing node i. This surface integral becomes a sum over
//! all channel cross-sections connected to the node.
//!
//! **Assumptions**:
//! 1. Steady-state flow (no time derivatives)
//! 2. Incompressible fluid (constant density)
//! 3. No mass sources/sinks at nodes
//! 4. One-dimensional flow in channels
//!
//! **Validity Conditions**: Applies to all microfluidic networks with incompressible fluids
//! under steady-state conditions.
//!
//! **Literature**: Kirchhoff, G. (1847). "Ueber die Auflösung der Gleichungen".
//! Annalen der Physik, 72, 497-508.
//!
//! ### Theorem: Pressure-Flow Relationship (Ohm's Law Analogy)
//!
//! **Statement**: The volumetric flow rate Q through a channel is proportional to the
//! pressure difference across the channel and inversely proportional to the hydraulic
//! resistance R of the channel.
//!
//! **Mathematical Formulation**:
//! ```text
//! Q = (ΔP) / R
//! ```
//!
//! where:
//! - Q is the volumetric flow rate [m³/s]
//! - ΔP is the pressure difference [Pa]
//! - R is the hydraulic resistance [Pa·s/m³]
//!
//! **Derivation**: Starting from the Hagen-Poiseuille equation for laminar flow:
//!
//! ΔP = (128μLQ)/(πD⁴)
//!
//! Rearranging gives:
//!
//! Q = (πD⁴ ΔP)/(128μL) = ΔP / R
//!
//! where R = (128μL)/(πD⁴)
//!
//! **Ohm's Law Analogy**: This relationship is analogous to Ohm's law V = IR in
//! electrical circuits, where pressure difference plays the role of voltage and
//! hydraulic resistance plays the role of electrical resistance.
//!
//! **Assumptions**:
//! 1. Laminar flow (Re < 2300)
//! 2. Fully developed flow (L/D > 10)
//! 3. Newtonian fluid
//! 4. Constant cross-section channel
//! 5. No entrance/exit losses (or included in R)
//!
//! **Literature**: Hagen, G. (1839). "Über die Bewegung der Flüssigkeiten".
//! Poggendorff's Annalen, 46, 423-442.
//!
//! ### Theorem: Circuit Analogy Validation
//!
//! **Statement**: The electrical circuit analogy provides a mathematically valid
//! representation of fluid networks under the assumptions of incompressible laminar flow.
//!
//! **Mathematical Foundation**: The system of equations for a fluid network:
//!
//! ```text
//! ∑_{j≠i} G_{ij} (P_i - P_j) = 0    for all interior nodes i
//! ```
//!
//! is isomorphic to the electrical circuit equations:
//!
//! ```text
//! ∑_{j≠i} G_{ij} (V_i - V_j) = 0    for all interior nodes i
//! ```
//!
//! where:
//! - P_i, V_i are pressures/voltages at nodes
//! - G_{ij} = 1/R_{ij} is conductance (inverse resistance)
//!
//! **Proof of Isomorphism**: Both systems are linear and can be written in matrix form:
//!
//! A x = b
//!
//! where A is the conductance matrix, x is the pressure/voltage vector, and b accounts
//! for boundary conditions and sources.
//!
//! **Assumptions and Limitations**:
//! 1. **Linear Relationship**: Pressure-flow relationship must be linear (valid for laminar flow)
//! 2. **Incompressibility**: Fluid density constant (no capacitive effects)
//! 3. **Steady State**: No time derivatives (no inductive effects)
//! 4. **Passive Elements**: No active pressure sources (pumps, etc.) in basic analogy
//!
//! ## Invariants and Boundary Handling
//!
//! - Units: `[A] = [Pa·s/m³]`, `[x] = [Pa]`, `[b] = [Pa]`
//! - Positivity: all edge conductances `G = 1/R` strictly positive
//! - Dirichlet BCs: enforced exactly via identity rows; preserves SPD
//! - Neumann BCs: applied as source terms in `b`; discrete conservation holds at junctions
//! - Nonlinear losses: if `ΔP = R·Q + k·Q²`, assemble with `R_eff = R + 2k|Q|` per iteration
//!
//! **Extensions**: Pumps can be modeled as voltage sources, valves as variable resistors.
//!
//! **Literature**: White, F.M. (2015). Fluid Mechanics (8th ed.). McGraw-Hill. §6.8.
//!
//! ### Theorem: Matrix Assembly and System Formulation
//!
//! **Statement**: The network equations can be assembled into a sparse linear system
//! Ax = b where the matrix A represents network topology and hydraulic conductances.
//!
//! **Mathematical Formulation**: For a network with N nodes and M edges:
//!
//! ```text
//! A ∈ ℝ^{N×N},    x ∈ ℝ^N,    b ∈ ℝ^N
//! ```
//!
//! **Matrix Assembly Rules**:
//! For each edge connecting nodes i and j with conductance G:
//!
//! ```text
//! A[i,i] += G    (diagonal contribution)
//! A[j,j] += G    (diagonal contribution)
//! A[i,j] -= G    (off-diagonal coupling)
//! A[j,i] -= G    (symmetric coupling)
//! ```
//!
//! **Physical Interpretation**: The matrix A represents the discrete Laplacian operator
//! for the network, where diagonal terms represent the total conductance connected to
//! each node, and off-diagonal terms represent inter-node coupling.
//!
//! **Boundary Conditions**: Dirichlet (fixed pressure) and Neumann (fixed flow) conditions
//! modify the system appropriately.
//!
//! **Existence and Uniqueness**: For connected networks with appropriate boundary conditions,
//! the system has a unique solution.
//!
//! **Proof**: The assembled matrix A is symmetric positive definite for connected networks
//! with proper boundary conditions, ensuring well-posedness.
//!
//! **Literature**: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow.
//! Hemisphere Publishing. §6.3-6.5.

mod convergence;
mod domain;
mod linear_system;
mod matrix_assembly;
mod problem;
mod state;

pub use convergence::ConvergenceChecker;
pub use domain::NetworkDomain;
pub use linear_system::{LinearSystemSolver, LinearSolverMethod};
pub use matrix_assembly::MatrixAssembler;
pub use problem::NetworkProblem;
pub use state::NetworkState;

use crate::network::Network;
use cfd_core::error::Result;
use cfd_core::compute::solver::{Configurable, Solver, Validatable};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig<T: RealField + Copy> {
    /// Convergence tolerance for solution accuracy
    pub tolerance: T,
    /// Maximum number of solver iterations before termination
    pub max_iterations: usize,
}

impl<T: RealField + Copy> cfd_core::compute::solver::SolverConfiguration<T> for SolverConfig<T> {
    fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }

    fn use_preconditioning(&self) -> bool {
        false // No preconditioning for network solver
    }
}

/// Main network solver implementing the core CFD suite trait system
pub struct NetworkSolver<T: RealField + Copy> {
    /// Solver configuration
    config: SolverConfig<T>,
    /// Matrix assembler for building the linear system
    assembler: MatrixAssembler<T>,
    /// Linear system solver
    #[allow(dead_code)]
    linear_solver: LinearSystemSolver<T>,
    /// Convergence checker
    convergence: ConvergenceChecker<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for NetworkSolver<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> NetworkSolver<T> {
    /// Create a new network solver with default configuration
    #[must_use]
    pub fn new() -> Self {
        let tolerance = T::from_f64(1e-6).expect(
            "Failed to represent the default tolerance value (1e-6) in the target numeric type T",
        );
        let config = SolverConfig {
            tolerance,
            max_iterations: 1000,
        };
        Self {
            config: config.clone(),
            assembler: MatrixAssembler::new(),
            linear_solver: LinearSystemSolver::new(),
            convergence: ConvergenceChecker::new(config.tolerance),
        }
    }

    /// Create with specific configuration
    pub fn with_config(config: SolverConfig<T>) -> Self {
        Self {
            assembler: MatrixAssembler::new(),
            linear_solver: LinearSystemSolver::new(),
            convergence: ConvergenceChecker::new(config.tolerance),
            config,
        }
    }

    /// Solve the network flow problem iteratively for non-linear systems
    pub fn solve_network(&self, problem: &NetworkProblem<T>) -> Result<Network<T>> {
        let mut network = problem.network.clone();
        let n = network.node_count();

        // Initialize solution vectors for convergence checking
        let mut last_solution = nalgebra::DVector::zeros(n);

        // Main iteration loop for non-linear problems
        for _iter in 0..self.config.max_iterations {
            // 1. Assemble the system based on the CURRENT network state
            //    This is the linearization step for non-linear problems
            let (matrix, rhs) = self.assembler.assemble(&network)?;

            // Heuristic SPD detection: positive diagonal, non-positive off-diagonals,
            // and diagonal dominance on non-Dirichlet rows (identity rows allowed)
            let mut is_spd = true;
            for i in 0..matrix.nrows() {
                let row = matrix.row(i);
                let mut diag = T::zero();
                let mut sum_off = T::zero();
                for (j, val) in row.col_indices().iter().zip(row.values()) {
                    if *j == i {
                        diag = *val;
                    } else {
                        // Laplacian structure: off-diagonals should be <= 0
                        if *val > T::zero() {
                            is_spd = false;
                            break;
                        }
                        sum_off += val.abs();
                    }
                }
                let is_identity_dirichlet = diag == T::one() && sum_off == T::zero();
                if diag <= T::zero() && !is_identity_dirichlet {
                    is_spd = false;
                    break;
                }
                if !is_identity_dirichlet && diag < sum_off {
                    is_spd = false;
                    break;
                }
            }

            // Choose solver method adaptively
            let mut adaptive_solver = LinearSystemSolver::new();
            let selected_method = if is_spd { LinearSolverMethod::ConjugateGradient } else { LinearSolverMethod::BiCGSTAB };
            adaptive_solver = adaptive_solver.with_method(selected_method);

            // 2. Solve the linearized system
            let solution = adaptive_solver.solve(&matrix, &rhs)?;

            // 3. Check for convergence using the convergence checker
            // Diagnostics: residual norm and solution change
            let residual_norm = {
                let mut norm = T::zero();
                for i in 0..n {
                    let row = matrix.row(i);
                    let mut ax_i = T::zero();
                    for (j, val) in row.col_indices().iter().zip(row.values()) {
                        ax_i += *val * solution[*j];
                    }
                    let r_i = ax_i - rhs[i];
                    norm += r_i * r_i;
                }
                norm.sqrt()
            };
            let rhs_norm = rhs.norm();
            if network.last_solver_method.is_some() {} else { network.last_solver_method = Some(selected_method); }
            network.residuals.push(residual_norm);

            // For non-linear systems, we must check solution change to ensure the non-linear
            // iteration (Picard/Newton) has stabilized. Linear residual check is insufficient
            // because the linear solver minimizes it within the current linearization step.
            // Using has_converged() ensures we check |x_new - x_old|.
            let converged = self
                .convergence
                .has_converged_dual(&solution, &last_solution, residual_norm, rhs_norm)?;

            // Check if solution has converged
            if converged {
                // Apply final solution and return
                self.update_network_solution(&mut network, &solution)?;
                return Ok(network);
            }

            // 4. Update network state for next iteration (zero-copy: swap instead of clone)
            self.update_network_solution(&mut network, &solution)?;
            last_solution = solution;

            // Optional: Add logging for iteration progress (residual)
            let _ = residual_norm;
        }

        // Failed to converge within max iterations
        Err(cfd_core::error::Error::Convergence(
            cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    fn update_network_solution(
        &self,
        network: &mut Network<T>,
        solution: &nalgebra::DVector<T>,
    ) -> Result<()> {
        // Update network pressures and flows from solution vector
        network.update_from_solution(solution)
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> Solver<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;
    type Solution = Network<T>;

    fn solve(&self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_network(problem)
    }

    fn name(&self) -> &'static str {
        "NetworkSolver"
    }
}

impl<T: RealField + Copy> Configurable<T> for NetworkSolver<T> {
    type Config = SolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn set_config(&mut self, config: Self::Config) {
        self.config = config;
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> Validatable<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;

    fn validate_problem(&self, problem: &Self::Problem) -> Result<()> {
        // Validate network has nodes
        if problem.network.node_count() == 0 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Network has no nodes".to_string(),
            ));
        }
        // Validate tolerance
        if self.config.tolerance <= T::zero() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Tolerance must be positive".to_string(),
            ));
        }
        problem.network.validate_coefficients()?;
        Ok(())
    }
}
