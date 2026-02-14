//! Modularized network solver for 1D CFD analysis
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies.

mod convergence;
mod geometry;
mod linear_system;
mod matrix_assembly;
mod problem;
mod state;
mod transient_composition;
mod transient_droplets;

pub use convergence::ConvergenceChecker;
pub use geometry::NetworkDomain;
pub use linear_system::{LinearSolverMethod, LinearSystemSolver};
pub use matrix_assembly::MatrixAssembler;
pub use problem::NetworkProblem;
pub use state::NetworkState;
pub use transient_composition::{
    CompositionState, EdgeFlowEvent, InletCompositionEvent, MixtureComposition,
    PressureBoundaryEvent, SimulationTimeConfig, TransientCompositionSimulator,
};
pub use transient_droplets::{
    ChannelOccupancy, DropletBoundary, DropletInjection, DropletPosition, DropletSnapshot,
    DropletSplitPolicy, DropletState, DropletTrackingState, SplitMode, TransientDropletSimulator,
};

use crate::network::Network;
use cfd_core::compute::solver::{Configurable, Solver, Validatable};
use cfd_core::error::Result;
use cfd_core::physics::fluid::{ConstantPropertyFluid, FluidTrait};
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
pub struct NetworkSolver<T: RealField + Copy, F: FluidTrait<T> = ConstantPropertyFluid<T>> {
    /// Solver configuration
    config: SolverConfig<T>,
    /// Matrix assembler for building the linear system
    assembler: MatrixAssembler<T>,
    /// Linear system solver
    #[allow(dead_code)]
    linear_solver: LinearSystemSolver<T>,
    /// Convergence checker
    convergence: ConvergenceChecker<T>,
    _phantom: std::marker::PhantomData<F>,
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> Default
    for NetworkSolver<T, F>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> NetworkSolver<T, F> {
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
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create with specific configuration
    pub fn with_config(config: SolverConfig<T>) -> Self {
        Self {
            assembler: MatrixAssembler::new(),
            linear_solver: LinearSystemSolver::new(),
            convergence: ConvergenceChecker::new(config.tolerance),
            config,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Solve the network flow problem iteratively for non-linear systems
    pub fn solve_network(&self, problem: &NetworkProblem<T, F>) -> Result<Network<T, F>> {
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
            let selected_method = if is_spd {
                LinearSolverMethod::ConjugateGradient
            } else {
                LinearSolverMethod::BiCGSTAB
            };
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
            if network.last_solver_method.is_some() {
            } else {
                network.last_solver_method = Some(selected_method);
            }
            network.residuals.push(residual_norm);

            // For non-linear systems, we must check solution change to ensure the non-linear
            // iteration (Picard/Newton) has stabilized. Linear residual check is insufficient
            // because the linear solver minimizes it within the current linearization step.
            // Using has_converged() ensures we check |x_new - x_old|.
            let converged = self.convergence.has_converged_dual(
                &solution,
                &last_solution,
                residual_norm,
                rhs_norm,
            )?;

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
        network: &mut Network<T, F>,
        solution: &nalgebra::DVector<T>,
    ) -> Result<()> {
        // Update network pressures and flows from solution vector
        network.update_from_solution(solution)
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> Solver<T> for NetworkSolver<T, F> {
    type Problem = NetworkProblem<T, F>;
    type Solution = Network<T, F>;

    fn solve(&self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_network(problem)
    }

    fn name(&self) -> &'static str {
        "NetworkSolver"
    }
}

impl<T: RealField + Copy, F: FluidTrait<T>> Configurable<T> for NetworkSolver<T, F> {
    type Config = SolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }

    fn set_config(&mut self, config: Self::Config) {
        self.config = config;
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy, F: FluidTrait<T> + Clone> Validatable<T> for NetworkSolver<T, F> {
    type Problem = NetworkProblem<T, F>;

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
