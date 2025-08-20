//! Modularized network solver for 1D CFD analysis
//!
//! This module provides a comprehensive solver for analyzing fluid flow in microfluidic
//! networks using sparse linear algebra and circuit analogies.

mod domain;
mod problem;
mod state;
mod matrix_assembly;
mod linear_system;
mod convergence;

pub use domain::NetworkDomain;
pub use problem::NetworkProblem;
pub use state::NetworkState;
pub use matrix_assembly::MatrixAssembler;
pub use linear_system::LinearSystemSolver;
pub use convergence::ConvergenceChecker;



use crate::network::Network;
use cfd_core::{Result, Solver, Configurable, Validatable};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Main network solver implementing the core CFD suite trait system
pub struct NetworkSolver<T: RealField> {
    /// Solver configuration
    config: SolverConfig<T>,
    /// Matrix assembler for building the linear system
    assembler: MatrixAssembler<T>,
    /// Linear system solver
    linear_solver: LinearSystemSolver<T>,
    /// Convergence checker
    convergence: ConvergenceChecker<T>,
}

impl<T: RealField + FromPrimitive> NetworkSolver<T> {
    /// Create a new network solver
    pub fn new(config: SolverConfig<T>) -> Self {
        Self {
            assembler: MatrixAssembler::new(),
            linear_solver: LinearSystemSolver::new(config.base.linear_solver.clone()),
            convergence: ConvergenceChecker::new(config.base.tolerance.clone()),
            config,
        }
    }

    /// Solve the network flow problem
    pub fn solve_network(&self, problem: &NetworkProblem<T>) -> Result<Network<T>> {
        // Build the linear system
        let (matrix, rhs) = self.assembler.assemble(&problem.network)?;
        
        // Solve the linear system
        let solution = self.linear_solver.solve(matrix, rhs)?;
        
        // Check convergence
        self.convergence.check(&solution)?;
        
        // Update network with solution
        let mut network = problem.network.clone();
        self.update_network_solution(&mut network, solution)?;
        
        Ok(network)
    }

    fn update_network_solution(&self, network: &mut Network<T>, solution: nalgebra::DVector<T>) -> Result<()> {
        // Update network pressures and flows from solution vector
        network.update_from_solution(solution)
    }
}

impl<T: RealField + FromPrimitive> Solver<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;
    type Solution = Network<T>;

    fn solve(&self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_network(problem)
    }
}

impl<T: RealField> Configurable for NetworkSolver<T> {
    type Config = SolverConfig<T>;

    fn configure(&mut self, config: Self::Config) -> Result<()> {
        self.config = config.clone();
        self.linear_solver.update_config(&config.base.linear_solver)?;
        self.convergence.update_tolerance(config.base.tolerance);
        Ok(())
    }

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + FromPrimitive> Validatable for NetworkSolver<T> {
    fn validate(&self) -> Result<()> {
        // Validate configuration
        if self.config.base.tolerance <= T::zero() {
            return Err(cfd_core::Error::InvalidConfiguration(
                "Tolerance must be positive".to_string()
            ));
        }
        Ok(())
    }
}