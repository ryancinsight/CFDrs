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
use serde::{Serialize, Deserialize};

/// Solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig<T: RealField> {
    pub tolerance: T,
    pub max_iterations: usize,
}

impl<T: RealField> cfd_core::solver::SolverConfiguration<T> for SolverConfig<T> {
    fn tolerance(&self) -> T {
        self.tolerance
    }
    
    fn max_iterations(&self) -> usize {
        self.max_iterations
    }
    
    fn verbosity(&self) -> u8 {
        0 // Default verbosity
    }
    
    fn parallel(&self) -> bool {
        true // Enable parallel by default
    }
}

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
    /// Create a new network solver with default configuration
    pub fn new() -> Self {
        let config = SolverConfig {
            tolerance: T::from_f64(1e-6).unwrap_or_else(T::one),
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

    fn solve(&mut self, problem: &Self::Problem) -> Result<Self::Solution> {
        self.solve_network(problem)
    }
    
    fn name(&self) -> &str {
        "NetworkSolver"
    }
}

impl<T: RealField> Configurable<T> for NetworkSolver<T> {
    type Config = SolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + FromPrimitive> Validatable<T> for NetworkSolver<T> {
    type Problem = NetworkProblem<T>;

    fn validate_problem(&self, problem: &Self::Problem) -> Result<()> {
        // Validate network has nodes
        if problem.network.node_count() == 0 {
            return Err(cfd_core::Error::InvalidConfiguration(
                "Network has no nodes".to_string()
            ));
        }
        // Validate tolerance
        if self.config.tolerance <= T::zero() {
            return Err(cfd_core::Error::InvalidConfiguration(
                "Tolerance must be positive".to_string()
            ));
        }
        Ok(())
    }

    fn can_handle(&self, _problem: &Self::Problem) -> bool {
        true // Can handle any network problem
    }
}