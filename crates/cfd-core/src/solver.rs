//! Solver trait and common solver functionality.
//!
//! This module provides the core solver abstractions following SOLID principles
//! with zero-copy operations and iterator-based processing.

use crate::{Problem, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use std::fmt::Debug;

/// Main solver trait for CFD simulations following Dependency Inversion Principle
pub trait Solver<T: RealField>: Send + Sync {
    /// Problem type this solver can handle
    type Problem: Problem<T>;
    /// Solution type produced by this solver
    type Solution;
    /// Configuration type for this solver
    type Config: SolverConfiguration<T>;

    /// Solve the given problem
    fn solve(&mut self, problem: &Self::Problem) -> Result<Self::Solution>;

    /// Get solver name
    fn name(&self) -> &str;

    /// Get solver configuration
    fn config(&self) -> &Self::Config;

    /// Set solver configuration
    fn set_config(&mut self, config: Self::Config);
}

/// Base solver configuration trait following Interface Segregation Principle
pub trait SolverConfiguration<T: RealField>: Clone + Send + Sync {
    /// Get convergence tolerance
    fn tolerance(&self) -> T;
    /// Get maximum iterations
    fn max_iterations(&self) -> usize;
    /// Get verbosity level
    fn verbosity(&self) -> u8;
    /// Check if parallel execution is enabled
    fn parallel(&self) -> bool;
}

/// Core solver configuration components using composition pattern
#[derive(Debug, Clone)]
pub struct ConvergenceConfig<T: RealField> {
    /// Convergence tolerance
    pub tolerance: T,
    /// Maximum number of iterations
    pub max_iterations: usize,
}

/// Execution configuration for parallel processing
#[derive(Debug, Clone)]
pub struct ExecutionConfig {
    /// Use parallel execution
    pub parallel: bool,
    /// Number of threads (None = use all available)
    pub num_threads: Option<usize>,
    /// Verbosity level (0 = silent, 1 = summary, 2 = detailed)
    pub verbosity: u8,
}

/// Numerical method configuration
#[derive(Debug, Clone)]
pub struct NumericalConfig<T: RealField> {
    /// Relaxation factor for iterative methods
    pub relaxation_factor: T,
    /// Under-relaxation factor (alternative name for compatibility)
    pub under_relaxation: Option<T>,
}

/// Unified solver configuration using composition over inheritance
#[derive(Debug, Clone)]
pub struct SolverConfig<T: RealField> {
    /// Convergence parameters
    pub convergence: ConvergenceConfig<T>,
    /// Execution parameters
    pub execution: ExecutionConfig,
    /// Numerical method parameters
    pub numerical: NumericalConfig<T>,
}

impl<T: RealField> SolverConfiguration<T> for SolverConfig<T> {
    fn tolerance(&self) -> T {
        self.convergence.tolerance.clone()
    }

    fn max_iterations(&self) -> usize {
        self.convergence.max_iterations
    }

    fn verbosity(&self) -> u8 {
        self.execution.verbosity
    }

    fn parallel(&self) -> bool {
        self.execution.parallel
    }
}

impl<T: RealField + FromPrimitive> Default for ConvergenceConfig<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-6).unwrap(),
            max_iterations: 1000,
        }
    }
}

impl Default for ExecutionConfig {
    fn default() -> Self {
        Self {
            parallel: true,
            num_threads: None,
            verbosity: 1,
        }
    }
}

impl<T: RealField> Default for NumericalConfig<T> {
    fn default() -> Self {
        Self {
            relaxation_factor: T::one(),
            under_relaxation: None,
        }
    }
}

impl<T: RealField + FromPrimitive> Default for SolverConfig<T> {
    fn default() -> Self {
        Self {
            convergence: ConvergenceConfig::default(),
            execution: ExecutionConfig::default(),
            numerical: NumericalConfig::default(),
        }
    }
}

/// Linear solver specific configuration using composition
#[derive(Debug, Clone)]
pub struct LinearSolverConfig<T: RealField> {
    /// Base solver configuration
    pub base: SolverConfig<T>,
    /// Restart parameter for GMRES
    pub restart: usize,
    /// Use preconditioning
    pub use_preconditioner: bool,
}

impl<T: RealField> SolverConfiguration<T> for LinearSolverConfig<T> {
    fn tolerance(&self) -> T {
        self.base.tolerance()
    }

    fn max_iterations(&self) -> usize {
        self.base.max_iterations()
    }

    fn verbosity(&self) -> u8 {
        self.base.verbosity()
    }

    fn parallel(&self) -> bool {
        self.base.parallel()
    }
}

impl<T: RealField> LinearSolverConfig<T> {
    /// Get restart parameter
    pub fn restart(&self) -> usize {
        self.restart
    }

    /// Check if preconditioning is enabled
    pub fn use_preconditioner(&self) -> bool {
        self.use_preconditioner
    }

    /// Get relaxation factor from base configuration
    pub fn relaxation_factor(&self) -> T {
        self.base.numerical.relaxation_factor.clone()
    }
}

impl<T: RealField + FromPrimitive> Default for LinearSolverConfig<T> {
    fn default() -> Self {
        Self {
            base: SolverConfig::default(),
            restart: 30,
            use_preconditioner: false,
        }
    }
}

/// Network solver specific configuration using composition
#[derive(Debug, Clone)]
pub struct NetworkSolverConfig<T: RealField> {
    /// Base solver configuration
    pub base: SolverConfig<T>,
    /// Enable verbose output (overrides base verbosity)
    pub verbose: bool,
}

impl<T: RealField> SolverConfiguration<T> for NetworkSolverConfig<T> {
    fn tolerance(&self) -> T {
        self.base.tolerance()
    }

    fn max_iterations(&self) -> usize {
        self.base.max_iterations()
    }

    fn verbosity(&self) -> u8 {
        if self.verbose { 2 } else { self.base.verbosity() }
    }

    fn parallel(&self) -> bool {
        self.base.parallel()
    }
}

impl<T: RealField> NetworkSolverConfig<T> {
    /// Get relaxation factor from base configuration
    pub fn relaxation_factor(&self) -> T {
        self.base.numerical.relaxation_factor.clone()
    }

    /// Get under-relaxation factor if set
    pub fn under_relaxation(&self) -> Option<T> {
        self.base.numerical.under_relaxation.clone()
    }
}

impl<T: RealField + FromPrimitive> Default for NetworkSolverConfig<T> {
    fn default() -> Self {
        Self {
            base: SolverConfig::default(),
            verbose: false,
        }
    }
}

impl<T: RealField> SolverConfig<T> {
    /// Builder pattern for configuration
    pub fn builder() -> SolverConfigBuilder<T> {
        SolverConfigBuilder::default()
    }

    /// Get relaxation factor
    pub fn relaxation_factor(&self) -> T {
        self.numerical.relaxation_factor.clone()
    }

    /// Get under-relaxation factor if set
    pub fn under_relaxation(&self) -> Option<T> {
        self.numerical.under_relaxation.clone()
    }

    /// Get number of threads
    pub fn num_threads(&self) -> Option<usize> {
        self.execution.num_threads
    }
}

/// Enhanced builder for unified solver configuration using fluent interface
#[derive(Debug, Clone)]
pub struct SolverConfigBuilder<T: RealField> {
    convergence: ConvergenceConfig<T>,
    execution: ExecutionConfig,
    numerical: NumericalConfig<T>,
}

impl<T: RealField + FromPrimitive> Default for SolverConfigBuilder<T> {
    fn default() -> Self {
        Self {
            convergence: ConvergenceConfig::default(),
            execution: ExecutionConfig::default(),
            numerical: NumericalConfig::default(),
        }
    }
}

impl<T: RealField> SolverConfigBuilder<T> {
    /// Set convergence tolerance
    pub fn tolerance(mut self, tolerance: T) -> Self {
        self.convergence.tolerance = tolerance;
        self
    }

    /// Set maximum iterations
    pub fn max_iterations(mut self, max_iterations: usize) -> Self {
        self.convergence.max_iterations = max_iterations;
        self
    }

    /// Set relaxation factor
    pub fn relaxation_factor(mut self, factor: T) -> Self {
        self.numerical.relaxation_factor = factor;
        self
    }

    /// Set under-relaxation factor
    pub fn under_relaxation(mut self, factor: T) -> Self {
        self.numerical.under_relaxation = Some(factor);
        self
    }

    /// Set verbosity level
    pub fn verbosity(mut self, verbosity: u8) -> Self {
        self.execution.verbosity = verbosity;
        self
    }

    /// Enable/disable parallel execution
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.execution.parallel = parallel;
        self
    }

    /// Set number of threads
    pub fn num_threads(mut self, threads: Option<usize>) -> Self {
        self.execution.num_threads = threads;
        self
    }

    /// Build the unified configuration
    pub fn build(self) -> SolverConfig<T> {
        SolverConfig {
            convergence: self.convergence,
            execution: self.execution,
            numerical: self.numerical,
        }
    }

    /// Build a linear solver configuration with additional parameters
    pub fn build_linear(self, restart: usize, use_preconditioner: bool) -> LinearSolverConfig<T> {
        LinearSolverConfig {
            base: self.build(),
            restart,
            use_preconditioner,
        }
    }

    /// Build a network solver configuration with verbose option
    pub fn build_network(self, verbose: bool) -> NetworkSolverConfig<T> {
        NetworkSolverConfig {
            base: self.build(),
            verbose,
        }
    }
}

/// Iterative solver trait with iterator-based convergence
pub trait IterativeSolver<T: RealField>: Solver<T> 
where
    Self::Solution: Clone,
{
    /// Iterator over solver iterations
    type IterationState: IterationState<T, Solution = Self::Solution>;

    /// Create an iterator over solver iterations
    fn iterations(
        &self,
        initial_state: Self::Solution,
    ) -> SolverIterator<Self::IterationState, Self::Solution, T> {
        SolverIterator::new(self.create_iteration_state(initial_state))
    }

    /// Create iteration state from initial solution
    fn create_iteration_state(&self, initial: Self::Solution) -> Self::IterationState;

    /// Enhanced iterative solve using iterator combinators and functional patterns
    fn iterative_solve(&mut self, initial: Self::Solution) -> Result<Self::Solution> {
        let config = self.config();
        let max_iterations = config.max_iterations();
        let tolerance = config.tolerance();
        let verbosity = config.verbosity();

        // Use iterator-based approach with scan for stateful iteration
        let result = (0..max_iterations)
            .scan(&mut initial, |current_solution, iter| {
                let mut state = self.create_iteration_state(current_solution);
                match state.iterate() {
                    Ok((solution, residual)) => {
                        if verbosity >= 2 {
                            tracing::debug!("Iteration {}: residual = {:?}", iter + 1, residual);
                        }

                        *current_solution = solution;
                        Some(Ok((iter + 1, current_solution.clone(), residual)))
                    }
                    Err(e) => Some(Err(e)),
                }
            })
            .find_map(|result| match result {
                Ok((iter, solution, residual)) => {
                    if residual < tolerance {
                        if verbosity >= 1 {
                            tracing::info!(
                                "Converged after {} iterations (residual: {:?})",
                                iter,
                                residual
                            );
                        }
                        Some(Ok(solution))
                    } else {
                        None
                    }
                }
                Err(e) => Some(Err(e)),
            });

        match result {
            Some(result) => result,
            None => Err(crate::Error::ConvergenceFailure(format!(
                "Failed to converge after {} iterations",
                max_iterations
            ))),
        }
    }
}

/// State for solver iterations
pub trait IterationState<T: RealField> {
    /// Solution type
    type Solution;

    /// Perform one iteration and return new solution with residual
    fn iterate(&mut self) -> Result<(Self::Solution, T)>;
}

/// Iterator over solver iterations
pub struct SolverIterator<S, Solution, T> {
    state: S,
    _phantom: std::marker::PhantomData<(Solution, T)>,
}

impl<S, Solution, T> SolverIterator<S, Solution, T> {
    /// Create new solver iterator
    pub fn new(state: S) -> Self {
        Self {
            state,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<S, Solution, T> Iterator for SolverIterator<S, Solution, T>
where
    S: IterationState<T, Solution = Solution>,
    T: RealField,
{
    type Item = Result<(Solution, T)>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.state.iterate())
    }
}

/// Direct solver trait with lazy evaluation
pub trait DirectSolver<T: RealField>: Solver<T> {
    /// Matrix assembly iterator
    type MatrixIterator: Iterator<Item = MatrixEntry<T>>;

    /// Create iterator for matrix assembly
    fn matrix_entries(&self, problem: &Self::Problem) -> Self::MatrixIterator;

    /// Assemble the system matrix using iterator
    fn assemble(&self, problem: &Self::Problem) -> Result<SystemMatrix<T>> {
        let entries: Vec<_> = self.matrix_entries(problem).collect();
        
        // Use parallel assembly if configured
        if self.config().parallel() {
            self.parallel_assemble(entries)
        } else {
            self.sequential_assemble(entries)
        }
    }

    /// Sequential matrix assembly
    fn sequential_assemble(&self, entries: Vec<MatrixEntry<T>>) -> Result<SystemMatrix<T>>;

    /// Parallel matrix assembly using rayon
    fn parallel_assemble(&self, entries: Vec<MatrixEntry<T>>) -> Result<SystemMatrix<T>> {
        use rayon::prelude::*;

        // Group entries by row for efficient parallel assembly
        let mut row_groups: std::collections::HashMap<usize, Vec<MatrixEntry<T>>> =
            std::collections::HashMap::new();

        entries.into_iter().for_each(|entry| {
            row_groups.entry(entry.row).or_default().push(entry);
        });

        // Process rows in parallel
        let processed_entries: Vec<_> = row_groups
            .into_par_iter()
            .flat_map(|(_, row_entries)| {
                // Sort entries within each row by column index
                let mut sorted_entries = row_entries;
                sorted_entries.sort_by_key(|entry| entry.col);
                sorted_entries
            })
            .collect();

        self.sequential_assemble(processed_entries)
    }

    /// Solve the linear system
    fn solve_system(&self, system: &SystemMatrix<T>) -> Result<Self::Solution>;
}

/// Matrix entry for sparse assembly
#[derive(Debug, Clone)]
pub struct MatrixEntry<T> {
    /// Row index
    pub row: usize,
    /// Column index
    pub col: usize,
    /// Value
    pub value: T,
}

/// System matrix representation
pub struct SystemMatrix<T: RealField> {
    /// Sparse matrix
    pub matrix: nalgebra_sparse::CsrMatrix<T>,
    /// Right-hand side vector
    pub rhs: nalgebra::DVector<T>,
}

/// Solution monitoring trait with iterator support
pub trait SolutionMonitor<T: RealField> {
    /// Solution type to monitor
    type Solution;

    /// Create an iterator that monitors solutions
    fn monitor_iter<I>(self, iter: I) -> MonitoredIterator<Self, I>
    where
        I: Iterator<Item = (usize, Self::Solution, T)>,
        Self: Sized,
    {
        MonitoredIterator::new(self, iter)
    }

    /// Called at each iteration
    fn on_iteration(&mut self, iteration: usize, solution: &Self::Solution, residual: T);

    /// Called when converged
    fn on_converged(&mut self, iterations: usize, solution: &Self::Solution, residual: T);

    /// Called on failure
    fn on_failed(&mut self, iterations: usize, solution: &Self::Solution, residual: T);
}

/// Iterator wrapper that monitors solutions
pub struct MonitoredIterator<M, I> {
    monitor: M,
    iter: I,
}

impl<M, I> MonitoredIterator<M, I> {
    /// Create new monitored iterator
    pub fn new(monitor: M, iter: I) -> Self {
        Self { monitor, iter }
    }
}

impl<M, I, S, T> Iterator for MonitoredIterator<M, I>
where
    M: SolutionMonitor<T, Solution = S>,
    I: Iterator<Item = (usize, S, T)>,
    T: RealField,
{
    type Item = (usize, S, T);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(iter, solution, residual)| {
            self.monitor.on_iteration(iter, &solution, residual.clone());
            (iter, solution, residual)
        })
    }
}

/// Default monitor that does nothing
pub struct NullMonitor;

impl<T: RealField> SolutionMonitor<T> for NullMonitor {
    type Solution = ();

    fn on_iteration(&mut self, _: usize, _: &Self::Solution, _: T) {}
    fn on_converged(&mut self, _: usize, _: &Self::Solution, _: T) {}
    fn on_failed(&mut self, _: usize, _: &Self::Solution, _: T) {}
}

/// Convergence criteria trait
pub trait ConvergenceCriteria<T: RealField> {
    /// Check if converged
    fn is_converged(&self, residual: T, iteration: usize) -> bool;

    /// Chain with another criteria (AND operation)
    fn and<C: ConvergenceCriteria<T>>(self, other: C) -> AndCriteria<Self, C>
    where
        Self: Sized,
    {
        AndCriteria::new(self, other)
    }

    /// Chain with another criteria (OR operation)
    fn or<C: ConvergenceCriteria<T>>(self, other: C) -> OrCriteria<Self, C>
    where
        Self: Sized,
    {
        OrCriteria::new(self, other)
    }
}

/// Tolerance-based convergence
pub struct ToleranceCriteria<T> {
    tolerance: T,
}

impl<T: RealField> ToleranceCriteria<T> {
    /// Create new tolerance criteria
    pub fn new(tolerance: T) -> Self {
        Self { tolerance }
    }
}

impl<T: RealField> ConvergenceCriteria<T> for ToleranceCriteria<T> {
    fn is_converged(&self, residual: T, _: usize) -> bool {
        residual < self.tolerance
    }
}

/// AND combination of criteria
pub struct AndCriteria<C1, C2> {
    criteria1: C1,
    criteria2: C2,
}

impl<C1, C2> AndCriteria<C1, C2> {
    /// Create new AND criteria
    pub fn new(criteria1: C1, criteria2: C2) -> Self {
        Self { criteria1, criteria2 }
    }
}

impl<T: RealField, C1, C2> ConvergenceCriteria<T> for AndCriteria<C1, C2>
where
    C1: ConvergenceCriteria<T>,
    C2: ConvergenceCriteria<T>,
{
    fn is_converged(&self, residual: T, iteration: usize) -> bool {
        self.criteria1.is_converged(residual.clone(), iteration)
            && self.criteria2.is_converged(residual, iteration)
    }
}

/// OR combination of criteria
pub struct OrCriteria<C1, C2> {
    criteria1: C1,
    criteria2: C2,
}

impl<C1, C2> OrCriteria<C1, C2> {
    /// Create new OR criteria
    pub fn new(criteria1: C1, criteria2: C2) -> Self {
        Self { criteria1, criteria2 }
    }
}

impl<T: RealField, C1, C2> ConvergenceCriteria<T> for OrCriteria<C1, C2>
where
    C1: ConvergenceCriteria<T>,
    C2: ConvergenceCriteria<T>,
{
    fn is_converged(&self, residual: T, iteration: usize) -> bool {
        self.criteria1.is_converged(residual.clone(), iteration)
            || self.criteria2.is_converged(residual, iteration)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solver_config_builder() {
        let config = SolverConfig::<f64>::builder()
            .tolerance(1e-8)
            .max_iterations(500)
            .relaxation_factor(0.9)
            .verbosity(2)
            .parallel(false)
            .build();

        assert_eq!(config.convergence.tolerance, 1e-8);
        assert_eq!(config.convergence.max_iterations, 500);
        assert_eq!(config.numerical.relaxation_factor, 0.9);
        assert_eq!(config.execution.verbosity, 2);
        assert!(!config.execution.parallel);
    }

    #[test]
    fn test_unified_config_composition() {
        let linear_config = SolverConfig::<f64>::builder()
            .tolerance(1e-10)
            .max_iterations(2000)
            .under_relaxation(0.8)
            .parallel(true)
            .build_linear(50, true);

        assert_eq!(linear_config.base.tolerance(), 1e-10);
        assert_eq!(linear_config.restart, 50);
        assert!(linear_config.use_preconditioner);

        let network_config = SolverConfig::<f64>::builder()
            .tolerance(1e-6)
            .relaxation_factor(0.9)
            .build_network(true);

        assert_eq!(network_config.verbosity(), 2); // verbose overrides
        assert_eq!(network_config.relaxation_factor(), 0.9);
    }

    #[test]
    fn test_convergence_criteria() {
        let tol_criteria = ToleranceCriteria::new(1e-6);
        assert!(tol_criteria.is_converged(1e-7, 10));
        assert!(!tol_criteria.is_converged(1e-5, 10));

        // Test AND criteria would require implementing ConvergenceCriteria for closures
        // For now, we test the basic tolerance criteria
        // let max_iter_criteria = |_: f64, iter: usize| iter < 100;
        // let combined = tol_criteria.and(max_iter_criteria);
    }
}