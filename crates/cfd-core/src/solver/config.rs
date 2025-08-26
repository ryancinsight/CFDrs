//! Solver configuration types following SOLID principles

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
/// Core solver configuration trait
pub trait SolverConfiguration<T: RealField + Copy>: Clone + Send + Sync {
    /// Get maximum iterations
    fn max_iterations(&self) -> usize;
    /// Get convergence tolerance
    fn tolerance(&self) -> T;
    /// Check if preconditioning is enabled
    fn use_preconditioning(&self) -> bool;
}
/// Convergence configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConvergenceConfig<T: RealField + Copy> {
    /// Maximum iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Relative tolerance
    pub relative_tolerance: T,
/// Execution configuration
pub struct ExecutionConfig {
    /// Enable parallel execution
    pub parallel: bool,
    /// Number of threads (None = auto)
    pub num_threads: Option<usize>,
    /// Enable verbose output
    pub verbose: bool,
    /// Save intermediate results
    pub save_intermediate: bool,
/// Numerical configuration
pub struct NumericalConfig<T: RealField + Copy> {
    /// Time step size
    pub dt: T,
    /// CFL number
    pub cfl: T,
    /// Relaxation factor
    pub relaxation: T,
/// Complete solver configuration
pub struct SolverConfig<T: RealField + Copy> {
    /// Convergence settings
    pub convergence: ConvergenceConfig<T>,
    /// Execution settings
    pub execution: ExecutionConfig,
    /// Numerical settings
    pub numerical: NumericalConfig<T>,
impl<T: RealField + FromPrimitive + Copy> Default for SolverConfig<T> {
    fn default() -> Self {
        Self {
            convergence: ConvergenceConfig {
                max_iterations: 1000,
                tolerance: T::from_f64(1e-6).unwrap_or_else(T::zero),
                relative_tolerance: T::from_f64(1e-4).unwrap_or_else(T::zero),
            },
            execution: ExecutionConfig {
                parallel: true,
                num_threads: None,
                verbose: false,
                save_intermediate: false,
            numerical: NumericalConfig {
                dt: T::from_f64(0.01).unwrap_or_else(T::zero),
                cfl: T::from_f64(0.5).unwrap_or_else(T::zero),
                relaxation: T::one(),
        }
    }
impl<T: RealField + Copy> SolverConfig<T> {
    /// Create a builder for configuration
    #[must_use]
    pub fn builder() -> SolverConfigBuilder<T>
    where
        T: FromPrimitive,
    {
        SolverConfigBuilder::new()
    /// Get relaxation factor
    pub fn relaxation_factor(&self) -> T
        T: Copy,
        self.numerical.relaxation
    /// Check if verbose output is enabled
    pub fn verbose(&self) -> bool {
        self.execution.verbose
impl<T: RealField + Copy> SolverConfiguration<T> for SolverConfig<T> {
    fn max_iterations(&self) -> usize {
        self.convergence.max_iterations
    fn tolerance(&self) -> T {
        self.convergence.tolerance
    fn use_preconditioning(&self) -> bool {
        false // Default implementation
/// Linear solver configuration
pub struct LinearSolverConfig<T: RealField + Copy> {
    /// Enable preconditioning
    pub preconditioning: bool,
impl<T: RealField + Copy + num_traits::FromPrimitive> Default for LinearSolverConfig<T> {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).unwrap_or_else(T::zero),
            preconditioning: false,
impl<T: RealField + Copy> LinearSolverConfig<T> {
    /// Get tolerance
    pub fn tolerance(&self) -> T
        self.tolerance
    pub fn builder() -> LinearSolverConfigBuilder<T>
        LinearSolverConfigBuilder::new()
/// Builder for linear solver configuration
pub struct LinearSolverConfigBuilder<T: RealField + Copy> {
    config: LinearSolverConfig<T>,
impl<T: RealField + Copy + FromPrimitive> LinearSolverConfigBuilder<T> {
    /// Create a new builder
    pub fn new() -> Self {
            config: LinearSolverConfig::default(),
    /// Set maximum iterations
    pub fn max_iterations(mut self, max_iter: usize) -> Self {
        self.config.max_iterations = max_iter;
        self
    /// Set tolerance
    pub fn tolerance(mut self, tol: T) -> Self {
        self.config.tolerance = tol;
    pub fn preconditioning(mut self, enable: bool) -> Self {
        self.config.preconditioning = enable;
    /// Build the configuration
    pub fn build(self) -> LinearSolverConfig<T> {
        self.config
/// Network solver configuration
pub struct NetworkSolverConfig<T: RealField + Copy> {
    /// Linear solver settings
    pub linear_solver: LinearSolverConfig<T>,
    /// Network-specific settings
    pub network: NetworkConfig,
/// Network configuration
pub struct NetworkConfig {
    /// Enable parallel assembly
    pub parallel_assembly: bool,
    /// Minimum conductance threshold
    pub min_conductance: f64,
/// Builder for solver configuration
pub struct SolverConfigBuilder<T: RealField + Copy> {
    config: SolverConfig<T>,
impl<T: RealField + Copy> Default for SolverConfigBuilder<T> {
        Self::new()
impl<T: RealField + Copy> SolverConfigBuilder<T> {
    /// Create new builder with defaults
            config: SolverConfig {
                convergence: ConvergenceConfig {
                    max_iterations: 1000,
                    tolerance: T::from_f64(1e-6).unwrap_or_else(T::zero),
                    relative_tolerance: T::from_f64(1e-4).unwrap_or_else(T::zero),
                },
                execution: ExecutionConfig {
                    parallel: true,
                    num_threads: None,
                    verbose: false,
                    save_intermediate: false,
                numerical: NumericalConfig {
                    dt: T::from_f64(0.01).unwrap_or_else(T::zero),
                    cfl: T::from_f64(0.5).unwrap_or_else(T::zero),
                    relaxation: T::one(),
        self.config.convergence.max_iterations = max_iter;
        self.config.convergence.tolerance = tol;
    /// Set time step
    pub fn time_step(mut self, dt: T) -> Self {
        self.config.numerical.dt = dt;
    /// Enable/disable parallel execution
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.config.execution.parallel = parallel;
    /// Set relaxation factor
    pub fn relaxation_factor(mut self, factor: T) -> Self {
        self.config.numerical.relaxation = factor;
    /// Set CFL number
    pub fn cfl(mut self, cfl: T) -> Self {
        self.config.numerical.cfl = cfl;
    /// Enable/disable verbose output
    pub fn verbose(mut self, verbose: bool) -> Self {
        self.config.execution.verbose = verbose;
    pub fn build(self) -> SolverConfig<T> {
