//! Core solver traits following Single Responsibility Principle

use crate::error::Result;
use crate::problem::Problem;
use nalgebra::RealField;

/// Core solver trait following Single Responsibility Principle
/// Focused solely on solving problems
pub trait Solver<T: RealField + Copy>: Send + Sync + Configurable<T> {
    /// Problem type this solver can handle
    type Problem: Problem<T>;
    /// Solution type produced by this solver
    type Solution;

    /// Solve the given problem
    /// Takes &self to enable concurrent execution
    ///
    /// # Errors
    /// Returns an error if the problem cannot be solved due to:
    /// - Invalid problem configuration
    /// - Convergence failure
    /// - Numerical instability
    fn solve(&self, problem: &Self::Problem) -> Result<Self::Solution>;

    /// Get solver name for identification
    fn name(&self) -> &str;
}

/// Configuration management trait following Interface Segregation Principle
/// Separated from core solving functionality
pub trait Configurable<T: RealField + Copy> {
    /// Configuration type for this solver
    type Config: super::config::SolverConfiguration<T>;

    /// Get solver configuration
    fn config(&self) -> &Self::Config;

    /// Set solver configuration
    fn set_config(&mut self, config: Self::Config);
}

/// Validation trait following Single Responsibility Principle
/// Separated validation concerns from solving
pub trait Validatable<T: RealField + Copy> {
    /// Problem type to validate
    type Problem: Problem<T>;

    /// Validate problem before solving
    ///
    /// # Errors
    /// Returns an error if the problem is invalid due to:
    /// - Missing required parameters
    /// - Invalid boundary conditions
    /// - Incompatible mesh and field dimensions
    fn validate_problem(&self, problem: &Self::Problem) -> Result<()>;
}
