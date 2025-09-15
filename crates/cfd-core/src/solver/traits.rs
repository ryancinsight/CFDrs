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
    /// Returns error if problem cannot be solved due to convergence failure,
    /// numerical instability, or invalid problem specification
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
    /// Returns error if problem specification is invalid, has incompatible boundary
    /// conditions, or contains physically impossible configurations
    fn validate_problem(&self, problem: &Self::Problem) -> Result<()>;
}
