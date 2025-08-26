//! Core solver traits following Single Responsibility Principle

use crate::error::Result;
use crate::problem::Problem;
use nalgebra::RealField;

/// Core solver trait following Single Responsibility Principle
/// Focused solely on solving problems
pub trait Solver<T: RealField + Copy>: Send + Sync {
    /// Problem type this solver can handle
    type Problem: Problem<T>;
    /// Solution type produced by this solver
    type Solution;

    /// Solve the given problem
    fn solve(&mut self, problem: &Self::Problem) -> Result<Self::Solution>;

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
    fn validate_problem(&self, problem: &Self::Problem) -> Result<()>;

    /// Check if solver can handle this problem type
    fn can_solve(&self, problem: &Self::Problem) -> bool;
}
