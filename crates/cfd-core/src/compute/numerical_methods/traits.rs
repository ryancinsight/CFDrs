//! Core traits for numerical methods

use nalgebra::{DMatrix, DVector, RealField};

/// Discretization scheme abstraction following Strategy pattern
pub trait DiscretizationScheme<T: RealField + Copy>: Send + Sync {
    /// Apply discretization to a field
    fn discretize(&self, field: &[T], grid_spacing: T) -> Vec<T>;

    /// Get scheme name
    fn name(&self) -> &str;

    /// Get scheme order of accuracy
    fn order(&self) -> usize;
}

/// Time integration scheme abstraction
pub trait TimeIntegrationScheme<T: RealField + Copy>: Send + Sync {
    /// Advance solution in time
    fn advance(&self, current: &[T], derivative: &[T], dt: T) -> Vec<T>;

    /// Get scheme name
    fn name(&self) -> &str;

    /// Get scheme order of accuracy
    fn order(&self) -> usize;

    /// Check if scheme is implicit
    fn is_implicit(&self) -> bool;
}

/// Linear system solver abstraction
pub trait LinearSystemSolver<T: RealField + Copy>: Send + Sync {
    /// Solve linear system Ax = b
    fn solve(&self, a: &DMatrix<T>, b: &DVector<T>) -> Option<DVector<T>>;

    /// Get solver name
    fn name(&self) -> &str;

    /// Check if solver is iterative
    fn is_iterative(&self) -> bool;
}
