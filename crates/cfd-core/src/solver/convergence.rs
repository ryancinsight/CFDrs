//! Convergence criteria for iterative solvers

use nalgebra::RealField;
use std::marker::PhantomData;
/// Convergence criteria trait
pub trait ConvergenceCriteria<T: RealField + Copy>: Send + Sync {
    /// Check if converged
    fn is_converged(&self, iteration: usize, residual: T, initial_residual: T) -> bool;
    /// Get description
    fn description(&self) -> String;
}
/// Tolerance-based convergence criteria
pub struct ToleranceCriteria<T> {
    /// Absolute tolerance
    pub absolute_tolerance: T,
    /// Relative tolerance
    pub relative_tolerance: T,
    /// Maximum iterations
    pub max_iterations: usize,
impl<T: RealField + Copy> ConvergenceCriteria<T> for ToleranceCriteria<T> {
    fn is_converged(&self, iteration: usize, residual: T, initial_residual: T) -> bool {
        iteration >= self.max_iterations
            || residual <= self.absolute_tolerance
            || residual <= initial_residual * self.relative_tolerance
    }
    fn description(&self) -> String {
        format!("Tolerance criteria: max_iter={}", self.max_iterations)
/// Combine two criteria with AND logic
pub struct AndCriteria<C1, C2, T> {
    criterion1: C1,
    criterion2: C2,
    _phantom: PhantomData<T>,
impl<C1, C2, T> AndCriteria<C1, C2, T>
where
    C1: ConvergenceCriteria<T>,
    C2: ConvergenceCriteria<T>,
    T: RealField + Copy,
{
    /// Create new AND criteria
    pub fn new(criterion1: C1, criterion2: C2) -> Self {
        Self {
            criterion1,
            criterion2,
            _phantom: PhantomData,
        }
impl<C1, C2, T> ConvergenceCriteria<T> for AndCriteria<C1, C2, T>
        self.criterion1
            .is_converged(iteration, residual, initial_residual)
            && self
                .criterion2
                .is_converged(iteration, residual, initial_residual)
        format!(
            "({}) AND ({})",
            self.criterion1.description(),
            self.criterion2.description()
        )
/// Combine two criteria with OR logic
pub struct OrCriteria<C1, C2, T> {
impl<C1, C2, T> OrCriteria<C1, C2, T>
    /// Create new OR criteria
impl<C1, C2, T> ConvergenceCriteria<T> for OrCriteria<C1, C2, T>
            || self
            "({}) OR ({})",
