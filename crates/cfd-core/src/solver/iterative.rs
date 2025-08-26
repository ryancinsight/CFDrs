//! Iterative solver traits and implementations

use super::traits::{Configurable, Solver};
use crate::error::Result;
use nalgebra::RealField;
use std::marker::PhantomData;
/// Iterative solver trait
pub trait IterativeSolver<T: RealField + Copy>: Solver<T> + Configurable<T> {
    /// Get current iteration count
    fn iteration_count(&self) -> usize;
    /// Get current residual
    fn residual(&self) -> T;
    /// Check if converged
    fn is_converged(&self) -> bool;
    /// Perform single iteration
    fn iterate(&mut self) -> Result<()>;
}
/// State for iterative solvers
pub trait IterationState<T: RealField + Copy> {
    /// Get iteration number
    fn iteration(&self) -> usize;
    /// Check convergence
    fn converged(&self) -> bool;
/// Iterator adapter for iterative solvers
pub struct SolverIterator<S, Solution, T> {
    solver: S,
    converged: bool,
    _phantom: PhantomData<(Solution, T)>,
impl<S, Solution, T> SolverIterator<S, Solution, T>
where
    S: IterativeSolver<T>,
    T: RealField + Copy,
{
    /// Create new solver iterator
    pub fn new(solver: S) -> Self {
        Self {
            solver,
            converged: false,
            _phantom: PhantomData,
        }
    }
impl<S, Solution, T> Iterator for SolverIterator<S, Solution, T>
    type Item = Result<T>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.converged || self.solver.is_converged() {
            return None;
        match self.solver.iterate() {
            Ok(()) => {
                self.converged = self.solver.is_converged();
                Some(Ok(self.solver.residual()))
            }
            Err(e) => {
                self.converged = true;
                Some(Err(e))
