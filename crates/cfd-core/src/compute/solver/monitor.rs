//! Solution monitoring for iterative solvers

use nalgebra::RealField;
use std::marker::PhantomData;

/// Solution monitor trait
pub trait SolutionMonitor<T: RealField + Copy>: Send + Sync {
    /// Called at each iteration
    fn on_iteration(&mut self, iteration: usize, residual: T);

    /// Called when converged
    fn on_convergence(&mut self, iteration: usize, residual: T);

    /// Called on error
    fn on_error(&mut self, error: &str);
}

/// Iterator wrapper with monitoring
pub struct MonitoredIterator<M, I> {
    monitor: M,
    iterator: I,
    iteration: usize,
}

impl<M, I> MonitoredIterator<M, I> {
    /// Create new monitored iterator
    pub fn new(monitor: M, iterator: I) -> Self {
        Self {
            monitor,
            iterator,
            iteration: 0,
        }
    }
}

impl<M, I, T> Iterator for MonitoredIterator<M, I>
where
    M: SolutionMonitor<T>,
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        match self.iterator.next() {
            Some(value) => {
                self.iteration += 1;
                self.monitor.on_iteration(self.iteration, value);
                Some(value)
            }
            None => None,
        }
    }
}

/// Null monitor (no-op)
pub struct NullMonitor<T> {
    _phantom: PhantomData<T>,
}

impl<T: RealField + Copy> Default for NullMonitor<T> {
    fn default() -> Self {
        Self {
            _phantom: PhantomData,
        }
    }
}

impl<T: RealField + Copy> SolutionMonitor<T> for NullMonitor<T> {
    fn on_iteration(&mut self, _: usize, _: T) {}
    fn on_convergence(&mut self, _: usize, _: T) {}
    fn on_error(&mut self, _: &str) {}
}
