//! Direct solver traits

use super::traits::{Solver, Configurable};
use nalgebra::RealField;

/// Direct solver trait
pub trait DirectSolver<T: RealField>: Solver<T> + Configurable<T> {
    /// Check if the problem size is within limits
    fn can_handle_size(&self, size: usize) -> bool;
    
    /// Get maximum problem size this solver can handle
    fn max_size(&self) -> usize;
}