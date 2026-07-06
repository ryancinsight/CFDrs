//! Configuration types for linear solvers
//!
//! These are self-contained within cfd-math to avoid dependency inversion.

use eunomia::{FloatElement, RealField};
use serde::{Deserialize, Serialize};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Configuration for iterative linear solvers
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct IterativeSolverConfig<T: RealField + Copy> {
    /// Maximum number of iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Whether to use preconditioning
    pub use_preconditioner: bool,
    /// Whether to use parallel SpMV operations (Moirai-backed)
    pub use_parallel_spmv: bool,
}

impl<T: RealField + Copy> IterativeSolverConfig<T> {
    /// Create a new configuration with default values
    pub fn new(tolerance: T) -> Self {
        Self {
            max_iterations: 1000,
            tolerance,
            use_preconditioner: false,
            use_parallel_spmv: false,
        }
    }

    /// Set maximum iterations
    pub fn with_max_iterations(mut self, max: usize) -> Self {
        self.max_iterations = max;
        self
    }

    /// Enable preconditioning
    pub fn with_preconditioner(mut self) -> Self {
        self.use_preconditioner = true;
        self
    }

    /// Enable parallel SpMV operations
    pub fn with_parallel_spmv(mut self) -> Self {
        self.use_parallel_spmv = true;
        self
    }
}

impl<T: RealField + Copy + FloatElement> Default for IterativeSolverConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: from_f64(1e-6),
            use_preconditioner: false,
            use_parallel_spmv: false,
        }
    }
}
