//! Configuration types for linear solvers
//!
//! These are self-contained within cfd-math to avoid dependency inversion.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Configuration for iterative linear solvers
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct IterativeSolverConfig<T: RealField + Copy> {
    /// Maximum number of iterations
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Whether to use preconditioning
    pub use_preconditioner: bool,
    /// Whether to use parallel SpMV operations (rayon-based)
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

impl<T: RealField + Copy + FromPrimitive> Default for IterativeSolverConfig<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            tolerance: T::from_f64(1e-6).expect("Failed to convert default tolerance"),
            use_preconditioner: false,
            use_parallel_spmv: false,
        }
    }
}
