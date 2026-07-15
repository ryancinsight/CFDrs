//! Configuration for FVM solver
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

// Named constants for FVM
const DEFAULT_CONVERGENCE_TOLERANCE: f64 = 1e-6;
const DEFAULT_MAX_ITERATIONS: usize = 1000;
const DEFAULT_CFL_NUMBER: f64 = 0.5;
const DEFAULT_RELAXATION_FACTOR: f64 = 0.7;
const DEFAULT_DIFFUSION_COEFFICIENT: f64 = 1e-3;

/// Configuration for FVM solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FvmConfig<T: FloatElement + Copy> {
    /// Grid size in x-direction
    pub nx: usize,
    /// Grid size in y-direction
    pub ny: usize,
    /// Grid spacing in x-direction
    pub dx: T,
    /// Grid spacing in y-direction
    pub dy: T,
    /// Time step
    pub dt: T,
    /// Convergence tolerance
    pub convergence_tolerance: T,
    /// Maximum iterations
    pub max_iterations: usize,
    /// CFL number for time stepping
    pub cfl_number: T,
    /// Under-relaxation factor
    pub relaxation_factor: T,
    /// Diffusion coefficient
    pub diffusion_coefficient: T,
}

impl<T: FloatElement + Copy> Default for FvmConfig<T> {
    fn default() -> Self {
        Self {
            nx: 100,
            ny: 100,
            dx: <T as FloatElement>::from_f64(0.01),
            dy: <T as FloatElement>::from_f64(0.01),
            dt: <T as FloatElement>::from_f64(0.001),
            convergence_tolerance: <T as FloatElement>::from_f64(DEFAULT_CONVERGENCE_TOLERANCE),
            max_iterations: DEFAULT_MAX_ITERATIONS,
            cfl_number: <T as FloatElement>::from_f64(DEFAULT_CFL_NUMBER),
            relaxation_factor: <T as FloatElement>::from_f64(DEFAULT_RELAXATION_FACTOR),
            diffusion_coefficient: <T as FloatElement>::from_f64(DEFAULT_DIFFUSION_COEFFICIENT),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_fvm_config_preserves_physical_constants() {
        let config = FvmConfig::<f64>::default();

        assert_eq!(config.nx, 100);
        assert_eq!(config.ny, 100);
        assert_eq!(config.dx, 0.01);
        assert_eq!(config.dy, 0.01);
        assert_eq!(config.dt, 0.001);
        assert_eq!(config.convergence_tolerance, DEFAULT_CONVERGENCE_TOLERANCE);
        assert_eq!(config.max_iterations, DEFAULT_MAX_ITERATIONS);
        assert_eq!(config.cfl_number, DEFAULT_CFL_NUMBER);
        assert_eq!(config.relaxation_factor, DEFAULT_RELAXATION_FACTOR);
        assert_eq!(config.diffusion_coefficient, DEFAULT_DIFFUSION_COEFFICIENT);
    }
}
