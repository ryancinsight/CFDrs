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

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

// Named constants for FVM
const DEFAULT_CONVERGENCE_TOLERANCE: f64 = 1e-6;
const DEFAULT_MAX_ITERATIONS: usize = 1000;
const DEFAULT_CFL_NUMBER: f64 = 0.5;
const DEFAULT_RELAXATION_FACTOR: f64 = 0.7;
const DEFAULT_DIFFUSION_COEFFICIENT: f64 = 1e-3;

/// Configuration for FVM solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FvmConfig<T: RealField + Copy> {
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

impl<T: RealField + Copy + num_traits::FromPrimitive> Default for FvmConfig<T> {
    fn default() -> Self {
        Self {
            nx: 100,
            ny: 100,
            dx: T::from_f64(0.01).expect("analytical constant conversion"),
            dy: T::from_f64(0.01).expect("analytical constant conversion"),
            dt: T::from_f64(0.001).expect("analytical constant conversion"),
            convergence_tolerance: T::from_f64(DEFAULT_CONVERGENCE_TOLERANCE)
                .expect("analytical constant conversion"),
            max_iterations: DEFAULT_MAX_ITERATIONS,
            cfl_number: T::from_f64(DEFAULT_CFL_NUMBER).expect("analytical constant conversion"),
            relaxation_factor: T::from_f64(DEFAULT_RELAXATION_FACTOR).expect("analytical constant conversion"),
            diffusion_coefficient: T::from_f64(DEFAULT_DIFFUSION_COEFFICIENT)
                .expect("analytical constant conversion"),
        }
    }
}
