//! Configuration for Finite Difference Method solvers.
//!
//! Provides unified configuration following SSOT principle.
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

use cfd_core::compute::solver::SolverConfig;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Finite Difference Method solver configuration
/// Uses unified `SolverConfig` as base to follow SSOT principle
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FdmConfig<T: RealField + Copy> {
    /// Base solver configuration (SSOT)
    pub base: SolverConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> Default for FdmConfig<T> {
    fn default() -> Self {
        Self {
            base: SolverConfig::default(),
        }
    }
}

impl<T: RealField + Copy> FdmConfig<T> {
    /// Get maximum iterations from base config
    pub fn max_iterations(&self) -> usize {
        self.base.convergence.max_iterations
    }

    /// Get tolerance from base config
    pub fn tolerance(&self) -> T {
        self.base.convergence.tolerance
    }

    /// Get relaxation factor from base config
    pub fn relaxation_factor(&self) -> T {
        self.base.numerical.relaxation
    }

    /// Get verbose flag from base config
    pub fn verbose(&self) -> bool {
        self.base.execution.verbose
    }
}
