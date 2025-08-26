//! Configuration for Finite Difference Method solvers.
//!
//! Provides unified configuration following SSOT principle.

use cfd_core::solver::SolverConfig;
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
impl<T: RealField + Copy + FromPrimitive + Copy> Default for FdmConfig<T> {
    fn default() -> Self {
        Self {
            base: SolverConfig::default(),
        }
    }
impl<T: RealField + Copy> FdmConfig<T> {
    /// Get maximum iterations from base config
    pub fn max_iterations(&self) -> usize {
        self.base.convergence.max_iterations
    /// Get tolerance from base config
    }

    pub fn tolerance(&self) -> T {
        self.base.convergence.tolerance
    /// Get relaxation factor from base config
    }

    pub fn relaxation_factor(&self) -> T {
        self.base.numerical.relaxation
    /// Get verbose flag from base config
    }

    pub fn verbose(&self) -> bool {
        self.base.execution.verbose


    }

}
}
}
}
