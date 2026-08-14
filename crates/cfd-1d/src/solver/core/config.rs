use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};

/// Solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SolverConfig<T: CfdScalar + Copy> {
    /// Convergence tolerance for solution accuracy
    pub tolerance: T,
    /// Maximum number of solver iterations before termination
    pub max_iterations: usize,
    /// When `true` (the default), the outer Picard loop also requires the flow
    /// rates to change by less than `tolerance` between consecutive iterates.
    ///
    /// Set to `false` for the reference-trace path where **pressure convergence
    /// plus a tiny linear-system residual is sufficient** — the flow rates are
    /// uniquely derived from the converged pressure and do not need an
    /// independent secondary criterion.  This is particularly important for
    /// networks with venturi channels whose quadratic loss terms create a
    /// slowly converging secondary fixed-point that the pressure already
    /// satisfies in a handful of Picard steps.
    #[serde(default = "default_require_flow_convergence")]
    pub require_flow_convergence: bool,
}

fn default_require_flow_convergence() -> bool {
    true
}

impl<T: CfdScalar + Copy> cfd_core::compute::solver::SolverConfiguration<T> for SolverConfig<T> {
    fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }

    fn use_preconditioning(&self) -> bool {
        false // No preconditioning for network solver
    }
}
