use super::LinearSolverMethod;
use cfd_core::error::Error;
use serde::{Deserialize, Serialize};
use std::error::Error as StdError;
use std::fmt;

/// Stable classification for why the primary nonlinear solve left the trusted path.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolveFailureReason {
    /// Picard / fixed-point iteration exhausted the configured iteration budget.
    MaxIterationsExceeded,
    /// A residual, solution vector, or assembled coefficient became non-finite.
    NonFiniteResidual,
    /// Matrix assembly or coefficient validation failed before a linear solve could proceed.
    MatrixAssemblyInvalid,
    /// The inner linear solve failed for the assembled system.
    LinearSolverFailure,
    /// Linearized recovery solved a normalized system but could not recover nonzero inlet flow.
    ZeroInletFlowAfterRecovery,
    /// Required geometry or coefficient contract was invalid before solving.
    InvalidGeometryContract,
}

/// Public solve-path status contract used by diagnostics, auditing, and reporting.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolvePathStatus {
    /// The primary nonlinear or linear-static path converged without degrading the model.
    PrimaryConverged,
    /// The primary path failed, but a linearized diagnostic recovery path produced a solution.
    RecoveredLinearized {
        /// Stable classification for why the trusted primary path had to be abandoned.
        reason: SolveFailureReason,
    },
    /// Neither the primary path nor any recovery path produced a trusted solution.
    Failed {
        /// Stable classification for why the solve could not produce a trusted result.
        reason: SolveFailureReason,
    },
}

/// Diagnostics captured from the primary `cfd-1d` solve attempt.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[derive(Default)]
pub struct PrimarySolveDiagnostics {
    /// Number of Picard iterations executed on the primary path.
    pub picard_iterations: usize,
    /// Selected inner linear solver method, when matrix assembly succeeded.
    pub linear_solver_method: Option<LinearSolverMethod>,
    /// Last computed residual norm `||Ax - b||₂`, when available.
    pub last_residual_norm: Option<f64>,
    /// Last computed solution-change norm `||x_k - x_(k-1)||₂`, when available.
    pub last_solution_change_norm: Option<f64>,
    /// True when the network was solved as a static linear system.
    pub matrix_treated_as_linear_static: bool,
    /// True when `cfd-optim` had to degrade geometry-dependent coefficients for recovery.
    pub degraded_geometry_for_recovery: bool,
    /// Human-readable failure detail from the underlying solver/model error.
    pub failure_detail: Option<String>,
}


/// Error emitted when the trusted primary solve path failed.
#[derive(Debug)]
pub struct PrimarySolveError {
    /// Stable failure classification for the primary solve path.
    pub reason: SolveFailureReason,
    /// Diagnostics captured before failure.
    pub diagnostics: PrimarySolveDiagnostics,
    source: Error,
}

impl PrimarySolveError {
    /// Construct a new typed primary-solve error.
    #[must_use]
    pub fn new(
        reason: SolveFailureReason,
        mut diagnostics: PrimarySolveDiagnostics,
        source: Error,
    ) -> Self {
        if diagnostics.failure_detail.is_none() {
            diagnostics.failure_detail = Some(source.to_string());
        }
        Self {
            reason,
            diagnostics,
            source,
        }
    }

    /// Convert the typed error back into the workspace core error.
    #[must_use]
    pub fn into_source(self) -> Error {
        self.source
    }
}

impl fmt::Display for PrimarySolveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.reason, self.source)
    }
}

impl StdError for PrimarySolveError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        Some(&self.source)
    }
}
