//! Athena Krylov entry points for CFD systems.
//!
//! Athena fixes the GMRES restart width at compile time through a const
//! generic, while CFDrs selects it at runtime — `LinearSolverChain` exposes a
//! builder, JFNK carries it in its config, and the FEM solvers clamp it to the
//! degree-of-freedom count. The bridge between the two is
//! [`athena_leto::KrylovWorkspace`], the one runtime bridge the stack carries;
//! the ladder of instantiations and its dispatch live there (Atlas ADR 0062).
//! This module holds the CFD-facing surface on top of it: the configuration
//! vocabulary translated into a validated Athena policy, the entry points that
//! wrap a caller-owned CSR matrix in Athena's operator seam, and the
//! CFD-domain interpretation of a solve outcome.

use athena_core::{
    BiCgStab, BiCgStabWorkspace, Cg, CgWorkspace, ConvergencePolicy, Identity, LinearOperator,
    Preconditioner, SolveError, SolveReport, Termination,
};
use athena_leto::{BorrowedCsrOperator, KrylovWorkspace, LetoBackend, LetoBackendError};
use eunomia::{FloatElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, RealScalar};

use super::IterativeSolverConfig;

/// Result of a CFD Krylov solve.
pub type KrylovResult<T> = Result<SolveReport<T>, SolveError<LetoBackendError>>;

/// Translate a CFD solver configuration into a validated Athena policy.
///
/// # Errors
///
/// Returns [`LetoBackendError::Leto`] when the configuration carries a
/// non-finite or negative tolerance, or a zero iteration budget, which the
/// Athena policy rejects at construction rather than mid-solve.
pub fn convergence_policy<T>(
    config: &IterativeSolverConfig<T>,
) -> Result<ConvergencePolicy<T>, SolveError<LetoBackendError>>
where
    T: RealField + Copy,
{
    ConvergencePolicy::new(
        config.tolerance,
        config.relative_tolerance,
        config.max_iterations,
    )
    .map_err(|reason| {
        SolveError::Backend(LetoBackendError::Leto(leto::LetoError::InvalidInput(
            format!("invalid iterative solver configuration: {reason}"),
        )))
    })
}

/// Solve `A·x = b` with restarted GMRES and an explicit preconditioner.
///
/// `restart` is rounded up the instantiation ladder. `solution` carries the
/// initial guess in and the result out.
///
/// # Errors
///
/// Returns a dimension, configuration, or backend failure. Numerical
/// termination is reported value-semantically in the [`SolveReport`].
pub fn gmres_preconditioned<T, P>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    preconditioner: &P,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
    restart: usize,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
    P: Preconditioner<LetoBackend<T>>,
{
    let dimension = matrix.shape().0;
    let mut workspace = KrylovWorkspace::new(restart, dimension).map_err(SolveError::Backend)?;
    gmres_preconditioned_with_workspace(
        matrix,
        right_hand_side,
        preconditioner,
        solution,
        config,
        &mut workspace,
    )
}

/// Solve `A·x = b` with restarted GMRES and a caller-owned workspace.
///
/// The workspace dimension must match the matrix dimension. A mismatch is
/// reported by Athena's value-semantic dimension validation.
///
/// # Errors
///
/// Returns a dimension, configuration, or backend failure. Numerical
/// termination is reported value-semantically in the [`SolveReport`].
pub fn gmres_preconditioned_with_workspace<T, P>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    preconditioner: &P,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
    workspace: &mut KrylovWorkspace<T>,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
    P: Preconditioner<LetoBackend<T>>,
{
    let operator = BorrowedCsrOperator::new(matrix).map_err(SolveError::Backend)?;
    let policy = convergence_policy(config)?;
    workspace.solve(&operator, preconditioner, right_hand_side, solution, policy)
}

/// Solve `A·x = b` with restarted GMRES and no preconditioner.
///
/// # Errors
///
/// See [`gmres_preconditioned`].
pub fn gmres<T>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
    restart: usize,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
{
    gmres_preconditioned(
        matrix,
        right_hand_side,
        &Identity,
        solution,
        config,
        restart,
    )
}

/// Solve `A·x = b` with BiCGSTAB and an explicit preconditioner.
///
/// # Errors
///
/// See [`gmres_preconditioned`].
pub fn bicgstab_preconditioned<T, P>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    preconditioner: &P,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
    P: Preconditioner<LetoBackend<T>>,
{
    let backend = LetoBackend::<T>::default();
    let operator = BorrowedCsrOperator::new(matrix).map_err(SolveError::Backend)?;
    let policy = convergence_policy(config)?;
    let mut workspace = BiCgStabWorkspace::new(&backend, LinearOperator::dimension(&operator))
        .map_err(SolveError::Backend)?;
    BiCgStab::<LetoBackend<T>>::solve_into(
        &backend,
        &operator,
        preconditioner,
        right_hand_side,
        solution,
        &mut workspace,
        policy,
    )
}

/// Solve `A·x = b` with BiCGSTAB and no preconditioner.
///
/// # Errors
///
/// See [`gmres_preconditioned`].
pub fn bicgstab<T>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
{
    let backend = LetoBackend::<T>::default();
    let operator = BorrowedCsrOperator::new(matrix).map_err(SolveError::Backend)?;
    let policy = convergence_policy(config)?;
    let mut workspace = BiCgStabWorkspace::new(&backend, LinearOperator::dimension(&operator))
        .map_err(SolveError::Backend)?;
    BiCgStab::<LetoBackend<T>>::solve_into(
        &backend,
        &operator,
        &Identity,
        right_hand_side,
        solution,
        &mut workspace,
        policy,
    )
}

/// Solve `A·x = b` with preconditioned conjugate gradients.
///
/// The operator must be symmetric positive definite; CG has no contract for
/// anything else and reports non-positive curvature rather than converging.
///
/// # Errors
///
/// See [`gmres_preconditioned`].
pub fn cg_preconditioned<T, P>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    preconditioner: &P,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
    P: Preconditioner<LetoBackend<T>>,
{
    let backend = LetoBackend::<T>::default();
    let operator = BorrowedCsrOperator::new(matrix).map_err(SolveError::Backend)?;
    let policy = convergence_policy(config)?;
    let mut workspace = CgWorkspace::new(&backend, LinearOperator::dimension(&operator))
        .map_err(SolveError::Backend)?;
    Cg::<LetoBackend<T>>::solve_into(
        &backend,
        &operator,
        preconditioner,
        right_hand_side,
        solution,
        &mut workspace,
        policy,
    )
}

/// Solve `A·x = b` with conjugate gradients and no preconditioner.
///
/// # Errors
///
/// See [`gmres_preconditioned`].
pub fn cg<T>(
    matrix: &CsrMatrix<T>,
    right_hand_side: &Array1<T>,
    solution: &mut Array1<T>,
    config: &IterativeSolverConfig<T>,
) -> KrylovResult<T>
where
    T: RealScalar + RealField + FloatElement,
{
    cg_preconditioned(matrix, right_hand_side, &Identity, solution, config)
}

/// Interpret an Athena solve outcome in CFD terms.
///
/// Athena reports numerical termination value-semantically in its
/// [`SolveReport`] rather than as an error, so an exhausted budget arrives as
/// `Ok`. CFD callers treat a stalled solve as a recoverable condition — the
/// last iterate is still a usable approximation inside an outer nonlinear
/// loop — while a genuine breakdown is not, so the two are separated here
/// instead of at each of the call sites.
///
/// # Errors
///
/// Returns [`cfd_core::error::Error::Solver`] for a dimension, backend, or
/// numerical-breakdown outcome, naming `context`.
pub fn interpret<T>(
    context: &str,
    outcome: KrylovResult<T>,
) -> Result<SolveOutcome<T>, cfd_core::error::Error>
where
    T: RealField + Copy,
{
    let report = outcome.map_err(|error| {
        cfd_core::error::Error::Solver(format!("{context}: linear solve failed: {error}"))
    })?;
    match report.termination {
        Termination::Converged | Termination::InitialResidual | Termination::NormalEquations => {
            Ok(SolveOutcome::Converged(report))
        }
        Termination::MaxIterations => Ok(SolveOutcome::Stalled(report)),
        // Breakdown leaves the last iterate intact. Whether that is usable is
        // the caller's judgement — an outer nonlinear loop can often continue
        // from it — so it is reported rather than raised.
        Termination::Breakdown | Termination::NonPositiveCurvature => {
            Ok(SolveOutcome::BrokenDown(report))
        }
        other => Err(cfd_core::error::Error::Solver(format!(
            "{context}: linear solve terminated as {other:?}"
        ))),
    }
}

/// Outcome of an interpreted CFD linear solve.
#[derive(Clone, Copy, Debug)]
pub enum SolveOutcome<T> {
    /// The solve met its convergence policy.
    Converged(SolveReport<T>),
    /// The iteration budget was exhausted. The last iterate remains a usable
    /// approximation for an outer nonlinear iteration.
    Stalled(SolveReport<T>),
    /// The recurrence broke down. The last iterate is still present; whether
    /// it is usable is the caller's judgement.
    BrokenDown(SolveReport<T>),
}

impl<T> SolveOutcome<T> {
    /// The report, whichever way the solve ended.
    pub const fn report(&self) -> &SolveReport<T> {
        match self {
            Self::Converged(report) | Self::Stalled(report) | Self::BrokenDown(report) => report,
        }
    }

    /// Whether the solve met its convergence policy.
    #[must_use]
    pub const fn converged(&self) -> bool {
        matches!(self, Self::Converged(_))
    }
}

/// Report only a converged solve, logging any other outcome.
///
/// A fallback chain treats every non-convergence the same way — move to the
/// next tier — so the distinction between a stalled solve, a breakdown, and a
/// hard failure belongs in the log rather than in each tier's control flow.
/// Athena reports the first two value-semantically, which is why this cannot
/// just be an `ok()`.
pub fn converged_or_none<T>(context: &str, outcome: KrylovResult<T>) -> Option<SolveReport<T>>
where
    T: RealField + Copy,
{
    match interpret(context, outcome) {
        Ok(SolveOutcome::Converged(report)) => Some(report),
        Ok(other) => {
            tracing::warn!(
                iterations = other.report().iterations,
                initial_residual = ?other.report().initial_residual_norm,
                final_residual = ?other.report().final_residual_norm,
                threshold = ?other.report().threshold,
                "{context}: did not converge ({:?})",
                other.report().termination
            );
            None
        }
        Err(error) => {
            tracing::warn!("{context}: {error}");
            None
        }
    }
}

/// The Krylov recurrences CFDrs selects between at runtime.
///
/// Athena's solvers are zero-sized markers carrying const generics and backend
/// GATs, so they are not object-safe and cannot be held behind `dyn`. The
/// closed set is small and fixed, which makes an exhaustively-matched enum the
/// right dispatch: it keeps runtime selection without a vtable, and adding a
/// recurrence becomes a compile error at every selection site rather than a
/// silent fallthrough.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SolverKind {
    /// Conjugate gradients. Requires a symmetric positive definite operator.
    ConjugateGradient,
    /// Stabilised biconjugate gradients, for general nonsymmetric operators.
    BiCgStab,
    /// Restarted GMRES with the given restart width.
    Gmres {
        /// Requested restart width, rounded up the instantiation ladder.
        restart: usize,
    },
}

impl SolverKind {
    /// Human-readable name, for validation reports and logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::ConjugateGradient => "ConjugateGradient",
            Self::BiCgStab => "BiCGSTAB",
            Self::Gmres { .. } => "GMRES",
        }
    }

    /// Solve with this recurrence and no preconditioner.
    ///
    /// # Errors
    ///
    /// See [`gmres_preconditioned`].
    pub fn solve<T>(
        self,
        matrix: &CsrMatrix<T>,
        right_hand_side: &Array1<T>,
        solution: &mut Array1<T>,
        config: &IterativeSolverConfig<T>,
    ) -> KrylovResult<T>
    where
        T: RealScalar + RealField + FloatElement,
    {
        match self {
            Self::ConjugateGradient => cg(matrix, right_hand_side, solution, config),
            Self::BiCgStab => bicgstab(matrix, right_hand_side, solution, config),
            Self::Gmres { restart } => gmres(matrix, right_hand_side, solution, config, restart),
        }
    }

    /// Solve with this recurrence and an explicit preconditioner.
    ///
    /// # Errors
    ///
    /// See [`gmres_preconditioned`].
    pub fn solve_preconditioned<T, P>(
        self,
        matrix: &CsrMatrix<T>,
        right_hand_side: &Array1<T>,
        preconditioner: &P,
        solution: &mut Array1<T>,
        config: &IterativeSolverConfig<T>,
    ) -> KrylovResult<T>
    where
        T: RealScalar + RealField + FloatElement,
        P: Preconditioner<LetoBackend<T>>,
    {
        match self {
            Self::ConjugateGradient => {
                cg_preconditioned(matrix, right_hand_side, preconditioner, solution, config)
            }
            Self::BiCgStab => {
                bicgstab_preconditioned(matrix, right_hand_side, preconditioner, solution, config)
            }
            Self::Gmres { restart } => gmres_preconditioned(
                matrix,
                right_hand_side,
                preconditioner,
                solution,
                config,
                restart,
            ),
        }
    }
}
