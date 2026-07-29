//! Athena Krylov entry points for CFD systems.
//!
//! Athena fixes the GMRES restart width at compile time through a const
//! generic, while CFDrs selects it at runtime — `LinearSolverChain` exposes a
//! builder, JFNK carries it in its config, and the FEM solvers clamp it to the
//! degree-of-freedom count. This module is that bridge: a fixed ladder of
//! instantiations and a dispatch that picks the smallest width covering the
//! request.
//!
//! Widening a restart never costs correctness — GMRES(m) with larger `m`
//! searches a superset of the same Krylov space — so rounding a request up the
//! ladder is safe, and only trades memory for a slightly deeper subspace.

use athena_core::{
    BiCgStab, BiCgStabWorkspace, ConvergencePolicy, Gmres, GmresWorkspace, Identity,
    LinearOperator, Preconditioner, SolveError, SolveReport,
};
use athena_leto::{BorrowedCsrOperator, LetoBackend, LetoBackendError};
use eunomia::{FloatElement, RealField};
use leto::Array1;
use leto_ops::{CsrMatrix, RealScalar};

use super::IterativeSolverConfig;

/// Result of a CFD Krylov solve.
pub type KrylovResult<T> = Result<SolveReport<T>, SolveError<LetoBackendError>>;

/// Restart widths Athena is instantiated at.
///
/// The ladder is geometric so a request of any size lands within a factor of
/// two of its width, bounding both the memory a solve reserves and the number
/// of monomorphisations the crate carries.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RestartWidth {
    W8,
    W16,
    W32,
    W64,
    W128,
    W256,
}

impl RestartWidth {
    /// Smallest ladder width covering `requested`, saturating at the largest.
    ///
    /// A request above the ceiling is served by the ceiling rather than
    /// rejected: the restart is a tuning parameter, and capping it costs
    /// convergence depth rather than correctness.
    const fn covering(requested: usize) -> Self {
        if requested <= 8 {
            Self::W8
        } else if requested <= 16 {
            Self::W16
        } else if requested <= 32 {
            Self::W32
        } else if requested <= 64 {
            Self::W64
        } else if requested <= 128 {
            Self::W128
        } else {
            Self::W256
        }
    }
}

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
    let backend = LetoBackend::<T>::default();
    let operator = BorrowedCsrOperator::new(matrix).map_err(SolveError::Backend)?;
    let policy = convergence_policy(config)?;
    let dimension = LinearOperator::dimension(&operator);

    macro_rules! run {
        ($width:literal) => {{
            let mut workspace = GmresWorkspace::<LetoBackend<T>, $width>::new(&backend, dimension)
                .map_err(SolveError::Backend)?;
            Gmres::<LetoBackend<T>, $width>::solve_into(
                &backend,
                &operator,
                preconditioner,
                right_hand_side,
                solution,
                &mut workspace,
                policy,
            )
        }};
    }

    match RestartWidth::covering(restart) {
        RestartWidth::W8 => run!(8),
        RestartWidth::W16 => run!(16),
        RestartWidth::W32 => run!(32),
        RestartWidth::W64 => run!(64),
        RestartWidth::W128 => run!(128),
        RestartWidth::W256 => run!(256),
    }
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

#[cfg(test)]
mod tests {
    use super::RestartWidth;

    #[test]
    fn the_ladder_covers_every_request() {
        // Each request must land on the smallest width that is at least as
        // large, so a widened restart never searches a smaller subspace than
        // the caller asked for.
        for (requested, expected) in [
            (1, RestartWidth::W8),
            (8, RestartWidth::W8),
            (9, RestartWidth::W16),
            (30, RestartWidth::W32),
            (100, RestartWidth::W128),
            (200, RestartWidth::W256),
        ] {
            assert_eq!(
                RestartWidth::covering(requested),
                expected,
                "request {requested}"
            );
        }
    }

    #[test]
    fn a_request_above_the_ceiling_saturates() {
        assert_eq!(RestartWidth::covering(10_000), RestartWidth::W256);
    }
}
