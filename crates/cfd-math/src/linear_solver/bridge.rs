//! Adapter bridge between the `cfd-math` operator/preconditioner traits and
//! the `leto-ops` equivalents, enabling the CFD-domain solver wrappers to
//! delegate algorithm execution to the Atlas SSOT implementations.
//!
//! ## Design
//!
//! `cfd-math` defines its own `LinearOperator` / `Preconditioner` traits that
//! return `cfd_core::error::Result` and carry CFD-specific methods
//! (`is_symmetric`, `norm_estimate`, …).  `leto-ops` defines the same
//! interfaces against `leto::Result`.  Two zero-sized-field adapter structs
//! bridge the gap at the call site:
//!
//! - [`CfdLinearOpBridge`] — wraps any `cfd-math::LinearOperator<T>` and
//!   presents it as a `leto-ops::LinearOperator<T>`.
//! - [`CfdPrecondBridge`] — wraps any `cfd-math::Preconditioner<T>` and
//!   presents it as a `leto-ops::Preconditioner<T>`.
//!
//! Neither adapter allocates; both convert errors via the `From<LetoError>`
//! impl that lives in `cfd-core`.

use std::marker::PhantomData;

use eunomia::RealField;
use leto::{Array1, LetoError, Result as LetoResult};
use leto_ops::Scalar as LetoScalar;

use super::config::IterativeSolverConfig as CfdSolverConfig;
use super::traits::{
    ConvergenceMonitor as CfdMonitor, LinearOperator as CfdLinearOperator,
    Preconditioner as CfdPreconditioner,
};

// ── Adapter: cfd-math LinearOperator → leto-ops LinearOperator ───────────────

/// Adapts a reference to a `cfd-math::LinearOperator<T>` so it implements
/// `leto-ops::LinearOperator<T>`.  Zero-copy; zero-allocation.
pub(super) struct CfdLinearOpBridge<'a, T, Op>
where
    T: RealField + Copy,
    Op: CfdLinearOperator<T> + ?Sized,
{
    op: &'a Op,
    _phantom: PhantomData<T>,
}

impl<'a, T, Op> CfdLinearOpBridge<'a, T, Op>
where
    T: RealField + Copy,
    Op: CfdLinearOperator<T> + ?Sized,
{
    pub(super) fn new(op: &'a Op) -> Self {
        Self {
            op,
            _phantom: PhantomData,
        }
    }
}

// SAFETY: The underlying `Op: Send + Sync` (required by `CfdLinearOperator`).
unsafe impl<T, Op> Send for CfdLinearOpBridge<'_, T, Op>
where
    T: RealField + Copy,
    Op: CfdLinearOperator<T> + ?Sized,
{
}
unsafe impl<T, Op> Sync for CfdLinearOpBridge<'_, T, Op>
where
    T: RealField + Copy,
    Op: CfdLinearOperator<T> + ?Sized,
{
}

impl<T, Op> leto_ops::LinearOperator<T> for CfdLinearOpBridge<'_, T, Op>
where
    T: RealField + Copy + LetoScalar,
    Op: CfdLinearOperator<T> + ?Sized,
{
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> LetoResult<()> {
        self.op
            .apply(x, y)
            .map_err(|e| LetoError::InvalidInput(e.to_string()))
    }

    fn size(&self) -> usize {
        self.op.size()
    }

    fn nrows(&self) -> usize {
        self.op.size()
    }

    fn ncols(&self) -> usize {
        self.op.size()
    }

    fn apply_transpose(&self, x: &Array1<T>, y: &mut Array1<T>) -> LetoResult<()> {
        self.op
            .apply_transpose(x, y)
            .map_err(|e| LetoError::InvalidInput(e.to_string()))
    }
}

// ── Adapter: cfd-math Preconditioner → leto-ops Preconditioner ───────────────

/// Adapts a reference to a `cfd-math::Preconditioner<T>` so it implements
/// `leto-ops::Preconditioner<T>`.  Zero-copy; zero-allocation.
pub(super) struct CfdPrecondBridge<'a, T, P>
where
    T: RealField + Copy,
    P: CfdPreconditioner<T>,
{
    p: &'a P,
    _phantom: PhantomData<T>,
}

impl<'a, T, P> CfdPrecondBridge<'a, T, P>
where
    T: RealField + Copy,
    P: CfdPreconditioner<T>,
{
    pub(super) fn new(p: &'a P) -> Self {
        Self {
            p,
            _phantom: PhantomData,
        }
    }
}

unsafe impl<T, P> Send for CfdPrecondBridge<'_, T, P>
where
    T: RealField + Copy,
    P: CfdPreconditioner<T>,
{
}
unsafe impl<T, P> Sync for CfdPrecondBridge<'_, T, P>
where
    T: RealField + Copy,
    P: CfdPreconditioner<T>,
{
}

impl<T, P> leto_ops::Preconditioner<T> for CfdPrecondBridge<'_, T, P>
where
    T: RealField + Copy,
    P: CfdPreconditioner<T>,
{
    fn apply_to(&self, r: &Array1<T>, z: &mut Array1<T>) -> LetoResult<()> {
        self.p
            .apply_to(r, z)
            .map_err(|e| LetoError::InvalidInput(e.to_string()))
    }
}

// ── Config conversion ─────────────────────────────────────────────────────────

/// Convert a `cfd-math` solver configuration into the equivalent `leto-ops`
/// configuration.  The extra `relative_tolerance` field introduced in
/// `leto-ops::IterativeSolverConfig` is zeroed (absolute tolerance only),
/// which preserves the existing `cfd-math` convergence semantics.
pub(super) fn to_leto_config<T: RealField + Copy>(
    cfg: &CfdSolverConfig<T>,
) -> leto_ops::IterativeSolverConfig<T> {
    leto_ops::IterativeSolverConfig::new(cfg.tolerance)
        .with_max_iterations(cfg.max_iterations)
}

// ── ConvergenceMonitor conversion ─────────────────────────────────────────────

/// Convert a `leto-ops::ConvergenceMonitor` into the `cfd-math` equivalent.
///
/// Both structs share the same field layout; this is a field-wise copy with no
/// semantic loss.
pub(super) fn monitor_from_leto<T: RealField + Copy>(
    m: leto_ops::ConvergenceMonitor<T>,
) -> CfdMonitor<T> {
    CfdMonitor {
        initial_residual: m.initial_residual,
        iteration: m.iteration,
        residual_history: m.residual_history,
        theoretical_bound: m.theoretical_bound,
        condition_number_estimate: m.condition_number_estimate,
    }
}
