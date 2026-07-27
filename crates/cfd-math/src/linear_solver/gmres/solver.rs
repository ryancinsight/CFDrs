//! GMRES solver.
//!
//! This is a compatibility wrapper over [`leto_ops::GMRES`], the SSOT
//! implementation.  All algorithmic logic lives in `leto-ops`; this module
//! provides the `cfd-math` trait surface and `cfd_core::Error` boundary.

use super::super::bridge::{monitor_from_leto, to_leto_config, CfdLinearOpBridge, CfdPrecondBridge};
use super::super::config::IterativeSolverConfig;
use super::super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{IterativeLinearSolver as LetoIterativeSolver, Scalar};
use std::fmt::Debug;

/// GMRES(m) solver for non-symmetric linear systems.
///
/// # Theorem — GMRES Krylov Optimality (Saad & Schultz 1986)
///
/// GMRES finds the iterate $x_m \in x_0 + K_m(A, r_0)$ that minimises the
/// 2-norm of the residual over the $m$-th Krylov subspace.
///
/// **Reference**: Saad & Schultz (1986), *SIAM J. Sci. Stat. Comput.* 7(3).
pub struct GMRES<T: RealField + Copy> {
    inner: leto_ops::GMRES<T>,
    config: IterativeSolverConfig<T>,
    restart_dim: usize,
}

impl<T: RealField + Copy + NumericElement + Scalar> GMRES<T> {
    /// Create a new GMRES(m) solver.
    ///
    /// # Panics
    ///
    /// Panics if `restart_dim` is zero.
    pub fn new(config: IterativeSolverConfig<T>, restart_dim: usize) -> Self {
        assert!(restart_dim > 0, "GMRES restart dimension must be positive");
        let leto_cfg = to_leto_config(&config);
        Self {
            inner: leto_ops::GMRES::new(leto_cfg, restart_dim),
            config,
            restart_dim,
        }
    }

    /// Restart dimension (maximum Krylov subspace size before restart).
    pub fn restart_dim(&self) -> usize {
        self.restart_dim
    }

    /// Convenience: solve with an explicit preconditioner.
    pub fn solve_preconditioned<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &Array1<T>,
        preconditioner: &P,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        <Self as IterativeLinearSolver<T>>::solve(self, a, b, x, Some(preconditioner))
    }

    /// Convenience: solve without preconditioning.
    pub fn solve_unpreconditioned<Op: LinearOperator<T> + ?Sized>(
        &self,
        a: &Op,
        b: &Array1<T>,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        <Self as IterativeLinearSolver<T>>::solve(
            self,
            a,
            b,
            x,
            None::<&super::super::preconditioners::IdentityPreconditioner>,
        )
    }
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> Configurable<T> for GMRES<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> IterativeLinearSolver<T> for GMRES<T> {
    fn solve<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &Array1<T>,
        x: &mut Array1<T>,
        preconditioner: Option<&P>,
    ) -> Result<ConvergenceMonitor<T>> {
        let bridge_op = CfdLinearOpBridge::new(a);
        let result = if let Some(p) = preconditioner {
            let bridge_p = CfdPrecondBridge::new(p);
            LetoIterativeSolver::solve(&self.inner, &bridge_op, b, x, Some(&bridge_p))
        } else {
            LetoIterativeSolver::solve(
                &self.inner,
                &bridge_op,
                b,
                x,
                None::<&leto_ops::IdentityPreconditioner>,
            )
        };
        result.map(monitor_from_leto).map_err(Error::from)
    }
}

impl<T: RealField + Copy + FloatElement + NumericElement + Scalar> super::super::traits::LinearSolver<T>
    for GMRES<T>
{
    fn solve_system(
        &self,
        a: &dyn LinearOperator<T>,
        b: &Array1<T>,
        x0: Option<&Array1<T>>,
    ) -> Result<Array1<T>> {
        let mut x = if let Some(initial) = x0 {
            initial.clone()
        } else {
            Array1::zeros(b.shape())
        };
        self.solve(
            a,
            b,
            &mut x,
            None::<&super::super::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}