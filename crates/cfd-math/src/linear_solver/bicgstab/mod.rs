//! BiCGSTAB solver.
//!
//! This is a compatibility wrapper over [`leto_ops::BiCGSTAB`], the SSOT
//! implementation.  All algorithmic logic lives in `leto-ops`; this module
//! provides the `cfd-math` trait surface and `cfd_core::Error` boundary.

use super::bridge::{monitor_from_leto, to_leto_config, CfdLinearOpBridge, CfdPrecondBridge};
use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{IterativeLinearSolver as LetoIterativeSolver, Scalar};
use std::fmt::Debug;

/// BiCGSTAB solver for non-symmetric linear systems.
///
/// # Convergence Theory
///
/// BiCGSTAB (BiConjugate Gradient Stabilized) stabilises Bi-CG by introducing
/// an additional GMRES-like step.  The residual satisfies
///
/// ```text
/// r_k = p_{2k}(A) r_0
/// ```
///
/// where $p_{2k}$ is a product of two polynomials, each of degree $k$.
/// Convergence is guaranteed when the field of values of $A$ lies in the right
/// half-plane.
///
/// **Reference**: Van der Vorst (1992), *SIAM J. Sci. Stat. Comput.* 13(2).
pub struct BiCGSTAB<T: RealField + Copy> {
    inner: leto_ops::BiCGSTAB<T>,
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy + NumericElement + Scalar> BiCGSTAB<T> {
    /// Create a new BiCGSTAB solver with the given configuration.
    pub fn new(config: IterativeSolverConfig<T>) -> Self {
        let leto_cfg = to_leto_config(&config);
        Self {
            inner: leto_ops::BiCGSTAB::new(leto_cfg),
            config,
        }
    }

    /// Create with default configuration.
    #[must_use]
    pub fn default() -> Self
    where
        T: FloatElement,
    {
        Self::new(IterativeSolverConfig::default())
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
            None::<&super::preconditioners::IdentityPreconditioner>,
        )
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
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> Configurable<T> for BiCGSTAB<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> IterativeLinearSolver<T>
    for BiCGSTAB<T>
{
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

impl<T: RealField + Copy + FloatElement + NumericElement + Scalar> super::traits::LinearSolver<T>
    for BiCGSTAB<T>
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
            None::<&super::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::preconditioners::IdentityPreconditioner;
    use super::super::traits::{Configurable, LinearSolver};
    use super::*;
    use eunomia::assert_relative_eq;
    use leto_ops::CsrMatrix;

    fn array(values: Vec<f64>) -> Array1<f64> {
        Array1::from_shape_vec([values.len()], values).expect("valid Leto vector shape")
    }

    fn assert_solves(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>, epsilon: f64) {
        let mut ax = Array1::zeros([b.shape()[0]]);
        a.apply(x, &mut ax).expect("operator application");
        for idx in 0..b.shape()[0] {
            assert_relative_eq!(ax[idx], b[idx], epsilon = epsilon);
        }
    }

    fn create_nonsymmetric_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![5.0, 1.0, 2.0, 4.0, 1.0, 2.0, 3.0];
        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
    }

    fn create_diagonal_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];
        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let _solver = BiCGSTAB::<f64>::new(IterativeSolverConfig::default());
    }

    #[test]
    fn test_solve_simple_system() {
        let a = create_nonsymmetric_matrix();
        let b = array(vec![6.0, 11.0, 8.0]);
        let mut x = Array1::zeros([3]);
        let solver = BiCGSTAB::new(IterativeSolverConfig::default());
        let precond = IdentityPreconditioner;
        let result = solver.solve(&a, &b, &mut x, Some(&precond));
        assert!(result.is_ok(), "BiCGSTAB should solve nonsymmetric CSR system");
        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let a = create_diagonal_matrix();
        let b = array(vec![2.0, 3.0, 4.0]);
        let mut x = Array1::zeros([3]);
        let solver = BiCGSTAB::new(IterativeSolverConfig::default());
        let precond = IdentityPreconditioner;
        let result = solver.solve(&a, &b, &mut x, Some(&precond));
        assert!(result.is_ok());
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_configurable_trait() {
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-8,
            max_iterations: 500,
            ..Default::default()
        };
        let solver = BiCGSTAB::new(config);
        let retrieved = solver.config();
        assert_relative_eq!(retrieved.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_nonsymmetric_matrix();
        let b = Array1::from_shape_vec([3], vec![6.0, 11.0, 8.0]).unwrap();
        let solver = BiCGSTAB::new(IterativeSolverConfig::default());
        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok(), "linear-solver facade should converge");
        assert_solves(&a, &result.unwrap(), &b, 1e-6);
    }
}