//! Preconditioned Conjugate Gradient (PCG) solver.
//!
//! This is a compatibility wrapper over [`leto_ops::ConjugateGradient`], the
//! SSOT implementation.  All algorithmic logic lives in `leto-ops`; this crate
//! provides the `cfd-math` trait surface and `cfd_core::Error` boundary.

use super::bridge::{monitor_from_leto, to_leto_config, CfdLinearOpBridge, CfdPrecondBridge};
use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::{
    IterativeLinearSolver as LetoIterativeSolver, Scalar,
};
use std::fmt::Debug;

/// Preconditioned Conjugate Gradient solver.
///
/// # Theorem (CG Convergence in the A-norm)
///
/// For an SPD matrix $A \in \mathbb{R}^{n \times n}$ with eigenvalues
/// $0 < \lambda_1 \le \cdots \le \lambda_n$, the CG iterates satisfy
///
/// $$\|x_k - x^*\|_A \le 2 \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|x_0 - x^*\|_A$$
///
/// where $\kappa = \lambda_n / \lambda_1$ is the spectral condition number.
///
/// **Reference**: Hestenes & Stiefel (1952); Saad (2003) §6.7.
pub struct ConjugateGradient<T: RealField + Copy> {
    inner: leto_ops::ConjugateGradient<T>,
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy + NumericElement + Scalar> ConjugateGradient<T> {
    /// Create a new CG solver with the given configuration.
    pub fn new(config: IterativeSolverConfig<T>) -> Self {
        let leto_cfg = to_leto_config(&config);
        Self {
            inner: leto_ops::ConjugateGradient::new(leto_cfg),
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
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> Configurable<T>
    for ConjugateGradient<T>
{
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + NumericElement + Scalar> IterativeLinearSolver<T>
    for ConjugateGradient<T>
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
    for ConjugateGradient<T>
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

    fn create_simple_spd_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0];
        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let config = IterativeSolverConfig::default();
        let _solver = ConjugateGradient::<f64>::new(config);
    }

    #[test]
    fn test_default_solver() {
        let _solver = ConjugateGradient::<f64>::default();
    }

    #[test]
    fn test_solve_simple_system() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        let result = solver.solve(&a, &b, &mut x, Some(&precond));
        assert!(result.is_ok(), "CG should converge");
        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = array(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        assert!(solver.solve(&a, &b, &mut x, Some(&precond)).is_ok());
    }

    #[test]
    fn test_mismatched_dimensions() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0]);
        let mut x = Array1::zeros([2]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        assert!(solver.solve(&a, &b, &mut x, Some(&precond)).is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;
        let result = solver.solve(&a, &b, &mut x, Some(&precond));
        assert!(result.is_ok(), "CG should converge with tight tolerance");
        assert_solves(&a, &x, &b, 1e-12);
    }

    #[test]
    fn test_configurable_trait() {
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-8,
            max_iterations: 500,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let retrieved = solver.config();
        assert_relative_eq!(retrieved.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());
        assert_solves(&a, &result.unwrap(), &b, 1e-6);
    }
}