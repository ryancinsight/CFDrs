//! Anderson acceleration for fixed-point iteration convergence.
//!
//! ## Algorithm: Anderson Acceleration (Walker & Ni 2011)
//!
//! Anderson acceleration applied to a fixed-point iteration x_{k+1} = G(x_k)
//! achieves superlinear convergence for contractive mappings, reducing the
//! Picard iteration count by a factor of 2–5× for typical microfluidic networks.
//!
//! The method maintains a sliding window of the last m iterates and residuals,
//! then solves a least-squares problem to find an optimal linear combination
//! that minimises the fixed-point residual in norm.

use cfd_math::nonlinear_solver::AndersonAccelerator;
use eunomia::{FloatElement, NumericElement};
use leto::Array1;

use super::vector_bridge::{array_from_dvector, dvector_from_array, dvector_from_owned_array};
use super::{NetworkSolveScalar, NetworkSolver};
use cfd_core::physics::fluid::FluidTrait;

impl<T: NetworkSolveScalar, F: FluidTrait<T> + Clone> NetworkSolver<T, F> {
    /// Apply Anderson acceleration (depth m=5) to the Picard iterate sequence.
    ///
    /// Delegates the least-squares subproblem to the Leto-backed `cfd-math`
    /// accelerator; nalgebra vectors remain only at the current 1D linear
    /// system boundary.
    pub(super) fn anderson_accelerate(
        n: usize,
        picard_solution: nalgebra::DVector<T>,
        last_solution: &Array1<T>,
        accelerator: &mut AndersonAccelerator<T>,
    ) -> nalgebra::DVector<T> {
        if n <= 1 {
            return picard_solution;
        }

        let image = array_from_dvector(&picard_solution);
        dvector_from_owned_array(accelerator.compute_next(last_solution, &image))
    }

    /// Select the best next iterate from among Anderson-accelerated, damped Picard,
    /// and raw Picard candidates.
    pub(super) fn select_next_iterate(
        iter: usize,
        n: usize,
        picard_solution: nalgebra::DVector<T>,
        workspace: &mut crate::solver::core::workspace::SolverWorkspace<T>,
        accelerator: &mut AndersonAccelerator<T>,
        matrix: &nalgebra_sparse::CsrMatrix<T>,
        picard_residual: T,
    ) -> nalgebra::DVector<T> {
        let last_solution = dvector_from_array(&workspace.last_solution);
        let picard_step_norm = (&picard_solution - &last_solution).norm();
        let accelerated = Self::anderson_accelerate(
            n,
            picard_solution.clone(),
            &workspace.last_solution,
            accelerator,
        );

        if Self::vector_is_finite(&accelerated) {
            let accelerated_step_norm = (&accelerated - &last_solution).norm();
            let accelerated_residual =
                Self::compute_residual_norm(matrix, &accelerated, &workspace.rhs, n);
            if <T as NumericElement>::is_finite(accelerated_residual)
                && accelerated_step_norm <= picard_step_norm
                && accelerated_residual <= picard_residual.max_scalar(T::default_epsilon())
            {
                return accelerated;
            }
        }

        let alpha = if iter < 20 {
            0.5
        } else if iter < 100 {
            0.35
        } else {
            0.2
        };
        let backup_damped = Self::damped_picard(&last_solution, &picard_solution, alpha);
        if Self::vector_is_finite(&backup_damped) {
            let damped_step_norm = (&backup_damped - &last_solution).norm();
            if damped_step_norm < picard_step_norm {
                return backup_damped;
            }
        }

        if Self::vector_is_finite(&accelerated) {
            accelerated
        } else {
            picard_solution
        }
    }

    /// Under-relaxed Picard step: x_{k+1} = x_k + α·(G(x_k) - x_k)
    pub(super) fn damped_picard(
        last_solution: &nalgebra::DVector<T>,
        picard_solution: &nalgebra::DVector<T>,
        alpha: f64,
    ) -> nalgebra::DVector<T> {
        let alpha = <T as FloatElement>::from_f64(alpha);
        last_solution + (picard_solution - last_solution) * alpha
    }
}
