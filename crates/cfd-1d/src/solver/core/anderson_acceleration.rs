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
use leto_ops::CsrMatrix as LetoCsrMatrix;

use super::vector_bridge::array_l2_norm;
use super::{NetworkSolveScalar, NetworkSolver};
use crate::solver::core::workspace::SolverWorkspace;
use cfd_core::physics::fluid::FluidTrait;

impl<T: NetworkSolveScalar, F: FluidTrait<T> + Clone> NetworkSolver<T, F> {
    /// Apply Anderson acceleration (depth m=5) to the Picard iterate sequence.
    ///
    /// Delegates the least-squares subproblem to the Leto-backed `cfd-math`
    /// accelerator. Both input and output are `Array1<T>`.
    pub(super) fn anderson_accelerate(
        n: usize,
        picard_solution: Array1<T>,
        last_solution: &Array1<T>,
        accelerator: &mut AndersonAccelerator<T>,
    ) -> Array1<T> {
        if n <= 1 {
            return picard_solution;
        }
        accelerator.compute_next(last_solution, &picard_solution)
    }

    /// Select the best next iterate from among Anderson-accelerated, damped Picard,
    /// and raw Picard candidates.
    pub(super) fn select_next_iterate(
        iter: usize,
        n: usize,
        picard_solution: Array1<T>,
        workspace: &mut SolverWorkspace<T>,
        accelerator: &mut AndersonAccelerator<T>,
        matrix: &LetoCsrMatrix<T>,
        picard_residual: T,
    ) -> Array1<T> {
        let picard_step_norm =
            array_l2_norm(&(&picard_solution - &workspace.last_solution));
        let accelerated = Self::anderson_accelerate(
            n,
            picard_solution.clone(),
            &workspace.last_solution,
            accelerator,
        );

        if Self::vector_is_finite(&accelerated) {
            let accelerated_step_norm =
                array_l2_norm(&(&accelerated - &workspace.last_solution));
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
        let backup_damped =
            Self::damped_picard(&workspace.last_solution, &picard_solution, alpha);
        if Self::vector_is_finite(&backup_damped) {
            let damped_step_norm =
                array_l2_norm(&(&backup_damped - &workspace.last_solution));
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
        last_solution: &Array1<T>,
        picard_solution: &Array1<T>,
        alpha: f64,
    ) -> Array1<T> {
        let alpha_t = <T as FloatElement>::from_f64(alpha);
        let n = last_solution.size();
        Array1::from_shape_fn([n], |[i]| {
            last_solution[[i]] + (picard_solution[[i]] - last_solution[[i]]) * alpha_t
        })
    }
}
