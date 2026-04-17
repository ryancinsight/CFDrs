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

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::NetworkSolver;
use cfd_core::physics::fluid::FluidTrait;

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + Float, F: FluidTrait<T> + Clone>
    NetworkSolver<T, F>
{
    /// Apply Anderson acceleration (depth m=5) to the Picard iterate sequence.
    ///
    /// Minimises the fixed-point residual in the least-squares sense over the
    /// last m iterates (Walker & Ni 2011).
    pub(super) fn anderson_accelerate(
        iter: usize,
        n: usize,
        picard_solution: nalgebra::DVector<T>,
        workspace: &mut crate::solver::core::workspace::SolverWorkspace<T>,
        depth: usize,
    ) -> nalgebra::DVector<T> {
        let residuals = &mut workspace.anderson_residuals;
        let iterates = &mut workspace.anderson_iterates;
        let last_solution = &workspace.last_solution;
        if iter == 0 || n <= 1 {
            if iter == 0 {
                let residual = &picard_solution - last_solution;
                residuals.push_back(residual);
                iterates.push_back(picard_solution.clone());
            }
            return picard_solution;
        }

        let residual = &picard_solution - last_solution;
        if residuals.len() >= depth {
            residuals.pop_front();
            iterates.pop_front();
        }
        residuals.push_back(residual);
        iterates.push_back(picard_solution.clone());

        let m = residuals.len();
        if m < 2 {
            return picard_solution;
        }

        let ncols = m - 1;
        let r_last = &residuals[m - 1];
        let mut gram = nalgebra::DMatrix::<T>::zeros(ncols, ncols);
        let mut rhs_ls = nalgebra::DVector::<T>::zeros(ncols);

        for j in 0..ncols {
            let dr_j = &residuals[j] - r_last;
            rhs_ls[j] = dr_j.dot(r_last);
            for k in j..ncols {
                let dr_k = &residuals[k] - r_last;
                let val = dr_j.dot(&dr_k);
                gram[(j, k)] = val;
                gram[(k, j)] = val;
            }
        }

        nalgebra::linalg::LU::new(gram)
            .solve(&rhs_ls)
            .map_or(picard_solution, |lu| {
                let alpha_sum: T = lu.iter().fold(T::zero(), |acc, &a| acc + a);
                let one_minus_sum = T::one() - alpha_sum;
                let x_last = &iterates[m - 1];
                let mut accelerated = x_last * one_minus_sum + r_last * one_minus_sum;
                for j in 0..ncols {
                    accelerated += (&iterates[j] + &residuals[j]) * lu[j];
                }
                accelerated
            })
    }

    /// Select the best next iterate from among Anderson-accelerated, damped Picard,
    /// and raw Picard candidates.
    pub(super) fn select_next_iterate(
        iter: usize,
        n: usize,
        picard_solution: nalgebra::DVector<T>,
        workspace: &mut crate::solver::core::workspace::SolverWorkspace<T>,
        depth: usize,
        matrix: &nalgebra_sparse::CsrMatrix<T>,
        picard_residual: T,
    ) -> nalgebra::DVector<T> {
        let picard_step_norm = (&picard_solution - &workspace.last_solution).norm();
        let accelerated =
            Self::anderson_accelerate(iter, n, picard_solution.clone(), workspace, depth);

        if Self::vector_is_finite(&accelerated) {
            let accelerated_step_norm = (&accelerated - &workspace.last_solution).norm();
            let accelerated_residual =
                Self::compute_residual_norm(matrix, &accelerated, &workspace.rhs, n);
            if accelerated_residual.is_finite()
                && accelerated_step_norm <= picard_step_norm
                && accelerated_residual <= <T as Float>::max(picard_residual, T::default_epsilon())
            {
                return accelerated;
            }
        }

        let backup_damped = Self::damped_picard(&workspace.last_solution, &picard_solution, 0.5);
        if Self::vector_is_finite(&backup_damped) {
            let damped_step_norm = (&backup_damped - &workspace.last_solution).norm();
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
        let alpha = T::from_f64(alpha).expect("Mathematical constant conversion compromised");
        last_solution + (picard_solution - last_solution) * alpha
    }
}
