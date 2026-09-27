use super::super::{DGSolution, matrix_cols};
use super::params::LimiterParams;
use super::traits::Limiter;
use crate::error::Result;

/// TVB (Total Variation Bounded) limiter
pub struct TVBLimiter;

impl Limiter for TVBLimiter {
    fn limit(
        &self,
        solution: &mut DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> Result<()> {
        if !params.adaptive || self.is_troubled_cell(solution, neighbors, params) {
            let h = 2.0; // Element size in reference space
            let m = params.tvb_m;

            // Get the cell average
            let u_avg = solution.average();

            // Get neighbor averages
            let u_left = if neighbors.is_empty() {
                u_avg.clone()
            } else {
                neighbors[0].average()
            };

            let u_right = if neighbors.len() > 1 {
                neighbors[1].average()
            } else {
                u_avg.clone()
            };

            // Compute limited slopes
            for i in 0..solution.num_components {
                let du_l = u_avg[i] - u_left[i];
                let du_r = u_right[i] - u_avg[i];
                let du = 0.5 * (du_l + du_r);

                // TVB modified minmod function
                let slope = if du.abs() <= m * h * h {
                    du
                } else if du_l * du_r <= 0.0 {
                    0.0
                } else if du_l.abs() <= du_r.abs() {
                    du_l
                } else {
                    du_r
                };

                // Update the solution coefficients
                for j in 1..matrix_cols(&solution.coefficients) {
                    solution.coefficients[[i, j]] = 0.0;
                }

                // Set the linear term (if any)
                if matrix_cols(&solution.coefficients) > 1 {
                    solution.coefficients[[i, 1]] = slope;
                }
            }
        }

        Ok(())
    }

    fn is_troubled_cell(
        &self,
        solution: &DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> bool {
        if neighbors.len() < 2 {
            return false;
        }

        let u_avg = solution.average();
        let u_left = neighbors[0].average();
        let u_right = neighbors[1].average();
        let h = 2.0; // Element size in reference space
        let m = params.tvb_m;

        // Check for local extrema that exceed the TVB bound
        for i in 0..solution.num_components {
            let du_l = u_avg[i] - u_left[i];
            let du_r = u_right[i] - u_avg[i];
            let du = 0.5 * (du_l + du_r);

            if du.abs() > m * h * h
                && (du_l * du_r <= 0.0 || du_l.abs() > m * h * h || du_r.abs() > m * h * h)
            {
                return true;
            }
        }

        false
    }
}
