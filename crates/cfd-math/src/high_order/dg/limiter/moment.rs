use super::super::{DGSolution, matrix_cols};
use super::params::LimiterParams;
use super::traits::Limiter;
use crate::error::Result;

/// Moment limiter for high-order methods
pub struct MomentLimiter;

impl Limiter for MomentLimiter {
    fn limit(
        &self,
        solution: &mut DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> Result<()> {
        if !params.adaptive || self.is_troubled_cell(solution, neighbors, params) {
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

            // Limit each component
            for i in 0..solution.num_components {
                // First, limit the highest order coefficients
                for j in (1..matrix_cols(&solution.coefficients)).rev() {
                    // Compute the limited coefficient
                    let c = solution.coefficients[[i, j]];

                    // Compute the minmod of the coefficient and its neighbors
                    let c_min = if j == 1 {
                        // For the first moment, use the minmod of the slopes
                        let du_l = u_avg[i] - u_left[i];
                        let du_r = u_right[i] - u_avg[i];

                        if du_l * du_r <= 0.0 {
                            0.0
                        } else if du_l.abs() <= du_r.abs() {
                            du_l
                        } else {
                            du_r
                        }
                    } else {
                        // For higher moments, use the minmod of the current coefficient
                        // and the same coefficient from the neighbors
                        let c_left = if !neighbors.is_empty()
                            && matrix_cols(&neighbors[0].coefficients) > j
                        {
                            neighbors[0].coefficients[[i, j]]
                        } else {
                            0.0
                        };

                        let c_right =
                            if neighbors.len() > 1 && matrix_cols(&neighbors[1].coefficients) > j {
                                neighbors[1].coefficients[[i, j]]
                            } else {
                                0.0
                            };

                        let mut c_min = c;

                        // Check left neighbor: if sign differs or neighbor is zero, result is zero
                        if c * c_left <= 0.0 {
                            c_min = 0.0;
                        } else if c_left.abs() < c_min.abs() {
                            c_min = c_left;
                        }

                        // Check right neighbor: if sign differs or neighbor is zero, result is zero
                        if c_min * c_right <= 0.0 {
                            c_min = 0.0;
                        } else if c_right.abs() < c_min.abs() {
                            c_min = c_right;
                        }

                        c_min
                    };

                    // Update the coefficient
                    solution.coefficients[[i, j]] = c_min;
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

        // Check if any coefficient is too large compared to the average
        let u_avg = solution.average();

        for i in 0..solution.num_components {
            for j in 1..matrix_cols(&solution.coefficients) {
                let c = solution.coefficients[[i, j]];

                if c.abs() > params.tolerance * u_avg[i].abs() {
                    return true;
                }
            }
        }

        false
    }
}
