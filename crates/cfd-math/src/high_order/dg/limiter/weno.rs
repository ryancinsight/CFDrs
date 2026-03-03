//! WENO (Weighted Essentially Non-Oscillatory) limiter for DG methods.

use super::{DGSolution, Limiter, LimiterParams, Result};

/// WENO (Weighted Essentially Non-Oscillatory) limiter
pub struct WENOLimiter {
    /// Number of stencils
    k: usize,
    /// WENO coefficients
    c: Vec<f64>,
    /// Smoothness indicators
    beta: Vec<f64>,
}

impl WENOLimiter {
    /// Create a new WENO limiter
    pub fn new(k: usize) -> Self {
        // Compute the WENO coefficients
        let c = match k {
            2 => vec![1.0 / 3.0, 2.0 / 3.0],
            3 => vec![0.1, 0.6, 0.3],
            _ => vec![1.0], // Default to first-order
        };

        Self {
            k: k.min(3), // Maximum k=3 supported
            c,
            beta: vec![0.0; k],
        }
    }

    /// Compute the smoothness indicators
    fn compute_smoothness_indicators(&mut self, u: &[f64]) -> Vec<f64> {
        match self.k {
            2 => {
                // Second-order WENO
                self.beta[0] = u[0] * u[0];
                self.beta[1] = u[1] * u[1];
            }
            3 => {
                // Third-order WENO
                self.beta[0] = (13.0 / 12.0) * (u[0] - 2.0 * u[1] + u[2]).powi(2)
                    + (1.0 / 4.0) * (3.0 * u[0] - 4.0 * u[1] + u[2]).powi(2);

                self.beta[1] = (13.0 / 12.0) * (u[1] - 2.0 * u[2] + u[3]).powi(2)
                    + (1.0 / 4.0) * (u[1] - u[3]).powi(2);

                self.beta[2] = (13.0 / 12.0) * (u[2] - 2.0 * u[3] + u[4]).powi(2)
                    + (1.0 / 4.0) * (u[2] - 4.0 * u[3] + 3.0 * u[4]).powi(2);
            }
            _ => {
                // First-order (no WENO)
                self.beta[0] = 0.0;
            }
        }

        self.beta.clone()
    }

    /// Compute the WENO weights
    fn compute_weights(&self, beta: &[f64], epsilon: f64) -> Vec<f64> {
        let mut alpha = vec![0.0; self.k];
        let mut sum_alpha = 0.0;

        for i in 0..self.k {
            alpha[i] = self.c[i] / (epsilon + beta[i]).powi(2);
            sum_alpha += alpha[i];
        }

        alpha.iter().map(|&a| a / sum_alpha).collect()
    }

    /// Reconstruct the solution using WENO
    fn weno_reconstruction(&mut self, u: &[f64], epsilon: f64) -> f64 {
        if self.k == 1 {
            return u[0];
        }

        let beta = self.compute_smoothness_indicators(u);
        let weights = self.compute_weights(&beta, epsilon);

        // Compute the reconstructed value

        match self.k {
            2 => {
                // Second-order WENO
                weights[0] * (1.5 * u[0] - 0.5 * u[1]) + weights[1] * (0.5 * u[0] + 0.5 * u[1])
            }
            3 => {
                // Third-order WENO
                weights[0] * (11.0 / 6.0 * u[0] - 7.0 / 6.0 * u[1] + 1.0 / 3.0 * u[2])
                    + weights[1] * (1.0 / 3.0 * u[1] + 5.0 / 6.0 * u[2] - 1.0 / 6.0 * u[3])
                    + weights[2] * (-1.0 / 6.0 * u[2] + 5.0 / 6.0 * u[3] + 1.0 / 3.0 * u[4])
            }
            _ => {
                // Should not happen
                u[0]
            }
        }
    }
}

impl Limiter for WENOLimiter {
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

            // Apply WENO limiting to each component
            for i in 0..solution.num_components {
                // Get the stencil values
                let mut u_stencil = vec![0.0; 5]; // Support up to 5-point stencil

                // Left neighbor of left neighbor (if available)
                u_stencil[0] = if neighbors.len() > 2 {
                    neighbors[0].average()[i]
                } else {
                    u_left[i]
                };

                // Left neighbor
                u_stencil[1] = u_left[i];

                // Current cell
                u_stencil[2] = u_avg[i];

                // Right neighbor
                u_stencil[3] = u_right[i];

                // Right neighbor of right neighbor (if available)
                u_stencil[4] = if neighbors.len() > 3 {
                    neighbors[2].average()[i]
                } else {
                    u_right[i]
                };

                // Apply WENO reconstruction
                let mut weno_limiter = WENOLimiter::new(3); // Use 3rd order WENO
                let u_rec = weno_limiter.weno_reconstruction(&u_stencil, params.weno_epsilon);

                // Update the solution coefficients
                for j in 1..solution.coefficients.ncols() {
                    solution.coefficients[(i, j)] = 0.0;
                }

                // Set the first moment (linear term)
                if solution.coefficients.ncols() > 1 {
                    solution.coefficients[(i, 1)] = u_rec - u_avg[i];
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
            for j in 1..solution.coefficients.ncols() {
                let c = solution.coefficients[(i, j)];

                if c.abs() > params.tolerance * u_avg[i].abs() {
                    return true;
                }
            }
        }

        false
    }
}
