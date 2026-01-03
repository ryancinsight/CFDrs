//! Slope limiters for Discontinuous Galerkin methods.
//!
//! This module provides various slope limiters that can be used to control
//! oscillations in DG solutions, especially in the presence of shocks or discontinuities.

use super::{DGSolution, Result};


/// Type of slope limiter
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LimiterType {
    /// No limiting (identity operator)
    None,
    /// Minmod limiter (most diffusive)
    Minmod,
    /// Monotonized Central (MC) limiter
    MC,
    /// Superbee limiter (least diffusive)
    Superbee,
    /// TVB (Total Variation Bounded) limiter
    TVB,
    /// Moment limiter (for high-order methods)
    Moment,
    /// WENO (Weighted Essentially Non-Oscillatory) limiter
    WENO,
}

/// Parameters for slope limiters
#[derive(Debug, Clone)]
pub struct LimiterParams {
    /// Type of limiter
    pub limiter_type: LimiterType,
    /// TVB parameter (for TVB limiter)
    pub tvb_m: f64,
    /// WENO parameters (epsilon and p)
    pub weno_epsilon: f64,
    /// Power parameter for WENO weights
    pub weno_p: f64,
    /// Whether to apply the limiter adaptively
    pub adaptive: bool,
    /// Tolerance for detecting troubled cells
    pub tolerance: f64,
}

impl Default for LimiterParams {
    fn default() -> Self {
        Self {
            limiter_type: LimiterType::Minmod,
            tvb_m: 1.0,
            weno_epsilon: 1e-6,
            weno_p: 2.0,
            adaptive: true,
            tolerance: 1e-4,
        }
    }
}

impl LimiterParams {
    /// Create a new set of limiter parameters
    pub fn new(limiter_type: LimiterType) -> Self {
        Self {
            limiter_type,
            ..Default::default()
        }
    }
    
    /// Set the TVB parameter
    pub fn with_tvb_m(mut self, m: f64) -> Self {
        self.tvb_m = m;
        self
    }
    
    /// Set the WENO parameters
    pub fn with_weno_params(mut self, epsilon: f64, p: f64) -> Self {
        self.weno_epsilon = epsilon;
        self.weno_p = p;
        self
    }
    
    /// Set the adaptive flag
    pub fn with_adaptive(mut self, adaptive: bool) -> Self {
        self.adaptive = adaptive;
        self
    }
    
    /// Set the tolerance for detecting troubled cells
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }
}

/// Trait for slope limiters
pub trait Limiter: Send + Sync {
    /// Apply the limiter to a DG solution
    fn limit(
        &self,
        solution: &mut DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> Result<()>;
    
    /// Check if a cell is troubled and needs limiting
    fn is_troubled_cell(
        &self,
        solution: &DGSolution,
        neighbors: &[DGSolution],
        params: &LimiterParams,
    ) -> bool;
}

/// No limiting (identity operator)
pub struct NoLimiter;

impl Limiter for NoLimiter {
    fn limit(
        &self,
        _solution: &mut DGSolution,
        _neighbors: &[DGSolution],
        _params: &LimiterParams,
    ) -> Result<()> {
        // Do nothing
        Ok(())
    }
    
    fn is_troubled_cell(
        &self,
        _solution: &DGSolution,
        _neighbors: &[DGSolution],
        _params: &LimiterParams,
    ) -> bool {
        false
    }
}

/// Minmod limiter (most diffusive)
pub struct MinmodLimiter;

impl Limiter for MinmodLimiter {
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
            
            // Compute minmod limited slopes
            for i in 0..solution.num_components {
                let du_l = u_avg[i] - u_left[i];
                let du_r = u_right[i] - u_avg[i];
                
                // Minmod function: minmod(a,b) = 0.5*(sgn(a)+sgn(b))*min(|a|,|b|)
                let slope = if du_l * du_r <= 0.0 {
                    0.0
                } else if du_l.abs() <= du_r.abs() {
                    du_l
                } else {
                    du_r
                };
                
                // Update the solution coefficients
                for j in 1..solution.coefficients.ncols() {
                    solution.coefficients[(i, j)] = 0.0;
                }
                
                // Set the linear term (if any)
                if solution.coefficients.ncols() > 1 {
                    solution.coefficients[(i, 1)] = slope;
                }
            }
        }
        
        Ok(())
    }
    
    fn is_troubled_cell(
        &self,
        solution: &DGSolution,
        neighbors: &[DGSolution],
        _params: &LimiterParams,
    ) -> bool {
        if neighbors.len() < 2 {
            return false;
        }
        
        let u_avg = solution.average();
        let u_left = neighbors[0].average();
        let u_right = neighbors[1].average();
        
        // Check for local extrema
        for i in 0..solution.num_components {
            if (u_avg[i] > u_left[i] && u_avg[i] > u_right[i]) ||
               (u_avg[i] < u_left[i] && u_avg[i] < u_right[i]) {
                return true;
            }
        }
        
        false
    }
}

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
                for j in 1..solution.coefficients.ncols() {
                    solution.coefficients[(i, j)] = 0.0;
                }
                
                // Set the linear term (if any)
                if solution.coefficients.ncols() > 1 {
                    solution.coefficients[(i, 1)] = slope;
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
            
            if du.abs() > m * h * h && (du_l * du_r <= 0.0 || du_l.abs() > m * h * h || du_r.abs() > m * h * h) {
                return true;
            }
        }
        
        false
    }
}

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
                for j in (1..solution.coefficients.ncols()).rev() {
                    // Compute the limited coefficient
                    let c = solution.coefficients[(i, j)];
                    
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
                        let c_left = if !neighbors.is_empty() && neighbors[0].coefficients.ncols() > j {
                            neighbors[0].coefficients[(i, j)]
                        } else {
                            0.0
                        };
                        
                        let c_right = if neighbors.len() > 1 && neighbors[1].coefficients.ncols() > j {
                            neighbors[1].coefficients[(i, j)]
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
                    solution.coefficients[(i, j)] = c_min;
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
                weights[0] * (1.5 * u[0] - 0.5 * u[1])
                    + weights[1] * (0.5 * u[0] + 0.5 * u[1])
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
            // For simplicity, we'll only implement the WENO limiter for the highest order term
            // A full implementation would apply WENO to all coefficients
            
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

/// Factory for creating limiter instances
pub struct LimiterFactory;

impl LimiterFactory {
    /// Create a new limiter
    pub fn create(limiter_type: LimiterType) -> Box<dyn Limiter> {
        match limiter_type {
            LimiterType::None => Box::new(NoLimiter),
            LimiterType::Minmod => Box::new(MinmodLimiter),
            LimiterType::TVB => Box::new(TVBLimiter),
            LimiterType::Moment => Box::new(MomentLimiter),
            LimiterType::WENO => Box::new(WENOLimiter::new(3)), // Default to 3rd order WENO
            _ => Box::new(MinmodLimiter), // Default to Minmod
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;
    
    #[test]
    fn test_minmod_limiter() {
        // Create a test solution with a discontinuity
        let mut solution = DGSolution::new(2, 1).unwrap();
        solution.coefficients = DMatrix::from_vec(1, 3, vec![1.0, 1.0, 0.5]);
        
        // Create left and right neighbors
        let mut left = DGSolution::new(2, 1).unwrap();
        left.coefficients = DMatrix::from_vec(1, 3, vec![0.0, 0.0, 0.0]);
        
        let mut right = DGSolution::new(2, 1).unwrap();
        right.coefficients = DMatrix::from_vec(1, 3, vec![2.0, 0.0, 0.0]);
        
        // Apply the limiter
        let limiter = MinmodLimiter;
        let mut params = LimiterParams::new(LimiterType::Minmod);
        params.adaptive = false; // Force limiting
        
        limiter.limit(&mut solution, &[left, right], &params).unwrap();
        
        // Check that the solution has been limited
        assert_eq!(solution.coefficients[(0, 0)], 1.0); // Average should be preserved
        assert_eq!(solution.coefficients[(0, 1)], 1.0); // Slope should be limited to 1.0
        assert_eq!(solution.coefficients[(0, 2)], 0.0); // Higher-order terms should be zeroed
    }
    
    #[test]
    fn test_tvb_limiter() {
        // Create a test solution with a smooth extremum
        let mut solution = DGSolution::new(2, 1).unwrap();
        solution.coefficients = DMatrix::from_vec(1, 3, vec![1.0, 0.1, 0.01]);
        
        // Create left and right neighbors
        let mut left = DGSolution::new(2, 1).unwrap();
        left.coefficients = DMatrix::from_vec(1, 3, vec![0.9, 0.1, 0.0]);
        
        let mut right = DGSolution::new(2, 1).unwrap();
        right.coefficients = DMatrix::from_vec(1, 3, vec![1.1, 0.1, 0.0]);
        
        // Apply the limiter with a large enough TVB parameter
        let limiter = TVBLimiter;
        let params = LimiterParams::new(LimiterType::TVB).with_tvb_m(1.0);
        
        limiter.limit(&mut solution, &[left, right], &params).unwrap();
        
        // Check that the solution has not been limited (smooth solution)
        assert_relative_eq!(solution.coefficients[(0, 0)], 1.0, epsilon = 1e-10);
        assert_relative_eq!(solution.coefficients[(0, 1)], 0.1, epsilon = 1e-10);
        assert_relative_eq!(solution.coefficients[(0, 2)], 0.01, epsilon = 1e-10);
    }
    
    #[test]
    fn test_moment_limiter() {
        // Create a test solution with high-order terms
        let mut solution = DGSolution::new(3, 1).unwrap();
        solution.coefficients = DMatrix::from_vec(1, 4, vec![1.0, 1.0, 0.5, 0.1]);
        
        // Create left and right neighbors
        let mut left = DGSolution::new(3, 1).unwrap();
        left.coefficients = DMatrix::from_vec(1, 4, vec![0.0, 0.0, 0.0, 0.0]);
        
        let mut right = DGSolution::new(3, 1).unwrap();
        right.coefficients = DMatrix::from_vec(1, 4, vec![2.0, 0.0, 0.0, 0.0]);
        
        // Apply the limiter
        let limiter = MomentLimiter;
        let mut params = LimiterParams::new(LimiterType::Moment);
        params.adaptive = false; // Force limiting
        
        limiter.limit(&mut solution, &[left, right], &params).unwrap();
        
        // Check that the solution has been limited
        assert_eq!(solution.coefficients[(0, 0)], 1.0); // Average should be preserved
        assert_eq!(solution.coefficients[(0, 1)], 1.0); // First moment should be limited
        assert_eq!(solution.coefficients[(0, 2)], 0.0); // Second moment should be zeroed
        assert_eq!(solution.coefficients[(0, 3)], 0.0); // Third moment should be zeroed
    }
}
