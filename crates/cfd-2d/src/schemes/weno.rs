//! Weighted Essentially Non-Oscillatory (WENO) schemes

use nalgebra::RealField;
use num_traits::FromPrimitive;
use super::{Grid2D, SpatialDiscretization, constants};

/// Fifth-order WENO scheme
pub struct WENO5<T: RealField + Copy> {
    epsilon: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + Copy> WENO5<T> {
    /// Create new WENO5 scheme
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(constants::WENO_EPSILON).unwrap_or_else(T::zero),
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Compute smoothness indicators
    fn smoothness_indicators(&self, v: &[T; 5]) -> [T; 3] {
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::zero);
        let three = T::from_f64(3.0).unwrap_or_else(T::zero);
        
        // Beta_0
        let beta0 = thirteen / T::from_f64(12.0).unwrap_or_else(T::zero) 
            * (v[0] - T::from_f64(2.0).unwrap_or_else(T::zero) * v[1] + v[2]).powi(2)
            + T::from_f64(0.25).unwrap_or_else(T::zero) 
            * (v[0] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[1] + three * v[2]).powi(2);
        
        // Beta_1
        let beta1 = thirteen / T::from_f64(12.0).unwrap_or_else(T::zero)
            * (v[1] - T::from_f64(2.0).unwrap_or_else(T::zero) * v[2] + v[3]).powi(2)
            + T::from_f64(0.25).unwrap_or_else(T::zero) * (v[1] - v[3]).powi(2);
        
        // Beta_2
        let beta2 = thirteen / T::from_f64(12.0).unwrap_or_else(T::zero)
            * (v[2] - T::from_f64(2.0).unwrap_or_else(T::zero) * v[3] + v[4]).powi(2)
            + T::from_f64(0.25).unwrap_or_else(T::zero)
            * (three * v[2] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3] + v[4]).powi(2);
        
        [beta0, beta1, beta2]
    }
    
    /// Compute WENO weights
    fn weno_weights(&self, beta: &[T; 3]) -> [T; 3] {
        let d0 = T::from_f64(constants::WENO5_WEIGHTS[0]).unwrap_or_else(T::zero);
        let d1 = T::from_f64(constants::WENO5_WEIGHTS[1]).unwrap_or_else(T::zero);
        let d2 = T::from_f64(constants::WENO5_WEIGHTS[2]).unwrap_or_else(T::zero);
        
        let alpha0 = d0 / (self.epsilon + beta[0]).powi(2);
        let alpha1 = d1 / (self.epsilon + beta[1]).powi(2);
        let alpha2 = d2 / (self.epsilon + beta[2]).powi(2);
        
        let sum = alpha0 + alpha1 + alpha2;
        
        [alpha0 / sum, alpha1 / sum, alpha2 / sum]
    }
}

impl<T: RealField + FromPrimitive + Copy> SpatialDiscretization<T> for WENO5<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // Extract stencil
        let v = [
            grid.data[(i - 2, j)],
            grid.data[(i - 1, j)],
            grid.data[(i, j)],
            grid.data[(i + 1, j)],
            grid.data[(i + 2, j)],
        ];
        
        // Compute smoothness indicators
        let beta = self.smoothness_indicators(&v);
        
        // Compute weights
        let w = self.weno_weights(&beta);
        
        // Compute flux
        let f0 = v[0] / T::from_f64(3.0).unwrap_or_else(T::zero)
            - T::from_f64(7.0).unwrap_or_else(T::zero) * v[1] / T::from_f64(6.0).unwrap_or_else(T::zero)
            + T::from_f64(11.0).unwrap_or_else(T::zero) * v[2] / T::from_f64(6.0).unwrap_or_else(T::zero);
        
        let f1 = -v[1] / T::from_f64(6.0).unwrap_or_else(T::zero)
            + T::from_f64(5.0).unwrap_or_else(T::zero) * v[2] / T::from_f64(6.0).unwrap_or_else(T::zero)
            + v[3] / T::from_f64(3.0).unwrap_or_else(T::zero);
        
        let f2 = v[2] / T::from_f64(3.0).unwrap_or_else(T::zero)
            + T::from_f64(5.0).unwrap_or_else(T::zero) * v[3] / T::from_f64(6.0).unwrap_or_else(T::zero)
            - v[4] / T::from_f64(6.0).unwrap_or_else(T::zero);
        
        (w[0] * f0 + w[1] * f1 + w[2] * f2) / grid.dx
    }
    
    fn order(&self) -> usize {
        5
    }
    
    fn is_conservative(&self) -> bool {
        true
    }
}