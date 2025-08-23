//! Fluid properties for FEM calculations

use nalgebra::{RealField, Matrix3};
use num_traits::FromPrimitive;

/// Fluid properties for FEM calculations
#[derive(Debug, Clone)]
pub struct FluidProperties<T: RealField + Copy> {
    /// Dynamic viscosity
    pub mu: T,
    /// Density
    pub rho: T,
    /// Kinematic viscosity (mu/rho)
    pub nu: T,
}

impl<T: RealField + FromPrimitive + Copy> FluidProperties<T> {
    /// Create new fluid properties
    pub fn new(mu: T, rho: T) -> Self {
        let nu = mu / rho;
        Self { mu, rho, nu }
    }
    
    /// Calculate stress tensor from strain rate
    /// σ = -pI + 2μ ε̇
    pub fn stress_tensor(&self, pressure: T, strain_rate: &Matrix3<T>) -> Matrix3<T> {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let mut stress = strain_rate * (two * self.mu);
        
        // Add pressure contribution to diagonal
        stress[(0, 0)] -= pressure;
        stress[(1, 1)] -= pressure;
        stress[(2, 2)] -= pressure;
        
        stress
    }
    
    /// Calculate strain rate tensor from velocity gradient
    /// ε̇ = 0.5 * (∇u + ∇u^T)
    pub fn strain_rate_tensor(&self, velocity_gradient: &Matrix3<T>) -> Matrix3<T> {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        (velocity_gradient + velocity_gradient.transpose()) * half
    }
}