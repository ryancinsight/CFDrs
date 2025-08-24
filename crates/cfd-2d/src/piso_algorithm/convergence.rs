//! Convergence criteria for PISO algorithm

use nalgebra::RealField;
use num_traits::FromPrimitive;
use crate::fields::SimulationFields;
use cfd_core::constants::*;

/// Convergence criteria for PISO iterations
#[derive(Debug, Clone)]
pub struct ConvergenceCriteria<T: RealField + Copy> {
    /// Tolerance for velocity residual
    pub velocity_tolerance: T,
    /// Tolerance for pressure residual
    pub pressure_tolerance: T,
    /// Tolerance for continuity residual
    pub continuity_tolerance: T,
    /// Maximum iterations
    pub max_iterations: usize,
}

impl<T: RealField + Copy + FromPrimitive> Default for ConvergenceCriteria<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            velocity_tolerance: T::from_f64(1e-5).unwrap_or_else(|| {
                T::from_f64(1e-6).unwrap_or_else(T::zero)
            }),
            pressure_tolerance: T::from_f64(1e-4).unwrap_or_else(|| {
                T::from_f64(1e-5).unwrap_or_else(T::zero)
            }),
            continuity_tolerance: T::from_f64(1e-5).unwrap_or_else(|| {
                T::from_f64(1e-6).unwrap_or_else(T::zero)
            }),
        }
    }
}

/// Convergence monitor for tracking residuals
pub struct ConvergenceMonitor<T: RealField + Copy> {
    /// Velocity residual history
    pub velocity_residuals: Vec<T>,
    /// Pressure residual history
    pub pressure_residuals: Vec<T>,
    /// Continuity residual history
    pub continuity_residuals: Vec<T>,
    /// Current iteration
    pub iteration: usize,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for ConvergenceMonitor<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> ConvergenceMonitor<T> {
    /// Create new convergence monitor
    #[must_use] pub fn new() -> Self {
        Self {
            velocity_residuals: Vec::new(),
            pressure_residuals: Vec::new(),
            continuity_residuals: Vec::new(),
            iteration: 0,
        }
    }

    /// Check if converged
    pub fn is_converged(&self) -> bool {
        if self.velocity_residuals.is_empty() || 
           self.pressure_residuals.is_empty() || 
           self.continuity_residuals.is_empty() {
            return false;
        }
        
        let vel_res = self.velocity_residuals.last().copied().unwrap_or(T::one());
        let pres_res = self.pressure_residuals.last().copied().unwrap_or(T::one());
        let cont_res = self.continuity_residuals.last().copied().unwrap_or(T::one());
        
        vel_res < self.criteria.velocity_tolerance &&
        pres_res < self.criteria.pressure_tolerance &&
        cont_res < self.criteria.continuity_tolerance
    }

    /// Update residuals
    pub fn update(
        &mut self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) {
        let vel_res = self.calculate_velocity_residual(fields_old, fields_new, nx, ny);
        let pres_res = self.calculate_pressure_residual(fields_old, fields_new, nx, ny);
        let cont_res = self.calculate_continuity_residual(fields_new, nx, ny);
        
        self.velocity_residuals.push(vel_res);
        self.pressure_residuals.push(pres_res);
        self.continuity_residuals.push(cont_res);
        self.iteration += 1;
    }

    /// Calculate velocity residual (L2 norm)
    fn calculate_velocity_residual(
        &self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) -> T {
        let mut sum = T::zero();
        let mut count = 0;
        
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let du = fields_new.u.at(i, j) - fields_old.u.at(i, j);
                let dv = fields_new.v.at(i, j) - fields_old.v.at(i, j);
                sum = sum + du * du + dv * dv;
                count += 2;
            }
        }
        
        (sum / T::from_usize(count).unwrap()).sqrt()
    }

    /// Calculate pressure residual (L2 norm)
    fn calculate_pressure_residual(
        &self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) -> T {
        let mut sum = T::zero();
        let mut count = 0;
        
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                let dp = fields_new.p.at(i, j) - fields_old.p.at(i, j);
                sum += dp * dp;
                count += 1;
            }
        }
        
        (sum / T::from_usize(count).unwrap()).sqrt()
    }

    /// Calculate continuity residual (mass imbalance)
    fn calculate_continuity_residual(
        &self,
        fields: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) -> T {
        let mut max_imbalance = T::zero();
        
        for i in 1..nx-1 {
            for j in 1..ny-1 {
                // Continuity equation check (du/dx + dv/dy = 0)
                let dudx = (fields.u.at(i+1, j) - fields.u.at(i-1, j)) / T::from_f64(2.0).unwrap();
                let dvdy = (fields.v.at(i, j+1) - fields.v.at(i, j-1)) / T::from_f64(2.0).unwrap();
                let imbalance = (dudx + dvdy).abs();
                
                if imbalance > max_imbalance {
                    max_imbalance = imbalance;
                }
            }
        }
        
        max_imbalance
    }
}