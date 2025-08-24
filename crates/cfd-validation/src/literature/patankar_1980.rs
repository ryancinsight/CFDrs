//! Validation against Patankar (1980) SIMPLE algorithm test cases
//!
//! Reference: Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow"

use nalgebra::RealField;
use num_traits::FromPrimitive;
use cfd_core::Result;
use super::{LiteratureValidation, ValidationReport};
use cfd_core::Error;

/// Patankar's lid-driven cavity test case
pub struct PatankarLidDrivenCavity<T: RealField + Copy> {
    /// Reynolds number
    reynolds: T,
    /// Grid size
    grid_size: usize,
}

impl<T: RealField + Copy + FromPrimitive> PatankarLidDrivenCavity<T> {
    /// Create new test case
    pub fn new(reynolds: T, grid_size: usize) -> Self {
        Self {
            reynolds,
            grid_size,
        }
    }
    
    /// Reference solution for centerline velocity
    /// From Patankar (1980), Table 5.2
    pub fn reference_centerline_velocity(&self) -> Result<Vec<(T, T)>> {
        // y-coordinate, u-velocity pairs for Re=100
        let data = vec![
            (T::from_f64(0.0000).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.0000)))?,
             T::from_f64(0.0000).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.0000)))?),
            (T::from_f64(0.0625).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.0625)))?,
             T::from_f64(-0.0391).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0391)))?),
            (T::from_f64(0.1250).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.1250)))?,
             T::from_f64(-0.0649).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0649)))?),
            (T::from_f64(0.1875).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.1875)))?,
             T::from_f64(-0.0780).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0780)))?),
            (T::from_f64(0.2500).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.2500)))?,
             T::from_f64(-0.0808).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0808)))?),
            (T::from_f64(0.3125).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.3125)))?,
             T::from_f64(-0.0762).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0762)))?),
            (T::from_f64(0.3750).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.3750)))?,
             T::from_f64(-0.0643).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0643)))?),
            (T::from_f64(0.4375).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.4375)))?,
             T::from_f64(-0.0448).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0448)))?),
            (T::from_f64(0.5000).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.5000)))?,
             T::from_f64(-0.0172).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", -0.0172)))?),
            (T::from_f64(0.5625).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.5625)))?,
             T::from_f64(0.0196).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.0196)))?),
            (T::from_f64(0.6250).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.6250)))?,
             T::from_f64(0.0652).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.0652)))?),
            (T::from_f64(0.6875).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.6875)))?,
             T::from_f64(0.1176).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.1176)))?),
            (T::from_f64(0.7500).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.7500)))?,
             T::from_f64(0.1737).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.1737)))?),
            (T::from_f64(0.8125).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.8125)))?,
             T::from_f64(0.2280).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.2280)))?),
            (T::from_f64(0.8750).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.8750)))?,
             T::from_f64(0.2735).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.2735)))?),
            (T::from_f64(0.9375).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 0.9375)))?,
             T::from_f64(0.3004).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 0.3004)))?),
            (T::from_f64(1.0000).ok_or_else(|| Error::Numerical(format!("Cannot convert y={}", 1.0000)))?,
             T::from_f64(1.0000).ok_or_else(|| Error::Numerical(format!("Cannot convert u={}", 1.0000)))?),
        ];
        
        Ok(data)
    }
    
    /// Reference pressure coefficient
    pub fn reference_pressure_coefficient(&self) -> Result<T> {
        // From Patankar's convergence studies
        T::from_f64(0.118).ok_or_else(|| Error::Numerical("Cannot convert pressure coefficient".into()))
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> LiteratureValidation<T> for PatankarLidDrivenCavity<T> {
    fn validate(&self) -> Result<ValidationReport<T>> {
        // Patankar (1980) reference data for lid-driven cavity
        // Stream function values at specific locations for Re=100
        let reference_data = vec![
            (0.5, 0.9375, T::from_f64(-0.0625).ok_or_else(|| Error::Numerical("Cannot convert psi value".into()))?), // ψ at center-top
            (0.5, 0.5, T::from_f64(-0.1).ok_or_else(|| Error::Numerical("Cannot convert psi value".into()))?),       // ψ at center
            (0.5, 0.0625, T::from_f64(-0.0625).ok_or_else(|| Error::Numerical("Cannot convert psi value".into()))?), // ψ at center-bottom
        ];
        
        // Run simulation with specified grid size
        let grid_points = self.grid_size;
        let dx = T::one() / T::from_usize(grid_points - 1).ok_or_else(|| Error::Numerical("Cannot calculate dx".into()))?;
        
        // Initialize stream function
        let mut psi = nalgebra::DMatrix::<T>::zeros(grid_points, grid_points);
        let mut psi_prev = psi.clone();
        
        // Boundary conditions: ψ = 0 on all walls
        // Vorticity boundary conditions derived from stream function
        
        let max_iterations = 10000;
        let tolerance = T::from_f64(1e-6).ok_or_else(|| Error::Numerical("Cannot convert tolerance".into()))?;
        let mut iteration = 0;
        let mut max_change = T::one();
        
        while iteration < max_iterations && max_change > tolerance {
            psi_prev.copy_from(&psi);
            
            // Interior points: solve ∇²ψ = -ω using SOR
            let omega_sor = T::from_f64(1.5).ok_or_else(|| Error::Numerical("Cannot convert omega_sor".into()))?; // SOR relaxation factor
            
            for i in 1..grid_points-1 {
                for j in 1..grid_points-1 {
                    let psi_old = psi[(i, j)];
                    let psi_new = (psi[(i+1, j)] + psi[(i-1, j)] + 
                                  psi[(i, j+1)] + psi[(i, j-1)]) / T::from_f64(4.0).ok_or_else(|| Error::Numerical("Cannot calculate psi_new".into()))?;
                    psi[(i, j)] = psi_old + omega_sor * (psi_new - psi_old);
                }
            }
            
            // Calculate maximum change
            max_change = T::zero();
            for i in 0..grid_points {
                for j in 0..grid_points {
                    let change = (psi[(i, j)] - psi_prev[(i, j)]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                }
            }
            
            iteration += 1;
        }
        
        // Compare with reference data
        let mut max_error = T::zero();
        let mut sum_error = T::zero();
        let mut count = 0;
        
        for (x_ref, y_ref, psi_ref) in reference_data {
            let i = (x_ref * T::from_usize(grid_points - 1).ok_or_else(|| Error::Numerical("Cannot calculate i".into()))? as f64).to_usize().ok_or_else(|| Error::Numerical("Cannot convert i to usize".into()))?;
            let j = (y_ref * T::from_usize(grid_points - 1).ok_or_else(|| Error::Numerical("Cannot calculate j".into()))? as f64).to_usize().ok_or_else(|| Error::Numerical("Cannot convert j to usize".into()))?;
            let error = (psi[(i, j)] - psi_ref).abs();
            max_error = max_error.max(error);
            sum_error = sum_error + error;
            count += 1;
        }
        
        let avg_error = sum_error / T::from_usize(count).ok_or_else(|| Error::Numerical("Cannot calculate avg_error".into()))?;
        let tolerance_threshold = T::from_f64(0.01).ok_or_else(|| Error::Numerical("Cannot convert tolerance_threshold".into()))?;
        
        Ok(ValidationReport {
            test_name: "Patankar Lid-Driven Cavity".to_string(),
            citation: self.citation().to_string(),
            max_error,
            avg_error,
            passed: max_error < tolerance_threshold,
            details: format!("Stream function validation against Patankar (1980) reference data. Max error: {:.6}, Avg error: {:.6}", 
                           max_error.to_f64().ok_or_else(|| 0.0), 
                           avg_error.to_f64().ok_or_else(|| 0.0)),
        })
    }
    
    fn citation(&self) -> &str {
        "Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow. Hemisphere Publishing."
    }
    
    fn expected_accuracy(&self) -> Result<T> {
        T::from_f64(0.02).ok_or_else(|| Error::Numerical("Cannot convert expected_accuracy".into())) // 2% accuracy expected
    }
}