//! Collision operators for LBM.
//!
//! This module implements various collision models including
//! BGK (Bhatnagar-Gross-Krook) single relaxation time.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use crate::solvers::lbm::lattice::{D2Q9, equilibrium};

/// Trait for collision operators
pub trait CollisionOperator<T: RealField + Copy> {
    /// Apply collision step to distribution functions
    fn collide(
        &self,
        f: &mut Vec<Vec<[T; 9]>>,
        density: &Vec<Vec<T>>,
        velocity: &Vec<Vec<[T; 2]>>,
    );
    
    /// Get relaxation time
    fn tau(&self) -> T;
}

/// BGK (Bhatnagar-Gross-Krook) collision operator
pub struct BgkCollision<T: RealField + Copy> {
    /// Relaxation time
    tau: T,
    /// Inverse relaxation time (omega = 1/tau)
    omega: T,
}

impl<T: RealField + Copy + FromPrimitive> BgkCollision<T> {
    /// Create new BGK collision operator
    pub fn new(tau: T) -> Self {
        Self {
            tau,
            omega: T::one() / tau,
        }
    }
    
    /// Create from kinematic viscosity
    pub fn from_viscosity(nu: T, dt: T, dx: T) -> Self {
        let cs2 = T::from_f64(1.0 / 3.0).unwrap_or_else(T::zero);
        let tau = T::from_f64(0.5).unwrap_or_else(T::zero) + nu * dt / (cs2 * dx * dx);
        Self::new(tau)
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for BgkCollision<T> {
    fn collide(
        &self,
        f: &mut Vec<Vec<[T; 9]>>,
        density: &Vec<Vec<T>>,
        velocity: &Vec<Vec<[T; 2]>>,
    ) {
        let ny = f.len();
        let nx = if ny > 0 { f[0].len() } else { 0 };
        
        for j in 0..ny {
            for i in 0..nx {
                let rho = density[j][i];
                let u = velocity[j][i];
                
                // Compute equilibrium distributions
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
                    let lattice_vel = &D2Q9::VELOCITIES[q];
                    let f_eq = equilibrium(rho, &u, q, weight, lattice_vel);
                    
                    // BGK collision: f = f - omega * (f - f_eq)
                    f[j][i][q] = f[j][i][q] - self.omega * (f[j][i][q] - f_eq);
                }
            }
        }
    }
    
    fn tau(&self) -> T {
        self.tau
    }
}

/// MRT (Multiple Relaxation Time) collision operator
pub struct MrtCollision<T: RealField + Copy> {
    /// Relaxation times for different moments
    relaxation_times: [T; 9],
}

impl<T: RealField + Copy + FromPrimitive> MrtCollision<T> {
    /// Create new MRT collision operator
    pub fn new(relaxation_times: [T; 9]) -> Self {
        Self { relaxation_times }
    }
    
    /// Create with default parameters
    pub fn default_parameters(tau: T) -> Self {
        let mut times = [T::one(); 9];
        times[7] = tau; // Shear viscosity
        times[8] = tau; // Shear viscosity
        Self::new(times)
    }
}

// Note: Full MRT implementation would require moment transformation matrices
// This is a placeholder for the interface

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for MrtCollision<T> {
    fn collide(
        &self,
        f: &mut Vec<Vec<[T; 9]>>,
        density: &Vec<Vec<T>>,
        velocity: &Vec<Vec<[T; 2]>>,
    ) {
        // MRT collision requires moment space transformation
        // For now, fallback to BGK-like behavior with primary relaxation time
        let bgk = BgkCollision::new(self.relaxation_times[7]);
        bgk.collide(f, density, velocity);
    }
    
    fn tau(&self) -> T {
        self.relaxation_times[7] // Return shear viscosity relaxation time
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bgk_creation() {
        let tau = 0.6_f64;
        let bgk = BgkCollision::new(tau);
        assert_eq!(bgk.tau(), tau);
        assert!((bgk.omega - 1.0 / tau).abs() < 1e-10);
    }
    
    #[test]
    fn test_bgk_from_viscosity() {
        let nu = 0.1_f64;
        let dt = 1.0;
        let dx = 1.0;
        let bgk = BgkCollision::from_viscosity(nu, dt, dx);
        
        let expected_tau = 0.5 + nu * dt / ((1.0 / 3.0) * dx * dx);
        assert!((bgk.tau() - expected_tau).abs() < 1e-10);
    }
}