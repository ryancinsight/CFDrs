//! Macroscopic quantity computations for LBM.
//!
//! This module handles the computation of macroscopic quantities
//! (density, velocity, pressure) from distribution functions.

use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;
use num_traits::FromPrimitive;
/// Container for macroscopic quantities
#[derive(Debug, Clone)]
pub struct MacroscopicQuantities<T: RealField + Copy> {
    /// Density field
    pub density: Vec<Vec<T>>,
    /// Velocity field
    pub velocity: Vec<Vec<[T; 2]>>,
    /// Pressure field (optional)
    pub pressure: Option<Vec<Vec<T>>>,
}
impl<T: RealField + Copy + FromPrimitive> MacroscopicQuantities<T> {
    /// Create new macroscopic quantities container
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            density: vec![vec![T::one(); nx]; ny],
            velocity: vec![vec![[T::zero(), T::zero()]; nx]; ny],
            pressure: None,
        }
    }
    /// Enable pressure computation
    pub fn with_pressure(mut self) -> Self {
        let ny = self.density.len();
        let nx = if ny > 0 { self.density[0].len() } else { 0 };
        self.pressure = Some(vec![vec![T::zero(); nx]; ny]);
        self
    /// Update all macroscopic quantities from distribution functions
    pub fn update_from_distributions(&mut self, f: &Vec<Vec<[T; 9]>>) {
        let ny = f.len();
        let nx = if ny > 0 { f[0].len() } else { 0 };
        for j in 0..ny {
            for i in 0..nx {
                // Compute density
                self.density[j][i] = compute_density(&f[j][i]);
                // Compute velocity
                self.velocity[j][i] = compute_velocity(&f[j][i], self.density[j][i]);
                // Compute pressure if enabled
                if let Some(ref mut pressure) = self.pressure {
                    pressure[j][i] = compute_pressure(self.density[j][i]);
                }
            }
/// Compute density from distribution functions
pub fn compute_density<T: RealField + Copy>(f_node: &[T; 9]) -> T {
    let mut rho = T::zero();
    for q in 0..9 {
        rho += f_node[q];
    rho
/// Compute velocity from distribution functions
    }

pub fn compute_velocity<T: RealField + Copy + FromPrimitive>(
    f_node: &[T; 9],
    density: T,
) -> [T; 2] {
    let mut ux = T::zero();
    let mut uy = T::zero();
        let (ex, ey) = D2Q9::VELOCITIES[q];
        let ex_t = T::from_i32(ex).unwrap_or_else(T::zero);
        let ey_t = T::from_i32(ey).unwrap_or_else(T::zero);
        ux += ex_t * f_node[q];
        uy += ey_t * f_node[q];
    [ux / density, uy / density]
/// Compute pressure from density (ideal gas equation of state)
    }

pub fn compute_pressure<T: RealField + Copy + FromPrimitive>(density: T) -> T {
    // For LBM, pressure = cs^2 * density where cs^2 = 1/3
    let cs2 = T::from_f64(1.0 / 3.0).unwrap_or_else(T::zero);
    cs2 * density
/// Compute momentum from distribution functions
    }

pub fn compute_momentum<T: RealField + Copy + FromPrimitive>(f_node: &[T; 9]) -> [T; 2] {
    let mut px = T::zero();
    let mut py = T::zero();
        px += ex_t * f_node[q];
        py += ey_t * f_node[q];
    [px, py]
/// Compute stress tensor from non-equilibrium distributions
    }

pub fn compute_stress_tensor<T: RealField + Copy + FromPrimitive>(
    f_eq: &[T; 9],
) -> [[T; 2]; 2] {
    let mut stress = [[T::zero(); 2]; 2];
        let f_neq = f_node[q] - f_eq[q];
        stress[0][0] += ex_t * ex_t * f_neq;
        stress[0][1] += ex_t * ey_t * f_neq;
        stress[1][0] += ey_t * ex_t * f_neq;
        stress[1][1] += ey_t * ey_t * f_neq;
    stress
/// Compute kinetic energy density
    }

pub fn compute_kinetic_energy<T: RealField + Copy>(density: T, velocity: &[T; 2]) -> T {
    let half = T::from_f64(0.5).unwrap_or_else(T::zero);
    half * density * (velocity[0] * velocity[0] + velocity[1] * velocity[1])
/// Compute vorticity (2D)
    }

pub fn compute_vorticity<T: RealField + Copy>(
    velocity: &Vec<Vec<[T; 2]>>,
    dx: T,
    dy: T,
    i: usize,
    j: usize,
) -> T {
    let ny = velocity.len();
    let nx = if ny > 0 { velocity[0].len() } else { 0 };
    if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
        return T::zero();
    // Central differences
    let dv_dx = (velocity[j][i + 1][1] - velocity[j][i - 1][1])
        / (T::from_f64(2.0).unwrap_or_else(T::zero) * dx);
    let du_dy = (velocity[j + 1][i][0] - velocity[j - 1][i][0])
        / (T::from_f64(2.0).unwrap_or_else(T::zero) * dy);
    dv_dx - du_dy
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_density_computation() {
        let mut f_node = [0.0_f64; 9];
        for q in 0..9 {
            f_node[q] = 0.1;
        let density = compute_density(&f_node);
        assert!((density - 0.9).abs() < 1e-10);
    }

    fn test_velocity_computation() {
        let mut f_node = [0.1_f64; 9];
        // Add some momentum in x-direction
        f_node[1] = 0.15; // East
        f_node[3] = 0.05; // West
        let velocity = compute_velocity(&f_node, density);
        // Should have positive x-velocity
        assert!(velocity[0] > 0.0);
        // Should have zero y-velocity (symmetric)
        assert!(velocity[1].abs() < 1e-10);
    }

    fn test_pressure_computation() {
        let density = 1.5_f64;
        let pressure = compute_pressure(density);
        let expected = (1.0 / 3.0) * density;
        assert!((pressure - expected).abs() < 1e-10);

    }


}
}
}
}
}
}
}
}
}
}
}
}
}
