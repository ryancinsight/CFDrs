//! Shan-Chen Pseudopotential Multiphase Model for LBM.
//!
//! Simulates multi-component or multiphase flows (liquid-gas or liquid-liquid)
//! by incorporating a cohesive/adhesive pseudo-force.
//!
//! # Theorem — Shan-Chen Cohesive Force (Shan & Chen 1993)
//!
//! Let $\psi(\vec{x})$ be the effective mass (pseudopotential) at site $\vec{x}$.
//! The nearest-neighbor interaction force is:
//!
//! ```text
//! \vec{F}(\vec{x}) = -G \cdot \psi(\vec{x}) \sum_i w_i \psi(\vec{x} + \vec{e}_i) \vec{e}_i
//! ```
//!
//! where $G$ is the interaction strength ($G < 0$ for attraction/cohesion).
//! To couple this force to the fluid without breaking momentum conservation:
//!
//! ```text
//! \vec{u}_{eq} = \vec{u} + \frac{\tau \vec{F}}{\rho}
//! ```
//!
//! **Proof of Phase Separation**:
//! A Taylor expansion of the force yields a non-ideal equation of state:
//! $P = \rho c_s^2 + \frac{G c_s^2}{2} \psi^2$. For $G < G_{critical}$, the EoS exhibits
//! a van der Waals-like loop with $\partial P / \partial \rho < 0$ (spinodal decomposition).
//! The fluid spontaneously separates into a high-density (liquid) and low-density (vapor)
//! phase, satisfying the Maxwell equal-area construction, which dictates the equilibrium
//! coexisting densities.
//!
//! # References
//! - Shan, X., & Chen, H. (1993). Lattice Boltzmann model for simulating flows with
//!   multiple phases and components. *Physical Review E*, 47(3), 1815.

use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Configuration for Shan-Chen pseudopotential multiphase simulation
#[derive(Debug, Clone, Copy)]
pub struct ShanChenMultiphase<T: RealField + Copy> {
    /// Interaction strength G ($G < 0$ for attraction). Usually around -5.0 to -6.0.
    pub g_c: T,
    /// Reference density $\rho_0$ for the exponential pseudopotential form
    pub rho_0: T,
}

impl<T: RealField + Copy + Float + FromPrimitive> ShanChenMultiphase<T> {
    /// Compute the exponential effective mass $\psi(\rho)$
    ///
    /// $\psi(\rho) = \rho_0 [1 - \exp(-\rho / \rho_0)]$
    #[inline]
    #[must_use]
    pub fn pseudopotential(&self, density: T) -> T {
        let one = T::one();
        self.rho_0 * (one - Float::exp(-density / self.rho_0))
    }

    /// Compute the local cohesive force field $\vec{F}(x,y)$ for the entire domain.
    ///
    /// Requires calculating the gradient-like sum of the pseudopotential over
    /// all neighboring lattice nodes, assuming periodic or bounce-back walls.
    pub fn compute_cohesive_force(
        &self,
        density: &[T],
        nx: usize,
        ny: usize,
    ) -> (Vec<T>, Vec<T>) {
        let mut f_x = vec![T::zero(); nx * ny];
        let mut f_y = vec![T::zero(); nx * ny];

        let zero = T::zero();

        // Precompute psi field
        let mut psi_field = vec![zero; nx * ny];
        for i in 0..(nx * ny) {
            psi_field[i] = self.pseudopotential(density[i]);
        }

        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;
                let psi = psi_field[cell];

                let mut sum_x = zero;
                let mut sum_y = zero;

                for q in 1..9 { // Skip q=0 (rest particle has e_i = 0)
                    let ex = D2Q9::VELOCITIES[q].0;
                    let ey = D2Q9::VELOCITIES[q].1;
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap();

                    // Calculate neighbor coordinates (assuming periodic boundaries for bulk force)
                    let nb_x = (i as i32 + ex).rem_euclid(nx as i32) as usize;
                    let nb_y = (j as i32 + ey).rem_euclid(ny as i32) as usize;
                    let nb_cell = nb_y * nx + nb_x;

                    let nb_psi = psi_field[nb_cell];
                    let e_x_t = T::from_i32(ex).unwrap();
                    let e_y_t = T::from_i32(ey).unwrap();

                    sum_x += weight * nb_psi * e_x_t;
                    sum_y += weight * nb_psi * e_y_t;
                }

                f_x[cell] = -self.g_c * psi * sum_x;
                f_y[cell] = -self.g_c * psi * sum_y;
            }
        }

        (f_x, f_y)
    }

    /// Shift the macroscopic velocity to incorporate the cohesive force
    /// before evaluating the equilibrium distribution in collision:
    ///
    /// $\vec{u}_{eq} = \vec{u} + \frac{\tau \vec{F}}{\rho}$
    #[inline]
    #[must_use]
    pub fn shift_equilibrium_velocity(
        u: [T; 2],
        force: [T; 2],
        tau: T,
        rho: T,
    ) -> [T; 2] {
        if rho > T::from_f64(1e-12).unwrap() {
            [
                u[0] + (tau * force[0]) / rho,
                u[1] + (tau * force[1]) / rho,
            ]
        } else {
            u
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Theorem: The cohesive interaction force for a uniform density field must be
    /// identically zero due to the isotropy of the lattice weights $\sum w_i \vec{e}_i = 0$.
    #[test]
    fn uniform_density_yields_zero_force() {
        let sc = ShanChenMultiphase::<f64> {
            g_c: -5.0,
            rho_0: 1.0,
        };
        let nx = 10;
        let ny = 10;
        let density = vec![1.5_f64; nx * ny];

        let (fx, fy) = sc.compute_cohesive_force(&density, nx, ny);

        for i in 0..(nx * ny) {
            assert!(fx[i].abs() < 1e-12, "Force X must be zero: {}", fx[i]);
            assert!(fy[i].abs() < 1e-12, "Force Y must be zero: {}", fy[i]);
        }
    }

    /// Theorem: A density gradient must produce a force directed toward the
    /// higher density region (for cohesive $G < 0$) to drive phase separation.
    #[test]
    fn density_gradient_yields_attractive_force() {
        let sc = ShanChenMultiphase::<f64> {
            g_c: -5.0, // Attractive
            rho_0: 1.0,
        };
        let nx = 5;
        let ny = 5;
        let mut density = vec![1.0_f64; nx * ny];
        
        // Create an artificial heavy droplet in the center
        density[2 * nx + 2] = 2.0;

        let (fx, _fy) = sc.compute_cohesive_force(&density, nx, ny);

        // The point directly to the left (1, 2) should experience a force pushing Right (+x)
        let left_idx = 2 * nx + 1;
        assert!(fx[left_idx] > 0.0, "Force should point towards droplet");

        // The point directly to the right (3, 2) should experience a force pushing Left (-x)
        let right_idx = 2 * nx + 3;
        assert!(fx[right_idx] < 0.0, "Force should point towards droplet");
    }
}
