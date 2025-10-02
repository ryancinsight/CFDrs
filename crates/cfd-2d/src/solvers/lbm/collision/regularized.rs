//! Regularized collision operator for improved stability

#![allow(dead_code)]

use super::traits::CollisionOperator;
use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Regularized collision operator
/// Based on Latt & Chopard (2006)
pub struct RegularizedCollision<T: RealField + Copy> {
    /// Relaxation time
    tau: T,
    /// Inverse relaxation time
    omega: T,
}

impl<T: RealField + Copy + FromPrimitive> RegularizedCollision<T> {
    /// Create new regularized collision operator
    pub fn new(tau: T) -> Self {
        Self {
            tau,
            omega: T::one() / tau,
        }
    }

    /// Compute non-equilibrium stress tensor
    fn compute_stress_tensor(f: &[T; 9], f_eq: &[T; 9]) -> [[T; 2]; 2] {
        let mut pi = [[T::zero(); 2]; 2];

        // Compute off-equilibrium part
        for q in 0..9 {
            let f_neq = f[q] - f_eq[q];
            let lattice_vel = D2Q9::VELOCITIES[q];
            let cx = T::from_i32(lattice_vel.0).unwrap_or_else(T::zero);
            let cy = T::from_i32(lattice_vel.1).unwrap_or_else(T::zero);

            // Pi_αβ = Σ_i c_iα c_iβ f_i^neq
            pi[0][0] += cx * cx * f_neq;
            pi[0][1] += cx * cy * f_neq;
            pi[1][0] += cy * cx * f_neq;
            pi[1][1] += cy * cy * f_neq;
        }

        pi
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for RegularizedCollision<T> {
    fn collide(
        &self,
        f: &mut Vec<Vec<[T; 9]>>,
        density: &[Vec<T>],
        velocity: &[Vec<[T; 2]>],
    ) {
        let ny = f.len();
        let nx = if ny > 0 { f[0].len() } else { 0 };

        for j in 0..ny {
            for i in 0..nx {
                let rho = density[j][i];
                let u = velocity[j][i];

                // Compute equilibrium
                let mut f_eq = [T::zero(); 9];
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
                    let lattice_vel = &D2Q9::VELOCITIES[q];
                    f_eq[q] = equilibrium(rho, &u, q, weight, lattice_vel);
                }

                // Compute non-equilibrium stress
                let pi = Self::compute_stress_tensor(&f[j][i], &f_eq);

                // Regularized collision
                for q in 0..9 {
                    // Standard BGK collision
                    let f_bgk = f[j][i][q] - self.omega * (f[j][i][q] - f_eq[q]);

                    // Regularization term
                    let lattice_vel = D2Q9::VELOCITIES[q];
                    let c = [
                        T::from_i32(lattice_vel.0).unwrap_or_else(T::zero),
                        T::from_i32(lattice_vel.1).unwrap_or_else(T::zero),
                    ];
                    let reg_term = self.regularization_term(c, pi);

                    f[j][i][q] = f_bgk + reg_term;
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.tau
    }

    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 =
            T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));
        cs2 * dx * dx * (self.tau - half) / dt
    }
}

impl<T: RealField + Copy + FromPrimitive> RegularizedCollision<T> {
    #[allow(dead_code)]
    fn regularization_term(&self, c: [T; 2], pi: [[T; 2]; 2]) -> T {
        // Compute regularization term Q_i^(1) = w_i * H_i : Pi^(1)
        // Based on Latt & Chopard (2006) Phys. Rev. E 74, 026701
        // H_i = (c_i ⊗ c_i - cs^2 I) where cs^2 = 1/3 for D2Q9
        let cs2 =
            T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));

        // Tensor contraction H_i : Pi^(1)
        let term = (c[0] * c[0] - cs2) * pi[0][0]
            + (c[1] * c[1] - cs2) * pi[1][1]
            + c[0] * c[1] * (pi[0][1] + pi[1][0]);

        // Apply relaxation factor omega/2
        term * self.omega / (T::one() + T::one())
    }
}
