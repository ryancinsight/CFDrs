//! Regularized collision operator for improved stability

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
            let c = D2Q9::velocity_vector::<T>(q);

            // Pi_αβ = Σ_i c_iα c_iβ f_i^neq
            pi[0][0] = pi[0][0] + c[0] * c[0] * f_neq;
            pi[0][1] = pi[0][1] + c[0] * c[1] * f_neq;
            pi[1][0] = pi[1][0] + c[1] * c[0] * f_neq;
            pi[1][1] = pi[1][1] + c[1] * c[1] * f_neq;
        }

        pi
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for RegularizedCollision<T> {
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

                // Compute equilibrium
                let mut f_eq = [T::zero(); 9];
                for q in 0..9 {
                    f_eq[q] = equilibrium::<T, D2Q9>(q, rho, u);
                }

                // Compute non-equilibrium stress
                let pi = Self::compute_stress_tensor(&f[j][i], &f_eq);

                // Regularized collision
                for q in 0..9 {
                    // Standard BGK collision
                    let f_bgk = f[j][i][q] - self.omega * (f[j][i][q] - f_eq[q]);

                    // Regularization term (simplified)
                    let c = D2Q9::velocity_vector::<T>(q);
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
    fn regularization_term(&self, c: [T; 2], pi: [[T; 2]; 2]) -> T {
        // Simplified regularization term
        // Full implementation would include proper tensor contractions
        let cs2 =
            T::from_f64(1.0 / 3.0).unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));

        let term = (c[0] * c[0] - cs2) * pi[0][0]
            + (c[1] * c[1] - cs2) * pi[1][1]
            + c[0] * c[1] * (pi[0][1] + pi[1][0]);

        term * self.omega / (T::one() + T::one())
    }
}
