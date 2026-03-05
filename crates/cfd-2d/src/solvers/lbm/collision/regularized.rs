//! Regularized collision operator for improved stability.
//!
//! # Theorem — Regularization (Latt & Chopard 2006)
//!
//! **Statement**: The regularized BGK operator reconstructs the off-equilibrium
//! distribution as
//! $f_q^{(1)\text{reg}} = \frac{w_q}{2 c_s^4} \mathbf{e}_q \mathbf{e}_q : \Pi^{(1)}$
//! from the non-equilibrium stress $\Pi^{(1)} = -(1/2)\tau \sum_q \mathbf{e}_q \mathbf{e}_q f_q^{neq}$,
//! and applies BGK relaxation only on this regularized part.
//!
//! **Properties**:
//! 1. The regularization step removes high-order ghost modes that cause instability.
//! 2. At steady state $f_q^{(1)} \equiv 0$ and the regularized scheme reduces to
//!    pure equilibrium BGK.
//! 3. Provides superior stability at high Reynolds numbers compared to standard BGK.
//!
//! **Reference**: Latt & Chopard (2006), *Phys. Rev. E* 74, 026701.

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
        for q in 0..9 {
            let f_neq = f[q] - f_eq[q];
            let lattice_vel = D2Q9::VELOCITIES[q];
            let cx = T::from_i32(lattice_vel.0)
                .expect("lattice velocity is a small integer");
            let cy = T::from_i32(lattice_vel.1)
                .expect("lattice velocity is a small integer");
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
        f: &mut Vec<T>,
        density: &[T],
        velocity: &[T],
        nx: usize,
        ny: usize,
    ) {
        use crate::solvers::lbm::streaming::f_idx;

        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;
                let rho = density[cell];
                let u   = [velocity[cell * 2], velocity[cell * 2 + 1]];

                // Compute equilibrium
                let mut f_eq = [T::zero(); 9];
                let mut f_node = [T::zero(); 9];
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q])
                        .expect("D2Q9 weights are exact f64 constants");
                    let lattice_vel = D2Q9::VELOCITIES[q];
                    f_eq[q]   = equilibrium(rho, &u, q, weight, lattice_vel);
                    f_node[q] = f[f_idx(j, i, q, nx)];
                }

                // Compute non-equilibrium stress Π^(1)
                let pi = Self::compute_stress_tensor(&f_node, &f_eq);

                // Regularized collision: BGK on equilibrium + regularized off-eq
                for q in 0..9 {
                    let f_bgk = f_node[q] - self.omega * (f_node[q] - f_eq[q]);
                    let lattice_vel = D2Q9::VELOCITIES[q];
                    let c = [
                        T::from_i32(lattice_vel.0).expect("lattice velocity is a small integer"),
                        T::from_i32(lattice_vel.1).expect("lattice velocity is a small integer"),
                    ];
                    let reg_term = self.regularization_term(c, pi);
                    f[f_idx(j, i, q, nx)] = f_bgk + reg_term;
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.tau
    }

    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 = T::from_f64(1.0 / 3.0)
            .expect("cs² = 1/3 is representable in IEEE 754");
        let half = T::from_f64(0.5)
            .expect("0.5 is representable in IEEE 754");
        cs2 * dx * dx * (self.tau - half) / dt
    }
}

impl<T: RealField + Copy + FromPrimitive> RegularizedCollision<T> {
    fn regularization_term(&self, c: [T; 2], pi: [[T; 2]; 2]) -> T {
        // Q_i^(1) = (ω/2) · [(c_x² - cs²)·Π_xx + (c_y² - cs²)·Π_yy + c_x c_y (Π_xy + Π_yx)]
        // (Latt & Chopard 2006, eq. 18)
        let cs2 = T::from_f64(1.0 / 3.0)
            .expect("cs² = 1/3 is representable in IEEE 754");
        let half = T::from_f64(0.5)
            .expect("0.5 is representable in IEEE 754");
        let term = (c[0] * c[0] - cs2) * pi[0][0]
            + (c[1] * c[1] - cs2) * pi[1][1]
            + c[0] * c[1] * (pi[0][1] + pi[1][0]);
        term * self.omega * half
    }
}
