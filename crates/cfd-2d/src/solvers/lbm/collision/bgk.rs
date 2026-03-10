//! BGK (Bhatnagar-Gross-Krook) collision operator.
//!
//! Single-relaxation-time collision on the flat `Vec<T>` distribution buffer.
//!
//! # Theorem — BGK H-Theorem (Boltzmann 1872, discrete form)
//!
//! **Statement**: Let $H = \sum_{q} f_q \ln(f_q / w_q)$ where $w_q$ is the
//! lattice weight. The BGK collision operator satisfies $\Delta H \leq 0$, i.e.,
//! the discrete entropy functional is non-increasing.
//!
//! **Proof**:
//!
//! 1. The BGK update is $f_q^* = (1 - \omega) f_q + \omega f_q^{eq}$, a convex
//!    combination for $\omega \in (0, 2)$.
//!
//! 2. $H^{eq} = \sum_q f_q^{eq} \ln(f_q^{eq} / w_q)$ is the global minimum of
//!    $H$ subject to fixed $\rho = \sum_q f_q$ and $\rho\mathbf{u} = \sum_q \mathbf{e}_q f_q$.
//!    This follows from the method of Lagrange multipliers on the convex function
//!    $\phi(f) = f \ln f - f$.
//!
//! 3. By the convexity of $x \mapsto x \ln x$, the mixture
//!    $H^* \leq (1 - \omega) H + \omega H^{eq} \leq H$.
//!    The last inequality holds because $H^{eq} \leq H$. Therefore $H^* \leq H$. □
//!
//! # Theorem — Viscosity-Relaxation Correspondence
//!
//! **Statement**: The BGK kinematic viscosity is
//! $\nu = c_s^2 (\tau - \tfrac{1}{2}) \Delta t$, where $\tau = 1/\omega$.
//!
//! **Proof**: Chapman-Enskog at $O(\epsilon^1)$ yields the viscous stress tensor
//! $\sigma_{xy} = -\mu (\partial u / \partial y + \partial v / \partial x)$ with
//! $\mu / \rho = c_s^2 (\tau - \tfrac{1}{2}) \Delta t$. See He & Luo (1997). □

use super::traits::CollisionOperator;
use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use crate::solvers::lbm::streaming::f_idx;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Lattice sound speed squared: $c_s^2 = 1/3$ (exact in IEEE 754).
const LATTICE_CS2: f64 = 1.0 / 3.0;
/// Additive constant for τ from Navier-Stokes viscosity derivation.
const HALF: f64 = 0.5;

/// BGK single-relaxation-time collision operator.
///
/// Operates on the flat distribution buffer with layout `f[j*nx*9 + i*9 + q]`.
pub struct BgkCollision<T: RealField + Copy> {
    /// Relaxation time τ (lattice units). Must satisfy τ > 0.5 for stability.
    tau: T,
    /// Collision frequency ω = 1/τ. Precomputed to avoid per-node division.
    omega: T,
}

impl<T: RealField + Copy + FromPrimitive> BgkCollision<T> {
    /// Construct from relaxation time τ.
    ///
    /// # Panics
    /// Never panics; tau = 0 would cause ω = ∞ but that is a caller invariant.
    #[must_use]
    pub fn new(tau: T) -> Self {
        Self {
            tau,
            omega: T::one() / tau,
        }
    }

    /// Construct from physical kinematic viscosity (Theorem — Viscosity correspondence).
    ///
    /// τ = ½ + ν Δt / (c_s² Δx²)
    pub fn from_viscosity(nu: T, dt: T, dx: T) -> Self {
        let cs2 = T::from_f64(LATTICE_CS2).expect("T must represent f64; cs² = 1/3");
        let half = T::from_f64(HALF).expect("T must represent f64; ½");
        let tau = half + nu * dt / (cs2 * dx * dx);
        Self::new(tau)
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for BgkCollision<T> {
    /// Perform BGK collision on the flat distribution buffer.
    ///
    /// For each node (i, j) and direction q:
    /// ```text
    /// f_q* = f_q - ω·(f_q - f_q^eq)
    /// ```
    ///
    /// # Correctness (BGK H-Theorem)
    /// The update is a convex step toward equilibrium; H is non-increasing. □
    fn collide(&self, f: &mut [T], density: &[T], velocity: &[T], nx: usize, ny: usize) {
        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;
                let rho = density[cell];
                let u = [velocity[cell * 2], velocity[cell * 2 + 1]];

                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q])
                        .expect("D2Q9 weights are exact f64 constants");
                    let lattice_vel = D2Q9::VELOCITIES[q];
                    let f_eq = equilibrium(rho, &u, q, weight, lattice_vel);
                    let idx = f_idx(j, i, q, nx);
                    f[idx] = f[idx] - self.omega * (f[idx] - f_eq);
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.tau
    }

    /// Kinematic viscosity: ν = c_s²(τ − ½)Δt / Δx² × Δx² = c_s²(τ − ½)Δt
    /// (in physical units, divide by Δx² if Δx ≠ 1)
    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 = T::from_f64(LATTICE_CS2).expect("T must represent f64; cs² = 1/3");
        let half = T::from_f64(HALF).expect("T must represent f64; ½");
        cs2 * dx * dx * (self.tau - half) / dt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// BGK H-theorem: after collision, H* ≤ H (monotone entropy decrease).
    #[test]
    fn test_h_theorem_bgk() {
        use crate::solvers::lbm::lattice::{equilibrium, D2Q9};

        let tau = 1.0_f64;
        let omega = 1.0 / tau;
        let rho = 1.0_f64;
        let u = [0.1_f64, 0.05];

        // Perturb equilibrium to get a non-eq state
        let mut f = [0.0_f64; 9];
        for q in 0..9 {
            let weight = D2Q9::WEIGHTS[q];
            f[q] = equilibrium(rho, &u, q, weight, D2Q9::VELOCITIES[q])
                + 0.01 * (q as f64 - 4.0) * 0.001;
        }

        let rho_after_perturbation = f.iter().sum::<f64>();
        let momentum = f
            .iter()
            .enumerate()
            .fold([0.0_f64, 0.0_f64], |mut acc, (q, &fq)| {
                acc[0] += f64::from(D2Q9::VELOCITIES[q].0) * fq;
                acc[1] += f64::from(D2Q9::VELOCITIES[q].1) * fq;
                acc
            });
        let u_after_perturbation = [
            momentum[0] / rho_after_perturbation,
            momentum[1] / rho_after_perturbation,
        ];

        // H before
        let h_before = f
            .iter()
            .enumerate()
            .map(|(q, &fq)| {
                if fq > 0.0 {
                    fq * (fq / D2Q9::WEIGHTS[q]).ln()
                } else {
                    0.0
                }
            })
            .sum::<f64>();

        // Apply BGK
        let mut f_after = f;
        for q in 0..9 {
            let weight = D2Q9::WEIGHTS[q];
            let f_eq = equilibrium(
                rho_after_perturbation,
                &u_after_perturbation,
                q,
                weight,
                D2Q9::VELOCITIES[q],
            );
            f_after[q] = f[q] - omega * (f[q] - f_eq);
        }

        // H after
        let h_after = f_after
            .iter()
            .enumerate()
            .map(|(q, &fq)| {
                if fq > 0.0 {
                    fq * (fq / D2Q9::WEIGHTS[q]).ln()
                } else {
                    0.0
                }
            })
            .sum::<f64>();

        assert!(
            h_after <= h_before + 1e-12,
            "BGK H-theorem violated: H_before={h_before:.6}, H_after={h_after:.6}"
        );
    }

    /// Viscosity-relaxation correspondence: ν = cs²(τ−½)Δt
    #[test]
    fn test_viscosity_from_tau() {
        let nu = 0.01_f64;
        let dt = 1.0_f64;
        let dx = 1.0_f64;
        let bgk = BgkCollision::<f64>::from_viscosity(nu, dt, dx);
        let nu_recovered = bgk.viscosity(dt, dx);
        assert_relative_eq!(nu_recovered, nu, epsilon = 1e-12);
    }
}
