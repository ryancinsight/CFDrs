//! BGK (Bhatnagar-Gross-Krook) collision operator

use super::traits::CollisionOperator;
use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use nalgebra::RealField;
use num_traits::FromPrimitive;

// Named constants for physics
const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;
const HALF: f64 = 0.5;

/// BGK collision operator with single relaxation time
pub struct BgkCollision<T: RealField + Copy> {
    /// Relaxation time
    tau: T,
    /// Inverse relaxation time (omega = 1/tau)
    omega: T,
}

impl<T: RealField + Copy + FromPrimitive> BgkCollision<T> {
    /// Create new BGK collision operator
    #[must_use]
    pub fn new(tau: T) -> Self {
        Self {
            tau,
            omega: T::one() / tau,
        }
    }

    /// Create from kinematic viscosity
    pub fn from_viscosity(nu: T, dt: T, dx: T) -> Self {
        let cs2 = T::from_f64(LATTICE_SOUND_SPEED_SQUARED)
            .unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));
        let half = T::from_f64(HALF).unwrap_or_else(|| T::one() / (T::one() + T::one()));
        let tau = half + nu * dt / (cs2 * dx * dx);
        Self::new(tau)
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for BgkCollision<T> {
    fn collide(&self, f: &mut Vec<Vec<[T; 9]>>, density: &[Vec<T>], velocity: &[Vec<[T; 2]>]) {
        let ny = f.len();
        let nx = if ny > 0 { f[0].len() } else { 0 };

        for j in 0..ny {
            for i in 0..nx {
                let rho = density[j][i];
                let u = velocity[j][i];

                // Compute equilibrium distributions
                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
                    let lattice_vel = D2Q9::VELOCITIES[q];
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

    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 = T::from_f64(LATTICE_SOUND_SPEED_SQUARED)
            .unwrap_or_else(|| T::one() / (T::one() + T::one() + T::one()));
        let half = T::from_f64(HALF).unwrap_or_else(|| T::one() / (T::one() + T::one()));
        cs2 * dx * dx * (self.tau - half) / dt
    }
}
