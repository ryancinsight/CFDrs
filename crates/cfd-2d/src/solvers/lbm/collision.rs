//! Collision operators for LBM.
//!
//! This module implements various collision models including
//! BGK (Bhatnagar-Gross-Krook) single relaxation time.

use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Trait for collision operators
pub trait CollisionOperator<T: RealField + Copy> {
    /// Apply collision step to distribution functions
    fn collide(&self, f: &mut Vec<Vec<[T; 9]>>, density: &Vec<Vec<T>>, velocity: &Vec<Vec<[T; 2]>>);

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

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for MrtCollision<T> {
    fn collide(
        &self,
        f: &mut Vec<Vec<[T; 9]>>,
        density: &Vec<Vec<T>>,
        velocity: &Vec<Vec<[T; 2]>>,
    ) {
        let nx = f.len();
        let ny = if nx > 0 { f[0].len() } else { 0 };

        // D2Q9 moment transformation matrix (orthogonal basis)
        let m = self.get_transformation_matrix();
        let m_inv = self.get_inverse_transformation_matrix();

        for i in 0..nx {
            for j in 0..ny {
                // Transform to moment space: m = M * f
                let mut moments = [T::zero(); 9];
                for k in 0..9 {
                    for l in 0..9 {
                        moments[k] = moments[k] + m[k][l] * f[i][j][l];
                    }
                }

                // Calculate equilibrium moments
                let rho = density[i][j];
                let ux = velocity[i][j][0];
                let uy = velocity[i][j][1];
                let u_sq = ux * ux + uy * uy;

                let mut eq_moments = [T::zero(); 9];
                eq_moments[0] = rho; // Density
                eq_moments[1] = -T::from_f64(2.0).unwrap_or_else(|| T::one()) * rho
                    + T::from_f64(3.0).unwrap_or_else(|| T::one()) * u_sq; // Energy
                eq_moments[2] = rho - T::from_f64(3.0).unwrap_or_else(|| T::one()) * u_sq; // Energy square
                eq_moments[3] = rho * ux; // Momentum x
                eq_moments[4] = -rho * ux; // Energy flux x
                eq_moments[5] = rho * uy; // Momentum y
                eq_moments[6] = -rho * uy; // Energy flux y
                eq_moments[7] = rho * (ux * ux - uy * uy); // Stress xx-yy
                eq_moments[8] = rho * ux * uy; // Stress xy

                // Collision in moment space
                for k in 0..9 {
                    moments[k] =
                        moments[k] - self.relaxation_times[k] * (moments[k] - eq_moments[k]);
                }

                // Transform back to velocity space: f = M^(-1) * m
                for k in 0..9 {
                    f[i][j][k] = T::zero();
                    for l in 0..9 {
                        f[i][j][k] = f[i][j][k] + m_inv[k][l] * moments[l];
                    }
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.relaxation_times[7] // Return shear viscosity relaxation time
    }
}

impl<T: RealField + Copy + FromPrimitive> MrtCollision<T> {
    /// Get the D2Q9 moment transformation matrix
    fn get_transformation_matrix(&self) -> [[T; 9]; 9] {
        let one = T::one();
        let two = T::from_f64(2.0).unwrap_or_else(|| one + one);
        let three = T::from_f64(3.0).unwrap_or_else(|| one + one + one);
        let four = T::from_f64(4.0).unwrap_or_else(|| two * two);
        let six = T::from_f64(6.0).unwrap_or_else(|| three * two);

        // Orthogonal moment basis for D2Q9
        [
            [one, one, one, one, one, one, one, one, one], // Density
            [-four, -one, -one, -one, -one, two, two, two, two], // Energy
            [four, -two, -two, -two, -two, one, one, one, one], // Energy square
            [
                T::zero(),
                one,
                T::zero(),
                -one,
                T::zero(),
                one,
                -one,
                -one,
                one,
            ], // Momentum x
            [
                T::zero(),
                -two,
                T::zero(),
                two,
                T::zero(),
                one,
                -one,
                -one,
                one,
            ], // Energy flux x
            [
                T::zero(),
                T::zero(),
                one,
                T::zero(),
                -one,
                one,
                one,
                -one,
                -one,
            ], // Momentum y
            [
                T::zero(),
                T::zero(),
                -two,
                T::zero(),
                two,
                one,
                one,
                -one,
                -one,
            ], // Energy flux y
            [
                T::zero(),
                one,
                -one,
                one,
                -one,
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
            ], // Stress xx-yy
            [
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
                one,
                -one,
                one,
                -one,
            ], // Stress xy
        ]
    }

    /// Get the inverse of the D2Q9 moment transformation matrix
    fn get_inverse_transformation_matrix(&self) -> [[T; 9]; 9] {
        // Pre-computed inverse for the orthogonal D2Q9 moment basis
        let inv9 = T::one() / T::from_f64(9.0).unwrap_or_else(|| T::one());
        let inv36 = T::one() / T::from_f64(36.0).unwrap_or_else(|| T::one());
        let inv18 = T::one() / T::from_f64(18.0).unwrap_or_else(|| T::one());
        let inv12 = T::one() / T::from_f64(12.0).unwrap_or_else(|| T::one());
        let inv6 = T::one() / T::from_f64(6.0).unwrap_or_else(|| T::one());
        let inv4 = T::one() / T::from_f64(4.0).unwrap_or_else(|| T::one());

        [
            [
                inv9,
                -inv9,
                inv9,
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
                T::zero(),
            ],
            [
                inv9,
                -inv36,
                -inv18,
                inv6,
                -inv6,
                T::zero(),
                T::zero(),
                inv4,
                T::zero(),
            ],
            [
                inv9,
                -inv36,
                -inv18,
                T::zero(),
                T::zero(),
                inv6,
                -inv6,
                -inv4,
                T::zero(),
            ],
            [
                inv9,
                -inv36,
                -inv18,
                -inv6,
                inv6,
                T::zero(),
                T::zero(),
                inv4,
                T::zero(),
            ],
            [
                inv9,
                -inv36,
                -inv18,
                T::zero(),
                T::zero(),
                -inv6,
                inv6,
                -inv4,
                T::zero(),
            ],
            [
                inv9,
                inv18,
                inv36,
                inv12,
                inv12,
                inv12,
                inv12,
                T::zero(),
                inv4,
            ],
            [
                inv9,
                inv18,
                inv36,
                -inv12,
                -inv12,
                inv12,
                inv12,
                T::zero(),
                -inv4,
            ],
            [
                inv9,
                inv18,
                inv36,
                -inv12,
                -inv12,
                -inv12,
                -inv12,
                T::zero(),
                inv4,
            ],
            [
                inv9,
                inv18,
                inv36,
                inv12,
                inv12,
                -inv12,
                -inv12,
                T::zero(),
                -inv4,
            ],
        ]
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
