//! Multiple Relaxation Time (MRT) collision operator
//!
//! Based on Lallemand & Luo (2000) - "Theory of the lattice Boltzmann method:
//! Dispersion, dissipation, isotropy, Galilean invariance, and stability"
//! Physical Review E, 61(6), 6546-6562.

#![allow(dead_code)]

use super::traits::CollisionOperator;
use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;
use num_traits::FromPrimitive;

// Physical constants from Lallemand & Luo (2000)
#[allow(dead_code)]
const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;
#[allow(dead_code)]
const RELAXATION_TIME_OFFSET: f64 = 0.5;

// MRT transformation matrix rows (D2Q9 moments)
// Row 0: Density ρ
#[allow(dead_code)]
const M_ROW_0: [f64; 9] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
// Row 1: Energy e
#[allow(dead_code)]
const M_ROW_1: [f64; 9] = [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0];
// Row 2: Energy squared ε
#[allow(dead_code)]
const M_ROW_2: [f64; 9] = [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0];
// Row 3: Momentum j_x
#[allow(dead_code)]
const M_ROW_3: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
// Row 4: Energy flux q_x
#[allow(dead_code)]
const M_ROW_4: [f64; 9] = [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0];
// Row 5: Momentum j_y
#[allow(dead_code)]
const M_ROW_5: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
// Row 6: Energy flux q_y
#[allow(dead_code)]
const M_ROW_6: [f64; 9] = [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0];
// Row 7: Stress tensor p_xx - p_yy
#[allow(dead_code)]
const M_ROW_7: [f64; 9] = [0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0];
// Row 8: Stress tensor p_xy
#[allow(dead_code)]
const M_ROW_8: [f64; 9] = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0];

/// Relaxation matrix for MRT
#[allow(dead_code)]
pub struct RelaxationMatrix<T: RealField + Copy> {
    /// Relaxation rates for each moment
    pub s: [T; 9],
}

impl<T: RealField + Copy + FromPrimitive> RelaxationMatrix<T> {
    /// Create default relaxation matrix
    #[allow(dead_code)]
    pub fn default_d2q9(tau: T) -> Self {
        let omega = T::one() / tau;

        // Standard D2Q9 relaxation rates optimized for stability and acoustic damping
        // per Lallemand & Luo (2000)
        let s_e = T::from_f64(1.64).unwrap_or(omega);
        let s_eps = T::from_f64(1.54).unwrap_or(omega);
        let s_q = T::from_f64(1.9).unwrap_or(omega);

        Self {
            s: [
                T::zero(), // s0: density (conserved)
                T::zero(), // s1: momentum x (conserved)
                T::zero(), // s2: momentum y (conserved)
                s_e,       // s3: energy (bulk viscosity control)
                s_eps,     // s4: energy squared
                s_q,       // s5: energy flux x
                s_q,       // s6: energy flux y
                omega,     // s7: stress xx (kinematic viscosity)
                omega,     // s8: stress xy (kinematic viscosity)
            ],
        }
    }
}

/// MRT collision operator
#[allow(dead_code)]
pub struct MrtCollision<T: RealField + Copy> {
    /// Transformation matrix to moment space (9x9)
    m: [[T; 9]; 9],
    /// Inverse transformation matrix (9x9)
    m_inv: [[T; 9]; 9],
    /// Relaxation matrix
    s: RelaxationMatrix<T>,
    /// Primary relaxation time
    tau: T,
}

impl<T: RealField + Copy + FromPrimitive> MrtCollision<T> {
    /// Create new MRT collision operator
    pub fn new(tau: T) -> Self {
        let s = RelaxationMatrix::default_d2q9(tau);
        let (m, m_inv) = Self::create_transform_matrices();

        Self { m, m_inv, s, tau }
    }

    fn create_transform_matrices() -> ([[T; 9]; 9], [[T; 9]; 9]) {
        // Proper D2Q9 MRT transformation matrices from Lallemand & Luo (2000)
        let mut m = [[T::zero(); 9]; 9];
        let mut m_inv = [[T::zero(); 9]; 9];

        // Build the transformation matrix M
        let rows = [
            M_ROW_0, M_ROW_1, M_ROW_2, M_ROW_3, M_ROW_4, M_ROW_5, M_ROW_6, M_ROW_7, M_ROW_8,
        ];
        for i in 0..9 {
            for j in 0..9 {
                m[i][j] = T::from_f64(rows[i][j]).unwrap_or_else(T::zero);
            }
        }

        // Build the inverse transformation matrix M^(-1)
        // Pre-computed values for D2Q9 from the literature
        let inv_9 = T::from_f64(1.0 / 9.0).unwrap_or_else(T::zero);
        let inv_36 = T::from_f64(1.0 / 36.0).unwrap_or_else(T::zero);
        let inv_6 = T::from_f64(1.0 / 6.0).unwrap_or_else(T::zero);
        let inv_12 = T::from_f64(1.0 / 12.0).unwrap_or_else(T::zero);
        let inv_4 = T::from_f64(1.0 / 4.0).unwrap_or_else(T::zero);

        for q in 0..9 {
            let vel = D2Q9::VELOCITIES[q];
            let cx = T::from_i32(vel.0).unwrap_or_else(T::zero);
            let cy = T::from_i32(vel.1).unwrap_or_else(T::zero);

            m_inv[q][0] = inv_9; // Density component

            // Energy components
            if q == 0 {
                m_inv[q][1] = -inv_9 * T::from_f64(4.0).unwrap_or_else(T::zero);
                m_inv[q][2] = inv_9 * T::from_f64(4.0).unwrap_or_else(T::zero);
            } else if q < 5 {
                m_inv[q][1] = -inv_36;
                m_inv[q][2] = -inv_36 * T::from_f64(2.0).unwrap_or_else(T::zero);
            } else {
                m_inv[q][1] = inv_36 * T::from_f64(2.0).unwrap_or_else(T::zero);
                m_inv[q][2] = inv_36;
            }

            // Momentum components
            m_inv[q][3] = inv_6 * cx;
            m_inv[q][4] = -inv_12 * cx;
            m_inv[q][5] = inv_6 * cy;
            m_inv[q][6] = -inv_12 * cy;

            // Stress tensor components
            m_inv[q][7] = inv_4 * (cx * cx - cy * cy);
            m_inv[q][8] = inv_4 * cx * cy;
        }

        (m, m_inv)
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for MrtCollision<T> {
    fn collide(&self, f: &mut Vec<Vec<[T; 9]>>, density: &[Vec<T>], velocity: &[Vec<[T; 2]>]) {
        let ny = f.len();
        let nx = if ny > 0 { f[0].len() } else { 0 };

        for j in 0..ny {
            for i in 0..nx {
                // Transform to moment space
                let mut moments = [T::zero(); 9];
                for q in 0..9 {
                    for p in 0..9 {
                        moments[p] += self.m[p][q] * f[j][i][q];
                    }
                }

                // Compute equilibrium moments
                let rho = density[j][i];
                let u = velocity[j][i];
                let m_eq = Self::equilibrium_moments(rho, u);

                // Apply relaxation
                for p in 0..9 {
                    moments[p] = moments[p] - self.s.s[p] * (moments[p] - m_eq[p]);
                }

                // Transform back to velocity space
                for q in 0..9 {
                    f[j][i][q] = T::zero();
                    for p in 0..9 {
                        f[j][i][q] += self.m_inv[q][p] * moments[p];
                    }
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.tau
    }

    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 = T::from_f64(LATTICE_SOUND_SPEED_SQUARED).unwrap_or_else(T::zero);
        let half = T::from_f64(RELAXATION_TIME_OFFSET).unwrap_or_else(T::zero);
        cs2 * dx * dx * (self.tau - half) / dt
    }
}

impl<T: RealField + Copy + FromPrimitive> MrtCollision<T> {
    fn equilibrium_moments(rho: T, u: [T; 2]) -> [T; 9] {
        // Proper equilibrium moments from Lallemand & Luo (2000)
        // These are the equilibrium values for the 9 moments in D2Q9
        let mut m_eq = [T::zero(); 9];
        let u_sq = u[0] * u[0] + u[1] * u[1];
        let two = T::one() + T::one();
        let three = two + T::one();

        // Conserved moments (these don't relax)
        m_eq[0] = rho; // ρ: density
        m_eq[3] = rho * u[0]; // j_x: x-momentum (ρu_x)
        m_eq[5] = rho * u[1]; // j_y: y-momentum (ρu_y)

        // Non-conserved moments (these relax to equilibrium)
        // Energy e
        m_eq[1] = -two * rho + three * rho * u_sq;

        // Energy squared ε
        m_eq[2] = rho - three * rho * u_sq;

        // Energy flux q_x
        m_eq[4] = -rho * u[0];

        // Energy flux q_y
        m_eq[6] = -rho * u[1];

        // Diagonal stress p_xx - p_yy
        m_eq[7] = rho * (u[0] * u[0] - u[1] * u[1]);

        // Off-diagonal stress p_xy
        m_eq[8] = rho * u[0] * u[1];

        m_eq
    }
}
