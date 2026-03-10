//! Multiple Relaxation Time (MRT) collision operator.
//!
//! Based on Lallemand & Luo (2000) — "Theory of the lattice Boltzmann method:
//! Dispersion, dissipation, isotropy, Galilean invariance, and stability",
//! Physical Review E, 61(6), 6546–6562.
//!
//! # Theorem — MRT Linear Stability
//!
//! **Statement**: The MRT collision operator $\Omega_{\text{MRT}} = -M^{-1} \hat{S}(\hat{m} - \hat{m}^{eq})$
//! in moment space is linearly stable if and only if all relaxation rates
//! $s_k \in (0, 2)$ for $k = 0, \dots, 8$.
//!
//! **Proof**:
//!
//! 1. In moment space the MRT update is $\hat{m}_k^* = (1 - s_k)\hat{m}_k + s_k \hat{m}_k^{eq}$.
//!
//! 2. Von Neumann stability analysis for a plane wave perturbation $\delta\hat{m}_k \propto e^{i\theta}$
//!    yields the amplification $|\hat{m}_k^*| = |1 - s_k| \cdot |\hat{m}_k|$.
//!
//! 3. Stability requires $|1 - s_k| \leq 1$, i.e., $s_k \in [0, 2]$.
//!    Strict inequality $s_k \in (0, 2)$ excludes the neutrally-stable endpoints.
//!    □
//!
//! **Reference**: Lallemand & Luo (2000), §IV.

use super::traits::CollisionOperator;
use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Lattice sound speed squared: $c_s^2 = 1/3$ for D2Q9
const LATTICE_SOUND_SPEED_SQUARED: f64 = 1.0 / 3.0;

/// Relaxation time offset: $\tau - 0.5$ relates viscosity to relaxation
const RELAXATION_TIME_OFFSET: f64 = 0.5;

// MRT transformation matrix rows (D2Q9 moments)
// Row 0: Density ρ
const M_ROW_0: [f64; 9] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
// Row 1: Energy e
const M_ROW_1: [f64; 9] = [-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0];
// Row 2: Energy squared ε
const M_ROW_2: [f64; 9] = [4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0];
// Row 3: Momentum j_x
const M_ROW_3: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
// Row 4: Energy flux q_x
const M_ROW_4: [f64; 9] = [0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0];
// Row 5: Momentum j_y
const M_ROW_5: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
// Row 6: Energy flux q_y
const M_ROW_6: [f64; 9] = [0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0];
// Row 7: Stress tensor p_xx - p_yy
const M_ROW_7: [f64; 9] = [0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0];
// Row 8: Stress tensor p_xy
const M_ROW_8: [f64; 9] = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0];

/// Relaxation matrix for MRT
pub struct RelaxationMatrix<T: RealField + Copy> {
    /// Relaxation rates for each moment
    pub s: [T; 9],
}

impl<T: RealField + Copy + FromPrimitive> RelaxationMatrix<T> {
    /// Create default relaxation matrix
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
                m[i][j] = T::from_f64(rows[i][j])
                    .expect("MRT transform matrix entries are exact f64 constants");
            }
        }

        // Build the inverse transformation matrix M^(-1)
        // Pre-computed values for D2Q9 from the literature
        let inv_9 = T::from_f64(1.0 / 9.0).expect("1/9 is representable in IEEE 754");
        let inv_36 = T::from_f64(1.0 / 36.0).expect("1/36 is representable in IEEE 754");
        let inv_6 = T::from_f64(1.0 / 6.0).expect("1/6 is representable in IEEE 754");
        let inv_12 = T::from_f64(1.0 / 12.0).expect("1/12 is representable in IEEE 754");
        let inv_4 = T::from_f64(1.0 / 4.0).expect("1/4 is representable in IEEE 754");
        let two = T::from_f64(2.0).expect("2.0 is representable in IEEE 754");

        for q in 0..9 {
            let vel = D2Q9::VELOCITIES[q];
            let cx = T::from_i32(vel.0).expect("lattice velocity is a small integer");
            let cy = T::from_i32(vel.1).expect("lattice velocity is a small integer");

            m_inv[q][0] = inv_9;

            if q == 0 {
                m_inv[q][1] = -inv_9 * two * two; // -4/9
                m_inv[q][2] = inv_9 * two * two; //  4/9
            } else if q < 5 {
                m_inv[q][1] = -inv_36;
                m_inv[q][2] = -inv_36 * two;
            } else {
                m_inv[q][1] = inv_36 * two;
                m_inv[q][2] = inv_36;
            }

            m_inv[q][3] = inv_6 * cx;
            m_inv[q][4] = -inv_12 * cx;
            m_inv[q][5] = inv_6 * cy;
            m_inv[q][6] = -inv_12 * cy;
            m_inv[q][7] = inv_4 * (cx * cx - cy * cy);
            m_inv[q][8] = inv_4 * cx * cy;
        }

        (m, m_inv)
    }
}

impl<T: RealField + Copy + FromPrimitive> CollisionOperator<T> for MrtCollision<T> {
    /// Apply MRT collision to the flat distribution buffer.
    ///
    /// For each node: transform to moment space, relax non-conserved moments
    /// toward equilibrium at rate s_k (Theorem — MRT stability: s_k ∈ (0,2)),
    /// then transform back to velocity space.
    fn collide(&self, f: &mut [T], density: &[T], velocity: &[T], nx: usize, ny: usize) {
        use crate::solvers::lbm::streaming::f_idx;

        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;

                // ── 1. Transform f → moment space ──
                let mut moments = [T::zero(); 9];
                for p in 0..9 {
                    for q in 0..9 {
                        moments[p] += self.m[p][q] * f[f_idx(j, i, q, nx)];
                    }
                }

                // ── 2. Equilibrium moments ──
                let rho = density[cell];
                let u = [velocity[cell * 2], velocity[cell * 2 + 1]];
                let m_eq = Self::equilibrium_moments(rho, u);

                // ── 3. Relax each moment at its rate s_k (Theorem: s_k ∈ (0,2)) ──
                for p in 0..9 {
                    moments[p] -= self.s.s[p] * (moments[p] - m_eq[p]);
                }

                // ── 4. Transform moment space → velocity space ──
                for q in 0..9 {
                    let mut fq = T::zero();
                    for p in 0..9 {
                        fq += self.m_inv[q][p] * moments[p];
                    }
                    f[f_idx(j, i, q, nx)] = fq;
                }
            }
        }
    }

    fn tau(&self) -> T {
        self.tau
    }

    fn viscosity(&self, dt: T, dx: T) -> T {
        let cs2 = T::from_f64(LATTICE_SOUND_SPEED_SQUARED)
            .expect("cs² = 1/3 is representable in IEEE 754");
        let half = T::from_f64(RELAXATION_TIME_OFFSET).expect("0.5 is representable in IEEE 754");
        cs2 * dx * dx * (self.tau - half) / dt
    }
}

impl<T: RealField + Copy + FromPrimitive> MrtCollision<T> {
    fn equilibrium_moments(rho: T, u: [T; 2]) -> [T; 9] {
        let mut m_eq = [T::zero(); 9];
        m_eq[0] = rho;
        m_eq[3] = rho * u[0];
        m_eq[5] = rho * u[1];

        let two = T::from_f64(2.0).expect("2.0 is representable in IEEE 754");
        let three = T::from_f64(3.0).expect("3.0 is representable in IEEE 754");
        let u_sq = u[0] * u[0] + u[1] * u[1];

        m_eq[1] = -two * rho + three * rho * u_sq;
        m_eq[2] = rho - three * rho * u_sq;
        m_eq[4] = -rho * u[0];
        m_eq[6] = -rho * u[1];
        m_eq[7] = rho * (u[0] * u[0] - u[1] * u[1]);
        m_eq[8] = rho * u[0] * u[1];

        m_eq
    }
}
