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
    /// Create default relaxation matrix.
    ///
    /// The M-matrix row ordering (Lallemand & Luo 2000) maps moment indices as:
    /// 0→ρ, 1→e, 2→ε, 3→j_x, 4→q_x, 5→j_y, 6→q_y, 7→p_xx, 8→p_xy.
    ///
    /// Conserved moments (ρ, j_x, j_y) at indices 0, 3, 5 must have s_k = 0.
    pub fn default_d2q9(tau: T) -> Self {
        let omega = T::one() / tau;

        // Standard D2Q9 relaxation rates optimized for stability and acoustic damping
        // per Lallemand & Luo (2000)
        let s_e = T::from_f64(1.64).unwrap_or(omega);
        let s_eps = T::from_f64(1.54).unwrap_or(omega);
        let s_q = T::from_f64(1.9).unwrap_or(omega);

        Self {
            s: [
                T::zero(), // s0: density ρ (conserved)
                s_e,       // s1: energy e (bulk viscosity control)
                s_eps,     // s2: energy squared ε
                T::zero(), // s3: momentum j_x (conserved)
                s_q,       // s4: energy flux q_x
                T::zero(), // s5: momentum j_y (conserved)
                s_q,       // s6: energy flux q_y
                omega,     // s7: stress p_xx − p_yy (kinematic viscosity)
                omega,     // s8: stress p_xy (kinematic viscosity)
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

        // Exact M⁻¹ derived from the orthogonality of M.
        //
        // The Lallemand-Luo M matrix is row-orthogonal: M·Mᵀ = diag(‖r_k‖²)
        // with ‖r_k‖² = [9, 36, 36, 6, 12, 6, 12, 4, 4] for k = 0..8.
        //
        // Therefore M⁻¹ = Mᵀ · diag(1/‖r_k‖²), i.e. M⁻¹[q][k] = M[k][q] / ‖r_k‖².
        //
        // Pre-computed rows (one per velocity direction q = 0..8):
        #[rustfmt::skip]
        const M_INV_ROWS: [[f64; 9]; 9] = [
            // q=0: (0,0) rest
            [ 1.0/9.0, -4.0/36.0,  4.0/36.0,  0.0/6.0,  0.0/12.0,  0.0/6.0,  0.0/12.0,  0.0/4.0,  0.0/4.0],
            // q=1: (1,0)
            [ 1.0/9.0, -1.0/36.0, -2.0/36.0,  1.0/6.0, -2.0/12.0,  0.0/6.0,  0.0/12.0,  1.0/4.0,  0.0/4.0],
            // q=2: (0,1)
            [ 1.0/9.0, -1.0/36.0, -2.0/36.0,  0.0/6.0,  0.0/12.0,  1.0/6.0, -2.0/12.0, -1.0/4.0,  0.0/4.0],
            // q=3: (-1,0)
            [ 1.0/9.0, -1.0/36.0, -2.0/36.0, -1.0/6.0,  2.0/12.0,  0.0/6.0,  0.0/12.0,  1.0/4.0,  0.0/4.0],
            // q=4: (0,-1)
            [ 1.0/9.0, -1.0/36.0, -2.0/36.0,  0.0/6.0,  0.0/12.0, -1.0/6.0,  2.0/12.0, -1.0/4.0,  0.0/4.0],
            // q=5: (1,1)
            [ 1.0/9.0,  2.0/36.0,  1.0/36.0,  1.0/6.0,  1.0/12.0,  1.0/6.0,  1.0/12.0,  0.0/4.0,  1.0/4.0],
            // q=6: (-1,1)
            [ 1.0/9.0,  2.0/36.0,  1.0/36.0, -1.0/6.0, -1.0/12.0,  1.0/6.0,  1.0/12.0,  0.0/4.0, -1.0/4.0],
            // q=7: (-1,-1)
            [ 1.0/9.0,  2.0/36.0,  1.0/36.0, -1.0/6.0, -1.0/12.0, -1.0/6.0, -1.0/12.0,  0.0/4.0,  1.0/4.0],
            // q=8: (1,-1)
            [ 1.0/9.0,  2.0/36.0,  1.0/36.0,  1.0/6.0,  1.0/12.0, -1.0/6.0, -1.0/12.0,  0.0/4.0, -1.0/4.0],
        ];

        for q in 0..9 {
            for k in 0..9 {
                m_inv[q][k] = T::from_f64(M_INV_ROWS[q][k])
                    .expect("MRT inverse matrix entries are exact f64 constants");
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Viscosity-relaxation round-trip: ν → τ → ν.
    /// MRT uses the same viscosity formula as BGK (stress moment relaxation).
    #[test]
    fn viscosity_round_trip() {
        let nu = 0.01_f64;
        let dt = 1.0;
        let dx = 1.0;
        let cs2 = 1.0 / 3.0;
        let tau = 0.5 + nu * dt / (cs2 * dx * dx);
        let mrt = MrtCollision::<f64>::new(tau);
        let nu_recovered = mrt.viscosity(dt, dx);
        assert_relative_eq!(nu_recovered, nu, epsilon = 1e-12);
    }

    /// Equilibrium moments: the equilibrium_moments function must produce
    /// moment 0 = ρ, moment 3 = ρ·u_x, moment 5 = ρ·u_y (conserved).
    #[test]
    fn equilibrium_moments_conserved_quantities() {
        let rho = 1.05_f64;
        let u = [0.08_f64, -0.03];
        let m_eq = MrtCollision::<f64>::equilibrium_moments(rho, u);

        // Moment 0 = density
        assert_relative_eq!(m_eq[0], rho, epsilon = 1e-14);
        // Moment 3 = j_x = ρ u_x
        assert_relative_eq!(m_eq[3], rho * u[0], epsilon = 1e-14);
        // Moment 5 = j_y = ρ u_y
        assert_relative_eq!(m_eq[5], rho * u[1], epsilon = 1e-14);
    }

    /// Relaxation matrix: conserved moments (indices 0, 3, 5) must have
    /// relaxation rate s_k = 0. Non-conserved moments must have s_k ∈ (0, 2).
    #[test]
    fn relaxation_rates_conservation_and_stability() {
        let tau = 0.8_f64;
        let s = RelaxationMatrix::<f64>::default_d2q9(tau);

        // Conserved moments: s = 0
        assert_relative_eq!(s.s[0], 0.0, epsilon = 1e-14); // density
        assert_relative_eq!(s.s[3], 0.0, epsilon = 1e-14); // j_x
        assert_relative_eq!(s.s[5], 0.0, epsilon = 1e-14); // j_y

        // Non-conserved moments: s ∈ (0, 2) (Theorem — MRT Linear Stability)
        for k in [1, 2, 4, 6, 7, 8] {
            assert!(
                s.s[k] > 0.0 && s.s[k] < 2.0,
                "s[{k}] = {} not in (0, 2)",
                s.s[k]
            );
        }
    }

    /// Transformation matrix M: row 0 sums to 9 (density moment = Σf_q),
    /// row 3 is the x-momentum (dot with e_x), row 5 is y-momentum.
    #[test]
    fn transformation_matrix_structure() {
        let mrt = MrtCollision::<f64>::new(1.0);

        // Row 0 (density): all ones → sum = 9
        let row0_sum: f64 = mrt.m[0].iter().sum();
        assert_relative_eq!(row0_sum, 9.0, epsilon = 1e-14);

        // Row 3 (j_x): should match e_x velocities [0,1,0,-1,0,1,-1,-1,1]
        let expected_row3 = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
        for q in 0..9 {
            assert_relative_eq!(mrt.m[3][q], expected_row3[q], epsilon = 1e-14);
        }

        // Row 5 (j_y): should match e_y velocities [0,0,1,0,-1,1,1,-1,-1]
        let expected_row5 = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
        for q in 0..9 {
            assert_relative_eq!(mrt.m[5][q], expected_row5[q], epsilon = 1e-14);
        }
    }

    /// M·M⁻¹ = I₉ₓ₉: the product of the transformation matrix and its inverse
    /// must be the identity matrix.
    #[test]
    fn m_times_m_inv_is_identity() {
        let mrt = MrtCollision::<f64>::new(1.0);

        for i in 0..9 {
            for j in 0..9 {
                let mut dot = 0.0_f64;
                for k in 0..9 {
                    dot += mrt.m[i][k] * mrt.m_inv[k][j];
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                assert_relative_eq!(dot, expected, epsilon = 1e-12,);
            }
        }
    }

    /// Conserved-moment invariance: density ρ and momenta j_x, j_y must be
    /// unchanged by MRT collision (s₀ = s₃ = s₅ = 0 in default_d2q9).
    #[test]
    fn conserved_moments_preserved() {
        use crate::solvers::lbm::lattice::{equilibrium, D2Q9};

        let nx = 1_usize;
        let ny = 1_usize;
        let tau = 0.8_f64;
        let mrt = MrtCollision::<f64>::new(tau);

        let rho = 1.05_f64;
        let u = [0.08_f64, -0.03];

        // Build equilibrium + perturbation
        let mut f = vec![0.0_f64; 9];
        for q in 0..9 {
            let w = D2Q9::WEIGHTS[q];
            f[q] = equilibrium(rho, &u, q, w, D2Q9::VELOCITIES[q]) + 0.001 * (q as f64 - 4.0);
        }

        // Recompute exact density and momentum from perturbed f
        let rho_before: f64 = f.iter().sum();
        let jx_before: f64 = f
            .iter()
            .enumerate()
            .map(|(q, &fq)| f64::from(D2Q9::VELOCITIES[q].0) * fq)
            .sum();
        let jy_before: f64 = f
            .iter()
            .enumerate()
            .map(|(q, &fq)| f64::from(D2Q9::VELOCITIES[q].1) * fq)
            .sum();

        let density = vec![rho_before];
        let velocity = vec![jx_before / rho_before, jy_before / rho_before];

        mrt.collide(&mut f, &density, &velocity, nx, ny);

        let rho_after: f64 = f.iter().sum();
        let jx_after: f64 = f
            .iter()
            .enumerate()
            .map(|(q, &fq)| f64::from(D2Q9::VELOCITIES[q].0) * fq)
            .sum();
        let jy_after: f64 = f
            .iter()
            .enumerate()
            .map(|(q, &fq)| f64::from(D2Q9::VELOCITIES[q].1) * fq)
            .sum();

        assert_relative_eq!(rho_after, rho_before, epsilon = 1e-12);
        assert_relative_eq!(jx_after, jx_before, epsilon = 1e-12);
        assert_relative_eq!(jy_after, jy_before, epsilon = 1e-12);
    }

    /// Equilibrium fixed point: applying MRT to f = f^eq must leave f unchanged.
    #[test]
    fn equilibrium_is_fixed_point() {
        use crate::solvers::lbm::lattice::{equilibrium, D2Q9};

        let nx = 1_usize;
        let ny = 1_usize;
        let tau = 1.0_f64;
        let mrt = MrtCollision::<f64>::new(tau);

        let rho = 1.0_f64;
        let u = [0.05_f64, 0.02];

        let mut f = vec![0.0_f64; 9];
        for q in 0..9 {
            let w = D2Q9::WEIGHTS[q];
            f[q] = equilibrium(rho, &u, q, w, D2Q9::VELOCITIES[q]);
        }
        let f_orig = f.clone();

        let density = vec![rho];
        let velocity = vec![u[0], u[1]];
        mrt.collide(&mut f, &density, &velocity, nx, ny);

        for q in 0..9 {
            assert_relative_eq!(f[q], f_orig[q], epsilon = 1e-12);
        }
    }
}
