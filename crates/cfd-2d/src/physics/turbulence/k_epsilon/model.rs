//! Core k-ε turbulence model structure and solver.
//!
//! This module contains the `KEpsilonModel` struct — the primary data holder
//! and time-stepper for both the standard and realizable k-ε variants.

use super::kato_launder;
use super::realizable;
use crate::physics::turbulence::constants::{
    C1_EPSILON, C2_EPSILON, C_MU, EPSILON_MIN, K_MIN, SIGMA_EPSILON, SIGMA_K,
};
use crate::physics::turbulence::traits::TurbulenceModel;
use cfd_core::{
    error::Result,
    physics::constants::mathematical::numeric::{ONE_HALF, TWO},
};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// k-ε turbulence model.
///
/// Supports both the standard k-ε formulation (Launder & Spalding 1974) and the
/// Realizable k-ε variant (Shih, Zhu & Lumley 1995).
///
/// When `use_realizable` is `true`, the fixed constant C_mu = 0.09 is replaced by
/// a strain-rate-dependent formulation that prevents unphysical over-prediction of
/// turbulent viscosity in stagnation regions:
///
/// ```text
/// C_mu = 1 / (A_0 + A_s · S̃ · k / ε)
/// ```
///
/// # Reference
///
/// Shih, T.-H., Zhu, J., & Lumley, J. L. (1995). A New Reynolds Stress Algebraic
/// Equation Model. *Computers & Fluids*, 24(3), 227–238.
pub struct KEpsilonModel<T: RealField + Copy + num_traits::ToPrimitive> {
    /// Grid X dimension.
    pub(crate) nx: usize,
    /// Grid Y dimension.
    pub(crate) ny: usize,
    /// Model coefficient C_μ.
    pub(crate) c_mu: T,
    /// Model coefficient C_ε1.
    pub(crate) c1_epsilon: T,
    /// Model coefficient C_ε2.
    pub(crate) c2_epsilon: T,
    /// Turbulent Prandtl number for k.
    pub(crate) sigma_k: T,
    /// Turbulent Prandtl number for ε.
    pub(crate) sigma_epsilon: T,
    /// When `true`, use the Realizable k-ε formulation (Shih et al. 1995)
    /// with a strain-rate-dependent C_mu instead of the fixed constant.
    pub(crate) use_realizable: bool,
    /// When `true`, use the Kato-Launder (1993) vorticity-strain production
    /// term P_k = ν_t · S · Ω instead of the standard P_k = ν_t · S².
    pub(crate) use_kato_launder: bool,
    /// Scratch buffer for k values at previous timestep (avoids per-update allocation).
    pub(crate) k_scratch: Vec<T>,
    /// Scratch buffer for epsilon values at previous timestep.
    pub(crate) eps_scratch: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive> KEpsilonModel<T> {
    /// Create a new standard k-ε model (C_mu = 0.09 fixed).
    pub fn new(nx: usize, ny: usize) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            c_mu: T::from_f64(C_MU).expect("analytical constant conversion"),
            c1_epsilon: T::from_f64(C1_EPSILON).expect("analytical constant conversion"),
            c2_epsilon: T::from_f64(C2_EPSILON).expect("analytical constant conversion"),
            sigma_k: T::from_f64(SIGMA_K).expect("analytical constant conversion"),
            sigma_epsilon: T::from_f64(SIGMA_EPSILON).expect("analytical constant conversion"),
            use_realizable: false,
            use_kato_launder: false,
            k_scratch: vec![T::zero(); n],
            eps_scratch: vec![T::zero(); n],
        }
    }

    /// Create a new Realizable k-ε model (Shih et al. 1995).
    ///
    /// The realizable variant computes C_mu from the local strain rate,
    /// bounding it above by `1/A_0 ≈ 0.247` and reducing it in regions of
    /// strong strain to prevent unphysical turbulent viscosity.
    pub fn new_realizable(nx: usize, ny: usize) -> Self {
        let mut model = Self::new(nx, ny);
        model.use_realizable = true;
        model
    }

    /// Enable or disable the Realizable k-ε formulation.
    pub fn set_realizable(&mut self, enabled: bool) {
        self.use_realizable = enabled;
    }

    /// Returns whether the Realizable formulation is active.
    pub fn is_realizable(&self) -> bool {
        self.use_realizable
    }

    /// Enable or disable the Kato-Launder (1993) vorticity-strain production.
    pub fn set_kato_launder(&mut self, enabled: bool) {
        self.use_kato_launder = enabled;
    }

    /// Returns whether the Kato-Launder formulation is active.
    pub fn is_kato_launder(&self) -> bool {
        self.use_kato_launder
    }

    /// Calculate strain rate tensor magnitude |S| = √(2 S_ij S_ij).
    ///
    /// ## Theorem — Strain Rate Magnitude
    ///
    /// For the symmetric part S_ij = ½(∂u_i/∂x_j + ∂u_j/∂x_i), the
    /// Frobenius norm |S| = √(2 S_ij S_ij) is always non-negative and
    /// vanishes if and only if the flow is in rigid-body rotation.
    pub(crate) fn strain_rate(&self, velocity_gradient: &[[T; 2]; 2]) -> T {
        let mut s_ij = [[T::zero(); 2]; 2];

        // S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
        for i in 0..2 {
            for j in 0..2 {
                s_ij[i][j] = (velocity_gradient[i][j] + velocity_gradient[j][i])
                    * T::from_f64(ONE_HALF).expect("analytical constant conversion");
            }
        }

        // Calculate magnitude: sqrt(2 * S_ij * S_ij)
        let mut s_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                s_squared += s_ij[i][j] * s_ij[i][j];
            }
        }

        (T::from_f64(TWO).expect("analytical constant conversion") * s_squared).sqrt()
    }

    /// Apply boundary conditions using the boundary condition system.
    pub(crate) fn apply_boundary_conditions(&self, k: &mut [T], epsilon: &mut [T]) {
        use crate::physics::turbulence::boundary_conditions::{
            TurbulenceBoundaryCondition, TurbulenceBoundaryManager,
        };
        use crate::physics::turbulence::wall_functions::{WallFunction, WallTreatment};

        let manager = TurbulenceBoundaryManager::new(self.nx, self.ny, T::one(), T::one());

        let boundaries = vec![
            (
                "west".to_string(),
                TurbulenceBoundaryCondition::Wall {
                    wall_treatment: WallTreatment::new(WallFunction::Standard),
                },
            ),
            (
                "east".to_string(),
                TurbulenceBoundaryCondition::Wall {
                    wall_treatment: WallTreatment::new(WallFunction::Standard),
                },
            ),
            (
                "south".to_string(),
                TurbulenceBoundaryCondition::Wall {
                    wall_treatment: WallTreatment::new(WallFunction::Standard),
                },
            ),
            (
                "north".to_string(),
                TurbulenceBoundaryCondition::Wall {
                    wall_treatment: WallTreatment::new(WallFunction::Standard),
                },
            ),
        ];

        manager.apply_k_epsilon_boundaries(k, epsilon, &boundaries);
    }
}

impl<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive> TurbulenceModel<T>
    for KEpsilonModel<T>
{
    fn turbulent_viscosity(&self, k: T, epsilon: T, density: T) -> T {
        let eps_min = T::from_f64(EPSILON_MIN).expect("analytical constant conversion");
        density * self.c_mu * k * k / epsilon.max(eps_min)
    }

    fn production_term(
        &self,
        velocity_gradient: &[[T; 2]; 2],
        turbulent_viscosity: T,
        _turbulence_variable: T,
        _wall_distance: T,
        _molecular_viscosity: T,
    ) -> T {
        // P_k = ν_t · |S|² where strain_rate() returns |S| = sqrt(2·S_ij·S_ij)
        // Reference: Launder & Spalding (1974), Eq. (2.5)
        let strain = self.strain_rate(velocity_gradient);
        turbulent_viscosity * strain * strain
    }

    fn dissipation_term(&self, _k: T, epsilon: T) -> T {
        epsilon
    }

    fn update(
        &mut self,
        k: &mut [T],
        epsilon: &mut [T],
        velocity: &[Vector2<T>],
        density: T,
        molecular_viscosity: T,
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;

        // Copy previous timestep values into persistent scratch buffers
        self.k_scratch.copy_from_slice(k);
        self.eps_scratch.copy_from_slice(epsilon);

        // Hoist T::from_f64 conversions out of the inner loop
        let two = T::from_f64(TWO).expect("analytical constant conversion");
        let two_f = T::from_f64(2.0).expect("analytical constant conversion");
        let k_min = T::from_f64(K_MIN).expect("analytical constant conversion");
        let eps_min = T::from_f64(EPSILON_MIN).expect("analytical constant conversion");

        let two_dx = two * dx;
        let two_dy = two * dy;
        let dx_sq = dx * dx;
        let dy_sq = dy * dy;

        // Update interior points
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                let k_prev = self.k_scratch[idx];
                let eps_prev = self.eps_scratch[idx];

                // Calculate velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x) / two_dx;
                let du_dy = (velocity[idx + nx].x - velocity[idx - nx].x) / two_dy;
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y) / two_dx;
                let dv_dy = (velocity[idx + nx].y - velocity[idx - nx].y) / two_dy;

                let grad = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate turbulent viscosity
                let nu_t = if self.use_realizable {
                    let c_mu_local = realizable::realizable_c_mu(self, &grad, k_prev, eps_prev);
                    density * c_mu_local * k_prev * k_prev / eps_prev.max(eps_min)
                } else {
                    self.turbulent_viscosity(k_prev, eps_prev, density)
                };

                // Production term
                let p_k = if self.use_kato_launder {
                    let grad_f64 = [
                        [
                            grad[0][0].to_f64().unwrap_or(0.0),
                            grad[0][1].to_f64().unwrap_or(0.0),
                        ],
                        [
                            grad[1][0].to_f64().unwrap_or(0.0),
                            grad[1][1].to_f64().unwrap_or(0.0),
                        ],
                    ];
                    let nu_t_f64 = nu_t.to_f64().unwrap_or(0.0);
                    T::from_f64(kato_launder::kato_launder_production(&grad_f64, nu_t_f64))
                        .expect("analytical constant conversion")
                } else {
                    self.production_term(&grad, nu_t, k_prev, T::zero(), molecular_viscosity)
                };

                // Diffusion terms
                let nu_eff_k = molecular_viscosity + nu_t / self.sigma_k;
                let nu_eff_eps = molecular_viscosity + nu_t / self.sigma_epsilon;

                // k equation diffusion
                let diff_k_x =
                    (self.k_scratch[idx + 1] - two_f * k_prev + self.k_scratch[idx - 1]) / dx_sq;
                let diff_k_y =
                    (self.k_scratch[idx + nx] - two_f * k_prev + self.k_scratch[idx - nx]) / dy_sq;
                let diff_k = nu_eff_k * (diff_k_x + diff_k_y);

                // epsilon equation diffusion
                let diff_eps_x = (self.eps_scratch[idx + 1] - two_f * eps_prev
                    + self.eps_scratch[idx - 1])
                    / dx_sq;
                let diff_eps_y = (self.eps_scratch[idx + nx] - two_f * eps_prev
                    + self.eps_scratch[idx - nx])
                    / dy_sq;
                let diff_eps = nu_eff_eps * (diff_eps_x + diff_eps_y);

                // Update k with realizability constraints
                let k_new = k_prev + dt * (p_k - eps_prev + diff_k);
                k[idx] = k_new.max(k_min);

                // Update epsilon with realizability constraints
                let k_denom = k_prev.max(eps_min);
                let eps_source = self.c1_epsilon * eps_prev / k_denom * p_k;
                let eps_sink = self.c2_epsilon * eps_prev * eps_prev / k_denom;
                let eps_new = eps_prev + dt * (eps_source - eps_sink + diff_eps);
                epsilon[idx] = eps_new.max(eps_min);
            }
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(k, epsilon);

        Ok(())
    }

    fn name(&self) -> &'static str {
        "k-epsilon"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        reynolds > T::from_f64(1e4).expect("analytical constant conversion")
    }
}
