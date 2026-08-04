use super::{
    super::constants::{
        EPSILON_MIN, SA_CB1, SA_CB2, SA_CT1, SA_CT2, SA_CT3, SA_CT4, SA_CV1, SA_CV2, SA_CW1,
        SA_CW2, SA_CW3, SA_KAPPA_SQ, SA_SIGMA,
    },
    wall_distance::{self, cbrt},
};
use cfd_core::physics::constants::mathematical::numeric::{ONE, ONE_HALF, TWO};
use eunomia::{FloatElement, NumericElement, RealField};
use tracing::instrument;

/// Spalart-Allmaras turbulence model
///
/// Solves transport equation for ν̃ (modified turbulent kinematic viscosity):
/// Dν̃/Dt = Cb1 * S̃ * ν̃ - Cw1 * fw * (ν̃/d)² + (1/σ) * ∇·[(ν + ν̃) ∇ν̃] + Cb2/σ * |∇ν̃|²
///
/// where:
/// - S̃ = modified vorticity magnitude
/// - d = distance to nearest wall
/// - fw = wall destruction function
#[derive(Debug)]
pub struct SpalartAllmaras<T: RealField> {
    /// Grid dimensions
    pub(super) nx: usize,
    pub(super) ny: usize,
    /// Model coefficients
    pub(super) cb1: T,
    pub(super) cb2: T,
    pub(super) cw1: T,
    pub(super) cw2: T,
    pub(super) cw3: T,
    pub(super) cv1: T,
    pub(super) cv2: T,
    pub(super) ct1: T,
    pub(super) ct2: T,
    pub(super) ct3: T,
    pub(super) ct4: T,
    pub(super) sigma: T,
    pub(super) kappa_sq: T,
}

impl<T: RealField> SpalartAllmaras<T> {
    /// Create new Spalart-Allmaras model
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            cb1: T::from_f64(SA_CB1),
            cb2: T::from_f64(SA_CB2),
            cw1: T::from_f64(SA_CW1),
            cw2: T::from_f64(SA_CW2),
            cw3: T::from_f64(SA_CW3),
            cv1: T::from_f64(SA_CV1),
            cv2: T::from_f64(SA_CV2),
            ct1: T::from_f64(SA_CT1),
            ct2: T::from_f64(SA_CT2),
            ct3: T::from_f64(SA_CT3),
            ct4: T::from_f64(SA_CT4),
            sigma: T::from_f64(SA_SIGMA),
            kappa_sq: T::from_f64(SA_KAPPA_SQ),
        }
    }

    /// Calculate eddy viscosity from modified viscosity
    ///
    /// νt = ν̃ * fv1
    /// where fv1 = χ³/(χ³ + Cv1³) and χ = ν̃/ν
    #[instrument(skip(self))]
    pub fn eddy_viscosity(&self, nu_tilde: T, molecular_viscosity: T) -> T {
        let chi = nu_tilde / molecular_viscosity.max_scalar(T::from_f64(EPSILON_MIN));
        let chi_cubed = chi * chi * chi;
        let cv1_cubed = self.cv1 * self.cv1 * self.cv1;
        let fv1 = chi_cubed / (chi_cubed + cv1_cubed);
        nu_tilde * fv1
    }

    /// Calculate vorticity magnitude
    ///
    /// Ω = √(2 * Ωij * Ωij)
    /// where Ωij = 0.5 * (∂ui/∂xj - ∂uj/∂xi) is the rotation tensor
    #[instrument(skip(self))]
    pub(super) fn vorticity_magnitude(&self, velocity_gradient: &[[T; 2]; 2]) -> T {
        // Rotation tensor Ωij = 0.5 * (dui/dxj - duj/dxi)
        // Ω12 = 0.5 * (du/dy - dv/dx)
        // Ω21 = 0.5 * (dv/dx - du/dy) = -Ω12
        let omega12 = (velocity_gradient[0][1] - velocity_gradient[1][0]) * T::from_f64(ONE_HALF);

        // For 2D: Ω = √(2 * (Ω12² + Ω21²)) = √(2 * 2 * Ω12²) = 2|Ω12|
        let two = T::from_f64(TWO);
        two * <T as NumericElement>::abs(omega12)
    }

    /// Calculate modified vorticity S̃
    ///
    /// S̃ = Ω + (ν̃/κ²d²) * fv2
    /// where fv2 = 1 - χ/(1 + χ * fv1)
    #[instrument(skip(self))]
    pub(super) fn modified_vorticity(
        &self,
        vorticity: T,
        nu_tilde: T,
        molecular_viscosity: T,
        wall_distance: T,
    ) -> T {
        let chi = nu_tilde / molecular_viscosity.max_scalar(T::from_f64(EPSILON_MIN));
        let chi_cubed = chi * chi * chi;
        let cv1_cubed = self.cv1 * self.cv1 * self.cv1;
        let fv1 = chi_cubed / (chi_cubed + cv1_cubed);

        // fv2 = 1 - χ/(1 + χ * fv1)
        let one = T::from_f64(ONE);
        let fv2 = one - chi / (one + chi * fv1);

        // S̃ = Ω + ν̃/(κ²d²) * fv2
        let d_sq = wall_distance * wall_distance;
        let modification =
            (nu_tilde / (self.kappa_sq * d_sq.max_scalar(T::from_f64(EPSILON_MIN)))) * fv2;

        vorticity + modification
    }

    /// Calculate production term
    ///
    /// P = Cb1 * S̃ * ν̃
    #[instrument(skip(self))]
    pub fn production(&self, nu_tilde: T, s_tilde: T) -> T {
        self.cb1 * s_tilde * nu_tilde
    }

    /// Calculate wall destruction function fw
    ///
    /// fw = g * [(1 + Cw3⁶)/(g⁶ + Cw3⁶)]^(1/6)
    /// where g = r + Cw2 * (r⁶ - r) and r = ν̃/(S̃ * κ²d²)
    #[instrument(skip(self))]
    pub(super) fn wall_destruction_function(&self, nu_tilde: T, s_tilde: T, wall_distance: T) -> T {
        let d_sq = wall_distance * wall_distance;
        let denominator = s_tilde * self.kappa_sq * d_sq;

        // Limit r to prevent division by zero
        let r = if <T as NumericElement>::abs(denominator) > T::from_f64(EPSILON_MIN) {
            nu_tilde / denominator
        } else {
            T::from_f64(EPSILON_MIN)
        };

        // Limit r to reasonable range (prevents numerical issues)
        let r = r.min_scalar(T::from_f64(10.0));

        // g = r + Cw2 * (r⁶ - r)
        let r_sq = r * r;
        let r_6 = r_sq * r_sq * r_sq;
        let g = r + self.cw2 * (r_6 - r);

        // fw = g * [(1 + Cw3⁶)/(g⁶ + Cw3⁶)]^(1/6)
        let cw3_sq = self.cw3 * self.cw3;
        let cw3_6 = cw3_sq * cw3_sq * cw3_sq;
        let g_sq = g * g;
        let g_6 = g_sq * g_sq * g_sq;

        let one = T::from_f64(ONE);
        let ratio = (one + cw3_6) / (g_6 + cw3_6);

        // Compute sixth root using helper: x^(1/6) = (x^(1/2))^(1/3)
        let sqrt_ratio = <T as NumericElement>::sqrt(ratio);
        let cbrt_sqrt = cbrt(sqrt_ratio);

        g * cbrt_sqrt
    }

    /// Calculate destruction term
    ///
    /// D = Cw1 * fw * (ν̃/d)²
    #[instrument(skip(self))]
    pub fn destruction(&self, nu_tilde: T, wall_distance: T, fw: T) -> T {
        let ratio = nu_tilde / wall_distance.max_scalar(T::from_f64(EPSILON_MIN));
        self.cw1 * fw * ratio * ratio
    }

    /// Calculate trip term for transition
    ///
    /// Typically zero for fully turbulent flows, non-zero near laminar-turbulent transition
    /// Ft = Ct1 * gt * exp(-Ct2 * ωt²/Δν² - Ct3 * √(ωt * Δν) * d²)
    #[instrument(skip(self))]
    pub fn trip_term(&self, _nu_tilde: T, _wall_distance: T) -> T {
        // For fully turbulent flows, trip term is zero
        // Can be activated for transition modeling if needed
        T::ZERO
    }

    /// Laminar suppression function $f_{t2} = C_{t3} \exp(-C_{t4} \chi^2)$
    ///
    /// Prevents growth of $\tilde{\nu}$ in laminar regions.
    /// Uses `cv2` as the negative-SA cutoff: ignored when $\chi < cv2$.
    #[instrument(skip(self))]
    pub fn ft2(&self, nu_tilde: T, molecular_viscosity: T) -> T {
        let chi = nu_tilde / molecular_viscosity.max_scalar(T::from_f64(EPSILON_MIN));
        // Skip for large chi (fully turbulent) or when chi < cv2 (negative SA guard)
        if chi < self.cv2 {
            return T::ZERO;
        }
        self.ct3 * <T as FloatElement>::exp(-self.ct4 * chi * chi)
    }

    /// Trip-term coefficient $C_{t1}$
    #[must_use]
    pub fn ct1(&self) -> T {
        self.ct1
    }

    /// Trip-term coefficient $C_{t2}$
    #[must_use]
    pub fn ct2(&self) -> T {
        self.ct2
    }

    /// Calculate the rectangular wall-distance field used by the turbulence
    /// boundary manager and the SA wall-damping terms.
    #[instrument(skip(self, dx, dy))]
    pub fn wall_distance_field(&self, dx: T, dy: T) -> Vec<T> {
        wall_distance::wall_distance_field_2d(self.nx, self.ny, dx, dy)
    }

    /// Apply boundary conditions
    #[instrument(skip(self, nu_tilde))]
    pub(super) fn apply_boundary_conditions(&self, nu_tilde: &mut [T]) {
        // Wall boundaries: ν̃ = 0
        for i in 0..self.nx {
            nu_tilde[i] = T::ZERO; // Bottom wall
            nu_tilde[(self.ny - 1) * self.nx + i] = T::ZERO; // Top wall
        }

        for j in 0..self.ny {
            nu_tilde[j * self.nx] = T::ZERO; // Left wall
            nu_tilde[j * self.nx + (self.nx - 1)] = T::ZERO; // Right wall
        }
    }
}
