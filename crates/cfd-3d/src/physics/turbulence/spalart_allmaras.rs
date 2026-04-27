//! Spalart-Allmaras one-equation turbulence model (Spalart & Allmaras 1992).
//!
//! # Theorem -- Spalart-Allmaras nu_tilde Transport (Spalart & Allmaras 1992)
//!
//! A single modified turbulent working variable nu_tilde is governed by:
//!
//! ```text
//! D(nu_tilde)/Dt = C_b1*(1-f_t2)*S_tilde*nu_tilde
//!                + (1/sigma)*[div((nu+nu_tilde)*grad(nu_tilde)) + C_b2*(grad nu_tilde)^2]
//!                - [C_w1*f_w - (C_b1/kappa^2)*f_t2]*(nu_tilde/d)^2
//! ```
//!
//! with eddy viscosity:
//! ```text
//! nu_t = nu_tilde * f_v1
//! f_v1 = chi^3 / (chi^3 + C_v1^3)
//! chi  = nu_tilde / nu
//! ```
//!
//! and the modified vorticity:
//! ```text
//! S_tilde = S + nu_tilde/(kappa*d)^2 * f_v2
//! f_v2    = 1 - chi / (1 + chi*f_v1)
//! ```
//!
//! ## Algorithm -- Steady-State Initialization
//!
//! ```text
//! 1. Set nu_tilde_0 = 3*nu  (standard free-stream initial condition)
//! 2. Evaluate f_v1 = chi^3/(chi^3 + C_v1^3)  (wall damping)
//! 3. nu_t = nu_tilde * f_v1  (bounded near walls by f_v1 -> 0)
//! ```
//!
//! ## References
//!
//! - Spalart, P.R. & Allmaras, S.R. (1992). A one-equation turbulence model
//!   for aerodynamic flows. AIAA Paper 92-0439.
//! - Allmaras, S.R., Johnson, F.T. & Spalart, P.R. (2012). Modifications and
//!   clarifications for the Spalart-Allmaras model. ICCFD7-1902.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::SA_CV1;
#[cfg(test)]
use super::constants::{SA_CW2, SA_CW3};
use super::sgs_energy::kinetic_energy_from_eddy_viscosity;

/// Spalart-Allmaras one-equation RANS model (Spalart & Allmaras 1992).
///
/// Uses a single modified turbulent viscosity nu_tilde as the working variable.
/// Suitable for aerodynamic flows with mild separation; inexpensive for 3D FEM.
#[derive(Debug, Clone)]
pub struct SpalartAllmarasModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Kinematic viscosity nu [m^2/s].
    pub nu: T,
    /// Modified turbulent working variable nu_tilde at each grid point.
    pub nu_tilde: Vec<T>,
    /// Per-point wall distance d [m].
    pub wall_distance: Vec<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> SpalartAllmarasModel<T> {
    /// Create SA model with uniform free-stream initialisation.
    ///
    /// Sets nu_tilde_0 = 3*nu (Spalart & Allmaras 1992, Section 4).
    ///
    /// # Arguments
    /// * `n_points`       -- number of grid points (nx*ny*nz)
    /// * `nu`             -- kinematic viscosity [m^2/s]
    /// * `wall_distances` -- per-point wall distances [m]
    pub fn new(n_points: usize, nu: T, wall_distances: Vec<T>) -> Self {
        let three = T::one() + T::one() + T::one();
        let nu_tilde = vec![three * nu; n_points];
        assert_eq!(
            wall_distances.len(),
            n_points,
            "SA wall_distances must match n_points"
        );
        Self {
            nu,
            nu_tilde,
            wall_distance: wall_distances,
        }
    }

    /// Initialise with a prescribed nu_tilde field.
    pub fn with_nu_tilde(nu: T, nu_tilde: Vec<T>, wall_distances: Vec<T>) -> Self {
        let n = nu_tilde.len();
        assert_eq!(
            wall_distances.len(),
            n,
            "SA wall_distances must match nu_tilde length"
        );
        Self {
            nu,
            nu_tilde,
            wall_distance: wall_distances,
        }
    }

    /// Compute the wall-damping function f_v1 = chi^3 / (chi^3 + C_v1^3).
    fn f_v1(chi: T) -> T {
        let cv1 = <T as FromPrimitive>::from_f64(SA_CV1)
            .expect("SA_CV1 is an IEEE 754 representable f64 constant");
        let cv1_3 = cv1 * cv1 * cv1;
        let chi_3 = chi * chi * chi;
        chi_3 / (chi_3 + cv1_3)
    }

    /// Compute chi = nu_tilde / nu.
    fn chi(nu_tilde: T, nu: T) -> T {
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");
        nu_tilde / (nu + eps)
    }

    /// Convert SA eddy viscosity to a diagnostic turbulence kinetic energy.
    ///
    /// # Theorem -- SA Eddy-Viscosity Energy Diagnostic
    ///
    /// The original Spalart-Allmaras model transports modified eddy viscosity
    /// rather than `k`. When a caller requires a kinetic-energy diagnostic, the
    /// wall distance supplies the local length scale and the Yoshizawa
    /// eddy-viscosity relation gives
    ///
    /// ```text
    /// k_d = (nu_t / (C_k d))^2.
    /// ```
    ///
    /// **Proof.** The SA closure produces `nu_t >= 0` after wall damping.
    /// Yoshizawa's algebraic relation `nu_t = C_k d sqrt(k_d)` is dimensionally
    /// consistent because `d sqrt(k_d)` has units `m^2/s`. Solving the relation
    /// for `k_d` gives the displayed expression. Squaring preserves
    /// non-negativity, and `d <= 0` maps to zero because no positive turbulent
    /// length scale exists at the wall point.
    #[inline]
    fn diagnostic_kinetic_energy(nu_t: T, wall_distance: T) -> T {
        kinetic_energy_from_eddy_viscosity(nu_t, wall_distance)
    }

    /// Compute the destruction function f_w (Spalart & Allmaras 1992 eq. 25).
    #[cfg(test)]
    fn f_w(r: T) -> T {
        let cw2 = <T as FromPrimitive>::from_f64(SA_CW2)
            .expect("SA_CW2 is an IEEE 754 representable f64 constant");
        let cw3 = <T as FromPrimitive>::from_f64(SA_CW3)
            .expect("SA_CW3 is an IEEE 754 representable f64 constant");
        let cw3_6 = num_traits::Float::powi(cw3, 6);
        let g = r + cw2 * (num_traits::Float::powi(r, 6) - r);
        let g_6 = num_traits::Float::powi(g, 6);
        g * num_traits::Float::powf(
            (T::one() + cw3_6) / (g_6 + cw3_6),
            T::one() / (T::one() + T::one() + T::one() + T::one() + T::one() + T::one()),
        )
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for SpalartAllmarasModel<T>
{
    /// Compute SA eddy viscosity nu_t = nu_tilde * f_v1.
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let n = flow_field.velocity.components.len();
        assert_eq!(
            self.nu_tilde.len(),
            n,
            "SA nu_tilde must match the flow-field size"
        );
        assert_eq!(
            self.wall_distance.len(),
            n,
            "SA wall_distance must match the flow-field size"
        );
        let mut viscosity = Vec::with_capacity(n);
        for idx in 0..n {
            let nu_t = self.nu_tilde[idx];
            let chi = Self::chi(nu_t, self.nu);
            let fv1 = Self::f_v1(chi);
            viscosity.push(num_traits::Float::max(T::zero(), nu_t * fv1));
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let eddy_viscosity = self.turbulent_viscosity(flow_field);
        eddy_viscosity
            .into_iter()
            .zip(self.wall_distance.iter().copied())
            .map(|(nu_t, wall_distance)| Self::diagnostic_kinetic_energy(nu_t, wall_distance))
            .collect()
    }

    fn name(&self) -> &'static str {
        "Spalart-Allmaras"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fv1_at_zero_chi() {
        // At chi = 0: f_v1 = 0/(0 + C_v1^3) = 0  (no turbulence at wall)
        let fv1 = SpalartAllmarasModel::<f64>::f_v1(0.0_f64);
        assert!(fv1.abs() < 1e-10, "f_v1(0) = {fv1}, expected 0");
    }

    #[test]
    fn test_fv1_large_chi() {
        // At chi -> inf: f_v1 -> 1  (fully turbulent)
        let fv1 = SpalartAllmarasModel::<f64>::f_v1(1000.0_f64);
        assert!((fv1 - 1.0).abs() < 1e-3, "f_v1(1000) = {fv1}, expected ~1");
    }

    #[test]
    fn test_viscosity_nonnegative() {
        // nu_t must be >= 0 everywhere
        let n = 8;
        let nu = 1e-6_f64;
        let model = SpalartAllmarasModel::<f64>::new(n, nu, vec![1.0; n]);
        let chi = SpalartAllmarasModel::<f64>::chi(model.nu_tilde[0], nu);
        let fv1 = SpalartAllmarasModel::<f64>::f_v1(chi);
        assert!(fv1 >= 0.0, "f_v1 must be non-negative: {fv1}");
    }

    #[test]
    fn test_f_w_positive() {
        // f_w must be positive for any r > 0
        let fw = SpalartAllmarasModel::<f64>::f_w(1.0_f64);
        assert!(fw > 0.0, "f_w(1.0) must be positive: {fw}");
    }

    #[test]
    fn test_tke_diagnostic_matches_yoshizawa_relation() {
        let nu = 1.0e-6_f64;
        let wall_distance = 2.0e-3_f64;
        let nu_tilde = 6.0e-6_f64;
        let model =
            SpalartAllmarasModel::<f64>::with_nu_tilde(nu, vec![nu_tilde], vec![wall_distance]);
        let flow = FlowField::<f64>::new(1, 1, 1);
        let nu_t = model.turbulent_viscosity(&flow)[0];
        let k = model.turbulent_kinetic_energy(&flow)[0];
        let expected = (nu_t / (0.094 * wall_distance)).powi(2);

        assert!((k - expected).abs() < 1.0e-18);
        assert!(k > 0.0);
    }

    #[test]
    fn test_tke_diagnostic_preserves_zero_wall_length_scale() {
        let model = SpalartAllmarasModel::<f64>::with_nu_tilde(1.0e-6, vec![3.0e-6], vec![0.0]);
        let flow = FlowField::<f64>::new(1, 1, 1);

        assert_eq!(model.turbulent_kinetic_energy(&flow)[0], 0.0);
    }
}
