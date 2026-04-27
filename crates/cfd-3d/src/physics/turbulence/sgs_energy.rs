//! Subgrid kinetic-energy closures shared by LES eddy-viscosity models.
//!
//! # Theorem — Yoshizawa SGS Energy Relation
//!
//! One-equation and algebraic LES closures commonly relate subgrid eddy
//! viscosity to SGS kinetic energy by
//!
//! ```text
//! ν_t = C_k Δ sqrt(k_sgs)
//! ```
//!
//! where `Δ` is the local LES length scale and `C_k` is the Yoshizawa
//! coefficient. Solving for `k_sgs` gives
//!
//! ```text
//! k_sgs = (ν_t / (C_k Δ))².
//! ```
//!
//! **Proof.** The eddy viscosity has dimensions `L²/T`, while
//! `Δ sqrt(k_sgs)` has dimensions `L · L/T = L²/T`; the coefficient is
//! dimensionless. Algebraic inversion preserves non-negativity because both
//! `ν_t` and `Δ` are non-negative in the validated LES closures. ∎

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Yoshizawa SGS kinetic-energy coefficient used in algebraic LES closures.
const YOSHIZAWA_C_K: f64 = 0.094;

/// Convert eddy viscosity to SGS turbulent kinetic energy.
#[inline]
pub(crate) fn kinetic_energy_from_eddy_viscosity<T>(eddy_viscosity: T, length_scale: T) -> T
where
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive,
{
    let eps = <T as FromPrimitive>::from_f64(1.0e-30)
        .expect("1e-30 is an IEEE 754 representable f64 constant");
    if eddy_viscosity <= T::zero() || length_scale <= eps {
        return T::zero();
    }

    let c_k = <T as FromPrimitive>::from_f64(YOSHIZAWA_C_K)
        .expect("YOSHIZAWA_C_K is an IEEE 754 representable f64 constant");
    let velocity_scale = eddy_viscosity / (c_k * length_scale);
    velocity_scale * velocity_scale
}

#[cfg(test)]
mod tests {
    use super::kinetic_energy_from_eddy_viscosity;

    #[test]
    fn yoshizawa_relation_inverts_eddy_viscosity_dimensions() {
        let nu_t = 2.0e-5_f64;
        let delta = 1.0e-3_f64;
        let k = kinetic_energy_from_eddy_viscosity(nu_t, delta);
        let c_k = 0.094_f64;

        assert!((k - (nu_t / (c_k * delta)).powi(2)).abs() < 1.0e-18);
        assert!(k > 0.0);
    }

    #[test]
    fn yoshizawa_relation_preserves_zero_viscosity_state() {
        assert_eq!(kinetic_energy_from_eddy_viscosity(0.0_f64, 1.0), 0.0);
    }
}
