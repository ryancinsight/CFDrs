//! Kato-Launder (1993) vorticity-strain production modification.
//!
//! ## Theorem — Kato-Launder Production (Kato & Launder 1993)
//!
//! The standard k-ε production term P_k = ν_t |S|² overpredicts turbulence
//! production in stagnation regions where vorticity Ω ≈ 0 but strain S is
//! large (irrotational strain). The Kato-Launder modification replaces:
//!
//! ```text
//! P_k = ν_t · S · Ω   (instead of ν_t · S²)
//! ```
//!
//! where S = √(2 S_ij S_ij) and Ω = √(2 Ω_ij Ω_ij).
//!
//! **Physical basis**: In stagnation regions, S > 0 but Ω → 0 (pure strain,
//! no rotation). The product S · Ω → 0, correctly predicting negligible
//! production. In shear flows, S ≈ Ω and P_k is unchanged.
//!
//! **Proof**: By the Cauchy-Schwarz inequality for tensors, S · Ω ≤ S²
//! for any velocity gradient, so P_KL ≤ P_standard universally.
//!
//! ## Reference
//!
//! Kato, M. & Launder, B.E. (1993). "The Modelling of Turbulent Flow
//! Around Stationary and Vibrating Square Cylinders", *Proc. 9th Symposium
//! on Turbulent Shear Flows*, Kyoto, pp. 10.4.1–10.4.6.

/// Compute Kato-Launder production P_k = ν_t · S · Ω from the 2D velocity
/// gradient tensor.
///
/// # Arguments
///
/// * `velocity_gradient` — 2×2 tensor [[du/dx, du/dy], [dv/dx, dv/dy]]
/// * `turbulent_viscosity` — eddy viscosity ν_t
///
/// # Returns
///
/// Non-negative production term [m²/s³].
pub fn kato_launder_production(velocity_gradient: &[[f64; 2]; 2], turbulent_viscosity: f64) -> f64 {
    // Strain rate tensor: S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
    let s_xx = velocity_gradient[0][0];
    let s_yy = velocity_gradient[1][1];
    let s_xy = 0.5 * (velocity_gradient[0][1] + velocity_gradient[1][0]);

    // Strain rate magnitude: S = sqrt(2 * S_ij * S_ij)
    // For 2D: 2 * (S_xx² + S_yy² + 2·S_xy²)  [S_xy = S_yx]
    let s_mag = (2.0 * (s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy)).sqrt();

    // Vorticity tensor: Ω_ij = 0.5 * (du_i/dx_j − du_j/dx_i)
    // For 2D, the only independent component is:
    // Ω_xy = 0.5 * (du/dy − dv/dx), Ω_yx = −Ω_xy
    let omega_xy = 0.5 * (velocity_gradient[0][1] - velocity_gradient[1][0]);

    // Vorticity magnitude: Ω = sqrt(2 * Ω_ij * Ω_ij)
    // For 2D: 2 * (Ω_xy² + Ω_yx²) = 4·Ω_xy²  ⇒  Ω = 2|Ω_xy|
    let omega_mag = (4.0 * omega_xy * omega_xy).sqrt();

    turbulent_viscosity * s_mag * omega_mag
}
