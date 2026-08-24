//! Viscous dissipation function Phi and Brinkman number for 2D incompressible flow.
//!
//! ## Theorem — Viscous Dissipation (Bejan 2013)
//!
//! The viscous dissipation function represents irreversible conversion
//! of kinetic energy to internal energy (heat):
//!
//! ```text
//! Phi = 2*mu*[(du/dx)^2 + (dv/dy)^2] + mu*(du/dy + dv/dx)^2
//! ```
//!
//! The Brinkman number Br = mu*U^2/(k*DeltaT) determines when Phi is significant:
//! - Br < 0.01: negligible (most millifluidic flows)
//! - Br ~ 0.1: moderate (high-shear venturi throats)
//! - Br > 1: dominant (polymer processing)
//!
//! **Reference**: Bejan, A. (2013). *Convection Heat Transfer* (4th ed.),
//! Wiley, Section 2.5.

use cfd_core::error::{Error, Result};

/// Viscous dissipation function Phi for 2D incompressible flow.
#[inline]
pub fn viscous_dissipation_2d(du_dx: f64, du_dy: f64, dv_dx: f64, dv_dy: f64, mu: f64) -> f64 {
    match viscous_dissipation_2d_validated(du_dx, du_dy, dv_dx, dv_dy, mu) {
        Ok(value) => value,
        Err(error) => panic!("viscous_dissipation_2d called with invalid inputs: {error}"),
    }
}

/// Viscous dissipation function Phi, validating the inputs.
///
/// # Errors
///
/// Returns [`Error::InvalidConfiguration`] when any input is non-finite or
/// when `mu` is negative (negative viscosity is unphysical; a viscous fluid
/// with negative `mu` releases kinetic energy instead of dissipating it).
#[inline]
pub fn viscous_dissipation_2d_validated(
    du_dx: f64,
    du_dy: f64,
    dv_dx: f64,
    dv_dy: f64,
    mu: f64,
) -> Result<f64> {
    if !du_dx.is_finite() || !du_dy.is_finite() || !dv_dx.is_finite() || !dv_dy.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "viscous_dissipation_2d_validated: velocity gradients must be finite, got ({du_dx}, {du_dy}, {dv_dx}, {dv_dy})"
        )));
    }
    if !mu.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "viscous_dissipation_2d_validated: dynamic viscosity mu must be finite, got {mu}"
        )));
    }
    if mu < 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "viscous_dissipation_2d_validated: dynamic viscosity mu must be non-negative (negative viscosity is unphysical), got {mu}"
        )));
    }
    Ok(2.0 * mu * (du_dx * du_dx + dv_dy * dv_dy) + mu * (du_dy + dv_dx) * (du_dy + dv_dx))
}

/// Brinkman number: ratio of viscous heating to conductive heat transfer.
///
/// Br = mu * U_ref^2 / (k_thermal * delta_T)
#[inline]
pub fn brinkman_number(mu: f64, u_ref: f64, k_thermal: f64, delta_t: f64) -> f64 {
    match brinkman_number_validated(mu, u_ref, k_thermal, delta_t) {
        Ok(value) => value,
        Err(error) => panic!("brinkman_number called with invalid inputs: {error}"),
    }
}

/// Brinkman number, validating the inputs.
///
/// # Errors
///
/// Returns [`Error::InvalidConfiguration`] when:
/// - `mu` is non-finite or negative (negative viscosity is unphysical),
/// - `u_ref`, `k_thermal`, or `delta_t` is non-finite (the ratio is
///   undefined for non-finite values),
/// - `k_thermal` is zero or negative (the conductive heat transfer
///   denominator `k_thermal * delta_T` is the physical carrier of the
///   Brinkman number and must be strictly positive),
/// - `delta_t` is zero or negative (the temperature difference
///   representing the conductive heat scale must be strictly positive;
///   the previous silent `delta_t.max(1e-30)` floor masked the invalid input
///   by silently inflating Br by 30 orders of magnitude).
#[inline]
pub fn brinkman_number_validated(mu: f64, u_ref: f64, k_thermal: f64, delta_t: f64) -> Result<f64> {
    if !mu.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: dynamic viscosity mu must be finite, got {mu}"
        )));
    }
    if mu < 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: dynamic viscosity mu must be non-negative, got {mu}"
        )));
    }
    if !u_ref.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: u_ref must be finite, got {u_ref}"
        )));
    }
    if !k_thermal.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: k_thermal must be finite, got {k_thermal}"
        )));
    }
    if k_thermal <= 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: k_thermal must be strictly positive (it is the conductive heat transfer scale), got {k_thermal}"
        )));
    }
    if !delta_t.is_finite() {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: delta_t must be finite, got {delta_t}"
        )));
    }
    if delta_t <= 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "brinkman_number_validated: delta_t must be strictly positive (the conductive heat scale), got {delta_t}"
        )));
    }
    Ok(mu * u_ref * u_ref / (k_thermal * delta_t))
}
