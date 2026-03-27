/// Viscous dissipation function Phi for 2D incompressible flow.
///
/// ## Theorem — Viscous Dissipation (Bejan 2013)
///
/// The viscous dissipation function represents irreversible conversion
/// of kinetic energy to internal energy (heat):
///
/// ```text
/// Phi = 2*mu*[(du/dx)^2 + (dv/dy)^2] + mu*(du/dy + dv/dx)^2
/// ```
///
/// The Brinkman number Br = mu*U^2/(k*DeltaT) determines when Phi is significant:
/// - Br < 0.01: negligible (most millifluidic flows)
/// - Br ~ 0.1: moderate (high-shear venturi throats)
/// - Br > 1: dominant (polymer processing)
///
/// **Reference**: Bejan, A. (2013). *Convection Heat Transfer* (4th ed.),
/// Wiley, Section 2.5.
#[inline]
pub fn viscous_dissipation_2d(
    du_dx: f64,
    du_dy: f64,
    dv_dx: f64,
    dv_dy: f64,
    mu: f64,
) -> f64 {
    2.0 * mu * (du_dx * du_dx + dv_dy * dv_dy)
        + mu * (du_dy + dv_dx) * (du_dy + dv_dx)
}

/// Brinkman number: ratio of viscous heating to conductive heat transfer.
///
/// Br = mu * U_ref^2 / (k_thermal * delta_T)
///
/// A floor of 1e-30 is applied to delta_T to avoid division by zero.
#[inline]
pub fn brinkman_number(mu: f64, u_ref: f64, k_thermal: f64, delta_t: f64) -> f64 {
    mu * u_ref * u_ref / (k_thermal * delta_t.max(1e-30))
}

