//! TVD scheme correction computations
//!
//! Implements Total Variation Diminishing (TVD) flux limiters for
//! deferred correction approaches in high-Peclet number flows.
//!
//! # References
//! * Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
//! * Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"

use super::super::tvd_limiters::TvdLimiter;
use crate::fields::SimulationFields;
use nalgebra::RealField;

use super::super::solver::MomentumComponent;

/// Compute TVD correction for X-direction convection
///
/// Returns the difference between TVD-limited flux and upwind flux.
/// Uses 3-point stencil: upwind, central, downwind.
///
/// # References
/// * Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
/// * Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
pub fn compute_tvd_correction_x<T, L>(
    i: usize,
    j: usize,
    u: T,
    rho: T,
    dy: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
    limiter: &L,
) -> T
where
    T: RealField + Copy,
    L: TvdLimiter<T>,
{
    // Get velocity values for 3-point stencil
    let (phi_u, phi_c, phi_d) = if u > T::zero() {
        // Positive velocity: upwind is i-1, central is i, downwind is i+1
        match component {
            MomentumComponent::U => (
                fields.u.at(i - 1, j),
                fields.u.at(i, j),
                fields.u.at(i + 1, j),
            ),
            MomentumComponent::V => (
                fields.v.at(i - 1, j),
                fields.v.at(i, j),
                fields.v.at(i + 1, j),
            ),
        }
    } else {
        // Negative velocity: upwind is i+1, central is i, downwind is i-1
        match component {
            MomentumComponent::U => (
                fields.u.at(i + 1, j),
                fields.u.at(i, j),
                fields.u.at(i - 1, j),
            ),
            MomentumComponent::V => (
                fields.v.at(i + 1, j),
                fields.v.at(i, j),
                fields.v.at(i - 1, j),
            ),
        }
    };

    // TVD-limited face value
    let phi_e_tvd = limiter.interpolate_face(phi_u, phi_c, phi_d);

    // Upwind face value
    let phi_e_upwind = phi_c;

    // Flux correction = mass_flux * (φ_TVD - φ_upwind)
    let mass_flux = rho * u * dy;
    mass_flux * (phi_e_tvd - phi_e_upwind)
}

/// Compute TVD correction for Y-direction convection
///
/// Returns the difference between TVD-limited flux and upwind flux.
/// Uses 3-point stencil: upwind, central, downwind.
///
/// # References
/// * Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
/// * Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
pub fn compute_tvd_correction_y<T, L>(
    i: usize,
    j: usize,
    v: T,
    rho: T,
    dx: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
    limiter: &L,
) -> T
where
    T: RealField + Copy,
    L: TvdLimiter<T>,
{
    // Get velocity values for 3-point stencil
    let (phi_u, phi_c, phi_d) = if v > T::zero() {
        // Positive velocity: upwind is j-1, central is j, downwind is j+1
        match component {
            MomentumComponent::U => (
                fields.u.at(i, j - 1),
                fields.u.at(i, j),
                fields.u.at(i, j + 1),
            ),
            MomentumComponent::V => (
                fields.v.at(i, j - 1),
                fields.v.at(i, j),
                fields.v.at(i, j + 1),
            ),
        }
    } else {
        // Negative velocity: upwind is j+1, central is j, downwind is j-1
        match component {
            MomentumComponent::U => (
                fields.u.at(i, j + 1),
                fields.u.at(i, j),
                fields.u.at(i, j - 1),
            ),
            MomentumComponent::V => (
                fields.v.at(i, j + 1),
                fields.v.at(i, j),
                fields.v.at(i, j - 1),
            ),
        }
    };

    // TVD-limited face value
    let phi_n_tvd = limiter.interpolate_face(phi_u, phi_c, phi_d);

    // Upwind face value
    let phi_n_upwind = phi_c;

    // Flux correction = mass_flux * (φ_TVD - φ_upwind)
    let mass_flux = rho * v * dx;
    mass_flux * (phi_n_tvd - phi_n_upwind)
}
