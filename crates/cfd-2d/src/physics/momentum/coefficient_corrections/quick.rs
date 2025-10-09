//! QUICK scheme correction computations
//!
//! Implements QUICK (Quadratic Upstream Interpolation for Convective Kinematics)
//! corrections for deferred correction approach.
//!
//! # References
//! * Leonard, B.P. (1979). "A stable and accurate convective modelling procedure"
//! * Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow", §5.4.3

use crate::discretization::extended_stencil::{ExtendedStencilScheme, QuickScheme};
use crate::fields::SimulationFields;
use nalgebra::RealField;

use super::super::solver::MomentumComponent;

/// Compute QUICK correction for X-direction convection
///
/// Returns the difference between QUICK flux and upwind flux
///
/// # References
/// * Leonard (1979), Patankar (1980) §5.4.3
pub fn compute_quick_correction_x<T: RealField + Copy>(
    i: usize,
    j: usize,
    u: T,
    rho: T,
    dy: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let quick = QuickScheme;

    // Get velocity values for 5-point stencil [i-2, i-1, i, i+1, i+2]
    let phi_values: [T; 5] = match component {
        MomentumComponent::U => [
            fields.u.at(i - 2, j),
            fields.u.at(i - 1, j),
            fields.u.at(i, j),
            fields.u.at(i + 1, j),
            fields.u.at(i + 2, j),
        ],
        MomentumComponent::V => [
            fields.v.at(i - 2, j),
            fields.v.at(i - 1, j),
            fields.v.at(i, j),
            fields.v.at(i + 1, j),
            fields.v.at(i + 2, j),
        ],
    };

    // QUICK face value at east face (between i and i+1)
    let phi_e_quick = quick.face_value(&phi_values, u, None);

    // Upwind face value at east face
    let phi_e_upwind = if u > T::zero() {
        phi_values[2] // From cell i (current)
    } else {
        phi_values[3] // From cell i+1 (downstream)
    };

    // Flux correction = mass_flux * (φ_QUICK - φ_upwind)
    let mass_flux = rho * u * dy;
    mass_flux * (phi_e_quick - phi_e_upwind)
}

/// Compute QUICK correction for Y-direction convection
///
/// Returns the difference between QUICK flux and upwind flux
///
/// # References
/// * Leonard (1979), Patankar (1980) §5.4.3
pub fn compute_quick_correction_y<T: RealField + Copy>(
    i: usize,
    j: usize,
    v: T,
    rho: T,
    dx: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let quick = QuickScheme;

    // Get velocity values for 5-point stencil [j-2, j-1, j, j+1, j+2]
    let phi_values: [T; 5] = match component {
        MomentumComponent::U => [
            fields.u.at(i, j - 2),
            fields.u.at(i, j - 1),
            fields.u.at(i, j),
            fields.u.at(i, j + 1),
            fields.u.at(i, j + 2),
        ],
        MomentumComponent::V => [
            fields.v.at(i, j - 2),
            fields.v.at(i, j - 1),
            fields.v.at(i, j),
            fields.v.at(i, j + 1),
            fields.v.at(i, j + 2),
        ],
    };

    // QUICK face value at north face (between j and j+1)
    let phi_n_quick = quick.face_value(&phi_values, v, None);

    // Upwind face value at north face
    let phi_n_upwind = if v > T::zero() {
        phi_values[2] // From cell j (current)
    } else {
        phi_values[3] // From cell j+1 (downstream)
    };

    // Flux correction = mass_flux * (φ_QUICK - φ_upwind)
    let mass_flux = rho * v * dx;
    mass_flux * (phi_n_quick - phi_n_upwind)
}
