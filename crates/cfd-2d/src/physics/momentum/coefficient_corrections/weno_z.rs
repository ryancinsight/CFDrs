//! WENO-Z convection corrections.
//!
//! This deferred-correction path uses the fifth-order WENO-Z reconstruction of
//! Borges et al. to compute a high-order face flux correction on the momentum
//! equations. The implementation reuses the shared WENO-5 stencil helpers from
//! the scheme module so the coefficient path and the scalar scheme remain
//! mathematically consistent.

use crate::fields::Field2D;
use crate::fields::SimulationFields;
use crate::schemes::constants::WENO_EPSILON;
use crate::schemes::weno_helpers::{
    weno5_candidate_fluxes, weno5_smoothness_indicators, weno5_z_weights,
};
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::super::solver::MomentumComponent;

#[inline]
fn weno5_face_value<T>(values: [T; 5], velocity: T, epsilon: T) -> T
where
    T: RealField + Copy + FromPrimitive,
{
    let samples = if velocity > T::zero() {
        values
    } else {
        [values[4], values[3], values[2], values[1], values[0]]
    };
    let beta = weno5_smoothness_indicators(&samples);
    let weights = weno5_z_weights(epsilon, &beta);
    let flux = weno5_candidate_fluxes(&samples);
    weights[0] * flux[0] + weights[1] * flux[1] + weights[2] * flux[2]
}

#[inline]
pub fn compute_weno_z_correction_x<T: RealField + Copy + FromPrimitive>(
    i: usize,
    j: usize,
    u: T,
    rho: T,
    dy: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let epsilon =
        T::from_f64(WENO_EPSILON).expect("analytical constant conversion for WENO epsilon");
    let values = match component {
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

    let phi_face = weno5_face_value(values, u, epsilon);
    let phi_upwind = if u > T::zero() { values[2] } else { values[3] };

    rho * u * dy * (phi_face - phi_upwind)
}

#[inline]
pub fn compute_weno_z_correction_y<T: RealField + Copy + FromPrimitive>(
    i: usize,
    j: usize,
    v: T,
    rho: T,
    dx: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let epsilon =
        T::from_f64(WENO_EPSILON).expect("analytical constant conversion for WENO epsilon");
    let values = match component {
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

    let phi_face = weno5_face_value(values, v, epsilon);
    let phi_upwind = if v > T::zero() { values[2] } else { values[3] };

    rho * v * dx * (phi_face - phi_upwind)
}

#[inline]
pub fn apply_weno_z_deferred_correction<T: RealField + Copy + FromPrimitive>(
    source: &mut Field2D<T>,
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    u: T,
    v: T,
    rho: T,
    dx: T,
    dy: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
    relaxation_factor: f64,
) {
    let alpha = T::from_f64(relaxation_factor)
        .unwrap_or_else(|| T::from_f64(0.7).expect("analytical constant conversion"));

    let correction_x = if i >= 2 && i < nx - 2 {
        compute_weno_z_correction_x(i, j, u, rho, dy, fields, component)
    } else {
        T::zero()
    };

    let correction_y = if j >= 2 && j < ny - 2 {
        compute_weno_z_correction_y(i, j, v, rho, dx, fields, component)
    } else {
        T::zero()
    };

    if let Some(source) = source.at_mut(i, j) {
        *source += alpha * (correction_x + correction_y);
    }
}
