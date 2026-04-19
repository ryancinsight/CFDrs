//! Second-order upwind convection corrections.
//!
//! These corrections implement the actual second-order upwind stencil used by
//! the SIMPLEC/PIMPLE default configuration. The deferred-correction form keeps
//! the solver structure unchanged while avoiding the silent fallback to
//! first-order upwind.

use crate::fields::Field2D;
use crate::fields::SimulationFields;
use nalgebra::RealField;

use super::super::solver::MomentumComponent;

#[inline]
pub fn compute_second_order_correction_x<T: RealField + Copy>(
    i: usize,
    j: usize,
    u: T,
    rho: T,
    dy: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let half = T::one() / (T::one() + T::one());
    let phi_face = if u > T::zero() {
        match component {
            MomentumComponent::U => {
                let phi_im1 = fields.u.at(i - 1, j);
                let phi_i = fields.u.at(i, j);
                (T::one() + half) * phi_i - half * phi_im1
            }
            MomentumComponent::V => {
                let phi_im1 = fields.v.at(i - 1, j);
                let phi_i = fields.v.at(i, j);
                (T::one() + half) * phi_i - half * phi_im1
            }
        }
    } else {
        match component {
            MomentumComponent::U => {
                let phi_ip1 = fields.u.at(i + 1, j);
                let phi_ip2 = fields.u.at(i + 2, j);
                (T::one() + half) * phi_ip1 - half * phi_ip2
            }
            MomentumComponent::V => {
                let phi_ip1 = fields.v.at(i + 1, j);
                let phi_ip2 = fields.v.at(i + 2, j);
                (T::one() + half) * phi_ip1 - half * phi_ip2
            }
        }
    };

    let phi_upwind = if u > T::zero() {
        match component {
            MomentumComponent::U => fields.u.at(i, j),
            MomentumComponent::V => fields.v.at(i, j),
        }
    } else {
        match component {
            MomentumComponent::U => fields.u.at(i + 1, j),
            MomentumComponent::V => fields.v.at(i + 1, j),
        }
    };

    rho * u * dy * (phi_face - phi_upwind)
}

#[inline]
pub fn compute_second_order_correction_y<T: RealField + Copy>(
    i: usize,
    j: usize,
    v: T,
    rho: T,
    dx: T,
    fields: &SimulationFields<T>,
    component: MomentumComponent,
) -> T {
    let half = T::one() / (T::one() + T::one());
    let phi_face = if v > T::zero() {
        match component {
            MomentumComponent::U => {
                let phi_jm1 = fields.u.at(i, j - 1);
                let phi_j = fields.u.at(i, j);
                (T::one() + half) * phi_j - half * phi_jm1
            }
            MomentumComponent::V => {
                let phi_jm1 = fields.v.at(i, j - 1);
                let phi_j = fields.v.at(i, j);
                (T::one() + half) * phi_j - half * phi_jm1
            }
        }
    } else {
        match component {
            MomentumComponent::U => {
                let phi_jp1 = fields.u.at(i, j + 1);
                let phi_jp2 = fields.u.at(i, j + 2);
                (T::one() + half) * phi_jp1 - half * phi_jp2
            }
            MomentumComponent::V => {
                let phi_jp1 = fields.v.at(i, j + 1);
                let phi_jp2 = fields.v.at(i, j + 2);
                (T::one() + half) * phi_jp1 - half * phi_jp2
            }
        }
    };

    let phi_upwind = if v > T::zero() {
        match component {
            MomentumComponent::U => fields.u.at(i, j),
            MomentumComponent::V => fields.v.at(i, j),
        }
    } else {
        match component {
            MomentumComponent::U => fields.u.at(i, j + 1),
            MomentumComponent::V => fields.v.at(i, j + 1),
        }
    };

    rho * v * dx * (phi_face - phi_upwind)
}

#[inline]
pub fn apply_second_order_deferred_correction<T: RealField + Copy>(
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

    let correction_x = if i >= 1 && i < nx - 2 {
        compute_second_order_correction_x(i, j, u, rho, dy, fields, component)
    } else {
        T::zero()
    };

    let correction_y = if j >= 1 && j < ny - 2 {
        compute_second_order_correction_y(i, j, v, rho, dx, fields, component)
    } else {
        T::zero()
    };

    if let Some(source) = source.at_mut(i, j) {
        *source += alpha * (correction_x + correction_y);
    }
}
