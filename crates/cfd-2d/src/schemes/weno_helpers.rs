//! Shared WENO-5 stencil helpers.
//!
//! The fifth-order WENO family in this crate uses the same 5-point candidate
//! stencils and Jiang-Shu smoothness indicators.  Keeping the stencil algebra
//! in one place avoids duplicating the reconstruction logic across WENO-JS and
//! WENO-Z variants.

use super::{constants, weno_constants};
use nalgebra::RealField;
use num_traits::FromPrimitive;

#[inline]
pub(crate) fn weno5_smoothness_indicators<T>(v: &[T; 5]) -> [T; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let coeff_13_12 = T::from_f64(weno_constants::WENO5_BETA_COEFF_13_12)
        .expect("analytical constant conversion");
    let coeff_quarter = T::from_f64(weno_constants::WENO5_BETA_COEFF_QUARTER)
        .expect("analytical constant conversion");
    let two = T::from_f64(2.0).expect("analytical constant conversion");
    let three = T::from_f64(3.0).expect("analytical constant conversion");
    let four = T::from_f64(4.0).expect("analytical constant conversion");

    let beta0 = coeff_13_12 * (v[0] - two * v[1] + v[2]).powi(2)
        + coeff_quarter * (v[0] - four * v[1] + three * v[2]).powi(2);

    let beta1 =
        coeff_13_12 * (v[1] - two * v[2] + v[3]).powi(2) + coeff_quarter * (v[1] - v[3]).powi(2);

    let beta2 = coeff_13_12 * (v[2] - two * v[3] + v[4]).powi(2)
        + coeff_quarter * (three * v[2] - four * v[3] + v[4]).powi(2);

    [beta0, beta1, beta2]
}

#[inline]
pub(crate) fn weno5_candidate_fluxes<T>(v: &[T; 5]) -> [T; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let one_third =
        T::from_f64(weno_constants::WENO5_FLUX_ONE_THIRD).expect("analytical constant conversion");
    let seven_sixth = T::from_f64(weno_constants::WENO5_FLUX_SEVEN_SIXTH)
        .expect("analytical constant conversion");
    let eleven_sixth = T::from_f64(weno_constants::WENO5_FLUX_ELEVEN_SIXTH)
        .expect("analytical constant conversion");
    let one_sixth =
        T::from_f64(weno_constants::WENO5_FLUX_ONE_SIXTH).expect("analytical constant conversion");
    let five_sixth =
        T::from_f64(weno_constants::WENO5_FLUX_FIVE_SIXTH).expect("analytical constant conversion");

    let f0 = v[0] * one_third - seven_sixth * v[1] + eleven_sixth * v[2];
    let f1 = -v[1] * one_sixth + v[2] * five_sixth + v[3] * one_third;
    let f2 = v[2] * one_third + v[3] * five_sixth - v[4] * one_sixth;

    [f0, f1, f2]
}

#[inline]
pub(crate) fn weno5_js_weights<T>(epsilon: T, beta: &[T; 3]) -> [T; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let d0 = T::from_f64(constants::WENO5_WEIGHTS[0]).expect("analytical constant conversion");
    let d1 = T::from_f64(constants::WENO5_WEIGHTS[1]).expect("analytical constant conversion");
    let d2 = T::from_f64(constants::WENO5_WEIGHTS[2]).expect("analytical constant conversion");

    let alpha0 = d0 / (epsilon + beta[0]).powi(2);
    let alpha1 = d1 / (epsilon + beta[1]).powi(2);
    let alpha2 = d2 / (epsilon + beta[2]).powi(2);

    let sum = alpha0 + alpha1 + alpha2;
    [alpha0 / sum, alpha1 / sum, alpha2 / sum]
}

#[inline]
pub(crate) fn weno5_z_weights<T>(epsilon: T, beta: &[T; 3]) -> [T; 3]
where
    T: RealField + Copy + FromPrimitive,
{
    let d0 = T::from_f64(constants::WENO5_WEIGHTS[0]).expect("analytical constant conversion");
    let d1 = T::from_f64(constants::WENO5_WEIGHTS[1]).expect("analytical constant conversion");
    let d2 = T::from_f64(constants::WENO5_WEIGHTS[2]).expect("analytical constant conversion");
    let tau5 = (beta[0] - beta[2]).abs();

    let alpha0 = d0 * (T::one() + (tau5 / (epsilon + beta[0])).powi(2));
    let alpha1 = d1 * (T::one() + (tau5 / (epsilon + beta[1])).powi(2));
    let alpha2 = d2 * (T::one() + (tau5 / (epsilon + beta[2])).powi(2));

    let sum = alpha0 + alpha1 + alpha2;
    [alpha0 / sum, alpha1 / sum, alpha2 / sum]
}
