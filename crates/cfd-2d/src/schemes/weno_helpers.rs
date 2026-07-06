//! Shared WENO-5 stencil helpers.
//!
//! The fifth-order WENO family in this crate uses the same 5-point candidate
//! stencils and Jiang-Shu smoothness indicators.  Keeping the stencil algebra
//! in one place avoids duplicating the reconstruction logic across WENO-JS and
//! WENO-Z variants.

use super::{constants, weno_constants};
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use eunomia::{FloatElement, NumericElement};

#[inline]
pub(crate) fn weno5_smoothness_indicators<T>(v: &[T; 5]) -> [T; 3]
where
    T: Cfd2dScalar + Copy + FloatElement,
{
    let coeff_13_12: T = scalar::from_f64(weno_constants::WENO5_BETA_COEFF_13_12);
    let coeff_quarter: T = scalar::from_f64(weno_constants::WENO5_BETA_COEFF_QUARTER);
    let two: T = scalar::from_f64(2.0);
    let three: T = scalar::from_f64(3.0);
    let four: T = scalar::from_f64(4.0);

    let beta0 = coeff_13_12 * FloatElement::powi(v[0] - two * v[1] + v[2], 2)
        + coeff_quarter * FloatElement::powi(v[0] - four * v[1] + three * v[2], 2);

    let beta1 = coeff_13_12 * FloatElement::powi(v[1] - two * v[2] + v[3], 2)
        + coeff_quarter * FloatElement::powi(v[1] - v[3], 2);

    let beta2 = coeff_13_12 * FloatElement::powi(v[2] - two * v[3] + v[4], 2)
        + coeff_quarter * FloatElement::powi(three * v[2] - four * v[3] + v[4], 2);

    [beta0, beta1, beta2]
}

#[inline]
pub(crate) fn weno5_candidate_fluxes<T>(v: &[T; 5]) -> [T; 3]
where
    T: Cfd2dScalar + Copy + FloatElement,
{
    let one_third: T = scalar::from_f64(weno_constants::WENO5_FLUX_ONE_THIRD);
    let seven_sixth: T = scalar::from_f64(weno_constants::WENO5_FLUX_SEVEN_SIXTH);
    let eleven_sixth: T = scalar::from_f64(weno_constants::WENO5_FLUX_ELEVEN_SIXTH);
    let one_sixth: T = scalar::from_f64(weno_constants::WENO5_FLUX_ONE_SIXTH);
    let five_sixth: T = scalar::from_f64(weno_constants::WENO5_FLUX_FIVE_SIXTH);

    let f0 = v[0] * one_third - seven_sixth * v[1] + eleven_sixth * v[2];
    let f1 = -v[1] * one_sixth + v[2] * five_sixth + v[3] * one_third;
    let f2 = v[2] * one_third + v[3] * five_sixth - v[4] * one_sixth;

    [f0, f1, f2]
}

#[inline]
pub(crate) fn weno5_js_weights<T>(epsilon: T, beta: &[T; 3]) -> [T; 3]
where
    T: Cfd2dScalar + Copy + FloatElement,
{
    let d0: T = scalar::from_f64(constants::WENO5_WEIGHTS[0]);
    let d1: T = scalar::from_f64(constants::WENO5_WEIGHTS[1]);
    let d2: T = scalar::from_f64(constants::WENO5_WEIGHTS[2]);

    let alpha0 = d0 / FloatElement::powi(epsilon + beta[0], 2);
    let alpha1 = d1 / FloatElement::powi(epsilon + beta[1], 2);
    let alpha2 = d2 / FloatElement::powi(epsilon + beta[2], 2);

    let sum = alpha0 + alpha1 + alpha2;
    [alpha0 / sum, alpha1 / sum, alpha2 / sum]
}

#[inline]
pub(crate) fn weno5_z_weights<T>(epsilon: T, beta: &[T; 3]) -> [T; 3]
where
    T: Cfd2dScalar + Copy + FloatElement,
{
    let d0: T = scalar::from_f64(constants::WENO5_WEIGHTS[0]);
    let d1: T = scalar::from_f64(constants::WENO5_WEIGHTS[1]);
    let d2: T = scalar::from_f64(constants::WENO5_WEIGHTS[2]);
    let tau5 = NumericElement::abs(beta[0] - beta[2]);

    let one: T = scalar::one();
    let alpha0 = d0 * (one + FloatElement::powi(tau5 / (epsilon + beta[0]), 2));
    let alpha1 = d1 * (one + FloatElement::powi(tau5 / (epsilon + beta[1]), 2));
    let alpha2 = d2 * (one + FloatElement::powi(tau5 / (epsilon + beta[2]), 2));

    let sum = alpha0 + alpha1 + alpha2;
    [alpha0 / sum, alpha1 / sum, alpha2 / sum]
}
