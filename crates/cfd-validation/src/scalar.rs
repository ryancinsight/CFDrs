//! Numeric helpers shared by the validation kernels.
//!
//! The scalar seam itself is [`cfd_core::CfdScalar`]; this module holds only
//! crate-local shorthands over the Eunomia numeric surface.

use eunomia::{FloatElement, NumericElement};

#[inline]
pub(crate) fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: grid index fits in u64");
    <T as FloatElement>::from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

#[inline]
pub(crate) fn zero<T: NumericElement>() -> T {
    <T as NumericElement>::ZERO
}

#[inline]
pub(crate) fn one<T: NumericElement>() -> T {
    <T as NumericElement>::ONE
}

#[inline]
pub(crate) fn abs<T: NumericElement>(value: T) -> T {
    <T as NumericElement>::abs(value)
}

#[inline]
pub(crate) fn sqrt<T: NumericElement>(value: T) -> T {
    <T as NumericElement>::sqrt(value)
}

#[inline]
pub(crate) fn sin<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::sin(value)
}

#[inline]
pub(crate) fn cos<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::cos(value)
}

#[inline]
pub(crate) fn cosh<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::cosh(value)
}

#[inline]
pub(crate) fn exp<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::exp(value)
}

#[inline]
pub(crate) fn ln<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::ln(value)
}

#[inline]
pub(crate) fn tanh<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::tanh(value)
}

#[inline]
pub(crate) fn atan2<T: FloatElement>(y: T, x: T) -> T {
    <T as FloatElement>::atan2(y, x)
}

#[inline]
pub(crate) fn powf<T: FloatElement>(value: T, exponent: T) -> T {
    <T as FloatElement>::powf(value, exponent)
}

#[inline]
pub(crate) fn powi<T: FloatElement>(value: T, exponent: i32) -> T {
    <T as FloatElement>::powi(value, exponent)
}

#[inline]
pub(crate) fn cbrt<T: FloatElement>(value: T) -> T {
    powf(value, <T as FloatElement>::from_f64(1.0 / 3.0))
}

#[inline]
pub(crate) fn min<T: NumericElement>(a: T, b: T) -> T {
    <T as NumericElement>::min_scalar(a, b)
}

#[inline]
pub(crate) fn max<T: NumericElement>(a: T, b: T) -> T {
    <T as NumericElement>::max_scalar(a, b)
}
