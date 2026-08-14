//! Numeric helpers shared by the 3D kernels.
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
pub(crate) fn min<T: NumericElement>(left: T, right: T) -> T {
    <T as NumericElement>::min_scalar(left, right)
}

#[inline]
pub(crate) fn max<T: NumericElement>(left: T, right: T) -> T {
    <T as NumericElement>::max_scalar(left, right)
}

#[inline]
pub(crate) fn sqrt<T: NumericElement>(value: T) -> T {
    <T as NumericElement>::sqrt(value)
}

#[inline]
pub(crate) fn cos<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::cos(value)
}

#[inline]
pub(crate) fn sin<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::sin(value)
}

#[inline]
pub(crate) fn floor<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::floor(value)
}

#[inline]
pub(crate) fn powi<T: FloatElement>(value: T, exponent: i32) -> T {
    <T as FloatElement>::powi(value, exponent)
}

#[inline]
pub(crate) fn powf<T: FloatElement>(value: T, exponent: T) -> T {
    <T as FloatElement>::powf(value, exponent)
}

#[inline]
pub(crate) fn ln<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::ln(value)
}

#[inline]
pub(crate) fn is_finite<T: NumericElement>(value: T) -> bool {
    <T as NumericElement>::is_finite(value)
}

#[inline]
pub(crate) fn tanh<T: FloatElement>(value: T) -> T {
    <T as FloatElement>::tanh(value)
}
