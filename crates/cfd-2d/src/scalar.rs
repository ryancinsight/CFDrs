//! Numeric helpers shared by the 2D kernels.
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
pub(crate) fn max<T: NumericElement>(left: T, right: T) -> T {
    <T as NumericElement>::max_scalar(left, right)
}
