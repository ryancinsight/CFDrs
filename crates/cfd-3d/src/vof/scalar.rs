//! VOF-local Eunomia scalar helpers.

use eunomia::{FloatElement, NumericElement, RealField};

/// Scalar contract for VOF kernels over Atlas numeric providers.
pub trait VofScalar: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy {}

impl<T> VofScalar for T where T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy {}

#[inline]
pub(super) fn constant<T: VofScalar>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
pub(super) fn from_usize<T: VofScalar>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: grid index fits in u64");
    constant(<u64 as NumericElement>::to_f64(value_u64))
}

#[inline]
pub(super) fn zero<T: VofScalar>() -> T {
    <T as NumericElement>::ZERO
}

#[inline]
pub(super) fn one<T: VofScalar>() -> T {
    <T as NumericElement>::ONE
}

#[inline]
pub(super) fn abs<T: VofScalar>(value: T) -> T {
    <T as NumericElement>::abs(value)
}

#[inline]
pub(super) fn min<T: VofScalar>(left: T, right: T) -> T {
    <T as NumericElement>::min_scalar(left, right)
}

#[inline]
pub(super) fn max<T: VofScalar>(left: T, right: T) -> T {
    <T as NumericElement>::max_scalar(left, right)
}

#[inline]
pub(super) fn sqrt<T: VofScalar>(value: T) -> T {
    <T as NumericElement>::sqrt(value)
}

#[inline]
pub(super) fn is_finite<T: VofScalar>(value: T) -> bool {
    <T as NumericElement>::is_finite(value)
}

#[inline]
pub(super) fn tanh<T: VofScalar>(value: T) -> T {
    <T as FloatElement>::tanh(value)
}
