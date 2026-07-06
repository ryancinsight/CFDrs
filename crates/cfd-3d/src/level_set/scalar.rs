//! Level-set-local Eunomia scalar helpers.

use eunomia::{FloatElement, NumericElement, RealField};

/// Scalar contract for level-set kernels over Atlas numeric providers.
pub trait LevelSetScalar: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy {}

impl<T> LevelSetScalar for T where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy
{
}

#[inline]
pub(super) fn from_f64<T: LevelSetScalar>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
pub(super) fn zero<T: LevelSetScalar>() -> T {
    <T as NumericElement>::ZERO
}

#[inline]
pub(super) fn one<T: LevelSetScalar>() -> T {
    <T as NumericElement>::ONE
}

#[inline]
pub(super) fn abs<T: LevelSetScalar>(value: T) -> T {
    <T as NumericElement>::abs(value)
}

#[inline]
pub(super) fn min<T: LevelSetScalar>(left: T, right: T) -> T {
    <T as NumericElement>::min_scalar(left, right)
}

#[inline]
pub(super) fn max<T: LevelSetScalar>(left: T, right: T) -> T {
    <T as NumericElement>::max_scalar(left, right)
}

#[inline]
pub(super) fn sqrt<T: LevelSetScalar>(value: T) -> T {
    <T as NumericElement>::sqrt(value)
}

#[inline]
pub(super) fn powi<T: LevelSetScalar>(value: T, exponent: i32) -> T {
    <T as FloatElement>::powi(value, exponent)
}

#[inline]
pub(super) fn is_finite<T: LevelSetScalar>(value: T) -> bool {
    <T as NumericElement>::is_finite(value)
}
