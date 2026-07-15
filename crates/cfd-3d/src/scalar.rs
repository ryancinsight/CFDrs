//! Crate-local Eunomia scalar helpers for provider migration seams.

use eunomia::{FloatElement, NumericElement};
use leto_ops::RealScalar as LetoRealScalar;

/// Real scalar supported by the 3D solver stack during Atlas provider migration.
///
/// This trait is the crate-level scalar seam for code that crosses those contracts.
pub trait Cfd3dScalar:
    cfd_mesh::domain::core::Scalar
    + eunomia::RealField
    + LetoRealScalar
    + FloatElement
    + NumericElement
    + Copy
    + std::fmt::Debug
    + std::fmt::Display
    + std::ops::DivAssign
    + Send
    + Sync
    + 'static
{
    /// Additive identity from the Eunomia numeric contract.
    #[inline]
    fn zero() -> Self {
        <Self as NumericElement>::ZERO
    }

    /// Multiplicative identity from the Eunomia numeric contract.
    #[inline]
    fn one() -> Self {
        <Self as NumericElement>::ONE
    }
}

impl<T> Cfd3dScalar for T where
    T: cfd_mesh::domain::core::Scalar
        + eunomia::RealField
        + LetoRealScalar
        + FloatElement
        + NumericElement
        + Copy
        + std::fmt::Debug
        + std::fmt::Display
        + std::ops::DivAssign
        + Send
        + Sync
        + 'static
{
}

#[inline]
pub(crate) fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
pub(crate) fn from_usize<T: FloatElement>(value: usize) -> T {
    let value_u64 = u64::try_from(value).expect("invariant: grid index fits in u64");
    <T as FloatElement>::from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

#[inline]
pub(crate) fn to_f64<T: NumericElement>(value: T) -> f64 {
    <T as NumericElement>::to_f64(value)
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
