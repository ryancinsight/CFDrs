//! Scalar contract and construction helpers for cfd-2d provider migrations.

use eunomia::{FloatElement, NumericElement};
use leto::ScalarOperand;
use leto_ops::RealScalar as LetoRealScalar;

/// Real scalar supported by the 2D solver stack.
///
/// Migrated cfd-core boundary/fluid contracts require Eunomia and the 1D
/// network coupling requires `cfd_1d::Cfd1dScalar`. This trait is the single
/// scalar seam for those provider contracts.
pub trait Cfd2dScalar:
    cfd_1d::Cfd1dScalar
    + eunomia::RealField
    + LetoRealScalar
    + ScalarOperand
    + FloatElement
    + NumericElement
    + Copy
    + std::fmt::Debug
    + Send
    + Sync
    + 'static
{
}

impl<T> Cfd2dScalar for T where
    T: cfd_1d::Cfd1dScalar
        + eunomia::RealField
        + LetoRealScalar
        + ScalarOperand
        + FloatElement
        + NumericElement
        + Copy
        + std::fmt::Debug
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
    from_f64(<u64 as NumericElement>::to_f64(value_u64))
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
