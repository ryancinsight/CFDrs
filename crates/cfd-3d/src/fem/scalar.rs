//! FEM-local Eunomia scalar construction helpers.

use eunomia::{FloatElement, NumericElement};

/// Convert an analytical constant through Eunomia's scalar provider surface.
#[inline]
pub(super) fn constant<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Additive identity through Eunomia's scalar provider surface.
#[inline]
pub(super) fn zero<T: NumericElement>() -> T {
    <T as NumericElement>::ZERO
}

/// Multiplicative identity through Eunomia's scalar provider surface.
#[inline]
pub(super) fn one<T: NumericElement>() -> T {
    <T as NumericElement>::ONE
}
