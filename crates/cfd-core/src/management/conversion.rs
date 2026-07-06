//! Safe numeric conversion utilities
//!
//! Provides safe conversion between numeric types without panics

use eunomia::{FloatElement, RealField};

/// Safe conversion from `f64` to generic Atlas floating-point types.
pub trait SafeFromF64: RealField + FloatElement {
    /// Convert from `f64`.
    #[inline]
    fn from_f64_or(value: f64, _fallback: Self) -> Self {
        <Self as FloatElement>::from_f64(value)
    }

    /// Convert from `f64`.
    #[inline]
    fn from_f64_or_zero(value: f64) -> Self {
        <Self as FloatElement>::from_f64(value)
    }

    /// Convert from `f64`.
    #[inline]
    fn from_f64_or_one(value: f64) -> Self {
        <Self as FloatElement>::from_f64(value)
    }

    /// Convert from `f64`, returning `Result` for API consistency.
    ///
    /// # Errors
    /// This implementation is infallible for `FloatElement` implementors.
    #[inline]
    fn try_from_f64(value: f64) -> crate::error::Result<Self> {
        Ok(<Self as FloatElement>::from_f64(value))
    }
}

impl<T> SafeFromF64 for T where T: RealField + FloatElement {}

/// Safe conversion from `i32` to generic Atlas floating-point types.
pub trait SafeFromI32: RealField + FloatElement {
    /// Convert from `i32`.
    #[inline]
    fn from_i32_or_zero(value: i32) -> Self {
        <Self as FloatElement>::from_f64(f64::from(value))
    }

    /// Convert from `i32`, returning `Result` for API consistency.
    ///
    /// # Errors
    /// This implementation is infallible for `FloatElement` implementors.
    #[inline]
    fn try_from_i32(value: i32) -> crate::error::Result<Self> {
        Ok(<Self as FloatElement>::from_f64(f64::from(value)))
    }
}

impl<T> SafeFromI32 for T where T: RealField + FloatElement {}

/// Safe conversion from `usize` to generic Atlas floating-point types.
pub trait SafeFromUsize: RealField + FloatElement {
    /// Convert from `usize`.
    #[inline]
    fn from_usize_or_one(value: usize) -> Self {
        <Self as FloatElement>::from_f64(value as f64)
    }

    /// Convert from `usize`, returning `Result` for API consistency.
    ///
    /// # Errors
    /// This implementation is infallible for `FloatElement` implementors.
    #[inline]
    fn try_from_usize(value: usize) -> crate::error::Result<Self> {
        Ok(<Self as FloatElement>::from_f64(value as f64))
    }
}

impl<T> SafeFromUsize for T where T: RealField + FloatElement {}

#[cfg(test)]
mod tests {
    use super::{SafeFromF64, SafeFromI32, SafeFromUsize};

    #[test]
    fn safe_conversion_uses_eunomia_float_construction() {
        assert_eq!(<f64 as SafeFromF64>::from_f64_or_zero(1.25), 1.25);
        assert_eq!(<f64 as SafeFromF64>::from_f64_or_one(2.5), 2.5);
        assert_eq!(<f64 as SafeFromF64>::try_from_f64(-3.75).unwrap(), -3.75);
        assert_eq!(<f64 as SafeFromI32>::from_i32_or_zero(-7), -7.0);
        assert_eq!(<f64 as SafeFromI32>::try_from_i32(11).unwrap(), 11.0);
        assert_eq!(<f64 as SafeFromUsize>::from_usize_or_one(13), 13.0);
        assert_eq!(<f64 as SafeFromUsize>::try_from_usize(17).unwrap(), 17.0);
    }
}
