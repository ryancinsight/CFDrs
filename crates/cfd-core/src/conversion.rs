//! Safe numeric conversion utilities
//!
//! Provides safe conversion between numeric types without panics

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Safe conversion from f64 to generic `RealField` types
pub trait SafeFromF64: RealField + FromPrimitive {
    /// Convert from f64 with fallback value
    #[inline]
    fn from_f64_or(value: f64, fallback: Self) -> Self {
        Self::from_f64(value).unwrap_or(fallback)
    }

    /// Convert from f64 with zero fallback
    #[inline]
    fn from_f64_or_zero(value: f64) -> Self {
        Self::from_f64(value).unwrap_or_else(Self::zero)
    }

    /// Convert from f64 with one fallback
    #[inline]
    fn from_f64_or_one(value: f64) -> Self {
        Self::from_f64(value).unwrap_or_else(Self::one)
    }

    /// Try to convert from f64, returning Result
    /// 
    /// # Errors
    /// Returns `ConversionError` if the f64 value cannot be represented in the target type
    #[inline]
    fn try_from_f64(value: f64) -> crate::error::Result<Self> {
        Self::from_f64(value).ok_or_else(|| {
            crate::error::Error::ConversionError(format!(
                "Failed to convert f64 value {value} to target type"
            ))
        })
    }
}

// Implement for all types that satisfy the bounds
impl<T> SafeFromF64 for T where T: RealField + FromPrimitive {}

/// Safe conversion from i32 to generic `RealField` types  
pub trait SafeFromI32: RealField + FromPrimitive {
    /// Convert from i32 with zero fallback
    #[inline]
    fn from_i32_or_zero(value: i32) -> Self {
        Self::from_i32(value).unwrap_or_else(Self::zero)
    }

    /// Try to convert from i32, returning Result
    /// 
    /// # Errors
    /// Returns `ConversionError` if the i32 value cannot be represented in the target type
    #[inline]
    fn try_from_i32(value: i32) -> crate::error::Result<Self> {
        Self::from_i32(value).ok_or_else(|| {
            crate::error::Error::ConversionError(format!(
                "Failed to convert i32 value {value} to target type"
            ))
        })
    }
}

impl<T> SafeFromI32 for T where T: RealField + FromPrimitive {}

/// Safe conversion from usize to generic `RealField` types
pub trait SafeFromUsize: RealField + FromPrimitive {
    /// Convert from usize with one fallback
    #[inline]
    fn from_usize_or_one(value: usize) -> Self {
        Self::from_usize(value).unwrap_or_else(Self::one)
    }

    /// Try to convert from usize, returning Result
    /// 
    /// # Errors
    /// Returns `ConversionError` if the usize value cannot be represented in the target type
    #[inline]
    fn try_from_usize(value: usize) -> crate::error::Result<Self> {
        Self::from_usize(value).ok_or_else(|| {
            crate::error::Error::ConversionError(format!(
                "Failed to convert usize value {value} to target type"
            ))
        })
    }
}

impl<T> SafeFromUsize for T where T: RealField + FromPrimitive {}
