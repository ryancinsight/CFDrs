//! Safe numeric conversions for CFD computations
//!
//! This module provides safe conversions between numeric types that properly
//! handle errors instead of silently falling back to incorrect values like zero.

use crate::error::{Error, NumericalErrorKind, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;
/// Safe conversion from f64 to generic RealField type
///
/// # Errors
/// Returns error if conversion fails instead of using dangerous fallback values
pub fn from_f64<T: RealField + FromPrimitive>(value: f64) -> Result<T> {
    T::from_f64(value).ok_or_else(|| {
        Error::Numerical(NumericalErrorKind::ConversionFailed {
            from_type: "f64",
            to_type: std::any::type_name::<T>(),
            value: value.to_string(),
        })
    })
}
/// Safe conversion from usize to generic RealField type
pub fn from_usize<T: RealField + FromPrimitive>(value: usize) -> Result<T> {
    T::from_usize(value).ok_or_else(|| {
            from_type: "usize",
/// Safe conversion from i32 to generic RealField type  
pub fn from_i32<T: RealField + FromPrimitive>(value: i32) -> Result<T> {
    T::from_i32(value).ok_or_else(|| {
            from_type: "i32",
/// Get mathematical constant PI
pub fn pi<T: RealField + FromPrimitive>() -> Result<T> {
    from_f64(std::f64::consts::PI)
/// Get mathematical constant E (Euler's number)
pub fn e<T: RealField + FromPrimitive>() -> Result<T> {
    from_f64(std::f64::consts::E)
/// Get mathematical constant TAU (2*PI)
pub fn tau<T: RealField + FromPrimitive>() -> Result<T> {
    from_f64(std::f64::consts::TAU)
}
}
}
}
}
}
}
