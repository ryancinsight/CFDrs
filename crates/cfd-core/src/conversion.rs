//! Safe numeric conversion utilities
//!
//! Provides safe conversion functions that return proper errors instead of panicking.

use crate::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;
/// Safe conversion from f64 to generic type T
pub fn from_f64<T: RealField + FromPrimitive>(value: f64) -> Result<T> {
    T::from_f64(value).ok_or_else(|| {
        Error::InvalidInput(format!("Cannot convert f64 value {} to target type", value))
    })
}
/// Safe conversion from usize to generic type T  
pub fn from_usize<T: RealField + FromPrimitive>(value: usize) -> Result<T> {
    T::from_usize(value).ok_or_else(|| {
        Error::InvalidInput(format!("Cannot convert usize value {} to target type", value))
/// Common numerical constants with safe conversion
pub struct NumericConstants;
impl NumericConstants {
    /// Get 0.5 in type T
    pub fn half<T: RealField + FromPrimitive>() -> Result<T> {
        from_f64(0.5)
    }
    
    /// Get 2.0 in type T
    pub fn two<T: RealField + FromPrimitive>() -> Result<T> {
        from_f64(2.0)
    /// Get 4.0 in type T
    pub fn four<T: RealField + FromPrimitive>() -> Result<T> {
        from_f64(4.0)
}
}
}
}
}
