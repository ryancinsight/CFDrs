//! Database of common fluids with validated properties
//!
//! All fluid properties are sourced from:
//! - NIST Chemistry `WebBook`
//! - Perry's Chemical Engineers' Handbook
//! - CRC Handbook of Chemistry and Physics

use super::newtonian::{ConstantPropertyFluid, IdealGas};
use crate::physics::constants::physics::fluid as fluid_constants;
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Create water at 20°C and 1 atm
///
/// Properties from NIST Chemistry WebBook
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn water_20c<T: RealField + FromPrimitive + Copy>() -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Water (20°C)".to_string(),
        T::from_f64(fluid_constants::WATER_DENSITY)
            .ok_or_else(|| Error::ConversionError("Failed to convert water density".to_string()))?,
        T::from_f64(fluid_constants::WATER_VISCOSITY).ok_or_else(|| {
            Error::ConversionError("Failed to convert water viscosity".to_string())
        })?,
        T::from_f64(fluid_constants::WATER_SPECIFIC_HEAT).ok_or_else(|| {
            Error::ConversionError("Failed to convert water specific heat".to_string())
        })?,
        T::from_f64(fluid_constants::WATER_THERMAL_CONDUCTIVITY).ok_or_else(|| {
            Error::ConversionError("Failed to convert water thermal conductivity".to_string())
        })?,
        T::from_f64(fluid_constants::WATER_SPEED_OF_SOUND).ok_or_else(|| {
            Error::ConversionError("Failed to convert water speed of sound".to_string())
        })?,
    ))
}

/// Create air at 20°C and 1 atm
///
/// Properties from NIST Chemistry WebBook
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn air_20c<T: RealField + FromPrimitive + Copy>() -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Air (20°C, 1 atm)".to_string(),
        T::from_f64(fluid_constants::AIR_DENSITY)
            .ok_or_else(|| Error::ConversionError("Failed to convert air density".to_string()))?,
        T::from_f64(fluid_constants::AIR_VISCOSITY)
            .ok_or_else(|| Error::ConversionError("Failed to convert air viscosity".to_string()))?,
        T::from_f64(fluid_constants::AIR_SPECIFIC_HEAT).ok_or_else(|| {
            Error::ConversionError("Failed to convert air specific heat".to_string())
        })?,
        T::from_f64(fluid_constants::AIR_THERMAL_CONDUCTIVITY).ok_or_else(|| {
            Error::ConversionError("Failed to convert air thermal conductivity".to_string())
        })?,
        T::from_f64(fluid_constants::AIR_SPEED_OF_SOUND).ok_or_else(|| {
            Error::ConversionError("Failed to convert air speed of sound".to_string())
        })?,
    ))
}

/// Create ideal air model for variable temperature/pressure
///
/// Uses Sutherland's law for viscosity variation
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn ideal_air<T: RealField + FromPrimitive + Copy>() -> Result<IdealGas<T>, Error> {
    Ok(IdealGas::new(
        "Air (Ideal Gas)".to_string(),
        T::from_f64(287.0)
            .ok_or_else(|| Error::ConversionError("Failed to convert gas constant".to_string()))?, // J/(kg·K)
        T::from_f64(1005.0)
            .ok_or_else(|| Error::ConversionError("Failed to convert cp".to_string()))?, // J/(kg·K)
        T::from_f64(1.716e-5).ok_or_else(|| {
            Error::ConversionError("Failed to convert reference viscosity".to_string())
        })?, // Pa·s at 273K
        T::from_f64(273.15).ok_or_else(|| {
            Error::ConversionError("Failed to convert reference temperature".to_string())
        })?, // K
        T::from_f64(110.4).ok_or_else(|| {
            Error::ConversionError("Failed to convert Sutherland constant".to_string())
        })?, // K
    ))
}

/// Create engine oil SAE 30 at 40°C
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn engine_oil_sae30<T: RealField + FromPrimitive + Copy>(
) -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Engine Oil SAE 30 (40°C)".to_string(),
        T::from_f64(870.0)
            .ok_or_else(|| Error::ConversionError("Failed to convert oil density".to_string()))?, // kg/m³
        T::from_f64(0.1)
            .ok_or_else(|| Error::ConversionError("Failed to convert oil viscosity".to_string()))?, // Pa·s
        T::from_f64(2000.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert oil specific heat".to_string())
        })?, // J/(kg·K)
        T::from_f64(0.145).ok_or_else(|| {
            Error::ConversionError("Failed to convert oil thermal conductivity".to_string())
        })?, // W/(m·K)
        T::from_f64(1740.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert oil speed of sound".to_string())
        })?, // m/s (approximate)
    ))
}

/// Create glycerin at 20°C
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn glycerin_20c<T: RealField + FromPrimitive + Copy>() -> Result<ConstantPropertyFluid<T>, Error>
{
    Ok(ConstantPropertyFluid::new(
        "Glycerin (20°C)".to_string(),
        T::from_f64(1260.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert glycerin density".to_string())
        })?, // kg/m³
        T::from_f64(1.49).ok_or_else(|| {
            Error::ConversionError("Failed to convert glycerin viscosity".to_string())
        })?, // Pa·s
        T::from_f64(2430.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert glycerin specific heat".to_string())
        })?, // J/(kg·K)
        T::from_f64(0.285).ok_or_else(|| {
            Error::ConversionError("Failed to convert glycerin thermal conductivity".to_string())
        })?, // W/(m·K)
        T::from_f64(1920.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert glycerin speed of sound".to_string())
        })?, // m/s
    ))
}

/// Create mercury at 20°C
///
/// # Errors
/// Returns an error if numeric conversion from f64 fails for the target type T
pub fn mercury_20c<T: RealField + FromPrimitive + Copy>() -> Result<ConstantPropertyFluid<T>, Error>
{
    Ok(ConstantPropertyFluid::new(
        "Mercury (20°C)".to_string(),
        T::from_f64(13534.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert mercury density".to_string())
        })?, // kg/m³
        T::from_f64(1.526e-3).ok_or_else(|| {
            Error::ConversionError("Failed to convert mercury viscosity".to_string())
        })?, // Pa·s
        T::from_f64(139.5).ok_or_else(|| {
            Error::ConversionError("Failed to convert mercury specific heat".to_string())
        })?, // J/(kg·K)
        T::from_f64(8.3).ok_or_else(|| {
            Error::ConversionError("Failed to convert mercury thermal conductivity".to_string())
        })?, // W/(m·K)
        T::from_f64(1450.0).ok_or_else(|| {
            Error::ConversionError("Failed to convert mercury speed of sound".to_string())
        })?, // m/s
    ))
}
