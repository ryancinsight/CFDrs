//! Database of common fluids with validated properties
//!
//! All fluid properties are sourced from:
//! - NIST Chemistry `WebBook`
//! - Perry's Chemical Engineers' Handbook
//! - CRC Handbook of Chemistry and Physics

use super::newtonian::{ConstantPropertyFluid, IdealGas};
use crate::error::Error;
use crate::physics::constants::physics::fluid as fluid_constants;
use eunomia::FloatElement;
use eunomia::RealField;

/// Create water at 20°C and 1 atm
///
/// Properties from NIST Chemistry WebBook
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn water_20c<T: RealField + FloatElement + Copy>() -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Water (20°C)".to_string(),
        <T as FloatElement>::from_f64(fluid_constants::WATER_DENSITY),
        <T as FloatElement>::from_f64(fluid_constants::WATER_VISCOSITY),
        <T as FloatElement>::from_f64(fluid_constants::WATER_SPECIFIC_HEAT),
        <T as FloatElement>::from_f64(fluid_constants::WATER_THERMAL_CONDUCTIVITY),
        <T as FloatElement>::from_f64(fluid_constants::WATER_SPEED_OF_SOUND),
    ))
}

/// Create air at 20°C and 1 atm
///
/// Properties from NIST Chemistry WebBook
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn air_20c<T: RealField + FloatElement + Copy>() -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Air (20°C, 1 atm)".to_string(),
        <T as FloatElement>::from_f64(fluid_constants::AIR_DENSITY),
        <T as FloatElement>::from_f64(fluid_constants::AIR_VISCOSITY),
        <T as FloatElement>::from_f64(fluid_constants::AIR_SPECIFIC_HEAT),
        <T as FloatElement>::from_f64(fluid_constants::AIR_THERMAL_CONDUCTIVITY),
        <T as FloatElement>::from_f64(fluid_constants::AIR_SPEED_OF_SOUND),
    ))
}

/// Create ideal air model for variable temperature/pressure
///
/// Uses Sutherland's law for viscosity variation
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn ideal_air<T: RealField + FloatElement + Copy>() -> Result<IdealGas<T>, Error> {
    Ok(IdealGas::new(
        "Air (Ideal Gas)".to_string(),
        <T as FloatElement>::from_f64(287.0),    // J/(kg·K)
        <T as FloatElement>::from_f64(1005.0),   // J/(kg·K)
        <T as FloatElement>::from_f64(1.716e-5), // Pa·s at 273K
        <T as FloatElement>::from_f64(273.15),   // K
        <T as FloatElement>::from_f64(110.4),    // K
    ))
}

/// Create engine oil SAE 30 at 40°C
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn engine_oil_sae30<T: RealField + FloatElement + Copy>(
) -> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Engine Oil SAE 30 (40°C)".to_string(),
        <T as FloatElement>::from_f64(870.0),  // kg/m³
        <T as FloatElement>::from_f64(0.1),    // Pa·s
        <T as FloatElement>::from_f64(2000.0), // J/(kg·K)
        <T as FloatElement>::from_f64(0.145),  // W/(m·K)
        <T as FloatElement>::from_f64(1740.0), // m/s (approximate)
    ))
}

/// Create glycerin at 20°C
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn glycerin_20c<T: RealField + FloatElement + Copy>() -> Result<ConstantPropertyFluid<T>, Error>
{
    Ok(ConstantPropertyFluid::new(
        "Glycerin (20°C)".to_string(),
        <T as FloatElement>::from_f64(1260.0), // kg/m³
        <T as FloatElement>::from_f64(1.49),   // Pa·s
        <T as FloatElement>::from_f64(2430.0), // J/(kg·K)
        <T as FloatElement>::from_f64(0.285),  // W/(m·K)
        <T as FloatElement>::from_f64(1920.0), // m/s
    ))
}

/// Create mercury at 20°C
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn mercury_20c<T: RealField + FloatElement + Copy>() -> Result<ConstantPropertyFluid<T>, Error>
{
    Ok(ConstantPropertyFluid::new(
        "Mercury (20°C)".to_string(),
        <T as FloatElement>::from_f64(13534.0),  // kg/m³
        <T as FloatElement>::from_f64(1.526e-3), // Pa·s
        <T as FloatElement>::from_f64(139.5),    // J/(kg·K)
        <T as FloatElement>::from_f64(8.3),      // W/(m·K)
        <T as FloatElement>::from_f64(1450.0),   // m/s
    ))
}
