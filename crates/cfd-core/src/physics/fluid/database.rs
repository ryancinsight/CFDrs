//! Database of common fluids with validated properties
//!
//! All fluid properties are sourced from:
//! - NIST Chemistry `WebBook`
//! - Perry's Chemical Engineers' Handbook
//! - CRC Handbook of Chemistry and Physics

use super::newtonian::{ConstantPropertyFluid, IdealGas};
use crate::error::Error;
use crate::physics::constants::physics::fluid as fluid_constants;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity, Velocity,
};
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
        MassDensity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::WATER_DENSITY,
        )),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::WATER_VISCOSITY,
        )),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::WATER_SPECIFIC_HEAT,
        )),
        ThermalConductivity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::WATER_THERMAL_CONDUCTIVITY,
        )),
        Velocity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::WATER_SPEED_OF_SOUND,
        )),
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
        MassDensity::from_base(<T as FloatElement>::from_f64(fluid_constants::AIR_DENSITY)),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::AIR_VISCOSITY,
        )),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::AIR_SPECIFIC_HEAT,
        )),
        ThermalConductivity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::AIR_THERMAL_CONDUCTIVITY,
        )),
        Velocity::from_base(<T as FloatElement>::from_f64(
            fluid_constants::AIR_SPEED_OF_SOUND,
        )),
    ))
}

/// Create ideal air model for variable temperature/pressure
///
/// Uses Sutherland's law for viscosity variation
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn ideal_air<T: RealField + FloatElement + Copy>() -> Result<IdealGas<T>, Error> {
    use aequitas::systems::si::quantities::{
        DynamicViscosity, SpecificHeatCapacity, TemperatureDifference, ThermodynamicTemperature,
    };

    Ok(IdealGas::new(
        "Air (Ideal Gas)".to_string(),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(287.0)),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(1005.0)),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(1.716e-5)),
        ThermodynamicTemperature::from_base(<T as FloatElement>::from_f64(273.15)),
        TemperatureDifference::from_base(<T as FloatElement>::from_f64(110.4)),
    ))
}

/// Create engine oil SAE 30 at 40°C
///
/// # Errors
/// This constructor currently has no input-dependent failure path.
pub fn engine_oil_sae30<T: RealField + FloatElement + Copy>()
-> Result<ConstantPropertyFluid<T>, Error> {
    Ok(ConstantPropertyFluid::new(
        "Engine Oil SAE 30 (40°C)".to_string(),
        MassDensity::from_base(<T as FloatElement>::from_f64(870.0)),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(0.1)),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(2000.0)),
        ThermalConductivity::from_base(<T as FloatElement>::from_f64(0.145)),
        Velocity::from_base(<T as FloatElement>::from_f64(1740.0)),
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
        MassDensity::from_base(<T as FloatElement>::from_f64(1260.0)),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(1.49)),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(2430.0)),
        ThermalConductivity::from_base(<T as FloatElement>::from_f64(0.285)),
        Velocity::from_base(<T as FloatElement>::from_f64(1920.0)),
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
        MassDensity::from_base(<T as FloatElement>::from_f64(13534.0)),
        DynamicViscosity::from_base(<T as FloatElement>::from_f64(1.526e-3)),
        SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(139.5)),
        ThermalConductivity::from_base(<T as FloatElement>::from_f64(8.3)),
        Velocity::from_base(<T as FloatElement>::from_f64(1450.0)),
    ))
}
