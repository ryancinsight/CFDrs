//! Constants for WENO schemes following SSOT principle

// WENO5 smoothness indicator coefficients
/// WENO5 smoothness indicator coefficient 13/12 from Jiang & Shu (1996)
pub const WENO5_BETA_COEFF_13_12: f64 = 13.0 / 12.0;
/// WENO5 smoothness indicator quarter coefficient
pub const WENO5_BETA_COEFF_QUARTER: f64 = 0.25;
/// WENO5 coefficient for value 2.0
pub const WENO5_BETA_COEFF_TWO: f64 = 2.0;
/// WENO5 coefficient for value 3.0  
pub const WENO5_BETA_COEFF_THREE: f64 = 3.0;
/// WENO5 coefficient for value 4.0
pub const WENO5_BETA_COEFF_FOUR: f64 = 4.0;

// WENO5 flux reconstruction coefficients
/// WENO5 flux reconstruction coefficient 1/3
pub const WENO5_FLUX_ONE_THIRD: f64 = 1.0 / 3.0;
/// WENO5 flux reconstruction coefficient 7/6
pub const WENO5_FLUX_SEVEN_SIXTH: f64 = 7.0 / 6.0;
/// WENO5 flux reconstruction coefficient 11/6
pub const WENO5_FLUX_ELEVEN_SIXTH: f64 = 11.0 / 6.0;
/// WENO5 flux reconstruction coefficient 1/6
pub const WENO5_FLUX_ONE_SIXTH: f64 = 1.0 / 6.0;
/// WENO5 flux reconstruction coefficient 5/6
pub const WENO5_FLUX_FIVE_SIXTH: f64 = 5.0 / 6.0;
