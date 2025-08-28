//! Constants for WENO schemes following SSOT principle

// WENO5 smoothness indicator coefficients
pub const WENO5_BETA_COEFF_13_12: f64 = 13.0 / 12.0;
pub const WENO5_BETA_COEFF_QUARTER: f64 = 0.25;
pub const WENO5_BETA_COEFF_TWO: f64 = 2.0;
pub const WENO5_BETA_COEFF_THREE: f64 = 3.0;
pub const WENO5_BETA_COEFF_FOUR: f64 = 4.0;

// WENO5 flux reconstruction coefficients
pub const WENO5_FLUX_ONE_THIRD: f64 = 1.0 / 3.0;
pub const WENO5_FLUX_SEVEN_SIXTH: f64 = 7.0 / 6.0;
pub const WENO5_FLUX_ELEVEN_SIXTH: f64 = 11.0 / 6.0;
pub const WENO5_FLUX_ONE_SIXTH: f64 = 1.0 / 6.0;
pub const WENO5_FLUX_FIVE_SIXTH: f64 = 5.0 / 6.0;

// Numerical constants
pub const TWO: f64 = 2.0;
pub const THREE: f64 = 3.0;
pub const FOUR: f64 = 4.0;
pub const FIVE: f64 = 5.0;
pub const SIX: f64 = 6.0;
pub const SEVEN: f64 = 7.0;
pub const ELEVEN: f64 = 11.0;
pub const THIRTEEN: f64 = 13.0;
