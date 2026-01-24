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

// WENO9 flux reconstruction coefficients
/// WENO9 linear weights for optimal 9th-order accuracy
pub const WENO9_LINEAR_WEIGHTS: [f64; 5] = [
    1.0 / 256.0,
    36.0 / 256.0,
    126.0 / 256.0,
    84.0 / 256.0,
    9.0 / 256.0,
];

/// WENO9 stencil coefficients for the 5 candidate stencils (5 points each).
/// Each row corresponds to a stencil k=0..4.
/// Values are divided by 128 in the reconstruction, so these are the numerators.
/// Stencil 0: 35/128, -180/128, 378/128, -420/128, 315/128
/// Stencil 1: -5/128, 28/128, -70/128, 140/128, 35/128
/// ...
pub const WENO9_STENCIL_COEFFS: [[f64; 5]; 5] = [
    [35.0, -180.0, 378.0, -420.0, 315.0],
    [-5.0, 28.0, -70.0, 140.0, 35.0],
    [3.0, -20.0, 90.0, 60.0, -5.0],
    [-5.0, 60.0, 90.0, -20.0, 3.0],
    [35.0, 140.0, -70.0, 28.0, -5.0],
];

/// WENO9 common denominator for stencil coefficients
pub const WENO9_STENCIL_DENOM: f64 = 128.0;
