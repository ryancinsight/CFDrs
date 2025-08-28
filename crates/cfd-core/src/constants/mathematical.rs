//! Mathematical constants for CFD calculations
//!
//! Provides commonly used mathematical constants to avoid magic numbers

use std::f64::consts;

/// Pi constant
pub const PI: f64 = consts::PI;

/// Two pi (full circle in radians)
pub const TWO_PI: f64 = 2.0 * consts::PI;

/// Pi over two (quarter circle)
pub const PI_OVER_TWO: f64 = consts::PI / 2.0;

/// Pi over four (eighth circle)
pub const PI_OVER_FOUR: f64 = consts::PI / 4.0;

/// Euler's number
pub const E: f64 = consts::E;

/// Square root of 2
pub const SQRT_TWO: f64 = consts::SQRT_2;

/// Natural logarithm of 2
pub const LN_TWO: f64 = consts::LN_2;

/// Common numeric constants
pub mod numeric {
    /// Zero
    pub const ZERO: f64 = 0.0;

    /// One
    pub const ONE: f64 = 1.0;

    /// Two
    pub const TWO: f64 = 2.0;

    /// Three
    pub const THREE: f64 = 3.0;

    /// Four  
    pub const FOUR: f64 = 4.0;

    /// Five
    pub const FIVE: f64 = 5.0;

    /// Six
    pub const SIX: f64 = 6.0;

    /// Eight
    pub const EIGHT: f64 = 8.0;

    /// Ten
    pub const TEN: f64 = 10.0;

    /// One half
    pub const ONE_HALF: f64 = 0.5;

    /// One third
    pub const ONE_THIRD: f64 = 1.0 / 3.0;

    /// Two thirds
    pub const TWO_THIRDS: f64 = 2.0 / 3.0;

    /// One quarter
    pub const ONE_QUARTER: f64 = 0.25;

    /// Three quarters
    pub const THREE_QUARTERS: f64 = 0.75;
}
