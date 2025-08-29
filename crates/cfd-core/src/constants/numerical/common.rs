//! Common numerical constants
//!
//! Frequently used numerical values to eliminate magic numbers

/// Zero value
pub const ZERO: f64 = 0.0;

/// Unity value
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

/// Seven
pub const SEVEN: f64 = 7.0;

/// Eight
pub const EIGHT: f64 = 8.0;

/// Nine
pub const NINE: f64 = 9.0;

/// Ten
pub const TEN: f64 = 10.0;

/// One half
pub const HALF: f64 = 0.5;

/// One third
pub const ONE_THIRD: f64 = 1.0 / 3.0;

/// Two thirds
pub const TWO_THIRDS: f64 = 2.0 / 3.0;

/// One quarter
pub const QUARTER: f64 = 0.25;

/// Three quarters
pub const THREE_QUARTERS: f64 = 0.75;

/// One sixth
pub const ONE_SIXTH: f64 = 1.0 / 6.0;

/// Five sixths
pub const FIVE_SIXTHS: f64 = 5.0 / 6.0;

/// One twelfth
pub const ONE_TWELFTH: f64 = 1.0 / 12.0;

/// Machine epsilon for f64
pub const EPSILON: f64 = f64::EPSILON;

/// Small value to prevent division by zero
pub const SMALL: f64 = 1e-10;

/// Large value for clamping
pub const LARGE: f64 = 1e10;
