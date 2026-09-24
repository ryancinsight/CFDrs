//! Raw configuration constants: bounds and canonical defaults.
//!
//! These are the single definition site for every constant in the
//! [`super`] groups; the group `DEFAULT` aggregates reference them.

/// Minimum allowed wall clearance (mm)
pub const MIN_WALL_CLEARANCE: f64 = 0.01;
/// Maximum allowed wall clearance (mm)
pub const MAX_WALL_CLEARANCE: f64 = 10.0;
/// Default wall clearance (mm)
pub const DEFAULT_WALL_CLEARANCE: f64 = 0.5;

// Rendering and geometry generation constants
/// Default number of points for serpentine path generation
pub const DEFAULT_SERPENTINE_POINTS: usize = 200;
/// Minimum number of points for serpentine path generation
pub const MIN_SERPENTINE_POINTS: usize = 10;
/// Maximum number of points for serpentine path generation
pub const MAX_SERPENTINE_POINTS: usize = 1000;

/// Default number of points for optimization path generation
pub const DEFAULT_OPTIMIZATION_POINTS: usize = 50;
/// Minimum number of points for optimization path generation
pub const MIN_OPTIMIZATION_POINTS: usize = 10;
/// Maximum number of points for optimization path generation
pub const MAX_OPTIMIZATION_POINTS: usize = 200;

/// Default number of middle points for smooth straight channels
pub const DEFAULT_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 10;
/// Minimum number of middle points for smooth straight channels
pub const MIN_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 2;
/// Maximum number of middle points for smooth straight channels
pub const MAX_SMOOTH_STRAIGHT_MIDDLE_POINTS: usize = 50;

/// Default wave multiplier for smooth transitions (2pi for one complete wave)
pub const DEFAULT_TRANSITION_WAVE_MULTIPLIER: f64 = 2.0;
/// Minimum wave multiplier for smooth transitions
pub const MIN_TRANSITION_WAVE_MULTIPLIER: f64 = 0.5;
/// Maximum wave multiplier for smooth transitions
pub const MAX_TRANSITION_WAVE_MULTIPLIER: f64 = 10.0;

// Adaptive serpentine control constants
/// Default distance normalization factor for node proximity effects
pub const DEFAULT_NODE_DISTANCE_NORMALIZATION: f64 = 10.0;
/// Minimum distance normalization factor
pub const MIN_NODE_DISTANCE_NORMALIZATION: f64 = 1.0;
/// Maximum distance normalization factor
pub const MAX_NODE_DISTANCE_NORMALIZATION: f64 = 50.0;

/// Default plateau width factor for horizontal channels (fraction of channel length)
pub const DEFAULT_PLATEAU_WIDTH_FACTOR: f64 = 0.4;
/// Minimum plateau width factor
pub const MIN_PLATEAU_WIDTH_FACTOR: f64 = 0.1;
/// Maximum plateau width factor
pub const MAX_PLATEAU_WIDTH_FACTOR: f64 = 0.8;

/// Default horizontal ratio threshold for middle section detection
/// This is a cosine-based ratio (`|dx| / distance`), bounded to [0, 1].
/// A value of 0.75 corresponds to channels within ±41.4° of horizontal.
pub const DEFAULT_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.75;
/// Minimum horizontal ratio threshold
pub const MIN_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.5;
/// Maximum horizontal ratio threshold
pub const MAX_HORIZONTAL_RATIO_THRESHOLD: f64 = 0.95;

/// Default middle section amplitude factor
pub const DEFAULT_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 0.7;
/// Minimum middle section amplitude factor
pub const MIN_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 0.1;
/// Maximum middle section amplitude factor
pub const MAX_MIDDLE_SECTION_AMPLITUDE_FACTOR: f64 = 1.0;

/// Default plateau amplitude factor
pub const DEFAULT_PLATEAU_AMPLITUDE_FACTOR: f64 = 0.8;
/// Minimum plateau amplitude factor
pub const MIN_PLATEAU_AMPLITUDE_FACTOR: f64 = 0.5;
/// Maximum plateau amplitude factor
pub const MAX_PLATEAU_AMPLITUDE_FACTOR: f64 = 1.0;

/// Minimum allowed channel width (mm)
pub const MIN_CHANNEL_WIDTH: f64 = 0.01;
/// Maximum allowed channel width (mm)
pub const MAX_CHANNEL_WIDTH: f64 = 50.0;
/// Default channel width (mm)
pub const DEFAULT_CHANNEL_WIDTH: f64 = 1.0;

/// Minimum allowed channel height (mm)
pub const MIN_CHANNEL_HEIGHT: f64 = 0.01;
/// Maximum allowed channel height (mm)
pub const MAX_CHANNEL_HEIGHT: f64 = 25.0;
/// Default channel height (mm)
pub const DEFAULT_CHANNEL_HEIGHT: f64 = 0.5;

/// Minimum fill factor for serpentine channels
pub const MIN_FILL_FACTOR: f64 = 0.1;
/// Maximum fill factor for serpentine channels
pub const MAX_FILL_FACTOR: f64 = 0.95;
/// Default fill factor for serpentine channels
pub const DEFAULT_FILL_FACTOR: f64 = 0.8;

/// Minimum wavelength factor for serpentine channels
pub const MIN_WAVELENGTH_FACTOR: f64 = 1.0;
/// Maximum wavelength factor for serpentine channels
pub const MAX_WAVELENGTH_FACTOR: f64 = 10.0;
/// Default wavelength factor for serpentine channels
pub const DEFAULT_WAVELENGTH_FACTOR: f64 = 6.0;

/// Minimum Gaussian width factor for serpentine channels
pub const MIN_GAUSSIAN_WIDTH_FACTOR: f64 = 2.0;
/// Maximum Gaussian width factor for serpentine channels
pub const MAX_GAUSSIAN_WIDTH_FACTOR: f64 = 20.0;
/// Default Gaussian width factor for serpentine channels
pub const DEFAULT_GAUSSIAN_WIDTH_FACTOR: f64 = 6.0;

/// Minimum wave density factor for serpentine channels
pub const MIN_WAVE_DENSITY_FACTOR: f64 = 0.5;
/// Maximum wave density factor for serpentine channels
pub const MAX_WAVE_DENSITY_FACTOR: f64 = 10.0;
/// Default wave density factor for serpentine channels
pub const DEFAULT_WAVE_DENSITY_FACTOR: f64 = 1.5;

/// Minimum curvature factor for arc channels
pub const MIN_CURVATURE_FACTOR: f64 = 0.0;
/// Maximum curvature factor for arc channels
pub const MAX_CURVATURE_FACTOR: f64 = 2.0;
/// Default curvature factor for arc channels
pub const DEFAULT_CURVATURE_FACTOR: f64 = 0.3;

/// Minimum smoothness for arc channels
pub const MIN_SMOOTHNESS: usize = 3;
/// Maximum smoothness for arc channels
pub const MAX_SMOOTHNESS: usize = 1000;
/// Default smoothness for arc channels
pub const DEFAULT_SMOOTHNESS: usize = 20;

/// Minimum separation distance between arc channels (mm)
pub const MIN_SEPARATION_DISTANCE: f64 = 0.1;
/// Maximum separation distance between arc channels (mm)
pub const MAX_SEPARATION_DISTANCE: f64 = 10.0;
/// Default minimum separation distance between arc channels (mm)
pub const DEFAULT_MIN_SEPARATION_DISTANCE: f64 = 1.0;

/// Minimum curvature reduction factor for collision prevention
pub const MIN_CURVATURE_REDUCTION: f64 = 0.1;
/// Maximum curvature reduction factor for collision prevention
pub const MAX_CURVATURE_REDUCTION_LIMIT: f64 = 1.0;
/// Default maximum curvature reduction factor
pub const DEFAULT_MAX_CURVATURE_REDUCTION: f64 = 0.5;

/// Strategy thresholds for smart channel type selection
pub mod strategy_thresholds {
    /// Threshold for long horizontal channels (fraction of box width)
    pub const LONG_HORIZONTAL_THRESHOLD: f64 = 0.6;
    /// Threshold for minimum arc length (fraction of box width)
    pub const MIN_ARC_LENGTH_THRESHOLD: f64 = 0.3;
    /// Threshold for horizontal vs angled channel detection
    pub const HORIZONTAL_ANGLE_THRESHOLD: f64 = 0.5;
    /// Threshold for angled channel detection (slope)
    pub const ANGLED_CHANNEL_SLOPE_THRESHOLD: f64 = 0.1;
    /// Default middle zone fraction for mixed by position
    pub const DEFAULT_MIDDLE_ZONE_FRACTION: f64 = 0.4;
}

// Adaptive collision constants.
/// Default maximum collision adjustment factor.
pub const DEFAULT_MAX_ADJUSTMENT_FACTOR: f64 = 0.5;
/// Default minimum distance between channels.
pub const DEFAULT_MIN_CHANNEL_DISTANCE: f64 = 0.5;
/// Default safety margin applied to collision checks.
pub const DEFAULT_SAFETY_MARGIN_FACTOR: f64 = 1.1;
/// Default maximum reduction factor for collision recovery.
pub const DEFAULT_MAX_REDUCTION_FACTOR: f64 = 0.8;
/// Default sensitivity used by collision detection.
pub const DEFAULT_DETECTION_SENSITIVITY: f64 = 0.1;
/// Default divisor for proximity normalization.
pub const DEFAULT_PROXIMITY_DIVISOR: f64 = 10.0;
/// Default minimum proximity factor.
pub const DEFAULT_MIN_PROXIMITY_FACTOR: f64 = 0.2;
/// Default maximum proximity factor.
pub const DEFAULT_MAX_PROXIMITY_FACTOR: f64 = 2.0;
/// Default divisor for branch adjustment scaling.
pub const DEFAULT_BRANCH_ADJUSTMENT_DIVISOR: f64 = 5.0;
/// Default maximum sensitivity multiplier.
pub const DEFAULT_MAX_SENSITIVITY_MULTIPLIER: f64 = 3.0;
/// Default multiplier for long-channel reduction.
pub const DEFAULT_LONG_CHANNEL_REDUCTION_MULTIPLIER: f64 = 0.9;
/// Default upper bound for collision reduction.
pub const DEFAULT_MAX_REDUCTION_LIMIT: f64 = 0.5;

// --- Strategy selection: frustum thresholds ---
/// Minimum length threshold (fraction of box width) for frustum selection
pub const FRUSTUM_MIN_LENGTH_THRESHOLD: f64 = 0.3;
/// Maximum length threshold (fraction of box width) for frustum selection
pub const FRUSTUM_MAX_LENGTH_THRESHOLD: f64 = 0.7;
/// Maximum angle threshold (dy/dx ratio) for frustum selection
pub const FRUSTUM_ANGLE_THRESHOLD: f64 = 0.5;

// --- Wave generation ---
/// Sharpness factor for square wave generation using tanh
pub const SQUARE_WAVE_SHARPNESS: f64 = 5.0;
/// Factor for smooth transition zones at wave endpoints
pub const TRANSITION_ZONE_FACTOR: f64 = 0.1;
/// Threshold for smooth endpoint transition start
pub const SMOOTH_ENDPOINT_START_THRESHOLD: f64 = 0.1;
/// Threshold for smooth endpoint transition end
pub const SMOOTH_ENDPOINT_END_THRESHOLD: f64 = 0.9;
/// Default length factor for smooth transitions
pub const DEFAULT_TRANSITION_LENGTH_FACTOR: f64 = 0.15;
/// Default amplitude factor for smooth transitions
pub const DEFAULT_TRANSITION_AMPLITUDE_FACTOR: f64 = 0.3;
/// Default number of points for transition smoothing
pub const DEFAULT_TRANSITION_SMOOTHNESS: usize = 20;
/// Scaling factor for neighbor avoidance calculations
pub const NEIGHBOR_AVOIDANCE_SCALING_FACTOR: f64 = 0.8;
/// Tolerance for geometric distance comparisons
pub const GEOMETRIC_TOLERANCE: f64 = 1e-6;

// --- Geometry generation ---
/// Multiplier for channel width to detect short channels
pub const SHORT_CHANNEL_WIDTH_MULTIPLIER: f64 = 2.0;
/// Floor applied to the *adaptive* curvature factor.
///
/// Distinct from [`MIN_CURVATURE_FACTOR`], which is the lower bound
/// accepted when *validating* a configured curvature factor. The two have
/// always held different values (0.1 here, 0.0 there) and govern different
/// decisions, so they are named apart rather than merged.
pub const MIN_ADAPTIVE_CURVATURE_FACTOR: f64 = 0.1;

// --- Optimization ---
/// Exponent for branch factor scaling
pub const BRANCH_FACTOR_EXPONENT: f64 = 0.75;
/// Maximum number of iterations for optimization algorithms
pub const MAX_OPTIMIZATION_ITERATIONS: usize = 100;
/// Tolerance for optimization convergence detection
pub const CONVERGENCE_TOLERANCE: f64 = 1e-6;
/// Wavelength factors for fast optimization
pub const FAST_WAVELENGTH_FACTORS: &[f64] = &[1.0, 2.0, 3.0, 4.0];
/// Wave density factors for fast optimization
pub const FAST_WAVE_DENSITY_FACTORS: &[f64] = &[1.0, 2.0, 3.0];
/// Fill factors for fast optimization
pub const FAST_FILL_FACTORS: &[f64] = &[0.7, 0.8, 0.9];

// --- Visualization ---
/// Default margin for chart rendering
pub const DEFAULT_CHART_MARGIN: u32 = 20;
/// Default right margin for chart rendering
pub const DEFAULT_CHART_RIGHT_MARGIN: u32 = 150;
/// Default label area size for the x-axis
pub const DEFAULT_X_LABEL_AREA_SIZE: u32 = 30;
/// Default label area size for the y-axis
pub const DEFAULT_Y_LABEL_AREA_SIZE: u32 = 30;
