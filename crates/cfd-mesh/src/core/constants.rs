//! Physical and geometric constants for millifluidic design.

use crate::core::scalar::Real;

/// π
pub const PI: Real = std::f64::consts::PI as Real;

/// 2π
pub const TAU: Real = std::f64::consts::TAU as Real;

/// π/2
pub const FRAC_PI_2: Real = std::f64::consts::FRAC_PI_2 as Real;

// ── Unit conversions (to meters) ──────────────────────────────

/// 1 mm in meters.
pub const MM: Real = 1e-3 as Real;

/// 1 μm in meters.
pub const UM: Real = 1e-6 as Real;

/// 1 cm in meters.
pub const CM: Real = 1e-2 as Real;

// ── Millifluidic defaults ─────────────────────────────────────

/// Default channel diameter for millifluidic devices (mm).
pub const DEFAULT_CHANNEL_DIAMETER_MM: Real = 1.0 as Real;

/// Default substrate height (mm).
pub const DEFAULT_SUBSTRATE_HEIGHT_MM: Real = 10.0 as Real;

/// Default wall thickness (mm).
pub const DEFAULT_WALL_THICKNESS_MM: Real = 2.0 as Real;

/// Minimum segment length before it is collapsed (mm).
pub const MIN_SEGMENT_LENGTH_MM: Real = 1e-3 as Real;

// ── Mesh quality defaults ─────────────────────────────────────

/// Minimum acceptable triangle quality score [0, 1].
pub const DEFAULT_MIN_QUALITY: Real = 0.3 as Real;

/// Maximum acceptable aspect ratio.
pub const DEFAULT_MAX_ASPECT_RATIO: Real = 10.0 as Real;

/// Minimum acceptable interior angle (degrees).
pub const DEFAULT_MIN_ANGLE_DEG: Real = 15.0 as Real;

/// Maximum acceptable interior angle (degrees).
pub const DEFAULT_MAX_ANGLE_DEG: Real = 150.0 as Real;

// ── Mesh quality constants needed by various modules ─────────

/// Minimum acceptable interior angle (radians) for quality checks.
pub const DEFAULT_MIN_ANGLE: Real = 15.0 * std::f64::consts::PI as Real / 180.0 as Real;

/// Maximum acceptable equiangle skewness [0, 1].
pub const DEFAULT_MAX_SKEWNESS: Real = 0.8 as Real;

/// Minimum acceptable edge-length ratio [0, 1].
pub const DEFAULT_MIN_EDGE_RATIO: Real = 0.1 as Real;

/// Default channel radius for millifluidic devices (m).
pub const DEFAULT_CHANNEL_RADIUS: Real = 0.5e-3 as Real;
