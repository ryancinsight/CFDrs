//! Physical and geometric constants for millifluidic design.

use crate::domain::core::scalar::Real;

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

// ── CSG / GWN numerical tolerances (SSOT) ────────────────────────────────────

/// GWN solid-angle denominator guard.
///
/// The van Oosterom–Strackee solid-angle formula uses `atan2(num, den)`.
/// When both `|num|` and `|den|` are below this threshold the face contributes
/// a near-zero solid angle and is skipped to avoid `atan2(0, 0) = NaN`.
///
/// For `f32` meshes use [`GWN_DENOMINATOR_GUARD_F32`] instead; this constant
/// is only safe for `f64` arithmetic.
pub const GWN_DENOMINATOR_GUARD: Real = 1e-30;

/// Per-vertex squared distance guard for the GWN near-vertex skip.
///
/// # Theorem — Guard Safety
///
/// `f64::MIN_POSITIVE` (≈ 2.2 × 10⁻³⁰⁸) is the smallest positive normal
/// `f64`.  Any `norm_squared()` below this value indicates that the query
/// point is within sub-ULP distance of a mesh vertex — geometrically
/// impossible for physical meshes.  Skipping such faces prevents
/// `atan2(0, 0) → NaN` without affecting any reachable query point.
///
/// For `T = f32`, the generic GWN path uses `T::min_positive_value()` at
/// runtime (≈ 1.2 × 10⁻³⁸), which is strictly greater than 0.0 in f32.
/// This constant is intentionally `f64`-only; do **not** cast it to `f32`.
pub const GWN_NEAR_VERTEX_GUARD_SQ: Real = 1e-40;

/// GWN threshold: `|wn| > GWN_INSIDE_THRESHOLD` → query is strictly inside.
pub const GWN_INSIDE_THRESHOLD: Real = 0.65;

/// GWN threshold: `|wn| < GWN_OUTSIDE_THRESHOLD` → query is strictly outside.
///
/// The band `[GWN_OUTSIDE_THRESHOLD, GWN_INSIDE_THRESHOLD]` triggers the
/// tiebreaker predicates in `classify_fragment`.
pub const GWN_OUTSIDE_THRESHOLD: Real = 0.35;

/// Sliver face exclusion ratio for Phase 4 fragment classification.
///
/// A fragment is considered a numerically degenerate sliver and excluded when
/// `area_sq < SLIVER_AREA_RATIO_SQ * max_edge_sq`.
///
/// # Theorem — Scale-Correct Threshold
///
/// `sqrt(1e-14)` = 1e-7.  A fragment is skipped only when its altitude-to-
/// edge ratio is below 1e-7.  For millifluidic meshes the minimum physically
/// meaningful ratio is ≈ 5e-4 (50 µm altitude on a 4 mm edge), safely above
/// the threshold.  Numerically degenerate slivers produced by near-parallel
/// face intersections have ratios of ~10⁻¹⁰ – 10⁻¹⁵, correctly below. ∎
///
/// **Previous value `1e-10`** (altitude ratio ~3 × 10⁻⁵) incorrectly skipped
/// valid 80:1 aspect-ratio millifluidic faces (50 µm / 4 mm edge).
pub const SLIVER_AREA_RATIO_SQ: Real = 1e-14;

/// CDT co-refinement weld tolerance squared (metres²).
///
/// A snap endpoint is classified as lying on an edge when its 3-D distance
/// to the edge's projection point is less than `2 * sqrt(COREFINE_WELD_TOL_SQ)`.
pub const COREFINE_WELD_TOL_SQ: Real = 1e-6;

/// CDT co-refinement edge-endpoint exclusion margin.
///
/// Snap endpoints within this normalised parameter distance of an edge corner
/// are treated as corner snaps, not interior edge Steiner insertions.
pub const COREFINE_EDGE_EPS: Real = 1e-6;

/// Seam propagation collinearity tolerance squared.
///
/// A point P is on edge [Va, Vb] if
/// `|cross(Vb − Va, P − Va)|² < SEAM_COLLINEAR_TOL_SQ × |Vb − Va|²`.
pub const SEAM_COLLINEAR_TOL_SQ: Real = 1e-6;

/// Maximum Steiner vertices per face during CDT co-refinement.
///
/// When the total count (edge Steiners + interior Steiners) exceeds this
/// bound, `corefine_face` falls back to `midpoint_subdivide` to prevent
/// O(s²) CDT blowup from complex multi-branch junction geometries.
pub const MAX_STEINER_PER_FACE: usize = 32768;
