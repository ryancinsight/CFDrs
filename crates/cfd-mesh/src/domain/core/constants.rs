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

/// Solid-angle clip margin for bounded GWN evaluation.
///
/// Each per-triangle solid angle `Ω` is clamped to `|Ω| ≤ 2π − δ` where
/// `δ = GWN_SOLID_ANGLE_CLIP`.  This prevents a single near-coincident
/// triangle from contributing a full half-winding (±0.5) to the total,
/// reducing numerical jitter for near-surface query points.
///
/// # Theorem — Clip Safety
///
/// For a query at distance `d` from the nearest mesh face, the dominant
/// face subtends `Ω ≈ 2π − O(d²/A)` where A is face area.  The clip
/// fires only when `O(d²/A) < δ = 1e-6`, i.e. `d < √(A × 1e-6)`.
/// For A = 1 mm² this is d < 1 nm — safely below any physical resolution.
/// Interior/exterior queries never trigger the clip. ∎
pub const GWN_SOLID_ANGLE_CLIP: Real = 1e-6;

/// GWN threshold: `|wn| > GWN_INSIDE_THRESHOLD` → query is strictly inside.
pub const GWN_INSIDE_THRESHOLD: Real = 0.65;

/// GWN threshold: `|wn| < GWN_OUTSIDE_THRESHOLD` → query is strictly outside.
///
/// The band `[GWN_OUTSIDE_THRESHOLD, GWN_INSIDE_THRESHOLD]` triggers the
/// tiebreaker predicates in `classify_fragment`.
pub const GWN_OUTSIDE_THRESHOLD: Real = 0.35;

/// Scale-relative tolerance for the nearest-face signed distance tiebreaker.
///
/// # Theorem — Scale Invariance
///
/// The signed distance `d = cp · n / ‖n‖` has unit [length].  The face
/// characteristic scale is `√(area) ≈ √(‖n‖/2)` where `n = ab × ac`.
/// A fragment is coplanar when `|d| < TIEBREAK_SIGN_REL_TOL × √(area)`.
///
/// For a 1 mm edge triangle (area ≈ 4.3 × 10⁻⁷ m², scale ≈ 6.6 × 10⁻⁴ m),
/// the threshold is ≈ 6.6 × 10⁻¹¹ m — well below any physical geometry.
/// For a 1 m edge triangle (scale ≈ 0.7 m), threshold ≈ 7 × 10⁻⁸ m.
///
/// **Previous absolute threshold 1e-9** breaks for meshes at scales ≫ 1 m
/// (threshold too tight) or ≪ 1 mm (threshold too loose). ∎
pub const TIEBREAK_SIGN_REL_TOL: Real = 1e-7;

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

/// Relative AABB expansion factor for broad-phase Boolean operations.
///
/// The mesh-level AABB is expanded by `AABB_RELATIVE_EXPANSION * diagonal`
/// to guard against floating-point precision misses on snapped vertices.
///
/// # Theorem — Scale Correctness
///
/// For a mesh with AABB diagonal `d`, the expansion is `1e-6 · d`.  This
/// ensures the relative guard is constant across scales: a 10 µm mesh
/// expands by ≈ 10 fm (sub-atomic, harmless), while a 1 m mesh expands by
/// ≈ 1 µm.  The previous absolute 1e-6 m expansion was 10 % of a 10 µm
/// mesh's diagonal — large enough to merge disjoint features. ∎
pub const AABB_RELATIVE_EXPANSION: Real = 1e-6;

/// Relative degenerate-normal threshold factor for CDT corefine and
/// fragment classification.
///
/// A cross-product normal `n = (B−A) × (C−A)` is degenerate when
/// `‖n‖² < DEGENERATE_NORMAL_REL_SQ · ‖B−A‖² · ‖C−A‖²`.
///
/// # Theorem — Dimensionless Bound
///
/// `‖n‖² = ‖B−A‖² ‖C−A‖² sin²θ` where θ is the included angle.  The
/// threshold `sin²θ < 1e-20` triggers only for θ < 1e-10 rad (≈ 6e-9°),
/// regardless of scale.  Previously the absolute check `‖n‖² < 1e-20`
/// triggered for any triangle with edge lengths below ≈ 3e-10, which
/// misclassified micro-scale millifluidic geometry. ∎
pub const DEGENERATE_NORMAL_REL_SQ: Real = 1e-20;

/// Relative degenerate-segment threshold factor.
///
/// A snap segment with `‖end − start‖² < DEGENERATE_SEGMENT_REL_SQ · diag²`
/// is collapsed, where `diag` is the face's maximum edge length.
///
/// # Theorem — Scale Correctness
///
/// At any scale, a segment shorter than `1e-12 · max_edge` is below
/// double-precision resolution for that geometry.  The previous absolute
/// threshold `1e-24` is `(1e-12)²` which is correct for unit-scale meshes
/// but too strict at macro scale (1 m mesh: segments up to 1e-12 m would
/// be kept, wasting CDT effort on sub-picometer features). ∎
pub const DEGENERATE_SEGMENT_REL_SQ: Real = 1e-24;

/// Relative 2-D projected-area threshold for sliver detection in corefine.
///
/// A boundary polygon is degenerate when `|area2D| < SLIVER_AREA2D_REL · ∑edge²`.
///
/// # Theorem — Scale Independence
///
/// The 2-D shoelace area scales as length², so comparing against ∑edge²
/// yields a dimensionless ratio that is independent of mesh scale. ∎
pub const SLIVER_AREA2D_REL: Real = 1e-10;

/// Relative interval-overlap tolerance for T-T intersection segment
/// computation.
///
/// Two 1-D intervals `[t₁_min, t₁_max]` and `[t₂_min, t₂_max]` are
/// considered non-overlapping when `t_enter > t_leave + EPS · span` where
/// `span = max(|t₁_max − t₁_min|, |t₂_max − t₂_min|, diag)`.
///
/// # Theorem — Scale Independence
///
/// The 1-D projections onto the intersection line scale linearly with
/// mesh dimensions.  Using a relative tolerance ensures the gap test
/// is uniform across scales. ∎
pub const INTERVAL_OVERLAP_REL: Real = 1e-12;
