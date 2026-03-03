//! `ShellCuboid` — a channelless millifluidic cuboid schematic.
//!
//! A `ShellCuboid` represents a millifluidic device modelled as a **shell**:
//! there are no internal channels, only an outer bounding wall and a concentric
//! inner cavity separated by a uniform wall thickness. A single inlet port and
//! a single outlet port pierce the shell on the left and right walls at the
//! horizontal midline, respectively.
//!
//! # 2D Schematic Geometry (units: mm)
//!
//! ```text
//!  (0, H) ──────────────────────────── (W, H)
//!    │  (t, H-t) ─────────── (W-t, H-t)  │
//!    │     │       (cavity)        │      │
//! ───┤──┤  │                       │  ├──├───
//! inlet  (t, t) ────────── (W-t, t)  outlet
//!    │                                    │
//!  (0, 0) ──────────────────────────── (W, 0)
//! ```
//!
//! - `W` = `outer_dims.0`, `H` = `outer_dims.1`, `t` = `shell_thickness_mm`.
//! - Inlet port: outer corner `(0, H/2)` ↔ inner corner `(t, H/2)`.
//! - Outlet port: inner corner `(W-t, H/2)` ↔ outer corner `(W, H/2)`.
//! - `box_outline` contains exactly 10 segments: 4 outer + 4 inner + 2 port stubs.

use crate::error::{GeometryError, GeometryResult};
use serde::{Deserialize, Serialize};

use super::tpms_fill::TpmsFillSpec;
use super::Point2D;

/// Schema/producer constants — mirror those in `ChannelSystem`.
const SCHEMA_VERSION: &str = "1.0.0";
const PRODUCER_PREFIX: &str = "scheme";
const LENGTH_UNITS: &str = "mm";

// ── Core type ──────────────────────────────────────────────────────────────────

/// A millifluidic cuboid modelled as a shell (no internal channels).
///
/// The device is a rectangular box with uniform wall thickness on all four
/// sides. There are no channels inside — the interior is an open cavity.
/// One inlet port (left wall) and one outlet port (right wall) pierce the
/// shell at the horizontal midline, represented as short line stubs between
/// the outer and inner rectangles.
///
/// Use [`ShellCuboid::new`] or the generator
/// [`create_shell_cuboid`](crate::geometry::generator::shell::create_shell_cuboid)
/// to construct this type.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ShellCuboid {
    /// Outer bounding-box dimensions `(width_mm, height_mm)`.
    pub outer_dims: (f64, f64),
    /// Uniform wall thickness applied on all four sides (mm).
    pub shell_thickness_mm: f64,
    /// Derived inner cavity dimensions `(width_mm, height_mm)`.
    ///
    /// Invariant: both components are strictly positive.
    pub inner_dims: (f64, f64),
    /// Inlet port midpoint on the **outer** left wall (mm).
    pub inlet_port: Point2D,
    /// Outlet port midpoint on the **outer** right wall (mm).
    pub outlet_port: Point2D,
    /// All 10 line segments defining the complete schematic outline.
    ///
    /// Ordering: 4 outer-rectangle edges, 4 inner-rectangle edges,
    /// 1 inlet stub, 1 outlet stub.
    pub box_outline: Vec<(Point2D, Point2D)>,
    /// Optional TPMS lattice fill specification for the inner cavity.
    ///
    /// When `Some`, downstream mesh pipelines will fill the cavity with
    /// the specified TPMS surface.  When `None`, the cavity is empty.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tpms_fill: Option<TpmsFillSpec>,
}

impl ShellCuboid {
    /// ANSI/SLAS 1-2004 96-well plate footprint: 127.76 × 85.47 mm.
    ///
    /// All millifluidic designs in this project target this standard
    /// microplate form factor.
    pub const WELL_PLATE_96_DIMS: (f64, f64) = (127.76, 85.47);

    /// Construct a `ShellCuboid` with the 96-well plate footprint.
    ///
    /// Equivalent to `ShellCuboid::new(Self::WELL_PLATE_96_DIMS, shell_thickness_mm)`.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError`] if `shell_thickness_mm` violates the
    /// geometric invariants (see [`new`](ShellCuboid::new)).
    pub fn well_plate_96(shell_thickness_mm: f64) -> GeometryResult<Self> {
        Self::new(Self::WELL_PLATE_96_DIMS, shell_thickness_mm)
    }

    /// Construct a `ShellCuboid` from outer dimensions and shell thickness.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError`] when:
    /// - any outer dimension is non-positive or non-finite, or
    /// - `shell_thickness_mm` is non-positive, non-finite, or
    ///   large enough that the inner cavity would be non-positive in either
    ///   dimension (i.e. `thickness ≥ outer_dim / 2`).
    pub fn new(outer_dims: (f64, f64), shell_thickness_mm: f64) -> GeometryResult<Self> {
        Self::validate_params(outer_dims, shell_thickness_mm)?;

        let (w, h) = outer_dims;
        let t = shell_thickness_mm;

        let inner_dims = (w - 2.0 * t, h - 2.0 * t);

        // Inlet: left-wall midline, outer→inner
        let inlet_port: Point2D = (0.0, h / 2.0);
        let inlet_inner: Point2D = (t, h / 2.0);

        // Outlet: right-wall midline, inner→outer
        let outlet_inner: Point2D = (w - t, h / 2.0);
        let outlet_port: Point2D = (w, h / 2.0);

        // Outer rectangle (CCW from bottom-left)
        let outer_rect: [(Point2D, Point2D); 4] = [
            ((0.0, 0.0), (w, 0.0)),
            ((w, 0.0), (w, h)),
            ((w, h), (0.0, h)),
            ((0.0, h), (0.0, 0.0)),
        ];

        // Inner rectangle (CCW from inner bottom-left)
        let (iw, ih) = inner_dims;
        let inner_rect: [(Point2D, Point2D); 4] = [
            ((t, t), (t + iw, t)),
            ((t + iw, t), (t + iw, t + ih)),
            ((t + iw, t + ih), (t, t + ih)),
            ((t, t + ih), (t, t)),
        ];

        let mut box_outline: Vec<(Point2D, Point2D)> = Vec::with_capacity(10);
        box_outline.extend_from_slice(&outer_rect);
        box_outline.extend_from_slice(&inner_rect);
        box_outline.push((inlet_port, inlet_inner)); // inlet stub
        box_outline.push((outlet_inner, outlet_port)); // outlet stub

        Ok(Self {
            outer_dims,
            shell_thickness_mm,
            inner_dims,
            inlet_port,
            outlet_port,
            box_outline,
            tpms_fill: None,
        })
    }

    /// Attach a TPMS fill specification to fill the inner cavity.
    ///
    /// Returns `self` for chaining.
    #[must_use]
    pub fn with_tpms_fill(mut self, fill: TpmsFillSpec) -> Self {
        self.tpms_fill = Some(fill);
        self
    }

    // ── Validation ─────────────────────────────────────────────────────────

    /// Validate outer dimensions and shell thickness independently of construction.
    ///
    /// Called by [`new`](ShellCuboid::new) and [`validate`](ShellCuboid::validate).
    fn validate_params(outer_dims: (f64, f64), shell_thickness_mm: f64) -> GeometryResult<()> {
        let (w, h) = outer_dims;
        if !w.is_finite() || w <= 0.0 || !h.is_finite() || h <= 0.0 {
            return Err(GeometryError::invalid_box_dimensions(w, h));
        }

        if !shell_thickness_mm.is_finite() || shell_thickness_mm <= 0.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "shell_thickness_mm must be a finite positive value, got {shell_thickness_mm}"
                ),
            });
        }

        let inner_w = w - 2.0 * shell_thickness_mm;
        let inner_h = h - 2.0 * shell_thickness_mm;
        if inner_w <= 0.0 || inner_h <= 0.0 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "shell_thickness_mm ({shell_thickness_mm}) is too large for outer_dims \
                         ({w}×{h}): inner cavity dimensions would be ({inner_w}×{inner_h})"
                ),
            });
        }

        Ok(())
    }

    /// Validate the invariants of an already-constructed [`ShellCuboid`].
    ///
    /// Useful for checking round-trip correctness after JSON deserialisation.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError`] on any invariant violation.
    pub fn validate(&self) -> GeometryResult<()> {
        Self::validate_params(self.outer_dims, self.shell_thickness_mm)?;

        // Derived fields must be consistent.
        let (w, h) = self.outer_dims;
        let t = self.shell_thickness_mm;
        let expected_inner = (w - 2.0 * t, h - 2.0 * t);
        if (self.inner_dims.0 - expected_inner.0).abs() > 1e-9
            || (self.inner_dims.1 - expected_inner.1).abs() > 1e-9
        {
            return Err(GeometryError::InvalidChannelPath {
                reason: "inner_dims are inconsistent with outer_dims and shell_thickness_mm"
                    .to_string(),
            });
        }

        if self.box_outline.len() != 10 {
            return Err(GeometryError::InvalidChannelPath {
                reason: format!(
                    "box_outline must contain exactly 10 segments, found {}",
                    self.box_outline.len()
                ),
            });
        }

        // Validate TPMS fill if present.
        if let Some(ref fill) = self.tpms_fill {
            fill.validate()?;
        }

        Ok(())
    }

    // ── JSON serialisation ──────────────────────────────────────────────────

    /// Serialise this `ShellCuboid` to a pretty-printed JSON string.
    ///
    /// # Errors
    ///
    /// Returns a [`serde_json::Error`] if serialisation fails (should not
    /// occur for well-formed instances).
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }

    /// Deserialise a `ShellCuboid` from a JSON string.
    ///
    /// Call [`validate`](ShellCuboid::validate) on the result to confirm
    /// all derived fields are consistent.
    ///
    /// # Errors
    ///
    /// Returns a [`serde_json::Error`] if the JSON is malformed.
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    // ── Interchange ─────────────────────────────────────────────────────────

    /// Schema version constant for interchange exports.
    pub const INTERCHANGE_SCHEMA_VERSION: &'static str = SCHEMA_VERSION;
    /// Length units used in all coordinate and dimension fields.
    pub const INTERCHANGE_LENGTH_UNITS: &'static str = LENGTH_UNITS;
    /// Producer identifier prefix for interchange exports.
    pub const INTERCHANGE_PRODUCER_PREFIX: &'static str = PRODUCER_PREFIX;
}

// ── Tests ──────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_80x40() -> ShellCuboid {
        ShellCuboid::new((80.0, 40.0), 2.0).expect("valid shell cuboid")
    }

    #[test]
    fn shell_cuboid_positive_cavity() {
        let sc = make_80x40();
        assert!(
            (sc.inner_dims.0 - 76.0).abs() < 1e-9,
            "inner width should be 76.0"
        );
        assert!(
            (sc.inner_dims.1 - 36.0).abs() < 1e-9,
            "inner height should be 36.0"
        );
    }

    #[test]
    fn shell_cuboid_inlet_outlet_positions() {
        let sc = make_80x40();
        // Inlet outer point on left wall at midline
        assert!((sc.inlet_port.0 - 0.0).abs() < 1e-9);
        assert!((sc.inlet_port.1 - 20.0).abs() < 1e-9);
        // Outlet outer point on right wall at midline
        assert!((sc.outlet_port.0 - 80.0).abs() < 1e-9);
        assert!((sc.outlet_port.1 - 20.0).abs() < 1e-9);
    }

    #[test]
    fn shell_cuboid_ten_segments() {
        let sc = make_80x40();
        assert_eq!(
            sc.box_outline.len(),
            10,
            "box_outline must have exactly 10 segments"
        );
    }

    #[test]
    fn shell_cuboid_validate_passes() {
        let sc = make_80x40();
        sc.validate()
            .expect("valid shell cuboid should pass validate()");
    }

    #[test]
    fn shell_cuboid_json_roundtrip() {
        let sc = make_80x40();
        let json = sc.to_json().expect("serialisation should succeed");
        let sc2 = ShellCuboid::from_json(&json).expect("deserialisation should succeed");
        sc2.validate()
            .expect("deserialized shell cuboid should be valid");
        assert_eq!(sc, sc2, "JSON round-trip must be lossless");
    }

    #[test]
    fn shell_cuboid_invalid_thickness_too_large() {
        // thickness ≥ half dimension on height axis (3.0 > 4.0/2 = 2.0)
        let result = ShellCuboid::new((10.0, 4.0), 3.0);
        assert!(
            result.is_err(),
            "thickness larger than half-height must fail"
        );
    }

    #[test]
    fn shell_cuboid_invalid_zero_outer_dim() {
        let result = ShellCuboid::new((0.0, 40.0), 2.0);
        assert!(result.is_err(), "zero outer width must fail");
    }

    #[test]
    fn shell_cuboid_invalid_negative_thickness() {
        let result = ShellCuboid::new((80.0, 40.0), -1.0);
        assert!(result.is_err(), "negative thickness must fail");
    }

    #[test]
    fn shell_cuboid_well_plate_96_dims() {
        let sc = ShellCuboid::well_plate_96(2.0).expect("valid plate shell");
        assert!(
            (sc.outer_dims.0 - 127.76).abs() < 1e-9,
            "width should be 127.76 mm"
        );
        assert!(
            (sc.outer_dims.1 - 85.47).abs() < 1e-9,
            "height should be 85.47 mm"
        );
        assert!(
            (sc.inner_dims.0 - 123.76).abs() < 1e-9,
            "inner width: 127.76 - 4"
        );
        assert!(
            (sc.inner_dims.1 - 81.47).abs() < 1e-9,
            "inner height: 85.47 - 4"
        );
        assert!((sc.inlet_port.1 - 42.735).abs() < 1e-9, "inlet at midline");
        sc.validate().expect("plate shell should pass validate()");
    }
}
