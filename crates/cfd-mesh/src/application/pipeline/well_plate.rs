//! ANSI SBS 96-well plate block dimensions.
//!
//! This module provides **only** the outer block dimensions from the SBS
//! microplate standard.  There are no wells — the constants are used purely
//! as the default substrate (chip body) block size.
//!
//! Channel centerlines should be placed at:
//! - `y = SbsWellPlate96::center_y()` — centered in the block width
//! - `z = chip_height_mm / 2.0`        — centered in the block height
//! - `x ∈ [0, WIDTH_MM]`               — inlet at left face, outlet at right

/// ANSI SBS 96-well plate outer dimensions used as the default substrate block size.
///
/// These are purely dimensional constants — the block has no wells.
pub struct SbsWellPlate96;

impl SbsWellPlate96 {
    /// Block length along X [mm].  Channels route from `x = 0` to `x = WIDTH_MM`.
    pub const WIDTH_MM: f64 = 127.76;

    /// Block width along Y [mm].  Channel centerline at `y = DEPTH_MM / 2`.
    pub const DEPTH_MM: f64 = 85.47;

    /// Standard center Y for all channels: `DEPTH_MM / 2 = 42.735 mm`.
    pub fn center_y() -> f64 {
        Self::DEPTH_MM / 2.0
    }

    /// `true` if a 2-D point `(x, y)` (in mm) lies strictly within the block
    /// minus `clearance_mm` on all sides.
    pub fn contains_point(x: f64, y: f64, clearance_mm: f64) -> bool {
        x >= clearance_mm
            && x <= Self::WIDTH_MM - clearance_mm
            && y >= clearance_mm
            && y <= Self::DEPTH_MM - clearance_mm
    }

    /// `true` if the straight segment from `(x0, y0)` to `(x1, y1)` (mm) lies
    /// within block bounds minus `clearance_mm` at **both** endpoints.
    pub fn segment_within_bounds(x0: f64, y0: f64, x1: f64, y1: f64, clearance_mm: f64) -> bool {
        Self::contains_point(x0, y0, clearance_mm)
            && Self::contains_point(x1, y1, clearance_mm)
    }

    /// `true` if the routing segment stays within the plate for channel layout.
    ///
    /// Inlet/outlet faces lie at `x = 0` and `x = WIDTH_MM` — segments are
    /// allowed to reach those faces (`x_clearance = 0`).  Side walls (Y
    /// direction) still require `side_clearance_mm` of keep-out margin.
    /// Any X coordinate outside `[0, WIDTH_MM]` is rejected.
    pub fn segment_within_routing_bounds(
        x0: f64,
        y0: f64,
        x1: f64,
        y1: f64,
        side_clearance_mm: f64,
    ) -> bool {
        let in_x = |x: f64| (0.0..=Self::WIDTH_MM).contains(&x);
        let in_y = |y: f64| (side_clearance_mm..=Self::DEPTH_MM - side_clearance_mm).contains(&y);
        in_x(x0) && in_y(y0) && in_x(x1) && in_y(y1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn center_y_is_half_depth() {
        let cy = SbsWellPlate96::center_y();
        assert!((cy - 42.735).abs() < 1e-6, "center_y = {cy}");
    }

    #[test]
    fn width_depth_constants() {
        assert!((SbsWellPlate96::WIDTH_MM - 127.76).abs() < 1e-9);
        assert!((SbsWellPlate96::DEPTH_MM - 85.47).abs() < 1e-9);
    }

    #[test]
    fn contains_point_inside() {
        assert!(SbsWellPlate96::contains_point(63.88, 42.735, 5.0));
    }

    #[test]
    fn contains_point_outside_left() {
        assert!(!SbsWellPlate96::contains_point(1.0, 42.7, 5.0));
    }

    #[test]
    fn segment_within_bounds_rejects_out_of_range() {
        // segment starting at x=0 must fail with 5 mm clearance
        assert!(!SbsWellPlate96::segment_within_bounds(
            0.0,
            42.735,
            127.76,
            42.735,
            5.0
        ));
    }
}
