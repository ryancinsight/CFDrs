//! ANSI SBS 96-well plate block dimensions.
//!
//! This module provides **only** the outer block dimensions from the SBS
//! microplate standard.  There are no wells — the constants are used purely
//! as the default substrate (chip body) block size.
//!
//! Channel centerlines should be placed at:
//! - `y = SbsWellPlate96::center_y()` — centered in the block width
//! - `z = chip_height_mm / 2.0`        — centered in the block height
//! - `x ∈ [0, WIDTH]`                  — inlet at left face, outlet at right

use aequitas::systems::si::{quantities::Length, units::Meter};

/// ANSI SBS 96-well plate outer dimensions used as the default substrate block size.
///
/// These are purely dimensional constants — the block has no wells.
pub struct SbsWellPlate96;

impl SbsWellPlate96 {
    /// Block length along X. Channels route from `x = 0` to `x = WIDTH`.
    pub const WIDTH: Length<f64> = Length::from_base(0.12776);

    /// Block width along Y. Channel centerline is at `y = DEPTH / 2`.
    pub const DEPTH: Length<f64> = Length::from_base(0.08547);

    /// Standard center Y for all channels: `DEPTH / 2 = 42.735 mm`.
    pub fn center_y() -> Length<f64> {
        Length::from_base(Self::DEPTH.as_base() / 2.0)
    }

    /// `true` if a 2-D point `(x, y)` lies within the block minus clearance
    /// on all sides.
    pub fn contains_point(x: Length<f64>, y: Length<f64>, clearance: Length<f64>) -> bool {
        let x = x.in_unit::<Meter>();
        let y = y.in_unit::<Meter>();
        let clearance = clearance.in_unit::<Meter>();
        x >= clearance
            && x <= Self::WIDTH.in_unit::<Meter>() - clearance
            && y >= clearance
            && y <= Self::DEPTH.in_unit::<Meter>() - clearance
    }

    /// `true` if the straight segment from `(x0, y0)` to `(x1, y1)` lies
    /// within block bounds minus clearance at **both** endpoints.
    pub fn segment_within_bounds(
        x0: Length<f64>,
        y0: Length<f64>,
        x1: Length<f64>,
        y1: Length<f64>,
        clearance: Length<f64>,
    ) -> bool {
        Self::contains_point(x0, y0, clearance) && Self::contains_point(x1, y1, clearance)
    }

    /// `true` if the routing segment stays within the plate for channel layout.
    ///
    /// Inlet/outlet faces lie at `x = 0` and `x = WIDTH` — segments are
    /// allowed to reach those faces (`x_clearance = 0`).  Side walls (Y
    /// direction) still require `side_clearance` of keep-out margin.
    /// Any X coordinate outside `[0, WIDTH]` is rejected.
    pub fn segment_within_routing_bounds(
        x0: Length<f64>,
        y0: Length<f64>,
        x1: Length<f64>,
        y1: Length<f64>,
        side_clearance: Length<f64>,
    ) -> bool {
        let in_x = |x: Length<f64>| {
            let x = x.in_unit::<Meter>();
            (0.0..=Self::WIDTH.in_unit::<Meter>()).contains(&x)
        };
        let in_y = |y: Length<f64>| {
            let y = y.in_unit::<Meter>();
            let clearance = side_clearance.in_unit::<Meter>();
            (clearance..=Self::DEPTH.in_unit::<Meter>() - clearance).contains(&y)
        };
        in_x(x0) && in_y(y0) && in_x(x1) && in_y(y1)
    }
}

#[cfg(test)]
mod tests {
    use aequitas::systems::si::units::Millimeter;

    use super::*;

    #[test]
    fn center_y_is_half_depth() {
        let cy = SbsWellPlate96::center_y();
        let cy_mm = cy.in_unit::<Millimeter>();
        assert!((cy_mm - 42.735).abs() < 1e-6, "center_y = {cy_mm}");
    }

    #[test]
    fn width_depth_constants() {
        assert!((SbsWellPlate96::WIDTH.in_unit::<Millimeter>() - 127.76).abs() < 1e-9);
        assert!((SbsWellPlate96::DEPTH.in_unit::<Millimeter>() - 85.47).abs() < 1e-9);
    }

    #[test]
    fn contains_point_inside() {
        assert!(SbsWellPlate96::contains_point(
            Length::from_unit::<Millimeter>(63.88),
            Length::from_unit::<Millimeter>(42.735),
            Length::from_unit::<Millimeter>(5.0),
        ));
    }

    #[test]
    fn contains_point_outside_left() {
        assert!(!SbsWellPlate96::contains_point(
            Length::from_unit::<Millimeter>(1.0),
            Length::from_unit::<Millimeter>(42.7),
            Length::from_unit::<Millimeter>(5.0),
        ));
    }

    #[test]
    fn segment_within_bounds_rejects_out_of_range() {
        // segment starting at x=0 must fail with 5 mm clearance
        assert!(!SbsWellPlate96::segment_within_bounds(
            Length::from_unit::<Millimeter>(0.0),
            Length::from_unit::<Millimeter>(42.735),
            Length::from_unit::<Millimeter>(127.76),
            Length::from_unit::<Millimeter>(42.735),
            Length::from_unit::<Millimeter>(5.0),
        ));
    }
}
