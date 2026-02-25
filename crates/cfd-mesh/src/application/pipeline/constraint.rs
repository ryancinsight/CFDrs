//! Inlet/outlet diameter and wall-clearance constraints.

use cfd_schematics::{NetworkBlueprint, NodeKind};
use thiserror::Error;

use super::well_plate::SbsWellPlate96;

/// Required hydraulic diameter for inlet and outlet channels [m].
pub const REQUIRED_HYD_DIAM_M: f64 = 4.0e-3;
/// Allowed tolerance around the required diameter [m].
pub const HYD_DIAM_TOLERANCE_M: f64 = 0.1e-3;
/// Default wall clearance [mm].
pub const DEFAULT_WALL_CLEARANCE: f64 = 5.0;

/// Error returned when a channel's hydraulic diameter does not meet the 4 mm spec.
#[derive(Debug, Error)]
pub enum DiameterConstraintError {
    /// Channel at an inlet/outlet node has the wrong hydraulic diameter.
    #[error(
        "channel '{channel_id}' at {node_kind} node '{node_id}' has hydraulic \
         diameter {actual_m:.4e} m (expected {expected_m:.4e} ± {tolerance_m:.4e} m)"
    )]
    WrongDiameter {
        /// Channel identifier.
        channel_id: String,
        /// Node identifier.
        node_id: String,
        /// Node kind label ("inlet" or "outlet").
        node_kind: &'static str,
        /// Actual measured hydraulic diameter.
        actual_m: f64,
        /// Required hydraulic diameter.
        expected_m: f64,
        /// Allowed tolerance.
        tolerance_m: f64,
    },
    /// A node has no adjacent channels.
    #[error("node '{node_id}' is isolated — no channels are adjacent")]
    IsolatedNode {
        /// Node identifier.
        node_id: String,
    },
}

/// Validates that all inlet and outlet nodes have channels with hydraulic
/// diameter `4.0 mm ± 0.1 mm`.
pub struct InletOutletConstraint;

impl InletOutletConstraint {
    /// Check all inlet and outlet channels in `bp`.
    pub fn check(bp: &NetworkBlueprint) -> Result<(), DiameterConstraintError> {
        for node in &bp.nodes {
            let node_kind = match node.kind {
                NodeKind::Inlet => "inlet",
                NodeKind::Outlet => "outlet",
                _ => continue,
            };

            // Find channels adjacent to this node
            let adjacent: Vec<_> = bp
                .channels
                .iter()
                .filter(|c| c.from.as_str() == node.id.as_str() || c.to.as_str() == node.id.as_str())
                .collect();

            if adjacent.is_empty() {
                return Err(DiameterConstraintError::IsolatedNode {
                    node_id: node.id.to_string(),
                });
            }

            for ch in adjacent {
                let actual = ch.cross_section.hydraulic_diameter();
                // 1e-12 m slop guards against floating-point rounding at the tolerance boundary.
                if (actual - REQUIRED_HYD_DIAM_M).abs() > HYD_DIAM_TOLERANCE_M + 1e-12 {
                    return Err(DiameterConstraintError::WrongDiameter {
                        channel_id: ch.id.as_str().to_string(),
                        node_id: node.id.to_string(),
                        node_kind,
                        actual_m: actual,
                        expected_m: REQUIRED_HYD_DIAM_M,
                        tolerance_m: HYD_DIAM_TOLERANCE_M,
                    });
                }
            }
        }
        Ok(())
    }
}

/// Wall clearance violation — a channel segment is too close to a plate edge.
#[derive(Debug, Error)]
#[error(
    "channel segment ({x0:.2}, {y0:.2}) → ({x1:.2}, {y1:.2}) mm violates \
     {clearance:.1} mm wall clearance on plate {plate_w:.2} × {plate_d:.2} mm"
)]
pub struct WallClearanceViolation {
    /// Start X of the violating segment [mm].
    pub x0: f64,
    /// Start Y of the violating segment [mm].
    pub y0: f64,
    /// End X of the violating segment [mm].
    pub x1: f64,
    /// End Y of the violating segment [mm].
    pub y1: f64,
    /// Required clearance [mm].
    pub clearance: f64,
    /// Plate width [mm].
    pub plate_w: f64,
    /// Plate depth [mm].
    pub plate_d: f64,
}

/// Validates that all channel segments remain within the SBS plate bounds.
pub struct WallClearanceConstraint;

impl WallClearanceConstraint {
    /// Check that all `(x0, y0) → (x1, y1)` segments stay inside the SBS plate
    /// with the given clearance (in mm).
    #[allow(clippy::type_complexity)]
    pub fn check(
        segments: &[((f64, f64), (f64, f64))],
        clearance_mm: f64,
    ) -> Result<(), WallClearanceViolation> {
        for &((x0, y0), (x1, y1)) in segments {
            if !SbsWellPlate96::segment_within_bounds(x0, y0, x1, y1, clearance_mm) {
                return Err(WallClearanceViolation {
                    x0,
                    y0,
                    x1,
                    y1,
                    clearance: clearance_mm,
                    plate_w: SbsWellPlate96::WIDTH_MM,
                    plate_d: SbsWellPlate96::DEPTH_MM,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use cfd_schematics::interface::presets::{serpentine_chain, venturi_chain};

    use super::*;

    #[test]
    fn four_mm_circular_passes() {
        let bp = venturi_chain("v", 0.030, 0.004, 0.002);
        assert!(InletOutletConstraint::check(&bp).is_ok());
    }

    #[test]
    fn two_mm_circular_fails() {
        let bp = serpentine_chain("s", 3, 0.010, 0.002);
        assert!(InletOutletConstraint::check(&bp).is_err());
    }

    #[test]
    fn four_mm_at_tolerance_boundary_passes() {
        // diameter_m = 4.1 mm is within ±0.1 mm
        use cfd_schematics::{ChannelSpec, NodeKind, NodeSpec};
        let mut bp = cfd_schematics::NetworkBlueprint::new("t");
        bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
        bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
        bp.add_channel(ChannelSpec::new_pipe("c", "inlet", "outlet", 0.01, 0.0041, 0.0, 0.0));
        assert!(InletOutletConstraint::check(&bp).is_ok());
    }

    #[test]
    fn isolated_node_fails() {
        use cfd_schematics::{NodeKind, NodeSpec};
        let mut bp = cfd_schematics::NetworkBlueprint::new("x");
        bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
        assert!(InletOutletConstraint::check(&bp).is_err());
    }

    #[test]
    fn wall_clearance_passes_for_center_segment() {
        let segments = vec![((10.0, 42.735), (117.76, 42.735))];
        assert!(WallClearanceConstraint::check(&segments, 5.0).is_ok());
    }

    #[test]
    fn wall_clearance_fails_for_out_of_bounds_segment() {
        let segments = vec![((0.0, 42.735), (127.76, 42.735))];
        assert!(WallClearanceConstraint::check(&segments, 5.0).is_err());
    }

    #[test]
    fn tolerance_boundary_near_4mm_passes() {
        // exactly REQUIRED - TOLERANCE should still pass
        use cfd_schematics::{ChannelSpec, NodeKind, NodeSpec};
        let mut bp = cfd_schematics::NetworkBlueprint::new("b");
        bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
        bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
        // 3.9 mm = 4.0 - 0.1 (borderline)
        bp.add_channel(ChannelSpec::new_pipe("c", "inlet", "outlet", 0.01, 0.0039, 0.0, 0.0));
        // 3.9 mm is exactly at the boundary (diff = 0.1e-3 = tolerance), should pass
        assert!(InletOutletConstraint::check(&bp).is_ok());
    }
}
