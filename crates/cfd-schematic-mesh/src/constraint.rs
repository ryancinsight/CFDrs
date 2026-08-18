//! Inlet/outlet diameter and wall-clearance constraints.

use std::fmt;

use aequitas::systems::si::{
    quantities::Length,
    units::{Meter, Millimeter},
};
use cfd_schematics::{NetworkBlueprint, NodeKind};

use super::well_plate::SbsWellPlate96;

/// Required hydraulic diameter for inlet and outlet channels.
pub const REQUIRED_HYDRAULIC_DIAMETER: Length<f64> = Length::from_base(4.0e-3);
/// Allowed tolerance around the required diameter.
pub const HYDRAULIC_DIAMETER_TOLERANCE: Length<f64> = Length::from_base(0.1e-3);
/// Default wall clearance.
pub const DEFAULT_WALL_CLEARANCE: Length<f64> = Length::from_base(5.0e-3);

/// Error returned when a channel's hydraulic diameter does not meet the 4 mm spec.
#[derive(Debug)]
pub enum DiameterConstraintError {
    /// Channel at an inlet/outlet node has the wrong hydraulic diameter.
    WrongDiameter {
        /// Channel identifier.
        channel_id: String,
        /// Node identifier.
        node_id: String,
        /// Node kind label ("inlet" or "outlet").
        node_kind: &'static str,
        /// Actual measured hydraulic diameter.
        actual: Length<f64>,
        /// Required hydraulic diameter.
        expected: Length<f64>,
        /// Allowed tolerance.
        tolerance: Length<f64>,
    },
    /// A node has no adjacent channels.
    IsolatedNode {
        /// Node identifier.
        node_id: String,
    },
}

impl fmt::Display for DiameterConstraintError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongDiameter {
                channel_id,
                node_id,
                node_kind,
                actual,
                expected,
                tolerance,
            } => write!(
                formatter,
                "channel '{channel_id}' at {node_kind} node '{node_id}' has hydraulic diameter {:.4e} m (expected {:.4e} ± {:.4e} m)",
                actual.in_unit::<Meter>(),
                expected.in_unit::<Meter>(),
                tolerance.in_unit::<Meter>(),
            ),
            Self::IsolatedNode { node_id } => {
                write!(
                    formatter,
                    "node '{node_id}' is isolated — no channels are adjacent"
                )
            }
        }
    }
}

impl std::error::Error for DiameterConstraintError {}

/// Validates that all inlet and outlet nodes have channels with hydraulic
/// diameter `4.0 mm ± 0.1 mm`.
pub struct InletOutletConstraint;

impl InletOutletConstraint {
    /// Check all inlet and outlet channels in `bp`.
    ///
    /// # Errors
    /// Returns `Err` if any inlet or outlet channel has a hydraulic diameter
    /// outside the required `4.0 mm ± 0.1 mm` range.
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
                .filter(|c| {
                    c.from.as_str() == node.id.as_str() || c.to.as_str() == node.id.as_str()
                })
                .collect();

            if adjacent.is_empty() {
                return Err(DiameterConstraintError::IsolatedNode {
                    node_id: node.id.to_string(),
                });
            }

            for ch in adjacent {
                let actual = ch.cross_section.hydraulic_diameter();
                // 1e-12 m slop guards against floating-point rounding at the tolerance boundary.
                if (actual.in_unit::<Meter>() - REQUIRED_HYDRAULIC_DIAMETER.in_unit::<Meter>())
                    .abs()
                    > HYDRAULIC_DIAMETER_TOLERANCE.in_unit::<Meter>() + 1e-12
                {
                    return Err(DiameterConstraintError::WrongDiameter {
                        channel_id: ch.id.as_str().to_string(),
                        node_id: node.id.to_string(),
                        node_kind,
                        actual,
                        expected: REQUIRED_HYDRAULIC_DIAMETER,
                        tolerance: HYDRAULIC_DIAMETER_TOLERANCE,
                    });
                }
            }
        }
        Ok(())
    }
}

/// Wall clearance violation — a channel segment is too close to a plate edge.
#[derive(Debug)]
pub struct WallClearanceViolation {
    /// Start X of the violating segment.
    pub x0: Length<f64>,
    /// Start Y of the violating segment.
    pub y0: Length<f64>,
    /// End X of the violating segment.
    pub x1: Length<f64>,
    /// End Y of the violating segment.
    pub y1: Length<f64>,
    /// Required clearance.
    pub clearance: Length<f64>,
    /// Plate width.
    pub plate_width: Length<f64>,
    /// Plate depth.
    pub plate_depth: Length<f64>,
}

impl fmt::Display for WallClearanceViolation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "channel segment ({:.2}, {:.2}) → ({:.2}, {:.2}) mm violates {:.1} mm wall clearance on plate {:.2} × {:.2} mm",
            self.x0.in_unit::<Millimeter>(),
            self.y0.in_unit::<Millimeter>(),
            self.x1.in_unit::<Millimeter>(),
            self.y1.in_unit::<Millimeter>(),
            self.clearance.in_unit::<Millimeter>(),
            self.plate_width.in_unit::<Millimeter>(),
            self.plate_depth.in_unit::<Millimeter>(),
        )
    }
}

impl std::error::Error for WallClearanceViolation {}

/// Validates that all channel segments remain within the SBS plate bounds.
pub struct WallClearanceConstraint;

impl WallClearanceConstraint {
    /// Check that all `(x0, y0) → (x1, y1)` segments stay inside the SBS plate
    /// with the given clearance.
    ///
    /// # Errors
    /// Returns `Err` if any segment falls outside the SBS plate boundary
    /// minus the required clearance margin.
    #[allow(clippy::type_complexity)]
    pub fn check(
        segments: &[((Length<f64>, Length<f64>), (Length<f64>, Length<f64>))],
        clearance: Length<f64>,
    ) -> Result<(), WallClearanceViolation> {
        for &((x0, y0), (x1, y1)) in segments {
            if !SbsWellPlate96::segment_within_bounds(x0, y0, x1, y1, clearance) {
                return Err(WallClearanceViolation {
                    x0,
                    y0,
                    x1,
                    y1,
                    clearance,
                    plate_width: SbsWellPlate96::WIDTH,
                    plate_depth: SbsWellPlate96::DEPTH,
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use aequitas::systems::si::units::Millimeter;
    use cfd_schematics::interface::presets::{serpentine_chain, venturi_chain};
    use cfd_schematics::NetworkBlueprint;

    use super::*;

    fn explicit_layout_blueprint(name: &str) -> NetworkBlueprint {
        NetworkBlueprint {
            name: name.to_string(),
            box_dims: (127.76, 85.47),
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        }
    }

    #[test]
    fn four_mm_circular_passes() {
        let bp = venturi_chain("v", 0.030, 0.004, 0.002);
        assert_eq!(
            InletOutletConstraint::check(&bp).map_err(|error| error.to_string()),
            Ok(())
        );
    }

    #[test]
    fn two_mm_circular_fails() {
        let bp = serpentine_chain("s", 3, 0.010, 0.002);
        let error = InletOutletConstraint::check(&bp)
            .expect_err("two millimetre inlet must violate the four millimetre contract");
        match error {
            DiameterConstraintError::WrongDiameter {
                actual,
                expected,
                tolerance,
                ..
            } => {
                assert!((actual.in_unit::<Millimeter>() - 2.0).abs() < f64::EPSILON);
                assert_eq!(expected, REQUIRED_HYDRAULIC_DIAMETER);
                assert_eq!(tolerance, HYDRAULIC_DIAMETER_TOLERANCE);
            }
            DiameterConstraintError::IsolatedNode { .. } => panic!("unexpected constraint error: {error}"),
        }
    }

    #[test]
    fn four_mm_at_tolerance_boundary_passes() {
        // diameter_m = 4.1 mm is within ±0.1 mm
        use cfd_schematics::{ChannelSpec, NodeKind, NodeSpec};
        let mut bp = explicit_layout_blueprint("t");
        bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
        bp.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (10.0, 0.0)));
        bp.add_channel(ChannelSpec::new_pipe(
            "c", "inlet", "outlet", 0.01, 0.0041, 0.0, 0.0,
        ));
        assert_eq!(
            InletOutletConstraint::check(&bp).map_err(|error| error.to_string()),
            Ok(())
        );
    }

    #[test]
    fn isolated_node_fails() {
        use cfd_schematics::{NodeKind, NodeSpec};
        let mut bp = explicit_layout_blueprint("x");
        bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
        assert!(matches!(
            InletOutletConstraint::check(&bp),
            Err(DiameterConstraintError::IsolatedNode { node_id }) if node_id == "inlet"
        ));
    }

    #[test]
    fn wall_clearance_passes_for_center_segment() {
        let segments = vec![(
            (
                Length::from_unit::<Millimeter>(10.0),
                Length::from_unit::<Millimeter>(42.735),
            ),
            (
                Length::from_unit::<Millimeter>(117.76),
                Length::from_unit::<Millimeter>(42.735),
            ),
        )];
        assert_eq!(
            WallClearanceConstraint::check(&segments, Length::from_unit::<Millimeter>(5.0),)
                .map_err(|error| error.to_string()),
            Ok(())
        );
    }

    #[test]
    fn wall_clearance_fails_for_out_of_bounds_segment() {
        let segments = vec![(
            (
                Length::from_unit::<Millimeter>(0.0),
                Length::from_unit::<Millimeter>(42.735),
            ),
            (
                Length::from_unit::<Millimeter>(127.76),
                Length::from_unit::<Millimeter>(42.735),
            ),
        )];
        let error = WallClearanceConstraint::check(&segments, Length::from_unit::<Millimeter>(5.0))
            .expect_err("a segment touching both plate faces must violate clearance");
        assert_eq!(error.x0, Length::from_unit::<Millimeter>(0.0));
        assert_eq!(error.x1, Length::from_unit::<Millimeter>(127.76));
        assert_eq!(error.clearance, Length::from_unit::<Millimeter>(5.0));
        assert_eq!(error.plate_width, SbsWellPlate96::WIDTH);
        assert_eq!(error.plate_depth, SbsWellPlate96::DEPTH);
    }

    #[test]
    fn tolerance_boundary_near_4mm_passes() {
        // exactly REQUIRED - TOLERANCE should still pass
        use cfd_schematics::{ChannelSpec, NodeKind, NodeSpec};
        let mut bp = explicit_layout_blueprint("b");
        bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
        bp.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (10.0, 0.0)));
        // 3.9 mm = 4.0 - 0.1 (borderline)
        bp.add_channel(ChannelSpec::new_pipe(
            "c", "inlet", "outlet", 0.01, 0.0039, 0.0, 0.0,
        ));
        // 3.9 mm is exactly at the boundary (diff = 0.1e-3 = tolerance), should pass
        assert_eq!(
            InletOutletConstraint::check(&bp).map_err(|error| error.to_string()),
            Ok(())
        );
    }
}
