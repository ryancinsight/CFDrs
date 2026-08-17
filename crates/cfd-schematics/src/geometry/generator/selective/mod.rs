mod builder;
mod path_geometry;
mod primitive;
mod routing;

use builder::SelectiveTreeBuilder;
pub use primitive::{
    create_primitive_selective_tree_geometry, create_primitive_selective_tree_geometry_from_spec,
    PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
};

use super::super::types::Point2D;
use crate::domain::model::NetworkBlueprint;
use aequitas::systems::si::quantities::Length;

#[derive(Debug, Clone, Copy, PartialEq)]
/// Describes the center serpentine overlay used by a selective tree.
pub struct CenterSerpentinePathSpec {
    /// Number of segments in the serpentine path.
    pub segments: usize,
    /// Bend radius of the path in metres.
    pub bend_radius_m: Length<f64>,
    /// Waveform type for the serpentine path.  Defaults to Sine when
    /// constructed from legacy topology specs that omit this field.
    pub wave_type: crate::topology::SerpentineWaveType,
}

#[derive(Debug, Clone, Copy)]
struct PendingVenturiPath {
    channel_idx: usize,
    start: Option<Point2D>,
    end: Option<Point2D>,
    preferred_y: f64,
    fallback_length_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
/// Selects the topology-specific construction parameters for a selective tree.
pub enum SelectiveTreeTopology {
    /// Builds a cascade of center trifurcations.
    CascadeCenterTrifurcation {
        /// Number of trifurcation levels in the cascade.
        n_levels: usize,
        /// Fraction of each split assigned to the center branch.
        center_frac: f64,
        /// Enables Venturi treatment geometry on treatment branches.
        venturi_treatment_enabled: bool,
        /// Optional center serpentine overlay.
        center_serpentine: Option<CenterSerpentinePathSpec>,
    },
    /// Builds incremental filtration stages with terminal tri-bifurcation.
    IncrementalFiltrationTriBi {
        /// Number of preliminary trifurcation stages.
        n_pretri: usize,
        /// Center fraction used by preliminary trifurcations.
        pretri_center_frac: f64,
        /// Center fraction used by the terminal trifurcation.
        terminal_tri_center_frac: f64,
        /// Fraction assigned to the treatment branch at the bifurcation.
        bi_treat_frac: f64,
        /// Enables Venturi treatment geometry on treatment branches.
        venturi_treatment_enabled: bool,
        /// Optional center serpentine overlay.
        center_serpentine: Option<CenterSerpentinePathSpec>,
        /// Length of the outlet tail in metres.
        outlet_tail_length_m: Length<f64>,
    },
    /// Builds a tri-bifurcation-bifurcation-trifurcation selective topology.
    TriBiTriSelective {
        /// Center fraction at the first trifurcation.
        first_center_frac: f64,
        /// Fraction assigned to the treatment branch at the bifurcation.
        bi_treat_frac: f64,
        /// Center fraction at the second trifurcation.
        second_center_frac: f64,
    },
    /// Builds a double-trifurcation CIF topology with a center Venturi train.
    DoubleTrifurcationCif {
        /// Center fraction at the first trifurcation.
        split1_center_frac: f64,
        /// Center fraction at the second trifurcation.
        split2_center_frac: f64,
        /// Number of serial center throats.
        center_throat_count: u8,
        /// Spacing between center throats in metres.
        inter_throat_spacing_m: Length<f64>,
    },
}

#[derive(Debug, Clone, PartialEq)]
/// Defines the physical dimensions and topology of a selective-tree blueprint.
pub struct SelectiveTreeRequest {
    /// Human-readable name assigned to the generated blueprint.
    pub name: String,
    /// Overall rectangular envelope dimensions in metres.
    pub box_dims_m: (Length<f64>, Length<f64>),
    /// Length of the inlet trunk in metres.
    pub trunk_length_m: Length<f64>,
    /// Length of ordinary branches in metres.
    pub branch_length_m: Length<f64>,
    /// Length of hybrid branches in metres.
    pub hybrid_branch_length_m: Length<f64>,
    /// Width of the main channels in metres.
    pub main_width_m: Length<f64>,
    /// Width of Venturi throats in metres.
    pub throat_width_m: Length<f64>,
    /// Length of Venturi throats in metres.
    pub throat_length_m: Length<f64>,
    /// Channel height in metres.
    pub channel_height_m: Length<f64>,
    /// Topology-specific construction parameters.
    pub topology: SelectiveTreeTopology,
}

impl SelectiveTreeRequest {
    fn box_dims_mm(&self) -> (f64, f64) {
        (
            self.box_dims_m.0.into_base() * 1.0e3,
            self.box_dims_m.1.into_base() * 1.0e3,
        )
    }

    fn geometry(&self) -> SelectiveTreeGeometry {
        SelectiveTreeGeometry {
            box_dims_mm: self.box_dims_mm(),
            trunk_length_m: self.trunk_length_m.into_base(),
            branch_length_m: self.branch_length_m.into_base(),
            hybrid_branch_length_m: self.hybrid_branch_length_m.into_base(),
            main_width_m: self.main_width_m.into_base(),
            throat_width_m: self.throat_width_m.into_base(),
            throat_length_m: self.throat_length_m.into_base(),
            channel_height_m: self.channel_height_m.into_base(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct SelectiveTreeGeometry {
    box_dims_mm: (f64, f64),
    trunk_length_m: f64,
    branch_length_m: f64,
    hybrid_branch_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    channel_height_m: f64,
}

/// Generates a selective-tree blueprint from its physical request and topology.
pub fn create_selective_tree_geometry(request: &SelectiveTreeRequest) -> NetworkBlueprint {
    let geometry = request.geometry();
    let mut builder = SelectiveTreeBuilder::new(request.name.clone(), geometry.box_dims_mm);
    match &request.topology {
        SelectiveTreeTopology::CascadeCenterTrifurcation {
            n_levels,
            center_frac,
            venturi_treatment_enabled,
            center_serpentine,
        } => builder.build_cct(
            *n_levels,
            *center_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            &geometry,
        ),
        SelectiveTreeTopology::IncrementalFiltrationTriBi {
            n_pretri,
            pretri_center_frac,
            terminal_tri_center_frac,
            bi_treat_frac,
            venturi_treatment_enabled,
            center_serpentine,
            outlet_tail_length_m,
        } => builder.build_cif(
            *n_pretri,
            *pretri_center_frac,
            *terminal_tri_center_frac,
            *bi_treat_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            outlet_tail_length_m.into_base(),
            &geometry,
        ),
        SelectiveTreeTopology::TriBiTriSelective {
            first_center_frac,
            bi_treat_frac,
            second_center_frac,
        } => builder.build_tbt(
            *first_center_frac,
            *bi_treat_frac,
            *second_center_frac,
            &geometry,
        ),
        SelectiveTreeTopology::DoubleTrifurcationCif {
            split1_center_frac,
            split2_center_frac,
            center_throat_count,
            inter_throat_spacing_m,
        } => builder.build_dtcv(
            *split1_center_frac,
            *split2_center_frac,
            *center_throat_count,
            inter_throat_spacing_m.into_base(),
            &geometry,
        ),
    }
    builder.finish()
}

#[cfg(test)]
mod tests {
    use super::path_geometry::path_intersects_any;
    use super::routing::route_monotone_treatment_path;
    use super::{
        create_primitive_selective_tree_geometry, CenterSerpentinePathSpec,
        PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
    };
    use aequitas::systems::si::quantities::Length;

    #[test]
    fn primitive_selective_tree_annotation_preserves_positive_channel_lengths() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-lengths".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(8.0e-3),
            throat_width_m: Length::from_base(55.0e-6),
            throat_length_m: Length::from_base(110.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        let invalid_lengths: Vec<String> = blueprint
            .channels
            .iter()
            .filter(|channel| {
                let length_m = channel.length_m.into_base();
                !(length_m.is_finite() && length_m > 0.0)
            })
            .map(|channel| format!("{}={}", channel.id.as_str(), channel.length_m.into_base()))
            .collect();

        assert!(
            invalid_lengths.is_empty(),
            "primitive selective annotation should preserve positive channel lengths: {}",
            invalid_lengths.join(", ")
        );
    }

    #[test]
    fn primitive_selective_venturi_paths_have_no_unresolved_crossings() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-no-crossings".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(8.0e-3),
            throat_width_m: Length::from_base(55.0e-6),
            throat_length_m: Length::from_base(110.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: true,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
        blueprint
            .validate()
            .expect("venturi treatment paths must remain planar");

        for channel in blueprint.venturi_channels() {
            if let (Some(start), Some(end)) = (channel.path.first(), channel.path.last()) {
                if (start.1 - end.1).abs() < 1e-9 {
                    assert!(
                        channel
                            .path
                            .iter()
                            .all(|point| (point.1 - start.1).abs() < 1e-9),
                        "equal-y treatment channel {} must stay on its own lane",
                        channel.id.as_str()
                    );
                }
            }
        }
    }

    #[test]
    fn primitive_selective_penta_quad_tri_is_geometry_authorable() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-penta-quad-tri".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Penta,
                PrimitiveSelectiveSplitKind::Quad,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(4.0e-3),
            throat_width_m: Length::from_base(35.0e-6),
            throat_length_m: Length::from_base(180.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.45,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: true,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        let envelope_bound = 8.0 * f64::EPSILON * 128.0;
        assert!((blueprint.box_dims.0 - 127.76).abs() <= envelope_bound);
        assert!((blueprint.box_dims.1 - 85.47).abs() <= envelope_bound);
        assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
        blueprint
            .validate()
            .expect("Penta->Quad->Tri primitive selective tree should remain planar");
    }

    #[test]
    fn primitive_selective_tri_tri_retains_multiple_treatment_window_lanes() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-tritri-lanes".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(8.0e-3),
            throat_width_m: Length::from_base(55.0e-6),
            throat_length_m: Length::from_base(110.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let mut lane_keys = std::collections::BTreeSet::new();
        for channel in &blueprint.channels {
            if channel.therapy_zone
                != Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
            {
                continue;
            }
            let points = if channel.path.is_empty() {
                let start = blueprint
                    .nodes
                    .iter()
                    .find(|node| node.id == channel.from)
                    .map(|node| node.point);
                let end = blueprint
                    .nodes
                    .iter()
                    .find(|node| node.id == channel.to)
                    .map(|node| node.point);
                [start, end].into_iter().flatten().collect::<Vec<_>>()
            } else {
                channel.path.clone()
            };
            if points.is_empty() {
                continue;
            }
            let min_x = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
            let max_x = points
                .iter()
                .map(|(x, _)| *x)
                .fold(f64::NEG_INFINITY, f64::max);
            if max_x <= mid_x + 1.0e-6 || min_x < mid_x - 1.0e-6 {
                continue;
            }
            let mean_y = points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64;
            lane_keys.insert((mean_y * 100.0).round() as i64);
        }

        assert!(
            lane_keys.len() >= 3,
            "Tri->Tri selective trees must preserve multiple treatment-window lanes instead of collapsing to a centerline surrogate, got {lane_keys:?}"
        );
    }

    #[test]
    fn monotone_treatment_routing_preserves_equal_y_lane() {
        let routed =
            route_monotone_treatment_path(Some((10.0, 24.0)), Some((40.0, 24.0)), 24.0, 42.0, &[]);
        assert_eq!(routed, vec![(10.0, 24.0), (40.0, 24.0)]);
    }

    #[test]
    fn monotone_treatment_routing_doglegs_around_existing_branch() {
        let existing_paths = vec![vec![(3.0, 0.0), (7.0, 4.0)]];
        let routed = route_monotone_treatment_path(
            Some((0.0, 2.0)),
            Some((10.0, 2.0)),
            0.0,
            2.0,
            &existing_paths,
        );

        assert!(
            routed.len() > 2,
            "crossing direct path should reroute with a dogleg"
        );
        assert!(
            !path_intersects_any(&routed, &existing_paths),
            "rerouted path must not cross the existing treatment lane"
        );
    }

    #[test]
    fn primitive_selective_serpentine_mirrors_on_both_sides_of_midline() {
        use crate::domain::model::ChannelShape;
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "serp-mirror-check".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(4.0e-3),
            throat_width_m: Length::from_base(55.0e-6),
            throat_length_m: Length::from_base(110.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: Some(CenterSerpentinePathSpec {
                segments: 5,
                bend_radius_m: Length::from_base(3.0e-3),
                wave_type: crate::SerpentineWaveType::default(),
            }),
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let serpentine_channels: Vec<_> = blueprint
            .channels
            .iter()
            .filter(|ch| matches!(ch.channel_shape, ChannelShape::Serpentine { .. }))
            .collect();

        assert!(
            serpentine_channels.len() >= 2,
            "at least one serpentine channel on each side expected, got {}",
            serpentine_channels.len()
        );

        let node_pts: std::collections::HashMap<_, _> = blueprint
            .nodes
            .iter()
            .map(|n| (n.id.clone(), n.point))
            .collect();
        let left_count = serpentine_channels
            .iter()
            .filter(|ch| {
                let max_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::NEG_INFINITY, f64::max);
                max_x <= mid_x + 1e-6
            })
            .count();
        let right_count = serpentine_channels
            .iter()
            .filter(|ch| {
                let min_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::INFINITY, f64::min);
                min_x >= mid_x - 1e-6
            })
            .count();

        assert!(
            left_count >= 1,
            "split-side (left) treatment channels must have serpentine overlays, got {left_count}"
        );
        assert!(
            right_count >= 1,
            "merge-side (right) treatment channels must have serpentine overlays, got {right_count}"
        );
    }

    #[test]
    fn treatment_channels_retain_center_treatment_role_on_merge_side() {
        use crate::geometry::metadata::ChannelVisualRole;
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "role-mirror-check".to_string(),
            box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: Length::from_base(4.0e-3),
            throat_width_m: Length::from_base(55.0e-6),
            throat_length_m: Length::from_base(110.0e-6),
            channel_height_m: Length::from_base(1.0e-3),
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: Some(CenterSerpentinePathSpec {
                segments: 5,
                bend_radius_m: Length::from_base(3.0e-3),
                wave_type: crate::SerpentineWaveType::default(),
            }),
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let node_pts: std::collections::HashMap<_, _> = blueprint
            .nodes
            .iter()
            .map(|n| (n.id.clone(), n.point))
            .collect();

        let merge_side_treatment: Vec<_> = blueprint
            .channels
            .iter()
            .filter(|ch| {
                ch.therapy_zone == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
            })
            .filter(|ch| {
                let max_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::NEG_INFINITY, f64::max);
                max_x > mid_x + 1e-6
            })
            .collect();

        assert!(
            !merge_side_treatment.is_empty(),
            "should have treatment channels on the merge side"
        );
        for ch in &merge_side_treatment {
            assert_eq!(
                ch.visual_role,
                Some(ChannelVisualRole::CenterTreatment),
                "merge-side treatment channel {} must retain CenterTreatment role",
                ch.id.as_str()
            );
        }
    }
}
