use std::collections::{HashMap, HashSet};

use super::super::super::super::types::Point2D;
use super::super::path_geometry::{infer_serpentine_shape, serpentine_overlay_path};
use super::super::routing::{
    channel_length_from_points_or_endpoints, channel_path_between, preferred_treatment_lane_y,
    primitive_treatment_leaf_indices, route_pending_venturi_paths,
};
use super::super::PendingVenturiPath;
use super::{PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest};
use crate::domain::model::{NetworkBlueprint, NodeId, NodeKind};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{
    ChannelVenturiSpec, ChannelVisualRole, JunctionFamily, JunctionGeometryMetadata,
    MetadataContainer, VenturiGeometryMetadata,
};

pub(super) fn annotate_primitive_tree(
    system: &mut NetworkBlueprint,
    request: &PrimitiveSelectiveTreeRequest,
) {
    if system.nodes.is_empty() || system.channels.is_empty() {
        return;
    }

    let mid_x = system.box_dims.0 * 0.5;
    let mid_y = system.box_dims.1 * 0.5;
    let inlet_node = system
        .nodes
        .iter()
        .min_by(|left, right| left.point.0.total_cmp(&right.point.0))
        .map(|node| node.id.clone());
    let outlet_node = system
        .nodes
        .iter()
        .max_by(|left, right| left.point.0.total_cmp(&right.point.0))
        .map(|node| node.id.clone());
    let (Some(inlet_node), Some(outlet_node)) = (inlet_node, outlet_node) else {
        return;
    };

    let mut adjacency: HashMap<NodeId, Vec<(NodeId, usize)>> = HashMap::new();
    for (idx, channel) in system.channels.iter().enumerate() {
        adjacency
            .entry(channel.from.clone())
            .or_default()
            .push((channel.to.clone(), idx));
        adjacency
            .entry(channel.to.clone())
            .or_default()
            .push((channel.from.clone(), idx));
    }

    let node_pts: HashMap<NodeId, Point2D> = system
        .nodes
        .iter()
        .map(|node| (node.id.clone(), node.point))
        .collect();

    let leaf_nodes: Vec<NodeId> = system
        .nodes
        .iter()
        .filter(|node| (node.point.0 - mid_x).abs() < 1e-6)
        .map(|node| node.id.clone())
        .collect();
    let mut ordered_leaves = leaf_nodes;
    ordered_leaves.sort_by(|left, right| node_pts[left].1.total_cmp(&node_pts[right].1));

    let treatment_leaf_indices = primitive_treatment_leaf_indices(
        ordered_leaves.len(),
        request.split_sequence.first().copied(),
    );
    let treatment_leaves: HashSet<NodeId> = treatment_leaf_indices
        .into_iter()
        .filter_map(|idx| ordered_leaves.get(idx).cloned())
        .collect();

    let mut treatment_channels = HashSet::new();
    for leaf in &treatment_leaves {
        if let Some(path) =
            channel_path_between(&inlet_node, leaf, &node_pts, &adjacency, Some(mid_x), true)
        {
            treatment_channels.extend(path);
        }
        if let Some(path) = channel_path_between(
            leaf,
            &outlet_node,
            &node_pts,
            &adjacency,
            Some(mid_x),
            false,
        ) {
            treatment_channels.extend(path);
        }
    }

    let first_tri = matches!(
        request.split_sequence.first(),
        Some(PrimitiveSelectiveSplitKind::Tri)
    );
    let bypass_width = if first_tri {
        (request.main_width_m * (1.0 - request.first_trifurcation_center_frac) * 0.5)
            .max(request.throat_width_m)
    } else {
        (request.main_width_m * (1.0 - request.bifurcation_treatment_frac))
            .max(request.throat_width_m)
    };
    let treatment_width = if first_tri {
        (request.main_width_m * request.first_trifurcation_center_frac).max(request.throat_width_m)
    } else {
        (request.main_width_m * request.bifurcation_treatment_frac).max(request.throat_width_m)
    };

    for node in &mut system.nodes {
        let degree = adjacency.get(&node.id).map_or(0, Vec::len);
        node.kind = if node.id == inlet_node {
            NodeKind::Inlet
        } else if node.id == outlet_node {
            NodeKind::Outlet
        } else {
            NodeKind::Junction
        };

        if degree >= 3 {
            let junction_family = if degree >= 4 {
                JunctionFamily::Trifurcation
            } else {
                JunctionFamily::Bifurcation
            };
            node.junction_geometry = Some(JunctionGeometryMetadata {
                junction_family,
                branch_angles_deg: Vec::new(),
                merge_angles_deg: Vec::new(),
            });
        }
    }

    let mut pending_venturi_paths = Vec::new();

    for (idx, channel) in system.channels.iter_mut().enumerate() {
        let points = channel.path.clone();
        let start_point = node_pts.get(&channel.from).copied();
        let end_point = node_pts.get(&channel.to).copied();
        let (min_x, max_x) = if points.is_empty() {
            let x0 = start_point.map_or(f64::INFINITY, |(x, _)| x);
            let x1 = end_point.map_or(f64::NEG_INFINITY, |(x, _)| x);
            (x0.min(x1), x0.max(x1))
        } else {
            let min_x = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
            let max_x = points
                .iter()
                .map(|(x, _)| *x)
                .fold(f64::NEG_INFINITY, f64::max);
            (min_x, max_x)
        };
        let avg_y = if points.is_empty() {
            mid_y
        } else {
            points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64
        };

        let is_treatment = treatment_channels.contains(&idx);
        let touches_inlet = channel.from == inlet_node || channel.to == inlet_node;
        let touches_outlet = channel.from == outlet_node || channel.to == outlet_node;
        let starts_at_treatment_leaf = treatment_leaves.contains(&channel.from);
        let is_treatment_window_channel =
            starts_at_treatment_leaf && max_x > mid_x + 1e-6 && min_x >= mid_x - 1e-6;
        let is_trunk = touches_inlet || touches_outlet;

        let physical_width = if touches_inlet || touches_outlet {
            request.main_width_m
        } else if is_treatment {
            treatment_width
        } else {
            bypass_width
        };

        channel.cross_section = crate::domain::model::CrossSectionSpec::Rectangular {
            width_m: physical_width,
            height_m: request.channel_height_m,
        };
        channel.length_m = channel_length_from_points_or_endpoints(
            &points,
            start_point,
            end_point,
            channel.length_m,
        );

        let mut role = if is_treatment {
            ChannelVisualRole::CenterTreatment
        } else {
            ChannelVisualRole::PeripheralBypass
        };
        if touches_inlet || touches_outlet {
            role = ChannelVisualRole::Trunk;
        } else if !is_treatment && max_x > mid_x {
            role = ChannelVisualRole::MergeCollector;
        }
        channel.visual_role = Some(role);
        channel.therapy_zone = Some(if is_trunk {
            TherapyZone::MixedFlow
        } else if is_treatment {
            TherapyZone::CancerTarget
        } else {
            TherapyZone::HealthyBypass
        });

        let should_overlay_serpentine = is_treatment
            && !is_trunk
            && request.center_serpentine.is_some()
            && (!request.treatment_branch_venturi_enabled || !is_treatment_window_channel);

        if should_overlay_serpentine {
            if let Some(spec) = request.center_serpentine {
                let source_points: Vec<Point2D> = if points.is_empty() {
                    [start_point, end_point].into_iter().flatten().collect()
                } else {
                    points.clone()
                };
                let serpentine_path = serpentine_overlay_path(
                    &source_points,
                    physical_width,
                    spec,
                    if request.treatment_branch_venturi_enabled {
                        1.1
                    } else {
                        1.5
                    },
                );
                channel.path = serpentine_path;
                channel.length_m = channel_length_from_points_or_endpoints(
                    &channel.path,
                    start_point,
                    end_point,
                    channel.length_m,
                );
                if let Some((start, end)) = channel.path.first().zip(channel.path.last()) {
                    channel.channel_shape =
                        infer_serpentine_shape(&channel.path, *start, *end, physical_width * 1.0e3);
                }
            }
        }

        if is_treatment && is_treatment_window_channel && request.treatment_branch_venturi_enabled {
            pending_venturi_paths.push(PendingVenturiPath {
                channel_idx: idx,
                start: points.first().copied().or(start_point),
                end: points.last().copied().or(end_point),
                preferred_y: preferred_treatment_lane_y(&points, start_point, end_point, avg_y),
                fallback_length_m: channel.length_m,
            });
            channel.visual_role = Some(ChannelVisualRole::VenturiThroat);
            channel.venturi_geometry = Some(VenturiGeometryMetadata {
                throat_width_m: request.throat_width_m,
                throat_height_m: request.channel_height_m,
                throat_length_m: request.throat_length_m,
                inlet_width_m: physical_width,
                outlet_width_m: physical_width,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 15.0,
                throat_position: 0.5,
            });
            let metadata = channel.metadata.get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: request.treatment_branch_throat_count.max(1),
                is_ctc_stream: true,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
            metadata.insert(
                channel
                    .venturi_geometry
                    .clone()
                    .expect("venturi geometry inserted"),
            );
        } else if !is_treatment {
            let metadata = channel.metadata.get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: 0,
                is_ctc_stream: false,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
        }
    }

    route_pending_venturi_paths(&mut system.channels, &pending_venturi_paths, mid_y);
}
