use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind};
use crate::geometry::metadata::ChannelPathMetadata;
use crate::geometry::{ChannelTypeCategory, Point2D};
use std::collections::HashMap;

use super::layout::{blueprint_node_positions, save_auto_layout_json, BlueprintNodeLayout};
use super::path_generation::{
    generated_parallel_path, generated_serpentine_path, path_has_visible_serpentine_curvature,
};
use super::path_simplification::render_path_for_display;

pub(crate) struct RenderChannelSystem {
    pub box_outline: Vec<(Point2D, Point2D)>,
    pub channel_paths: Vec<Vec<Point2D>>,
    pub channel_categories: Vec<ChannelTypeCategory>,
}

pub(crate) fn channel_system_from_blueprint(
    blueprint: &NetworkBlueprint,
    box_dims_hint: Option<(f64, f64)>,
    output_path: Option<&str>,
) -> RenderChannelSystem {
    let box_dims = box_dims_hint.unwrap_or(blueprint.box_dims);
    let node_layout = blueprint_node_positions(blueprint, box_dims);
    if !node_layout.auto_layout_indices().is_empty() {
        if let Some(path) = output_path {
            save_auto_layout_json(&node_layout, box_dims, path);
        }
    }
    let mut channel_paths = Vec::with_capacity(blueprint.channels.len());
    let mut channel_categories = Vec::with_capacity(blueprint.channels.len());
    let resolved_paths = resolved_channel_paths(blueprint, &node_layout, box_dims);

    for (channel_spec, raw_path) in blueprint.channels.iter().zip(resolved_paths) {
        let category = channel_category_from_blueprint(channel_spec);
        channel_paths.push(render_path_for_display(&raw_path, category, box_dims));
        channel_categories.push(category);
    }

    RenderChannelSystem {
        box_outline: blueprint_box_outline(blueprint, box_dims),
        channel_paths,
        channel_categories,
    }
}

pub(super) fn resolved_channel_paths(
    blueprint: &NetworkBlueprint,
    node_layout: &BlueprintNodeLayout<'_>,
    box_dims: (f64, f64),
) -> Vec<Vec<Point2D>> {
    let parallel_groups = parallel_channel_groups(blueprint, node_layout);
    let default_point = (box_dims.0 * 0.5, box_dims.1 * 0.5);

    blueprint
        .channels
        .iter()
        .enumerate()
        .map(|(channel_idx, channel_spec)| {
            let from = node_layout
                .position(channel_spec.from.as_str())
                .unwrap_or(default_point);
            let to = node_layout
                .position(channel_spec.to.as_str())
                .unwrap_or(default_point);
            let key = match (
                node_layout.index_of(channel_spec.from.as_str()),
                node_layout.index_of(channel_spec.to.as_str()),
            ) {
                (Some(from_idx), Some(to_idx)) => Some((from_idx, to_idx)),
                _ => None,
            };
            let parallel_slot = key.and_then(|key| {
                parallel_groups.get(&key).and_then(|indices| {
                    indices
                        .iter()
                        .position(|idx| *idx == channel_idx)
                        .map(|slot| (slot, indices.len()))
                })
            });
            explicit_or_generated_path(channel_spec, from, to, parallel_slot, box_dims)
        })
        .collect()
}

fn parallel_channel_groups(
    blueprint: &NetworkBlueprint,
    node_layout: &BlueprintNodeLayout<'_>,
) -> HashMap<(usize, usize), Vec<usize>> {
    let mut groups: HashMap<(usize, usize), Vec<usize>> =
        HashMap::with_capacity(blueprint.channels.len());
    for (channel_idx, channel) in blueprint.channels.iter().enumerate() {
        let Some(from_idx) = node_layout.index_of(channel.from.as_str()) else {
            continue;
        };
        let Some(to_idx) = node_layout.index_of(channel.to.as_str()) else {
            continue;
        };
        groups
            .entry((from_idx, to_idx))
            .or_default()
            .push(channel_idx);
    }
    groups
}

fn blueprint_box_outline(
    blueprint: &NetworkBlueprint,
    box_dims: (f64, f64),
) -> Vec<(Point2D, Point2D)> {
    let mut outline = vec![
        ((0.0, 0.0), (box_dims.0, 0.0)),
        ((box_dims.0, 0.0), (box_dims.0, box_dims.1)),
        ((box_dims.0, box_dims.1), (0.0, box_dims.1)),
        ((0.0, box_dims.1), (0.0, 0.0)),
    ];

    for node_spec in &blueprint.nodes {
        match node_spec.kind {
            NodeKind::Inlet => {
                outline.push(((0.0, node_spec.point.1), node_spec.point));
            }
            NodeKind::Outlet => {
                outline.push((node_spec.point, (box_dims.0, node_spec.point.1)));
            }
            NodeKind::Reservoir | NodeKind::Junction => {}
        }
    }

    outline
}

fn explicit_or_generated_path(
    channel_spec: &ChannelSpec,
    from: Point2D,
    to: Point2D,
    parallel_slot: Option<(usize, usize)>,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if let ChannelShape::Serpentine {
        segments,
        bend_radius_m,
        ..
    } = channel_spec.channel_shape
    {
        if !channel_spec.path.is_empty() {
            return if path_has_visible_serpentine_curvature(
                &channel_spec.path,
                bend_radius_m * 1.0e3,
            ) {
                channel_spec.path.clone()
            } else {
                generated_serpentine_path(from, to, segments, bend_radius_m * 1.0e3)
            };
        }
    }

    if !channel_spec.path.is_empty() {
        return channel_spec.path.clone();
    }

    if let Some(path) = channel_spec
        .metadata
        .as_ref()
        .and_then(|meta| meta.get::<ChannelPathMetadata>())
        .map(|path| path.polyline_mm.clone())
    {
        return path;
    }

    if let Some((slot, count)) = parallel_slot.filter(|(_, count)| *count > 1) {
        return generated_parallel_path(from, to, slot, count, box_dims);
    }

    match channel_spec.channel_shape {
        ChannelShape::Serpentine {
            segments,
            bend_radius_m,
            ..
        } => generated_serpentine_path(from, to, segments, bend_radius_m * 1.0e3),
        ChannelShape::Straight => vec![from, to],
    }
}

pub(super) fn channel_category_from_blueprint(channel_spec: &ChannelSpec) -> ChannelTypeCategory {
    if channel_spec.venturi_geometry.is_some() {
        return ChannelTypeCategory::Tapered;
    }
    match channel_spec.channel_shape {
        ChannelShape::Serpentine { .. } => ChannelTypeCategory::Curved,
        ChannelShape::Straight => ChannelTypeCategory::Straight,
    }
}
