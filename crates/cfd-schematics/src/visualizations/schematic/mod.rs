use crate::domain::model::NetworkBlueprint;
use crate::error::VisualizationResult;
use crate::geometry::metadata::NodeLayoutMetadata;
use crate::geometry::Point2D;
use crate::visualizations::annotations::SchematicAnnotations;
use crate::visualizations::plotters_backend::PlottersRenderer;
use crate::visualizations::traits::{RenderConfig, SchematicRenderer};

mod auto_annotations;
mod channel_system;
mod layout;
mod path_generation;
mod path_simplification;

use auto_annotations::build_auto_annotations;
use channel_system::resolved_channel_paths;
use layout::blueprint_node_positions;
pub(crate) use channel_system::channel_system_from_blueprint;

pub fn plot_geometry(blueprint: &NetworkBlueprint, output_path: &str) -> VisualizationResult<()> {
    plot_blueprint_auto_annotated(blueprint, output_path, &RenderConfig::default())
}

pub fn plot_geometry_with_config(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    renderer.render_system(blueprint, output_path, config)
}

pub fn plot_geometry_with_annotations(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let mut annotated_config = config.clone();
    annotated_config.annotations = Some(annotations.clone());
    renderer.render_system(blueprint, output_path, &annotated_config)
}

pub fn plot_blueprint(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    plot_geometry_with_config(blueprint, output_path, config)
}

pub fn plot_blueprint_with_annotations(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    plot_geometry_with_annotations(blueprint, output_path, config, annotations)
}

pub fn plot_blueprint_auto_annotated(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let annotations = build_auto_annotations(blueprint);
    plot_blueprint_with_annotations(blueprint, output_path, config, &annotations)
}

pub(crate) fn materialize_blueprint_layout(blueprint: &mut NetworkBlueprint) {
    let box_dims = blueprint.box_dims;
    let (node_points, auto_layout_used) = blueprint_node_positions(blueprint, box_dims);
    let channel_paths = resolved_channel_paths(blueprint, &node_points, box_dims);

    for node in &mut blueprint.nodes {
        let Some(&(x_mm, y_mm)) = node_points.get(node.id.as_str()) else {
            continue;
        };
        if auto_layout_used.contains_key(node.id.as_str()) {
            let layout = NodeLayoutMetadata { x_mm, y_mm };
            node.point = (x_mm, y_mm);
            node.layout = Some(layout);
            node.metadata
                .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
                .insert(layout);
        }
    }

    for (channel, path) in blueprint.channels.iter_mut().zip(channel_paths) {
        if channel.path.is_empty() && path.len() >= 2 {
            channel.path = path;
        }
    }
}

pub fn plot_geometry_auto_annotated(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    plot_blueprint_auto_annotated(blueprint, output_path, config)
}

pub fn plot_geometry_with_renderer<R: SchematicRenderer>(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    renderer: &R,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    renderer.render_system(blueprint, output_path, config)
}

/// Collect all centerline vertices from a rendered system, including routed
/// polyline bend points. Useful for report annotations that should expose the
/// actual routed path, not just graph-node endpoints.
#[must_use]
pub fn centerline_vertices(blueprint: &NetworkBlueprint) -> Vec<Point2D> {
    channel_system_from_blueprint(blueprint, Some(blueprint.box_dims), None)
        .channel_paths
        .into_iter()
        .flatten()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::channel_system::channel_category_from_blueprint;
    use super::path_generation::generated_serpentine_path;
    use super::path_simplification::render_path_for_display;
    use crate::domain::model::{ChannelShape, ChannelSpec};
    use crate::geometry::ChannelTypeCategory;

    #[test]
    fn generated_serpentine_path_has_mirrored_two_lobe_offsets() {
        let path = generated_serpentine_path((0.0, 0.0), (20.0, 0.0), 2, 2.0);
        let signed_offsets: Vec<f64> = path
            .iter()
            .skip(1)
            .take(path.len().saturating_sub(2))
            .map(|point| point.1)
            .filter(|offset| offset.abs() > 1.0e-6)
            .collect();
        assert!(
            signed_offsets.iter().any(|offset| *offset > 0.0),
            "serpentine path must lobe above the centerline"
        );
        assert!(
            signed_offsets.iter().any(|offset| *offset < 0.0),
            "serpentine path must lobe below the centerline"
        );
    }

    #[test]
    fn serpentine_channels_render_as_curved_paths() {
        let path = generated_serpentine_path((0.0, 0.0), (20.0, 0.0), 5, 2.5);
        let rendered = render_path_for_display(&path, ChannelTypeCategory::Curved, (40.0, 20.0));
        assert!(
            rendered.len() > path.len(),
            "curved serpentine display path should be spline-resampled"
        );
    }

    #[test]
    fn serpentine_channel_category_is_curved() {
        let mut channel = ChannelSpec::new_pipe_rect(
            "serp",
            "a".to_string(),
            "b".to_string(),
            1.0,
            1.0e-3,
            1.0e-3,
            0.0,
            0.0,
        );
        channel.channel_shape = ChannelShape::Serpentine {
            segments: 5,
            bend_radius_m: 2.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
        assert_eq!(
            channel_category_from_blueprint(&channel),
            ChannelTypeCategory::Curved
        );
    }
}