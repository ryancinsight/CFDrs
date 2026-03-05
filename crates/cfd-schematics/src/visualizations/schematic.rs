use crate::error::VisualizationResult;
use crate::geometry::ChannelSystem;
use crate::visualizations::annotations::SchematicAnnotations;
use crate::visualizations::plotters_backend::PlottersRenderer;
use crate::visualizations::traits::{RenderConfig, SchematicRenderer};

/// Plot a channel system using the default renderer and configuration
///
/// This function provides backward compatibility while using the new
/// abstracted visualization system internally.
pub fn plot_geometry(system: &ChannelSystem, output_path: &str) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let config = RenderConfig::default();
    renderer.render_system(system, output_path, &config)
}

/// Plot a channel system with custom configuration
///
/// This function allows for more control over the rendering process
/// by accepting a custom configuration.
pub fn plot_geometry_with_config(
    system: &ChannelSystem,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    renderer.render_system(system, output_path, config)
}

/// Plot a channel system with explicit annotation payload.
///
/// This helper leaves the base config unchanged and overlays caller-provided
/// markers/labels through `RenderConfig.annotations` for this render only.
pub fn plot_geometry_with_annotations(
    system: &ChannelSystem,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let mut annotated_config = config.clone();
    annotated_config.annotations = Some(annotations.clone());
    renderer.render_system(system, output_path, &annotated_config)
}

/// Plot a channel system using a custom renderer
///
/// This function demonstrates the flexibility of the abstracted system
/// by allowing any renderer that implements `SchematicRenderer` to be used.
pub fn plot_geometry_with_renderer<R: SchematicRenderer>(
    system: &ChannelSystem,
    output_path: &str,
    renderer: &R,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    renderer.render_system(system, output_path, config)
}
