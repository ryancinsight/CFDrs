use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::ChannelSystem;
use plotters::coord::{types::RangedCoordf64, Shift};
use plotters::element::PathElement;
use plotters::prelude::*;

/// Legacy visualization function for backward compatibility
///
/// This function provides the original visualization interface using plotters directly.
/// For new code, prefer using the abstracted visualization system through the
/// `SchematicRenderer` trait.
///
/// # Arguments
///
/// * `system` - The channel system to visualize
/// * `output_path` - Path where the output image will be saved
/// * `title` - Title for the visualization
///
/// # Returns
///
/// Returns a tuple containing the drawing area, chart context, and scale factor
/// for further customization of the visualization.
///
/// # Examples
///
/// ```rust,no_run
/// use cfd_schematics::geometry::generator::create_geometry;
/// use cfd_schematics::geometry::SplitType;
/// use cfd_schematics::config::{GeometryConfig, ChannelTypeConfig};
/// use cfd_schematics::visualizations::shared_utilities::visualize;
///
/// let system = create_geometry(
///     (200.0, 100.0),
///     &[SplitType::Bifurcation],
///     &GeometryConfig::default(),
///     &ChannelTypeConfig::AllStraight,
/// );
///
/// let result = visualize(&system, "output.png", "My Schematic");
/// ```
#[allow(clippy::type_complexity)]
pub fn visualize<'a, 'b>(
    system: &'a ChannelSystem,
    output_path: &'b str,
    title: &str,
) -> VisualizationResult<(
    DrawingArea<BitMapBackend<'b>, Shift>,
    ChartContext<'a, BitMapBackend<'b>, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    f64,
)> {
    let root = BitMapBackend::new(output_path, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    let (length, width) = system.box_dims;
    let x_buffer = length * 0.05;
    let y_buffer = width * 0.05;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 30).into_font())
        .margin(20)
        .margin_right(150)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-x_buffer..length + x_buffer, -y_buffer..width + y_buffer)
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    chart
        .configure_mesh()
        .x_desc("X (mm)")
        .y_desc("Y (mm)")
        .draw()
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    let plotting_area = chart.plotting_area();
    let y_range = chart.y_range();
    let pixel_range = plotting_area.get_pixel_range();
    let y_pixel_height = pixel_range.1.end - pixel_range.1.start;
    let y_data_height = y_range.end - y_range.start;
    let y_scale_factor = f64::from(y_pixel_height) / y_data_height;

    chart
        .draw_series(
            system
                .box_outline
                .iter()
                .map(|&(p1, p2)| PathElement::new(vec![p1, p2], BLACK.stroke_width(2))),
        )
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    Ok((root, chart, y_scale_factor))
}
