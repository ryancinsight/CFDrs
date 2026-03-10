use crate::config::ConstantsRegistry;
use crate::domain::model::{ChannelShape, NetworkBlueprint};
use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::{ChannelTypeCategory, Point2D};
use crate::visualizations::analysis_field::{
    colorize, AnalysisField, AnalysisOverlay, ColormapKind,
};
use crate::visualizations::traits::{
    Color as CfdColor, GeometricDrawer, LineStyle, OutputFormat, RenderConfig, SchematicRenderer,
    TextStyle, VisualizationEngine,
};
use plotters::coord::{types::RangedCoordf64, Shift};
use plotters::prelude::*;
use plotters::style::Color as PlottersColor;
use std::path::Path;

use super::convert_color;
use super::render_annotations::draw_annotation_overlay;
use crate::visualizations::schematic::channel_system_from_blueprint;

/// Plotters-based implementation of the schematic renderer.
pub struct PlottersRenderer;

impl SchematicRenderer for PlottersRenderer {
    fn render_system(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.render_analysis(
            system,
            output_path,
            config,
            &AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::BlueRed),
        )
    }

    fn render_analysis(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        if system.channels.is_empty() && system.nodes.is_empty() {
            return Err(VisualizationError::EmptyChannelSystem);
        }

        match self.detect_output_format(output_path)? {
            OutputFormat::PNG | OutputFormat::JPEG => {
                self.render_bitmap(system, output_path, config, overlay)
            }
            OutputFormat::SVG => self.render_svg(system, output_path, config, overlay),
            OutputFormat::PDF => Err(VisualizationError::UnsupportedFormat {
                format: "PDF".to_string(),
                message: "PDF output is not supported by the plotters backend".to_string(),
            }),
        }
    }

    fn supported_formats(&self) -> Vec<OutputFormat> {
        vec![OutputFormat::PNG, OutputFormat::JPEG, OutputFormat::SVG]
    }
}

impl PlottersRenderer {
    pub(super) fn detect_output_format(
        &self,
        output_path: &str,
    ) -> VisualizationResult<OutputFormat> {
        let path = Path::new(output_path);
        let extension = path
            .extension()
            .and_then(|ext| ext.to_str())
            .map(str::to_lowercase)
            .ok_or_else(|| {
                VisualizationError::invalid_output_path(
                    output_path,
                    "File must have a valid extension",
                )
            })?;

        match extension.as_str() {
            "png" => Ok(OutputFormat::PNG),
            "jpg" | "jpeg" => Ok(OutputFormat::JPEG),
            "svg" => Ok(OutputFormat::SVG),
            "pdf" => Ok(OutputFormat::PDF),
            _ => Err(VisualizationError::invalid_output_path(
                output_path,
                &format!("Unsupported file extension: .{extension}"),
            )),
        }
    }

    fn render_bitmap(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        let root =
            BitMapBackend::new(output_path, (config.width, config.height)).into_drawing_area();
        root.fill(&convert_color(&config.background_color))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        self.render_with_backend(system, config, root, output_path, overlay)
    }

    fn render_svg(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        let root = SVGBackend::new(output_path, (config.width, config.height)).into_drawing_area();
        root.fill(&convert_color(&config.background_color))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        self.render_with_backend(system, config, root, output_path, overlay)
    }

    fn render_with_backend<DB: DrawingBackend>(
        &self,
        system: &NetworkBlueprint,
        config: &RenderConfig,
        root: DrawingArea<DB, Shift>,
        output_path: &str,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        let renderable =
            channel_system_from_blueprint(system, Some(system.box_dims), Some(output_path))?;
        let (length, width) = system.box_dims;
        let x_buffer = length * config.margin_fraction;
        let y_buffer = width * config.margin_fraction;

        let constants = ConstantsRegistry::new();
        let mut chart = ChartBuilder::on(&root)
            .caption(
                &config.title,
                (
                    config.title_style.font_family.as_str(),
                    config.title_style.font_size as i32,
                ),
            )
            .margin(constants.get_default_chart_margin())
            .margin_right(constants.get_default_chart_right_margin())
            .x_label_area_size(constants.get_default_x_label_area_size())
            .y_label_area_size(constants.get_default_y_label_area_size())
            .build_cartesian_2d(-x_buffer..length + x_buffer, -y_buffer..width + y_buffer)
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        if config.show_axes {
            chart
                .configure_mesh()
                .x_desc("X (mm)")
                .y_desc("Y (mm)")
                .draw()
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }

        chart
            .draw_series(renderable.box_outline.iter().map(|(p1, p2)| {
                PathElement::new(
                    vec![*p1, *p2],
                    convert_color(&config.boundary_style.color)
                        .stroke_width(config.boundary_style.width as u32),
                )
            }))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        let (min_width, max_width) = {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for channel in &system.channels {
                let w = (channel.cross_section.dims().0 * 1e3 * 1e6).round() / 1e6; // to mm
                lo = lo.min(w);
                hi = hi.max(w);
            }
            (lo, hi)
        };

        const MAX_STROKE_MULT: u32 = 3;
        let width_multiplier = |w: f64| -> u32 {
            if max_width <= min_width {
                return 1;
            }
            let rounded = (w * 1e6).round() / 1e6;
            let ratio = (rounded - min_width) / (max_width - min_width);
            let mult = ratio.mul_add((MAX_STROKE_MULT - 1) as f64, 1.0);
            (mult.round() as u32).clamp(1, MAX_STROKE_MULT)
        };

        let has_edge_data = !overlay.edge_data.is_empty();
        for (i, channel) in system.channels.iter().enumerate() {
            let base_style = if has_edge_data {
                let color = overlay
                    .edge_color(i)
                    .unwrap_or_else(|| CfdColor::rgb(128, 128, 128));
                LineStyle::solid(color, config.channel_style.width)
            } else if let (Some(role), Some(role_styles)) =
                (channel.visual_role, config.role_styles.as_ref())
            {
                // Role-based styling: treatment vs. bypass vs. trunk get distinct
                // colors and base widths.  The width_multiplier (driven by the
                // physical cross-section) is applied on top so wider bypass
                // channels also render with a thicker stroke.
                role_styles.get_style(role).clone()
            } else {
                let category = renderable
                    .channel_categories
                    .get(i)
                    .copied()
                    .unwrap_or_else(|| match channel.channel_shape {
                        ChannelShape::Straight => ChannelTypeCategory::Straight,
                        ChannelShape::Serpentine { .. } => ChannelTypeCategory::Straight,
                    });
                config.channel_type_styles.get_style(category).clone()
            };

            let channel_width_mm = channel.cross_section.dims().0 * 1e3;
            let stroke = if max_width > min_width {
                (base_style.width * width_multiplier(channel_width_mm) as f64).round() as u32
            } else {
                base_style.width.round() as u32
            }
            .max(1);

            let path = renderable.channel_paths.get(i).cloned().unwrap_or_default();
            if path.len() < 2 {
                continue;
            }

            chart
                .draw_series(std::iter::once(PathElement::new(
                    path,
                    convert_color(&base_style.color).stroke_width(stroke),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }

        let has_node_data = !overlay.node_data.is_empty();
        if has_node_data {
            let node_radius = (length.min(width) * 0.02).max(0.5);
            for (i, node) in system.nodes.iter().enumerate() {
                let color = overlay
                    .node_color(i)
                    .unwrap_or_else(|| CfdColor::rgb(200, 200, 200));
                let (cx, cy) = node.point;
                let steps = 16usize;
                let circle_pts: Vec<Point2D> = (0..=steps)
                    .map(|i| {
                        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (steps as f64);
                        (
                            cx + node_radius * angle.cos(),
                            cy + node_radius * angle.sin(),
                        )
                    })
                    .collect();
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        circle_pts,
                        convert_color(&color).stroke_width(3),
                    )))
                    .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
            }
        }

        if let Some(annotations) = config.annotations.as_ref() {
            if let Err(message) = annotations.validate() {
                return Err(VisualizationError::InvalidParameters {
                    parameter: "annotations".to_string(),
                    value: "invalid annotation payload".to_string(),
                    constraint: message,
                });
            }
            draw_annotation_overlay(&mut chart, system, annotations)?;
        }

        if has_edge_data || has_node_data {
            draw_colorbar(&mut chart, system, x_buffer, y_buffer, overlay)?;
        }

        root.present()
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        ::tracing::info!("Schematic plot saved to {output_path}");
        Ok(())
    }
}

fn draw_colorbar<DB: DrawingBackend>(
    chart: &mut ChartContext<'_, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    system: &NetworkBlueprint,
    x_buffer: f64,
    y_buffer: f64,
    overlay: &AnalysisOverlay,
) -> VisualizationResult<()> {
    let (length, width) = system.box_dims;
    let has_edge_data = !overlay.edge_data.is_empty();

    let (min_val, max_val) = if has_edge_data {
        overlay.edge_range()
    } else {
        overlay.node_range()
    };

    let bar_steps = 20usize;
    let bar_x_start = length + x_buffer * 0.2;
    let bar_x_end = length + x_buffer * 0.7;
    let bar_y_start = y_buffer;
    let bar_y_end = width - y_buffer;
    let step_height = (bar_y_end - bar_y_start) / bar_steps as f64;

    for idx in 0..bar_steps {
        let t = 1.0 - (idx as f64) / (bar_steps as f64);
        let color = colorize(t, overlay.colormap);
        let y_top = bar_y_start + (idx as f64) * step_height;
        let y_bot = y_top + step_height;
        let rect_pts = vec![
            (bar_x_start, y_top),
            (bar_x_end, y_top),
            (bar_x_end, y_bot),
            (bar_x_start, y_bot),
            (bar_x_start, y_top),
        ];
        chart
            .draw_series(std::iter::once(PathElement::new(
                rect_pts,
                convert_color(&color).stroke_width(4),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
    }

    chart
        .draw_series(std::iter::once(Text::new(
            format!("{max_val:.2e}"),
            (bar_x_start, bar_y_start),
            ("sans-serif", 10).into_font(),
        )))
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    chart
        .draw_series(std::iter::once(Text::new(
            format!("{min_val:.2e}"),
            (bar_x_start, bar_y_end),
            ("sans-serif", 10).into_font(),
        )))
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    chart
        .draw_series(std::iter::once(Text::new(
            overlay.field.label().to_string(),
            (bar_x_start, f64::midpoint(bar_y_start, bar_y_end)),
            ("sans-serif", 9).into_font(),
        )))
        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

    Ok(())
}

/// Plotters-based implementation of geometric drawing.
pub struct PlottersDrawer<'a, DB: DrawingBackend> {
    drawing_area: &'a DrawingArea<DB, Shift>,
    chart: Option<&'a mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>>,
}

impl<'a, DB: DrawingBackend> PlottersDrawer<'a, DB> {
    #[must_use]
    pub const fn new(
        drawing_area: &'a DrawingArea<DB, Shift>,
        chart: Option<&'a mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>>,
    ) -> Self {
        Self {
            drawing_area,
            chart,
        }
    }

    #[must_use]
    pub const fn drawing_area(&self) -> &DrawingArea<DB, Shift> {
        self.drawing_area
    }
}

impl<DB: DrawingBackend> GeometricDrawer for PlottersDrawer<'_, DB> {
    fn draw_line(
        &mut self,
        from: Point2D,
        to: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        if let Some(chart) = &mut self.chart {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    vec![from, to],
                    convert_color(&style.color).stroke_width(style.width as u32),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }
        Ok(())
    }

    fn draw_path(&mut self, points: &[Point2D], style: &LineStyle) -> VisualizationResult<()> {
        if points.len() < 2 {
            return Err(VisualizationError::InvalidParameters {
                parameter: "points".to_string(),
                value: format!("{} points", points.len()),
                constraint: "Path must have at least 2 points".to_string(),
            });
        }

        if let Some(chart) = &mut self.chart {
            chart
                .draw_series(std::iter::once(PathElement::new(
                    points,
                    convert_color(&style.color).stroke_width(style.width as u32),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }
        Ok(())
    }

    fn draw_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        let points = vec![
            top_left,
            (bottom_right.0, top_left.1),
            bottom_right,
            (top_left.0, bottom_right.1),
            top_left,
        ];
        self.draw_path(&points, style)
    }

    fn fill_rectangle(
        &mut self,
        _top_left: Point2D,
        _bottom_right: Point2D,
        _color: &CfdColor,
    ) -> VisualizationResult<()> {
        Ok(())
    }

    fn draw_text(
        &mut self,
        _position: Point2D,
        _text: &str,
        _style: &TextStyle,
    ) -> VisualizationResult<()> {
        Ok(())
    }
}

/// Plotters-based visualization engine.
pub struct PlottersVisualizationEngine<'a, DB: DrawingBackend> {
    drawer: PlottersDrawer<'a, DB>,
}

impl<'a, DB: DrawingBackend> PlottersVisualizationEngine<'a, DB> {
    #[must_use]
    pub const fn new(drawer: PlottersDrawer<'a, DB>) -> Self {
        Self { drawer }
    }
}

impl<DB: DrawingBackend> VisualizationEngine for PlottersVisualizationEngine<'_, DB> {
    fn visualize_system(
        &mut self,
        system: &NetworkBlueprint,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.visualize_boundary(system, &config.boundary_style)?;
        self.visualize_channels(system, &config.channel_style)?;
        Ok(())
    }

    fn visualize_channels(
        &mut self,
        system: &NetworkBlueprint,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        for channel in &system.channels {
            if channel.path.len() >= 2 {
                self.drawer.draw_path(&channel.path, style)?;
            }
        }
        Ok(())
    }

    fn visualize_boundary(
        &mut self,
        system: &NetworkBlueprint,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        for &(p1, p2) in &system.box_outline {
            self.drawer.draw_line(p1, p2, style)?;
        }
        Ok(())
    }

    fn add_axes(
        &mut self,
        _system: &NetworkBlueprint,
        _config: &RenderConfig,
    ) -> VisualizationResult<()> {
        Ok(())
    }

    fn add_title(&mut self, _title: &str, _style: &TextStyle) -> VisualizationResult<()> {
        Ok(())
    }
}

/// Convenience function to create a plotters renderer.
#[must_use]
pub const fn create_plotters_renderer() -> PlottersRenderer {
    PlottersRenderer
}
