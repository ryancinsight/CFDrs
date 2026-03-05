//! `visualizations/plotters_backend.rs` - Plotters Implementation
//!
//! This module provides a concrete implementation of the visualization traits
//! using the plotters library. This demonstrates how the abstraction allows
//! for different rendering backends while maintaining the same interface.

use crate::config::ConstantsRegistry;
use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::{ChannelSystem, Point2D, ShellCuboid};
use crate::visualizations::analysis_field::{colorize, AnalysisOverlay, ColormapKind};
use crate::visualizations::traits::{
    Color, GeometricDrawer, LineStyle, OutputFormat, RenderConfig, SchematicRenderer, TextStyle,
    VisualizationEngine,
};
use plotters::coord::{types::RangedCoordf64, Shift};
use plotters::prelude::*;
use plotters::style::Color as PlottersColor;
use std::path::Path;

/// Plotters-based implementation of the schematic renderer
pub struct PlottersRenderer;

impl SchematicRenderer for PlottersRenderer {
    fn render_system(
        &self,
        system: &ChannelSystem,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.render_analysis(
            system,
            output_path,
            config,
            &AnalysisOverlay::new(
                crate::visualizations::analysis_field::AnalysisField::FlowRate,
                ColormapKind::BlueRed,
            ),
        )
    }

    fn render_analysis(
        &self,
        system: &ChannelSystem,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        // Validate input
        if system.channels.is_empty() && system.nodes.is_empty() {
            return Err(VisualizationError::EmptyChannelSystem);
        }

        self.validate_output_path(output_path)?;

        // Determine output format from file extension
        let output_format = self.detect_output_format(output_path)?;

        match output_format {
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
    /// Detect output format from file extension
    fn detect_output_format(&self, output_path: &str) -> VisualizationResult<OutputFormat> {
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

    /// Render to bitmap formats (PNG, JPEG)
    fn render_bitmap(
        &self,
        system: &ChannelSystem,
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

    /// Render to SVG format
    fn render_svg(
        &self,
        system: &ChannelSystem,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
        let root = SVGBackend::new(output_path, (config.width, config.height)).into_drawing_area();

        root.fill(&convert_color(&config.background_color))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        self.render_with_backend(system, config, root, output_path, overlay)
    }

    /// Common rendering logic for both bitmap and SVG backends
    ///
    /// Renders channels colored by `overlay.edge_data`, nodes colored by
    /// `overlay.node_data`, and a colorbar legend strip on the right margin.
    fn render_with_backend<DB: DrawingBackend>(
        &self,
        system: &ChannelSystem,
        config: &RenderConfig,
        root: DrawingArea<DB, Shift>,
        output_path: &str,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()> {
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

        // --- Boundary outline ---
        chart
            .draw_series(system.box_outline.iter().map(|(p1, p2)| {
                PathElement::new(
                    vec![*p1, *p2],
                    convert_color(&config.boundary_style.color)
                        .stroke_width(config.boundary_style.width as u32),
                )
            }))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        // --- Channels colored by overlay edge data ---
        // Collect distinct channel widths for proportional stroke rendering.
        // Channel physical width is mapped linearly to a stroke multiplier
        // in [1, MAX_STROKE_MULT], making wider channels visually thicker.
        let (min_width, max_width) = {
            let mut lo = f64::INFINITY;
            let mut hi = f64::NEG_INFINITY;
            for ch in &system.channels {
                let w = (ch.width * 1e6).round() / 1e6;
                if w < lo {
                    lo = w;
                }
                if w > hi {
                    hi = w;
                }
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
            // Linear map [0, 1] → [1, MAX_STROKE_MULT]
            let mult = ratio.mul_add((MAX_STROKE_MULT - 1) as f64, 1.0);
            (mult.round() as u32).clamp(1, MAX_STROKE_MULT)
        };

        let has_edge_data = !overlay.edge_data.is_empty();
        for channel in &system.channels {
            let base_style = if has_edge_data {
                // Use analysis color if data is present for this channel
                let color = overlay
                    .edge_color(channel.id)
                    .unwrap_or_else(|| Color::rgb(128, 128, 128));
                LineStyle::solid(color, config.channel_style.width)
            } else {
                // Fall back to type-based default style
                let category = crate::geometry::ChannelTypeCategory::from(&channel.channel_type);
                config.channel_type_styles.get_style(category).clone()
            };

            // Scale stroke width by proportional size when width range exists
            let stroke = if max_width > min_width {
                base_style.width as u32 * width_multiplier(channel.width)
            } else {
                base_style.width as u32
            };

            let path = crate::geometry::types::centerline_for_channel(channel, &system.nodes);
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

        // --- Nodes colored by overlay node data ---
        let has_node_data = !overlay.node_data.is_empty();
        if has_node_data {
            // Node radius in data-space units (2% of the shorter dimension)
            let node_radius = (length.min(width) * 0.02).max(0.5);
            for node in &system.nodes {
                let color = overlay
                    .node_color(node.id)
                    .unwrap_or_else(|| Color::rgb(200, 200, 200));
                let (cx, cy) = node.point;
                // Draw filled circle approximated as a small polygon
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

        // --- Colorbar legend (right margin, 10 gradient steps) ---
        if has_edge_data || has_node_data {
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

            for i in 0..bar_steps {
                // t=1 at top (high), t=0 at bottom (low)
                let t = 1.0 - (i as f64) / (bar_steps as f64);
                let color = colorize(t, overlay.colormap);
                let y_top = bar_y_start + (i as f64) * step_height;
                let y_bot = y_top + step_height;
                // Draw as a filled rectangle approximated by a closed path
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

            // Label min/max on colorbar
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
            // Field label
            chart
                .draw_series(std::iter::once(Text::new(
                    overlay.field.label().to_string(),
                    (bar_x_start, f64::midpoint(bar_y_start, bar_y_end)),
                    ("sans-serif", 9).into_font(),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }

        root.present()
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        println!("Schematic plot saved to {output_path}");
        Ok(())
    }
}

/// Plotters-based implementation of geometric drawing
pub struct PlottersDrawer<'a, DB: DrawingBackend> {
    drawing_area: &'a DrawingArea<DB, Shift>,
    chart: Option<&'a mut ChartContext<'a, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>>,
}

impl<'a, DB: DrawingBackend> PlottersDrawer<'a, DB> {
    /// Create a new `PlottersDrawer` with the given drawing area and optional chart context
    ///
    /// # Arguments
    ///
    /// * `drawing_area` - The plotters drawing area to draw on
    /// * `chart` - Optional chart context for coordinate transformations
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

    /// Access the underlying drawing area
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
            // Use slice directly instead of cloning to Vec - zero-copy optimization
            chart
                .draw_series(std::iter::once(PathElement::new(
                    points, // Pass slice directly instead of points.to_vec()
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
            top_left, // Close the rectangle
        ];
        self.draw_path(&points, style)
    }

    fn fill_rectangle(
        &mut self,
        _top_left: Point2D,
        _bottom_right: Point2D,
        _color: &Color,
    ) -> VisualizationResult<()> {
        // For now, just return Ok - filled rectangles would require more complex plotters usage
        Ok(())
    }

    fn draw_text(
        &mut self,
        _position: Point2D,
        _text: &str,
        _style: &TextStyle,
    ) -> VisualizationResult<()> {
        // For now, just return Ok - text drawing would require more complex plotters usage
        Ok(())
    }
}

/// Plotters-based implementation of the visualization engine
pub struct PlottersVisualizationEngine<'a, DB: DrawingBackend> {
    drawer: PlottersDrawer<'a, DB>,
}

impl<'a, DB: DrawingBackend> PlottersVisualizationEngine<'a, DB> {
    /// Create a new `PlottersVisualizationEngine` with the given drawer
    ///
    /// # Arguments
    ///
    /// * `drawer` - The `PlottersDrawer` to use for rendering operations
    #[must_use]
    pub const fn new(drawer: PlottersDrawer<'a, DB>) -> Self {
        Self { drawer }
    }
}

impl<DB: DrawingBackend> VisualizationEngine for PlottersVisualizationEngine<'_, DB> {
    fn visualize_system(
        &mut self,
        system: &ChannelSystem,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.visualize_boundary(system, &config.boundary_style)?;
        self.visualize_channels(system, &config.channel_style)?;
        Ok(())
    }

    fn visualize_channels(
        &mut self,
        system: &ChannelSystem,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        let lines = system.get_lines();
        for (p1, p2) in lines {
            self.drawer.draw_line(p1, p2, style)?;
        }
        Ok(())
    }

    fn visualize_boundary(
        &mut self,
        system: &ChannelSystem,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        for &(p1, p2) in &system.box_outline {
            self.drawer.draw_line(p1, p2, style)?;
        }
        Ok(())
    }

    fn add_axes(
        &mut self,
        _system: &ChannelSystem,
        _config: &RenderConfig,
    ) -> VisualizationResult<()> {
        // Axes are handled by the chart configuration in plotters
        Ok(())
    }

    fn add_title(&mut self, _title: &str, _style: &TextStyle) -> VisualizationResult<()> {
        // Title is handled by the chart configuration in plotters
        Ok(())
    }
}

/// Convert our Color type to plotters `RGBColor`
const fn convert_color(color: &Color) -> RGBColor {
    RGBColor(color.r, color.g, color.b)
}

/// Convenience function to create a plotters renderer
#[must_use]
pub const fn create_plotters_renderer() -> PlottersRenderer {
    PlottersRenderer
}

/// Convenience function for backward compatibility
pub fn plot_geometry_with_plotters(
    system: &ChannelSystem,
    output_path: &str,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let config = RenderConfig::default();
    renderer.render_system(system, output_path, &config)
}

impl PlottersRenderer {
    /// Render a [`ShellCuboid`] schematic to a file.
    ///
    /// Produces a 2D cross-section showing:
    /// - **Outer rectangle** (full device boundary)
    /// - **Inner rectangle** (cavity boundary, inset by the shell thickness)
    /// - **Inlet stub** — short horizontal line on the left wall (teal colour)
    /// - **Outlet stub** — short horizontal line on the right wall (teal colour)
    ///
    /// # Arguments
    ///
    /// * `cuboid`      — the shell geometry to render.
    /// * `output_path` — destination file path (`.svg`, `.png`, or `.jpg`).
    /// * `config`      — render configuration (title, colours, dimensions).
    ///
    /// # Errors
    ///
    /// Returns [`VisualizationError`] if the output path is invalid,
    /// the format is unsupported, or rendering fails.
    pub fn render_shell_cuboid(
        &self,
        cuboid: &ShellCuboid,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.validate_output_path(output_path)?;
        let output_format = self.detect_output_format(output_path)?;

        match output_format {
            OutputFormat::PNG | OutputFormat::JPEG => {
                self.render_shell_bitmap(cuboid, output_path, config)
            }
            OutputFormat::SVG => self.render_shell_svg(cuboid, output_path, config),
            OutputFormat::PDF => Err(VisualizationError::UnsupportedFormat {
                format: "PDF".to_string(),
                message: "PDF output is not supported by the plotters backend".to_string(),
            }),
        }
    }

    fn render_shell_bitmap(
        &self,
        cuboid: &ShellCuboid,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        let root =
            BitMapBackend::new(output_path, (config.width, config.height)).into_drawing_area();
        root.fill(&convert_color(&config.background_color))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        self.draw_shell_on_backend(cuboid, config, root, output_path)
    }

    fn render_shell_svg(
        &self,
        cuboid: &ShellCuboid,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        let root = SVGBackend::new(output_path, (config.width, config.height)).into_drawing_area();
        root.fill(&convert_color(&config.background_color))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        self.draw_shell_on_backend(cuboid, config, root, output_path)
    }

    fn draw_shell_on_backend<DB: DrawingBackend>(
        &self,
        cuboid: &ShellCuboid,
        config: &RenderConfig,
        root: DrawingArea<DB, Shift>,
        output_path: &str,
    ) -> VisualizationResult<()> {
        let (w, h) = cuboid.outer_dims;
        let x_buffer = w * config.margin_fraction;
        let y_buffer = h * config.margin_fraction;

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
            .build_cartesian_2d(-x_buffer..w + x_buffer, -y_buffer..h + y_buffer)
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        if config.show_axes {
            chart
                .configure_mesh()
                .x_desc("X (mm)")
                .y_desc("Y (mm)")
                .draw()
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }

        // Port stub colour: teal, visually distinct from wall outline.
        let port_color = RGBColor(0, 170, 136);
        let wall_style = convert_color(&config.boundary_style.color)
            .stroke_width(config.boundary_style.width as u32);

        // The first 8 segments are outer rect (4) + inner rect (4) → wall colour.
        // The last 2 segments are port stubs → teal.
        for (idx, (p1, p2)) in cuboid.box_outline.iter().enumerate() {
            let is_port_stub = idx >= 8;
            if is_port_stub {
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        vec![*p1, *p2],
                        port_color.stroke_width(config.boundary_style.width.max(3.0) as u32),
                    )))
                    .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
            } else {
                chart
                    .draw_series(std::iter::once(PathElement::new(
                        vec![*p1, *p2],
                        wall_style.clone(),
                    )))
                    .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
            }
        }

        // Label port midpoints.
        let (_, h_mid) = (0.0_f64, h / 2.0);
        let label_offset = x_buffer * 0.15;
        chart
            .draw_series(std::iter::once(Text::new(
                "IN",
                (0.0 - label_offset, h_mid),
                ("sans-serif", 11).into_font(),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        chart
            .draw_series(std::iter::once(Text::new(
                "OUT",
                (w + label_offset * 0.1, h_mid),
                ("sans-serif", 11).into_font(),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        // Draw TPMS fill hatch pattern inside inner rectangle if specified.
        if let Some(ref fill) = cuboid.tpms_fill {
            let t = cuboid.shell_thickness_mm;
            let (iw, ih) = cuboid.inner_dims;
            let hatch_color = RGBColor(100, 130, 200);
            let hatch_stroke = hatch_color.stroke_width(1);
            let spacing = fill.period_mm;
            let diag = iw + ih;
            let mut offset = -ih;
            while offset < diag {
                let x_start = (t + offset).max(t);
                let y_start = (t + (x_start - t - offset)).max(t);
                let x_end = (t + offset + ih).min(t + iw);
                let y_end = (t + (x_end - t - offset)).min(t + ih);
                if x_start < t + iw && x_end > t && y_start < t + ih && y_end > t {
                    chart
                        .draw_series(std::iter::once(PathElement::new(
                            vec![(x_start, y_start), (x_end, y_end)],
                            hatch_stroke.clone(),
                        )))
                        .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
                }
                offset += spacing;
            }
            let label = format!("{} λ={:.1}mm", fill.surface.label(), fill.period_mm);
            let cx = t + iw / 2.0;
            let cy = t + ih / 2.0;
            chart
                .draw_series(std::iter::once(Text::new(
                    label,
                    (cx - iw * 0.15, cy),
                    ("sans-serif", 13).into_font().color(&RGBColor(40, 60, 120)),
                )))
                .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        }

        root.present()
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;

        println!("Shell cuboid schematic saved to {output_path}");
        Ok(())
    }
}

/// Convenience function: render a [`ShellCuboid`] with default config.
pub fn plot_shell_cuboid(cuboid: &ShellCuboid, output_path: &str) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let config = RenderConfig::default();
    renderer.render_shell_cuboid(cuboid, output_path, &config)
}
