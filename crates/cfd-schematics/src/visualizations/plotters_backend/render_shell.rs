use crate::config::ConstantsRegistry;
use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::ShellCuboid;
use crate::visualizations::traits::{OutputFormat, RenderConfig, SchematicRenderer};
use plotters::coord::Shift;
use plotters::prelude::*;
use plotters::style::Color as PlottersColor;

use super::convert_color;
use super::render_core::PlottersRenderer;

impl PlottersRenderer {
    /// Render a shell-cuboid schematic to file.
    pub fn render_shell_cuboid(
        &self,
        cuboid: &ShellCuboid,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        self.validate_output_path(output_path)?;
        match self.detect_output_format(output_path)? {
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

        let port_color = RGBColor(0, 170, 136);
        let wall_style = convert_color(&config.boundary_style.color)
            .stroke_width(config.boundary_style.width as u32);

        for (idx, (p1, p2)) in cuboid.box_outline.iter().enumerate() {
            if idx >= 8 {
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
                        wall_style,
                    )))
                    .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
            }
        }

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
                            hatch_stroke,
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

        ::tracing::info!("Shell cuboid schematic saved to {output_path}");
        Ok(())
    }
}

/// Convenience function: render a shell-cuboid with default config.
pub fn plot_shell_cuboid(cuboid: &ShellCuboid, output_path: &str) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let config = RenderConfig::default();
    renderer.render_shell_cuboid(cuboid, output_path, &config)
}
