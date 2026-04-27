use plotters::coord::{types::RangedCoordf64, Shift};
use plotters::prelude::*;

use crate::domain::model::NetworkBlueprint;
use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::Point2D;
use crate::visualizations::traits::{
    Color as CfdColor, GeometricDrawer, LineStyle, RenderConfig, TextStyle, VisualizationEngine,
};

use super::convert_color;
use super::render_core::PlottersRenderer;

/// Plotters-based implementation of geometric drawing.
pub struct PlottersDrawer<'area, 'chart, DB: DrawingBackend> {
    drawing_area: &'area DrawingArea<DB, Shift>,
    chart: Option<&'chart mut ChartContext<'area, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>>,
}

impl<'area, 'chart, DB: DrawingBackend> PlottersDrawer<'area, 'chart, DB> {
    #[must_use]
    pub const fn new(
        drawing_area: &'area DrawingArea<DB, Shift>,
        chart: Option<
            &'chart mut ChartContext<'area, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
        >,
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

    fn chart_mut(
        &mut self,
    ) -> VisualizationResult<
        &mut ChartContext<'area, DB, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    > {
        self.chart
            .as_deref_mut()
            .ok_or_else(|| VisualizationError::CoordinateTransformError {
                message: "plotters drawer requires a chart context for schematic coordinates"
                    .to_string(),
            })
    }
}

impl<DB: DrawingBackend> GeometricDrawer for PlottersDrawer<'_, '_, DB> {
    fn draw_line(
        &mut self,
        from: Point2D,
        to: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        let stroke_width = style.width.max(1.0).round() as u32;
        self.chart_mut()?
            .draw_series(std::iter::once(PathElement::new(
                [from, to],
                convert_color(&style.color).stroke_width(stroke_width),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
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

        let stroke_width = style.width.max(1.0).round() as u32;
        self.chart_mut()?
            .draw_series(std::iter::once(PathElement::new(
                points,
                convert_color(&style.color).stroke_width(stroke_width),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        Ok(())
    }

    fn draw_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()> {
        let points = normalized_rectangle_outline(top_left, bottom_right);
        self.draw_path(&points, style)
    }

    fn fill_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        color: &CfdColor,
    ) -> VisualizationResult<()> {
        let (x0, y0, x1, y1) = normalized_rectangle_bounds(top_left, bottom_right);
        self.chart_mut()?
            .draw_series(std::iter::once(Rectangle::new(
                [(x0, y0), (x1, y1)],
                convert_color(color).filled(),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        Ok(())
    }

    fn draw_text(
        &mut self,
        position: Point2D,
        text: &str,
        style: &TextStyle,
    ) -> VisualizationResult<()> {
        let font_size = style.font_size.max(1.0).round() as i32;
        self.chart_mut()?
            .draw_series(std::iter::once(Text::new(
                text.to_string(),
                position,
                (style.font_family.as_str(), font_size)
                    .into_font()
                    .color(&convert_color(&style.color)),
            )))
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        Ok(())
    }
}

/// Plotters-based visualization engine.
pub struct PlottersVisualizationEngine<'area, 'chart, DB: DrawingBackend> {
    drawer: PlottersDrawer<'area, 'chart, DB>,
}

impl<'area, 'chart, DB: DrawingBackend> PlottersVisualizationEngine<'area, 'chart, DB> {
    #[must_use]
    pub const fn new(drawer: PlottersDrawer<'area, 'chart, DB>) -> Self {
        Self { drawer }
    }
}

impl<DB: DrawingBackend> VisualizationEngine for PlottersVisualizationEngine<'_, '_, DB> {
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
        config: &RenderConfig,
    ) -> VisualizationResult<()> {
        let mut mesh = self.drawer.chart_mut()?.configure_mesh();
        mesh.x_desc("X (mm)").y_desc("Y (mm)");
        if !config.show_grid {
            mesh.disable_mesh();
        }
        mesh.draw()
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))
    }

    fn add_title(&mut self, title: &str, style: &TextStyle) -> VisualizationResult<()> {
        let font_size = style.font_size.max(1.0).round() as i32;
        let title_style = (style.font_family.as_str(), font_size)
            .into_font()
            .color(&convert_color(&style.color));
        self.drawer
            .drawing_area()
            .titled(title, title_style)
            .map_err(|e| VisualizationError::rendering_error(&e.to_string()))?;
        Ok(())
    }
}

/// Convenience function to create a plotters renderer.
#[must_use]
pub const fn create_plotters_renderer() -> PlottersRenderer {
    PlottersRenderer
}

fn normalized_rectangle_bounds(top_left: Point2D, bottom_right: Point2D) -> (f64, f64, f64, f64) {
    let x0 = top_left.0.min(bottom_right.0);
    let y0 = top_left.1.min(bottom_right.1);
    let x1 = top_left.0.max(bottom_right.0);
    let y1 = top_left.1.max(bottom_right.1);
    (x0, y0, x1, y1)
}

fn normalized_rectangle_outline(top_left: Point2D, bottom_right: Point2D) -> [Point2D; 5] {
    let (x0, y0, x1, y1) = normalized_rectangle_bounds(top_left, bottom_right);
    [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
}
