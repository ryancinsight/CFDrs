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
