//! Plotters backend split into core, annotation, and shell renderers.

mod render_annotations;
mod render_core;
mod render_shell;

use crate::visualizations::traits::Color;
use plotters::prelude::RGBColor;

pub use render_core::{
    create_plotters_renderer, PlottersDrawer, PlottersRenderer, PlottersVisualizationEngine,
};
pub use render_shell::plot_shell_cuboid;

pub(super) const fn convert_color(color: &Color) -> RGBColor {
    RGBColor(color.r, color.g, color.b)
}
