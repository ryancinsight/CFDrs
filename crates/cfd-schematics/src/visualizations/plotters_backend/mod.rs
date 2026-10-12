//! Plotters backend split into core, annotation, colorbar, drawer, and shell
//! renderers, with text laid out by an embedded pure-Rust typeface.

mod colorbar;
mod drawer;
mod render_annotations;
mod render_core;
mod render_shell;
mod typeface;

use crate::visualizations::traits::Color;
use plotters::prelude::RGBColor;

pub use drawer::{PlottersDrawer, PlottersVisualizationEngine, create_plotters_renderer};
pub use render_core::PlottersRenderer;
pub use render_shell::plot_shell_cuboid;

pub(super) const fn convert_color(color: &Color) -> RGBColor {
    RGBColor(color.r, color.g, color.b)
}
