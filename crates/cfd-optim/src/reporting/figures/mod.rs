//! Dynamic figure generation, SVG builders, and rendering primitives for M12 reports.

mod manifest;
mod primitives;
pub(in crate::reporting) mod process;
mod svg;

pub use manifest::generate_m12_report_figures;
pub use manifest::{FigureGenerationInput, NarrativeFigureSpec};
