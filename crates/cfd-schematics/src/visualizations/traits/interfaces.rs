use crate::domain::model::NetworkBlueprint;
use crate::error::VisualizationResult;
use crate::geometry::Point2D;
use crate::visualizations::analysis_field::AnalysisOverlay;

use super::styles::{Color, LineStyle, OutputFormat, RenderConfig, TextStyle};

/// Trait for rendering 2D microfluidic schematics.
pub trait SchematicRenderer {
    /// Render a complete channel system to the specified output.
    fn render_system(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    /// Render a channel system with CFD simulation results overlaid as colors.
    fn render_analysis(
        &self,
        system: &NetworkBlueprint,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()>;

    /// Get the supported output formats for this renderer.
    fn supported_formats(&self) -> Vec<OutputFormat>;

    /// Validate that the output path has a supported format.
    fn validate_output_path(&self, path: &str) -> VisualizationResult<()> {
        let formats = self.supported_formats();
        let path_lower = path.to_lowercase();

        for format in &formats {
            if path_lower.ends_with(&format.extension()) {
                return Ok(());
            }
        }

        let supported_extensions: Vec<String> = formats
            .iter()
            .map(|format| format!(".{}", format.extension()))
            .collect();

        Err(crate::error::VisualizationError::invalid_output_path(
            path,
            &format!(
                "Unsupported format. Supported formats: {}",
                supported_extensions.join(", ")
            ),
        ))
    }
}

/// Trait for drawing basic geometric primitives.
pub trait GeometricDrawer {
    fn draw_line(
        &mut self,
        from: Point2D,
        to: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    fn draw_path(&mut self, points: &[Point2D], style: &LineStyle) -> VisualizationResult<()>;

    fn draw_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    fn fill_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        color: &Color,
    ) -> VisualizationResult<()>;

    fn draw_text(
        &mut self,
        position: Point2D,
        text: &str,
        style: &TextStyle,
    ) -> VisualizationResult<()>;
}

/// High-level visualization operations.
pub trait VisualizationEngine {
    fn visualize_system(
        &mut self,
        system: &NetworkBlueprint,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    fn visualize_channels(
        &mut self,
        system: &NetworkBlueprint,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    fn visualize_boundary(
        &mut self,
        system: &NetworkBlueprint,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    fn add_axes(
        &mut self,
        system: &NetworkBlueprint,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    fn add_title(&mut self, title: &str, style: &TextStyle) -> VisualizationResult<()>;
}