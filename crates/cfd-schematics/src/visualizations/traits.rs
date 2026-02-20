//! visualizations/traits.rs - Visualization Abstraction Traits
//!
//! This module defines traits for visualization that abstract away the specific
//! plotting library implementation. This follows the Dependency Inversion Principle
//! by allowing the visualization logic to depend on abstractions rather than
//! concrete implementations.

use crate::error::VisualizationResult;
use crate::geometry::{ChannelSystem, Point2D};
use crate::visualizations::analysis_field::AnalysisOverlay;

/// Trait for rendering 2D microfluidic schematics
///
/// This trait abstracts the rendering backend, allowing different
/// plotting libraries to be used without changing the core visualization logic.
pub trait SchematicRenderer {
    /// Render a complete channel system to the specified output
    ///
    /// # Arguments
    /// * `system` - The channel system to render
    /// * `output_path` - Path where the rendered output should be saved
    /// * `config` - Configuration for the rendering
    ///
    /// # Returns
    /// Result indicating success or failure of the rendering operation
    fn render_system(
        &self,
        system: &ChannelSystem,
        output_path: &str,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    /// Render a channel system with CFD simulation results overlaid as colors
    ///
    /// Channels are colored by the field quantity in `overlay.edge_data` and
    /// nodes are colored by `overlay.node_data`. Normalization and colormap
    /// application are handled internally by the renderer.
    ///
    /// # Arguments
    /// * `system` - The channel system to render
    /// * `output_path` - Path where the rendered output should be saved
    /// * `config` - Configuration for the rendering
    /// * `overlay` - Typed CFD analysis overlay with field data and colormap
    fn render_analysis(
        &self,
        system: &ChannelSystem,
        output_path: &str,
        config: &RenderConfig,
        overlay: &AnalysisOverlay,
    ) -> VisualizationResult<()>;

    /// Get the supported output formats for this renderer
    fn supported_formats(&self) -> Vec<OutputFormat>;

    /// Validate that the output path has a supported format
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
            .map(|f| format!(".{}", f.extension()))
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

/// Trait for drawing basic geometric primitives
///
/// This trait provides a low-level interface for drawing operations
/// that can be implemented by different rendering backends.
pub trait GeometricDrawer {
    /// Draw a line between two points
    fn draw_line(
        &mut self,
        from: Point2D,
        to: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    /// Draw a series of connected line segments
    fn draw_path(&mut self, points: &[Point2D], style: &LineStyle) -> VisualizationResult<()>;

    /// Draw a rectangle outline
    fn draw_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    /// Fill a rectangle
    fn fill_rectangle(
        &mut self,
        top_left: Point2D,
        bottom_right: Point2D,
        color: &Color,
    ) -> VisualizationResult<()>;

    /// Draw text at a specific position
    fn draw_text(
        &mut self,
        position: Point2D,
        text: &str,
        style: &TextStyle,
    ) -> VisualizationResult<()>;
}

/// Configuration for rendering operations
#[derive(Debug, Clone)]
pub struct RenderConfig {
    /// Width of the output image in pixels
    pub width: u32,
    /// Height of the output image in pixels
    pub height: u32,
    /// Background color
    pub background_color: Color,
    /// Title to display on the schematic
    pub title: String,
    /// Whether to show coordinate axes
    pub show_axes: bool,
    /// Whether to show grid lines
    pub show_grid: bool,
    /// Margin around the content as a fraction of the total size
    pub margin_fraction: f64,
    /// Style for channel lines (default for all channel types)
    pub channel_style: LineStyle,
    /// Style for boundary box lines
    pub boundary_style: LineStyle,
    /// Style for axis labels
    pub axis_label_style: TextStyle,
    /// Style for the title
    pub title_style: TextStyle,
    /// Channel type-specific styling
    pub channel_type_styles: ChannelTypeStyles,
}

/// Channel type-specific styling configuration
///
/// This struct allows different visual styles for different channel types,
/// enabling easy identification of channel types in visualizations.
#[derive(Debug, Clone)]
pub struct ChannelTypeStyles {
    /// Style for straight channels (Straight, `SmoothStraight`)
    pub straight_style: LineStyle,
    /// Style for curved channels (Serpentine, Arc)
    pub curved_style: LineStyle,
    /// Style for tapered channels (Frustum)
    pub tapered_style: LineStyle,
}

impl Default for ChannelTypeStyles {
    fn default() -> Self {
        Self {
            straight_style: LineStyle::solid(Color::rgb(0, 0, 0), 1.0), // Black
            curved_style: LineStyle::solid(Color::rgb(0, 100, 200), 1.5), // Blue
            tapered_style: LineStyle::solid(Color::rgb(200, 50, 50), 2.0), // Red - distinctive for frustum
        }
    }
}

impl ChannelTypeStyles {
    /// Get the appropriate line style for a given channel type category
    #[must_use]
    pub const fn get_style(&self, category: crate::geometry::ChannelTypeCategory) -> &LineStyle {
        use crate::geometry::ChannelTypeCategory;
        match category {
            ChannelTypeCategory::Straight => &self.straight_style,
            ChannelTypeCategory::Curved => &self.curved_style,
            ChannelTypeCategory::Tapered => &self.tapered_style,
        }
    }
}

impl Default for RenderConfig {
    fn default() -> Self {
        Self {
            width: 1024,
            height: 768,
            background_color: Color::rgb(255, 255, 255), // White
            title: "Channel Schematic".to_string(),
            show_axes: true,
            show_grid: false,
            margin_fraction: 0.05,
            channel_style: LineStyle {
                color: Color::rgb(0, 0, 0), // Black
                width: 1.0,
                dash_pattern: None,
            },
            boundary_style: LineStyle {
                color: Color::rgb(0, 0, 0), // Black
                width: 2.0,
                dash_pattern: None,
            },
            axis_label_style: TextStyle {
                color: Color::rgb(0, 0, 0), // Black
                font_size: 12.0,
                font_family: "sans-serif".to_string(),
            },
            title_style: TextStyle {
                color: Color::rgb(0, 0, 0), // Black
                font_size: 24.0,
                font_family: "sans-serif".to_string(),
            },
            channel_type_styles: ChannelTypeStyles::default(),
        }
    }
}

/// Supported output formats for schematic rendering
///
/// Different formats have different characteristics:
/// - PNG: Raster format, good for web display and documentation
/// - SVG: Vector format, scalable and editable
/// - PDF: Vector format, good for publications and printing
/// - JPEG: Compressed raster format, smaller file sizes
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputFormat {
    /// Portable Network Graphics - raster format with lossless compression
    PNG,
    /// Scalable Vector Graphics - vector format for web and editing
    SVG,
    /// Portable Document Format - vector format for publications
    PDF,
    /// Joint Photographic Experts Group - compressed raster format
    JPEG,
}

impl OutputFormat {
    /// Get the file extension for this format
    #[must_use]
    pub const fn extension(&self) -> &'static str {
        match self {
            Self::PNG => "png",
            Self::SVG => "svg",
            Self::PDF => "pdf",
            Self::JPEG => "jpg",
        }
    }
}

/// Color representation using RGBA values
///
/// Each component is represented as a value from 0 to 255.
///
/// # Examples
///
/// ```rust
/// use scheme::visualizations::Color;
///
/// // Create a red color
/// let red = Color::rgba(255, 0, 0, 255);
///
/// // Use predefined colors
/// let white = Color::WHITE;
/// let black = Color::BLACK;
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Color {
    /// Red component (0-255)
    pub r: u8,
    /// Green component (0-255)
    pub g: u8,
    /// Blue component (0-255)
    pub b: u8,
    /// Alpha (transparency) component (0-255, where 255 is opaque)
    pub a: u8,
}

impl Color {
    /// Pure white color (255, 255, 255, 255)
    pub const WHITE: Self = Self {
        r: 255,
        g: 255,
        b: 255,
        a: 255,
    };
    /// Pure black color (0, 0, 0, 255)
    pub const BLACK: Self = Self {
        r: 0,
        g: 0,
        b: 0,
        a: 255,
    };
    /// Pure red color (255, 0, 0, 255)
    pub const RED: Self = Self {
        r: 255,
        g: 0,
        b: 0,
        a: 255,
    };
    /// Pure green color (0, 255, 0, 255)
    pub const GREEN: Self = Self {
        r: 0,
        g: 255,
        b: 0,
        a: 255,
    };
    /// Pure blue color (0, 0, 255, 255)
    pub const BLUE: Self = Self {
        r: 0,
        g: 0,
        b: 255,
        a: 255,
    };

    /// Create a new color with RGB values
    #[must_use]
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: 255 }
    }

    /// Create a new color with RGBA values
    #[must_use]
    pub const fn rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }
}

/// Style configuration for drawing lines
///
/// Defines the visual appearance of lines including color, width, and dash patterns.
///
/// # Examples
///
/// ```rust
/// use scheme::visualizations::{LineStyle, Color};
///
/// // Solid black line
/// let solid_line = LineStyle::solid(Color::BLACK, 1.0);
///
/// // Dashed red line
/// let dashed_line = LineStyle::dashed(Color::RED, 2.0, vec![5.0, 3.0]);
/// ```
#[derive(Debug, Clone)]
pub struct LineStyle {
    /// Color of the line
    pub color: Color,
    /// Width of the line in pixels
    pub width: f64,
    /// Optional dash pattern - alternating lengths of dashes and gaps
    pub dash_pattern: Option<Vec<f64>>,
}

impl LineStyle {
    /// Create a solid line style
    #[must_use]
    pub const fn solid(color: Color, width: f64) -> Self {
        Self {
            color,
            width,
            dash_pattern: None,
        }
    }

    /// Create a dashed line style
    #[must_use]
    pub const fn dashed(color: Color, width: f64, dash_pattern: Vec<f64>) -> Self {
        Self {
            color,
            width,
            dash_pattern: Some(dash_pattern),
        }
    }
}

/// Style configuration for drawing text
///
/// Defines the visual appearance of text including color, size, and font family.
///
/// # Examples
///
/// ```rust
/// use scheme::visualizations::{TextStyle, Color};
///
/// // Default text style
/// let text_style = TextStyle::new(Color::BLACK, 12.0, "Arial");
/// ```
#[derive(Debug, Clone)]
pub struct TextStyle {
    /// Color of the text
    pub color: Color,
    /// Font size in points
    pub font_size: f64,
    /// Font family name
    pub font_family: String,
}

impl TextStyle {
    /// Create a new text style
    #[must_use]
    pub fn new(color: Color, font_size: f64, font_family: &str) -> Self {
        Self {
            color,
            font_size,
            font_family: font_family.to_string(),
        }
    }
}

/// High-level visualization operations
///
/// This trait provides higher-level operations for visualizing
/// microfluidic systems, built on top of the basic drawing primitives.
pub trait VisualizationEngine {
    /// Visualize a complete channel system
    fn visualize_system(
        &mut self,
        system: &ChannelSystem,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    /// Visualize just the channels without the boundary box
    fn visualize_channels(
        &mut self,
        system: &ChannelSystem,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    /// Visualize the boundary box
    fn visualize_boundary(
        &mut self,
        system: &ChannelSystem,
        style: &LineStyle,
    ) -> VisualizationResult<()>;

    /// Add coordinate axes to the visualization
    fn add_axes(
        &mut self,
        system: &ChannelSystem,
        config: &RenderConfig,
    ) -> VisualizationResult<()>;

    /// Add a title to the visualization
    fn add_title(&mut self, title: &str, style: &TextStyle) -> VisualizationResult<()>;
}
