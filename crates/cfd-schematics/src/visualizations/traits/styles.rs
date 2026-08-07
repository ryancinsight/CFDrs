use crate::visualizations::annotations::SchematicAnnotations;

/// Configuration for rendering operations.
#[derive(Debug, Clone)]
pub struct RenderConfig {
    /// Output width in pixels.
    pub width: u32,
    /// Output height in pixels.
    pub height: u32,
    /// Background color of the rendered image.
    pub background_color: Color,
    /// Title rendered above the schematic.
    pub title: String,
    /// Whether coordinate axes are drawn.
    pub show_axes: bool,
    /// Whether a background grid is drawn.
    pub show_grid: bool,
    /// Fraction of the canvas reserved as margin on each side.
    pub margin_fraction: f64,
    /// Style for channel paths.
    pub channel_style: LineStyle,
    /// Style for boundary paths.
    pub boundary_style: LineStyle,
    /// Style for axis tick labels.
    pub axis_label_style: TextStyle,
    /// Style for the schematic title.
    pub title_style: TextStyle,
    /// Style per channel geometry category.
    pub channel_type_styles: ChannelTypeStyles,
    /// Style per visual role when rendering selective trees.
    pub role_styles: Option<VisualRoleStyles>,
    /// Annotations overlaid on the schematic.
    pub annotations: Option<SchematicAnnotations>,
}

/// Channel type-specific styling configuration.
#[derive(Debug, Clone)]
pub struct ChannelTypeStyles {
    /// Style for straight channels.
    pub straight_style: LineStyle,
    /// Style for curved channels.
    pub curved_style: LineStyle,
    /// Style for tapered channels.
    pub tapered_style: LineStyle,
}

/// Role-based styling for selective-tree schematics.
#[derive(Debug, Clone)]
pub struct VisualRoleStyles {
    /// Style for trunk segments.
    pub trunk: LineStyle,
    /// Style for center-treatment segments.
    pub center_treatment: LineStyle,
    /// Style for peripheral-bypass segments.
    pub peripheral_bypass: LineStyle,
    /// Style for merge-collector segments.
    pub merge_collector: LineStyle,
    /// Style for venturi-throat segments.
    pub venturi_throat: LineStyle,
}

impl Default for VisualRoleStyles {
    fn default() -> Self {
        Self {
            trunk: LineStyle::solid(Color::rgb(60, 60, 60), 2.5),
            center_treatment: LineStyle::solid(Color::rgb(140, 50, 160), 1.5),
            peripheral_bypass: LineStyle::solid(Color::rgb(55, 120, 185), 1.5),
            merge_collector: LineStyle::solid(Color::rgb(40, 130, 130), 1.2),
            venturi_throat: LineStyle::solid(Color::rgb(140, 50, 160), 1.5),
        }
    }
}

impl VisualRoleStyles {
    /// Return the line style for the given visual role.
    #[must_use]
    pub fn get_style(&self, role: crate::geometry::metadata::ChannelVisualRole) -> &LineStyle {
        use crate::geometry::metadata::ChannelVisualRole;
        match role {
            ChannelVisualRole::Trunk
            | ChannelVisualRole::Diffuser
            | ChannelVisualRole::InternalLink => &self.trunk,
            ChannelVisualRole::CenterTreatment => &self.center_treatment,
            ChannelVisualRole::PeripheralBypass => &self.peripheral_bypass,
            ChannelVisualRole::MergeCollector => &self.merge_collector,
            ChannelVisualRole::VenturiThroat => &self.venturi_throat,
        }
    }
}

impl Default for ChannelTypeStyles {
    fn default() -> Self {
        Self {
            straight_style: LineStyle::solid(Color::rgb(0, 0, 0), 1.0),
            curved_style: LineStyle::solid(Color::rgb(0, 100, 200), 1.5),
            tapered_style: LineStyle::solid(Color::rgb(0, 0, 0), 1.0),
        }
    }
}

impl ChannelTypeStyles {
    /// Return the line style for the given channel category.
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
            background_color: Color::rgb(255, 255, 255),
            title: "Channel Schematic".to_string(),
            show_axes: true,
            show_grid: false,
            margin_fraction: 0.05,
            channel_style: LineStyle {
                color: Color::rgb(0, 0, 0),
                width: 1.0,
                dash_pattern: None,
            },
            boundary_style: LineStyle {
                color: Color::rgb(0, 0, 0),
                width: 2.0,
                dash_pattern: None,
            },
            axis_label_style: TextStyle {
                color: Color::rgb(0, 0, 0),
                font_size: 12.0,
                font_family: "sans-serif".to_string(),
            },
            title_style: TextStyle {
                color: Color::rgb(0, 0, 0),
                font_size: 24.0,
                font_family: "sans-serif".to_string(),
            },
            channel_type_styles: ChannelTypeStyles::default(),
            role_styles: None,
            annotations: None,
        }
    }
}

impl RenderConfig {
    /// Return a configuration sized for a 96-well plate schematic.
    #[must_use]
    pub fn well_plate_96() -> Self {
        Self {
            width: 1278,
            height: 855,
            title: "96-Well Plate Schematic (127.76 × 85.47 mm)".to_string(),
            ..Self::default()
        }
    }

    /// Return a 96-well plate configuration with report annotations and
    /// default role styling applied.
    #[must_use]
    pub fn well_plate_96_report_annotated() -> Self {
        let mut config = Self::well_plate_96();
        config.annotations = Some(SchematicAnnotations::report_default());
        config.role_styles = Some(VisualRoleStyles::default());
        config
    }
}

/// Supported output formats for schematic rendering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutputFormat {
    /// Portable Network Graphics.
    PNG,
    /// Scalable Vector Graphics.
    SVG,
    /// Portable Document Format.
    PDF,
    /// Joint Photographic Experts Group.
    JPEG,
}

impl OutputFormat {
    /// Return the file extension for this format, without a leading dot.
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

/// Color representation using RGBA values.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Color {
    /// Red channel, 0-255.
    pub r: u8,
    /// Green channel, 0-255.
    pub g: u8,
    /// Blue channel, 0-255.
    pub b: u8,
    /// Alpha channel, 0-255.
    pub a: u8,
}

impl Color {
    /// Opaque white.
    pub const WHITE: Self = Self {
        r: 255,
        g: 255,
        b: 255,
        a: 255,
    };
    /// Opaque black.
    pub const BLACK: Self = Self {
        r: 0,
        g: 0,
        b: 0,
        a: 255,
    };
    /// Opaque red.
    pub const RED: Self = Self {
        r: 255,
        g: 0,
        b: 0,
        a: 255,
    };
    /// Opaque green.
    pub const GREEN: Self = Self {
        r: 0,
        g: 255,
        b: 0,
        a: 255,
    };
    /// Opaque blue.
    pub const BLUE: Self = Self {
        r: 0,
        g: 0,
        b: 255,
        a: 255,
    };

    /// Construct an opaque color from RGB channels.
    #[must_use]
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: 255 }
    }

    /// Construct a color from RGBA channels.
    #[must_use]
    pub const fn rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }
}

impl From<iris::color::Rgba> for Color {
    fn from(color: iris::color::Rgba) -> Self {
        let [red, green, blue, alpha] = color.to_rgba8();
        Self::rgba(red, green, blue, alpha)
    }
}

/// Style configuration for drawing lines.
#[derive(Debug, Clone)]
pub struct LineStyle {
    /// Stroke color.
    pub color: Color,
    /// Stroke width in points.
    pub width: f64,
    /// Optional dash pattern: alternating on/off segment lengths.
    pub dash_pattern: Option<Vec<f64>>,
}

impl LineStyle {
    /// Construct a solid line style.
    #[must_use]
    pub const fn solid(color: Color, width: f64) -> Self {
        Self {
            color,
            width,
            dash_pattern: None,
        }
    }

    /// Construct a dashed line style with the given dash pattern.
    #[must_use]
    pub const fn dashed(color: Color, width: f64, dash_pattern: Vec<f64>) -> Self {
        Self {
            color,
            width,
            dash_pattern: Some(dash_pattern),
        }
    }
}

/// Style configuration for drawing text.
#[derive(Debug, Clone)]
pub struct TextStyle {
    /// Text color.
    pub color: Color,
    /// Font size in points.
    pub font_size: f64,
    /// Font family name.
    pub font_family: String,
}

impl TextStyle {
    /// Construct a text style.
    #[must_use]
    pub fn new(color: Color, font_size: f64, font_family: &str) -> Self {
        Self {
            color,
            font_size,
            font_family: font_family.to_string(),
        }
    }
}
