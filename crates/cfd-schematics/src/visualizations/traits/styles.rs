use crate::visualizations::annotations::SchematicAnnotations;

/// Configuration for rendering operations.
#[derive(Debug, Clone)]
pub struct RenderConfig {
    pub width: u32,
    pub height: u32,
    pub background_color: Color,
    pub title: String,
    pub show_axes: bool,
    pub show_grid: bool,
    pub margin_fraction: f64,
    pub channel_style: LineStyle,
    pub boundary_style: LineStyle,
    pub axis_label_style: TextStyle,
    pub title_style: TextStyle,
    pub channel_type_styles: ChannelTypeStyles,
    pub role_styles: Option<VisualRoleStyles>,
    pub annotations: Option<SchematicAnnotations>,
}

/// Channel type-specific styling configuration.
#[derive(Debug, Clone)]
pub struct ChannelTypeStyles {
    pub straight_style: LineStyle,
    pub curved_style: LineStyle,
    pub tapered_style: LineStyle,
}

/// Role-based styling for selective-tree schematics.
#[derive(Debug, Clone)]
pub struct VisualRoleStyles {
    pub trunk: LineStyle,
    pub center_treatment: LineStyle,
    pub peripheral_bypass: LineStyle,
    pub merge_collector: LineStyle,
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
    #[must_use]
    pub fn well_plate_96() -> Self {
        Self {
            width: 1278,
            height: 855,
            title: "96-Well Plate Schematic (127.76 × 85.47 mm)".to_string(),
            ..Self::default()
        }
    }

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
    PNG,
    SVG,
    PDF,
    JPEG,
}

impl OutputFormat {
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
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

impl Color {
    pub const WHITE: Self = Self {
        r: 255,
        g: 255,
        b: 255,
        a: 255,
    };
    pub const BLACK: Self = Self {
        r: 0,
        g: 0,
        b: 0,
        a: 255,
    };
    pub const RED: Self = Self {
        r: 255,
        g: 0,
        b: 0,
        a: 255,
    };
    pub const GREEN: Self = Self {
        r: 0,
        g: 255,
        b: 0,
        a: 255,
    };
    pub const BLUE: Self = Self {
        r: 0,
        g: 0,
        b: 255,
        a: 255,
    };

    #[must_use]
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: 255 }
    }

    #[must_use]
    pub const fn rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self { r, g, b, a }
    }
}

/// Style configuration for drawing lines.
#[derive(Debug, Clone)]
pub struct LineStyle {
    pub color: Color,
    pub width: f64,
    pub dash_pattern: Option<Vec<f64>>,
}

impl LineStyle {
    #[must_use]
    pub const fn solid(color: Color, width: f64) -> Self {
        Self {
            color,
            width,
            dash_pattern: None,
        }
    }

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
    pub color: Color,
    pub font_size: f64,
    pub font_family: String,
}

impl TextStyle {
    #[must_use]
    pub fn new(color: Color, font_size: f64, font_family: &str) -> Self {
        Self {
            color,
            font_size,
            font_family: font_family.to_string(),
        }
    }
}