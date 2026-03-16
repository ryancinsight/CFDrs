//! Shared schematic annotation types and geometry utilities.

mod geometry;

use crate::geometry::Point2D;
use crate::visualizations::traits::Color;
use std::collections::BTreeMap;

pub use geometry::{
    center_biased_main_path, classify_node_roles, infer_terminal_nodes_by_x, node_degrees,
    project_markers_along_path, therapy_zone_presence, throat_count_from_blueprint_metadata,
};

/// Marker role for node/throat annotation points.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MarkerRole {
    Inlet,
    Outlet,
    Split,
    Merge,
    VenturiThroat,
    TherapyTarget,
    Bypass,
    Internal,
}

impl MarkerRole {
    #[must_use]
    pub const fn legend_label(self) -> &'static str {
        match self {
            Self::Inlet => "Inlet",
            Self::Outlet => "Outlet",
            Self::Split => "Split",
            Self::Merge => "Merge",
            Self::VenturiThroat => "Venturi throat",
            Self::TherapyTarget => "Therapy target",
            Self::Bypass => "Bypass",
            Self::Internal => "Internal",
        }
    }
}

/// Text label density control for annotations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelDensity {
    None,
    Critical,
    All,
}

/// Visual style controls for marker rendering.
#[derive(Debug, Clone)]
pub struct AnnotationStyle {
    /// Dot radius in pixels.
    pub dot_radius_px: i32,
    /// Marker label font size in points.
    pub label_font_size_pt: i32,
    /// Whether to render a role legend.
    pub show_legend: bool,
    /// Per-role marker colors.
    pub role_colors: BTreeMap<MarkerRole, Color>,
}

impl Default for AnnotationStyle {
    fn default() -> Self {
        let role_colors = BTreeMap::from([
            (MarkerRole::Inlet, Color::rgb(0x0B, 0x72, 0x85)),
            (MarkerRole::Outlet, Color::rgb(0xA6, 0x1E, 0x4D)),
            (MarkerRole::Split, Color::rgb(0x1D, 0x4E, 0xD8)),
            (MarkerRole::Merge, Color::rgb(0xD9, 0x77, 0x06)),
            (MarkerRole::VenturiThroat, Color::rgb(0xD2, 0x28, 0x28)),
            (MarkerRole::TherapyTarget, Color::rgb(0x0F, 0x76, 0x6E)),
            (MarkerRole::Bypass, Color::rgb(0x6B, 0x72, 0x80)),
            (MarkerRole::Internal, Color::rgb(0x37, 0x41, 0x51)),
        ]);
        Self {
            dot_radius_px: 5,
            label_font_size_pt: 11,
            show_legend: true,
            role_colors,
        }
    }
}

impl AnnotationStyle {
    /// Report-friendly style preset.
    #[must_use]
    pub fn report_default() -> Self {
        Self {
            dot_radius_px: 6,
            label_font_size_pt: 12,
            show_legend: true,
            ..Self::default()
        }
    }

    /// Resolve marker color for a role.
    #[must_use]
    pub fn color_for_role(&self, role: MarkerRole) -> Color {
        self.role_colors
            .get(&role)
            .cloned()
            .unwrap_or_else(|| Color::rgb(64, 64, 64))
    }
}

/// A single schematic marker.
#[derive(Debug, Clone)]
pub struct AnnotationMarker {
    pub point: Point2D,
    pub role: MarkerRole,
    pub label: Option<String>,
    /// True when this label should still be shown under `LabelDensity::Critical`.
    pub critical_label: bool,
}

impl AnnotationMarker {
    #[must_use]
    pub fn new(point: Point2D, role: MarkerRole) -> Self {
        Self {
            point,
            role,
            label: None,
            critical_label: false,
        }
    }

    #[must_use]
    pub fn with_label(mut self, label: impl Into<String>, critical_label: bool) -> Self {
        self.label = Some(label.into());
        self.critical_label = critical_label;
        self
    }
}

/// Full annotation payload for one render.
#[derive(Debug, Clone)]
pub struct SchematicAnnotations {
    pub markers: Vec<AnnotationMarker>,
    pub label_density: LabelDensity,
    pub style: AnnotationStyle,
    /// Optional free-text legend note (e.g. throat count summary).
    pub legend_note: Option<String>,
}

impl Default for SchematicAnnotations {
    fn default() -> Self {
        Self {
            markers: Vec::new(),
            label_density: LabelDensity::Critical,
            style: AnnotationStyle::default(),
            legend_note: None,
        }
    }
}

impl SchematicAnnotations {
    /// Preset aligned to milestone report schematics.
    #[must_use]
    pub fn report_default() -> Self {
        Self {
            style: AnnotationStyle::report_default(),
            ..Self::default()
        }
    }

    /// Validate annotation coordinates.
    pub fn validate(&self) -> Result<(), String> {
        if self.style.dot_radius_px <= 0 {
            return Err("dot radius must be positive".to_string());
        }
        if self.style.label_font_size_pt <= 0 {
            return Err("label font size must be positive".to_string());
        }
        for marker in &self.markers {
            if !marker.point.0.is_finite() || !marker.point.1.is_finite() {
                return Err("annotation marker coordinates must be finite".to_string());
            }
        }
        Ok(())
    }
}

/// Whether a label should be rendered for a marker under the selected density.
#[must_use]
pub fn should_render_label(marker: &AnnotationMarker, density: LabelDensity) -> bool {
    marker.label.as_ref().is_some_and(|_| match density {
        LabelDensity::None => false,
        LabelDensity::Critical => marker.critical_label,
        LabelDensity::All => true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label_density_critical_filters_dense_labels() {
        let critical = AnnotationMarker::new((0.0, 0.0), MarkerRole::Inlet).with_label("IN", true);
        let dense = AnnotationMarker::new((0.0, 0.0), MarkerRole::Internal).with_label("N1", false);

        assert!(should_render_label(&critical, LabelDensity::Critical));
        assert!(!should_render_label(&dense, LabelDensity::Critical));
        assert!(should_render_label(&dense, LabelDensity::All));
        assert!(!should_render_label(&critical, LabelDensity::None));
    }
}
