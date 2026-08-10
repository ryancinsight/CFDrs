use serde::{Deserialize, Serialize};

/// Render hints consumed by visualizers and renderers.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlueprintRenderHints {
    /// Stage sequence label.
    pub stage_sequence: String,
    /// Number of visible split layers.
    pub split_layers: usize,
    /// Venturi throat count hint.
    pub throat_count_hint: usize,
    /// Treatment label, e.g. "venturi" or "ultrasound".
    pub treatment_label: String,
    /// Mirror the geometry about the vertical axis.
    #[serde(default)]
    pub mirror_x: bool,
    /// Mirror the geometry about the horizontal axis.
    #[serde(default)]
    pub mirror_y: bool,
}

crate::impl_metadata!(BlueprintRenderHints, "BlueprintRenderHints");

/// Provenance marking a blueprint as authored by the canonical geometry pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeometryAuthoringProvenance {
    /// Provenance source identifier.
    pub source: String,
}

impl GeometryAuthoringProvenance {
    /// Provenance for the canonical `create_geometry` pipeline.
    #[must_use]
    pub fn create_geometry() -> Self {
        Self {
            source: "create_geometry".to_string(),
        }
    }

    /// Provenance for the selective-tree wrapper around `create_geometry`.
    #[must_use]
    pub fn selective_wrapper() -> Self {
        Self {
            source: "create_geometry::selective_wrapper".to_string(),
        }
    }
}

crate::impl_metadata!(GeometryAuthoringProvenance, "GeometryAuthoringProvenance");
