use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlueprintRenderHints {
    pub stage_sequence: String,
    pub split_layers: usize,
    pub throat_count_hint: usize,
    pub treatment_label: String,
    #[serde(default)]
    pub mirror_x: bool,
    #[serde(default)]
    pub mirror_y: bool,
}

crate::impl_metadata!(BlueprintRenderHints, "BlueprintRenderHints");

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GeometryAuthoringProvenance {
    pub source: String,
}

impl GeometryAuthoringProvenance {
    #[must_use]
    pub fn create_geometry() -> Self {
        Self {
            source: "create_geometry".to_string(),
        }
    }

    #[must_use]
    pub fn selective_wrapper() -> Self {
        Self {
            source: "create_geometry::selective_wrapper".to_string(),
        }
    }
}

crate::impl_metadata!(GeometryAuthoringProvenance, "GeometryAuthoringProvenance");
