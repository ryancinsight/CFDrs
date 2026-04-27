use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NodeLayoutMetadata {
    pub x_mm: f64,
    pub y_mm: f64,
}

crate::impl_metadata!(NodeLayoutMetadata, "NodeLayoutMetadata");

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChannelVisualRole {
    Trunk,
    CenterTreatment,
    PeripheralBypass,
    MergeCollector,
    VenturiThroat,
    Diffuser,
    InternalLink,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChannelPathMetadata {
    pub polyline_mm: Vec<(f64, f64)>,
    pub visual_role: ChannelVisualRole,
}

crate::impl_metadata!(ChannelPathMetadata, "ChannelPathMetadata");

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionFamily {
    Bifurcation,
    Trifurcation,
    Tee,
    Cross,
    Merge,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JunctionGeometryMetadata {
    pub junction_family: JunctionFamily,
    pub branch_angles_deg: Vec<f64>,
    pub merge_angles_deg: Vec<f64>,
}

crate::impl_metadata!(JunctionGeometryMetadata, "JunctionGeometryMetadata");
