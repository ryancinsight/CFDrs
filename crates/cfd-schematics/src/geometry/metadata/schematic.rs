use serde::{Deserialize, Serialize};

/// Layout metadata positioning a node in schematic coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NodeLayoutMetadata {
    /// Horizontal position in millimetres.
    pub x_mm: f64,
    /// Vertical position in millimetres.
    pub y_mm: f64,
}

crate::impl_metadata!(NodeLayoutMetadata, "NodeLayoutMetadata");

/// Visual role of a channel for rendering and analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ChannelVisualRole {
    /// Inlet trunk segment.
    Trunk,
    /// Center treatment-path segment.
    CenterTreatment,
    /// Peripheral RBC-bypass segment.
    PeripheralBypass,
    /// Outlet merge collector segment.
    MergeCollector,
    /// Venturi throat segment.
    VenturiThroat,
    /// Diffuser segment adjacent to a throat.
    Diffuser,
    /// Internal link segment between stages.
    InternalLink,
}

/// Authored channel path metadata for rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChannelPathMetadata {
    /// Centerline polyline in millimetres.
    pub polyline_mm: Vec<(f64, f64)>,
    /// Visual role of the channel.
    pub visual_role: ChannelVisualRole,
}

crate::impl_metadata!(ChannelPathMetadata, "ChannelPathMetadata");

/// Family of a junction in the schematic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum JunctionFamily {
    /// Two-way split junction.
    Bifurcation,
    /// Three-way split junction.
    Trifurcation,
    /// Offset tee junction.
    Tee,
    /// Four-way crossing junction.
    Cross,
    /// Convergence junction.
    Merge,
}

/// Junction geometry metadata describing branch angles.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JunctionGeometryMetadata {
    /// Family of the junction.
    pub junction_family: JunctionFamily,
    /// Outgoing branch angles in degrees.
    pub branch_angles_deg: Vec<f64>,
    /// Incoming merge angles in degrees.
    pub merge_angles_deg: Vec<f64>,
}

crate::impl_metadata!(JunctionGeometryMetadata, "JunctionGeometryMetadata");
