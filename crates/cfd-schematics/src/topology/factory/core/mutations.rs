use crate::topology::model::{
    SerpentineSpec, SplitKind, ThroatGeometrySpec, TreatmentActuationMode, VenturiPlacementMode,
    VenturiPlacementSpec,
};

/// High-level mutations available to the GA optimization engine.
#[derive(Debug, Clone, PartialEq)]
pub enum BlueprintTopologyMutation {
    /// Update the branch width at a stage.
    UpdateBranchWidth {
        /// Identifier of the topology stage.
        stage_id: String,
        /// Label of the branch being resized.
        branch_label: String,
        /// New branch width in meters.
        new_width_m: f64,
    },
    /// Replace the split kind of a stage.
    ReplaceSplitKind {
        /// Identifier of the topology stage.
        stage_id: String,
        /// Replacement split kind.
        split_kind: SplitKind,
    },
    /// Insert a new stage at the given index.
    InsertStage {
        /// Position at which the stage is inserted.
        stage_index: usize,
        /// Split kind of the inserted stage.
        split_kind: SplitKind,
    },
    /// Update the venturi configuration.
    UpdateVenturiConfiguration {
        /// Venturi placement specifications.
        placements: Vec<VenturiPlacementSpec>,
        /// Treatment actuation mode.
        treatment_mode: TreatmentActuationMode,
    },
    /// Set the treatment serpentine of a branch.
    SetTreatmentSerpentine {
        /// Identifier of the topology stage.
        stage_id: String,
        /// Label of the branch being modified.
        branch_label: String,
        /// Serpentine specification, or `None` to clear it.
        serpentine: Option<SerpentineSpec>,
    },
    /// Configure a treatment channel with serial venturi throats.
    SetTreatmentChannelVenturi {
        /// Identifier of the target channel.
        target_channel_id: String,
        /// Number of serial throat segments.
        serial_throat_count: u8,
        /// Throat geometry specification.
        throat_geometry: ThroatGeometrySpec,
        /// Placement mode for the throats.
        placement_mode: VenturiPlacementMode,
    },
    /// Set the treatment serpentine of a channel.
    SetTreatmentChannelSerpentine {
        /// Identifier of the target channel.
        target_channel_id: String,
        /// Serpentine specification, or `None` to clear it.
        serpentine: Option<SerpentineSpec>,
    },
    /// Insert a split-merge treatment structure into a channel.
    InsertTreatmentSplitMerge {
        /// Identifier of the target channel.
        target_channel_id: String,
        /// Split kind of the inserted structure.
        split_kind: SplitKind,
        /// Serpentine applied to the treatment structure, if any.
        treatment_serpentine: Option<SerpentineSpec>,
        /// Number of serial venturi throats, if any.
        venturi_serial_throat_count: Option<u8>,
        /// Throat geometry specification, if any.
        venturi_throat_geometry: Option<ThroatGeometrySpec>,
        /// Placement mode for the venturi throats.
        venturi_placement_mode: VenturiPlacementMode,
    },
}
