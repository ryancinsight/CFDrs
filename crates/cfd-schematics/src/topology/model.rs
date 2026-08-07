//! Typed topology specifications shared by schematic builders and queries.

use crate::domain::therapy_metadata::TherapyZone;
use aequitas::systems::si::quantities::{Angle, Dimensionless, Length};
use serde::{Deserialize, Serialize};

/// Describes how a topology stage branches a channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SplitKind {
    /// Split one channel into the requested number of branches.
    NFurcation(usize),
}

impl SplitKind {
    /// Returns the number of branches described by this split.
    #[must_use]
    pub const fn branch_count(self) -> usize {
        match self {
            Self::NFurcation(n) => n,
        }
    }
}

/// Physiological role assigned to a branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BranchRole {
    /// Branch carrying the treatment flow.
    Treatment,
    /// Branch collecting white blood cells.
    WbcCollection,
    /// Branch bypassing red blood cells.
    RbcBypass,
    /// Branch with no specialized treatment role.
    Neutral,
}

/// Treatment actuation used by a topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TreatmentActuationMode {
    /// Apply ultrasound without a Venturi stage.
    UltrasoundOnly,
    /// Apply ultrasound together with Venturi cavitation.
    VenturiCavitation,
}

/// Placement rule for a Venturi stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VenturiPlacementMode {
    /// Place the Venturi in a straight channel segment.
    StraightSegment,
    /// Place the Venturi at the peak Dean-number site of a curve.
    CurvaturePeakDeanNumber,
    /// Place the Venturi at a diffuser shoulder.
    DiffuserShoulder,
}

/// Serpentine waveform type for the 1D physics model.
///
/// Controls both the 2D rendering path shape and the 1D solver's bend loss
/// model.  Different waveforms produce distinct curvature profiles that
/// affect Dean secondary flow intensity and minor loss coefficients.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SerpentineWaveType {
    /// Smooth sinusoidal bends.  Curvature is distributed continuously
    /// along each half-period.  Produces moderate Dean numbers with
    /// gradual onset/decay.
    #[default]
    Sine,
    /// Near-square-wave U-turns (tanh-smoothed).  Curvature is
    /// concentrated at the direction reversals with near-constant-velocity
    /// straight segments between them.  Highest K-factor minor losses.
    Square,
    /// Triangular (zigzag) path.  Linear ramps between sharp apices.
    /// Curvature is zero along the ramps and singular at the apex
    /// points, producing very localized Dean vortex peaks.
    Triangular,
}

/// Geometric parameters for a serpentine route.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SerpentineSpec {
    /// Number of repeated channel segments.
    pub segments: usize,
    /// Radius of each bend.
    pub bend_radius_m: Length<f64>,
    /// Length of each straight segment.
    pub segment_length_m: Length<f64>,
    /// Waveform type controlling bend geometry and 1D loss model.
    /// Defaults to `Sine` for backward compatibility with existing
    /// topology specs that omit this field.
    #[serde(default)]
    pub wave_type: SerpentineWaveType,
}

/// Geometric and therapy metadata for one channel route.
/// Authoring specification for a flow channel between two endpoints.
///
/// Composes a rectangular duct cross-section with an optional serpentine
/// waveform and the therapy-zone semantic tag driving mesh refinement and
/// physics coupling downstream of topology authoring.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChannelRouteSpec {
    /// Total route length.
    pub length_m: Length<f64>,
    /// Route width.
    pub width_m: Length<f64>,
    /// Route height.
    pub height_m: Length<f64>,
    /// Optional serpentine geometry applied to the route.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub serpentine: Option<SerpentineSpec>,
    /// Therapy zone served by the route.
    #[serde(default = "default_therapy_zone")]
    pub therapy_zone: TherapyZone,
}

fn default_therapy_zone() -> TherapyZone {
    TherapyZone::MixedFlow
}

/// Geometry and label for one recovery branch.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SubBranchSpec {
    /// Stable label used to derive the channel identifier.
    pub label: String,
    /// Sub-branch width.
    pub width_m: Length<f64>,
    /// Sub-branch height.
    pub height_m: Length<f64>,
}

/// Secondary split used to recover a branch of the flow.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RecoverySubSplit {
    /// Branches created by the recovery split.
    pub sub_branches: Vec<SubBranchSpec>,
    /// Index of the branch carrying the recovery flow.
    pub recovery_arm_index: usize,
}

/// Complete specification of one topology branch.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BranchSpec {
    /// Stable branch label used to derive the channel identifier.
    pub label: String,
    /// Physiological role of the branch.
    pub role: BranchRole,
    /// Whether this branch is part of the treatment path.
    pub treatment_path: bool,
    /// Geometry and therapy metadata for the branch route.
    pub route: ChannelRouteSpec,
    /// Optional second split used to recover the bypass flow.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recovery_sub_split: Option<RecoverySubSplit>,
}

/// One ordered branching stage in a topology.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SplitStageSpec {
    /// Stable identifier for the split stage.
    pub stage_id: String,
    /// Branching operation performed at this stage.
    pub split_kind: SplitKind,
    /// Branch definitions emitted by this stage.
    pub branches: Vec<BranchSpec>,
}

/// Identifier and route specification for one channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyChannelSpec {
    /// Stable channel identifier.
    pub channel_id: String,
    /// Geometry and therapy metadata for the channel.
    pub route: ChannelRouteSpec,
}

/// Channel specification used in a series segment.
pub type SeriesChannelSpec = TopologyChannelSpec;
/// Channel specification used in a parallel segment.
pub type ParallelChannelSpec = TopologyChannelSpec;

/// Cross-sectional and taper geometry for a Venturi throat.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThroatGeometrySpec {
    /// Venturi throat width.
    pub throat_width_m: Length<f64>,
    /// Venturi throat height.
    pub throat_height_m: Length<f64>,
    /// Venturi throat length.
    pub throat_length_m: Length<f64>,
    /// Inlet channel width before contraction.
    pub inlet_width_m: Length<f64>,
    /// Outlet channel width after expansion.
    pub outlet_width_m: Length<f64>,
    /// Half-angle of the convergent section.
    pub convergent_half_angle: Angle<f64>,
    /// Half-angle of the divergent section.
    pub divergent_half_angle: Angle<f64>,
}

/// Shared Venturi configuration for one or more channels.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiConfig {
    /// Channel identifiers targeted by the configuration.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub target_channel_ids: Vec<String>,
    /// Number of serial throats.
    pub serial_throat_count: u8,
    /// Geometry shared by the configured throats.
    pub throat_geometry: ThroatGeometrySpec,
    /// Rule used to choose the placement site.
    pub placement_mode: VenturiPlacementMode,
}

/// Explicit Venturi placement on one channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiPlacementSpec {
    /// Stable identifier for the placement.
    pub placement_id: String,
    /// Channel receiving the Venturi stage.
    pub target_channel_id: String,
    /// Number of serial throats at this placement.
    pub serial_throat_count: u8,
    /// Geometry of the placed throats.
    pub throat_geometry: ThroatGeometrySpec,
    /// Rule used to choose the placement site.
    pub placement_mode: VenturiPlacementMode,
}

/// Estimated Dean-number site used for geometry placement.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct DeanSiteEstimate {
    /// Dean number estimated at the selected site.
    pub dean_number: Dimensionless<f64>,
    /// Local curvature radius.
    pub curvature_radius_m: Length<f64>,
    /// Local arc length used for the estimate.
    pub arc_length_m: Length<f64>,
}

/// Complete authored topology input for schematic generation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlueprintTopologySpec {
    /// Stable blueprint identifier.
    pub topology_id: String,
    /// Human-readable design name.
    pub design_name: String,
    /// Width and height of the enclosing plate.
    pub box_dims_m: (Length<f64>, Length<f64>),
    /// Inlet width.
    pub inlet_width_m: Length<f64>,
    /// Outlet width.
    pub outlet_width_m: Length<f64>,
    /// Length of the main trunk.
    pub trunk_length_m: Length<f64>,
    /// Length of the outlet tail.
    pub outlet_tail_length_m: Length<f64>,
    /// Channels arranged in series.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub series_channels: Vec<SeriesChannelSpec>,
    /// Channels arranged in parallel.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parallel_channels: Vec<ParallelChannelSpec>,
    /// Ordered split stages in the topology.
    pub split_stages: Vec<SplitStageSpec>,
    /// Explicit Venturi placements.
    #[serde(default)]
    pub venturi_placements: Vec<VenturiPlacementSpec>,
    /// Physical actuation mode for the blueprint.
    pub treatment_mode: TreatmentActuationMode,
}

/// Stages in the topology optimization lineage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TopologyOptimizationStage {
    /// Optimize residence-time separation at an asymmetric split.
    AsymmetricSplitResidenceSeparation,
    /// Optimize Venturi cavitation selectivity at an asymmetric split.
    AsymmetricSplitVenturiCavitationSelectivity,
    /// Refine an in-place Dean-serpentine route.
    InPlaceDeanSerpentineRefinement,
}

/// One mutation recorded in a topology optimization lineage.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyLineageEvent {
    /// Optimization stage that produced the mutation.
    pub stage: TopologyOptimizationStage,
    /// Human-readable description of the mutation.
    pub mutation: String,
    /// Blueprint from which the mutation was derived, when known.
    pub source_blueprint: Option<String>,
}

/// Lineage metadata for an optimized topology blueprint.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyLineageMetadata {
    /// Name of the root blueprint in the lineage.
    pub root_blueprint_name: String,
    /// Current optimization stage.
    pub current_stage: TopologyOptimizationStage,
    /// Optional source blueprint for the first optimization option.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub option1_source_blueprint: Option<String>,
    /// Optional source blueprint for the second optimization option.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub option2_source_blueprint: Option<String>,
    /// Optional blueprint used to seed the genetic algorithm.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ga_seed_blueprint: Option<String>,
    /// Ordered mutations applied to the lineage.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub mutations: Vec<TopologyLineageEvent>,
}

impl Default for TopologyLineageMetadata {
    fn default() -> Self {
        Self {
            root_blueprint_name: String::new(),
            current_stage: TopologyOptimizationStage::AsymmetricSplitResidenceSeparation,
            option1_source_blueprint: None,
            option2_source_blueprint: None,
            ga_seed_blueprint: None,
            mutations: Vec::new(),
        }
    }
}

impl BlueprintTopologySpec {
    /// Returns the authored plate envelope in millimetres for mesh and layout
    /// coordinates.
    #[must_use]
    pub fn box_dims_mm(&self) -> (f64, f64) {
        (
            self.box_dims_m.0.into_base() * 1.0e3,
            self.box_dims_m.1.into_base() * 1.0e3,
        )
    }

    /// Finds a route by its series, parallel, or derived split-channel ID.
    #[must_use]
    pub fn channel_route(&self, channel_id: &str) -> Option<&ChannelRouteSpec> {
        self.series_channels
            .iter()
            .find(|channel| channel.channel_id == channel_id)
            .map(|channel| &channel.route)
            .or_else(|| {
                self.parallel_channels
                    .iter()
                    .find(|channel| channel.channel_id == channel_id)
                    .map(|channel| &channel.route)
            })
            .or_else(|| {
                self.split_stages.iter().find_map(|stage| {
                    stage.branches.iter().find_map(|branch| {
                        (Self::branch_channel_id(&stage.stage_id, &branch.label) == channel_id)
                            .then_some(&branch.route)
                    })
                })
            })
    }

    /// Returns the stable IDs of channels on the treatment path.
    #[must_use]
    pub fn treatment_channel_ids(&self) -> Vec<String> {
        if !self.split_stages.is_empty() {
            return self
                .split_stages
                .iter()
                .flat_map(|stage| {
                    stage
                        .branches
                        .iter()
                        .filter(|branch| branch.treatment_path)
                        .map(|branch| format!("{}_{}", stage.stage_id, branch.label))
                })
                .collect();
        }

        let mut ids = self
            .series_channels
            .iter()
            .chain(self.parallel_channels.iter())
            .filter(|channel| channel.route.therapy_zone == TherapyZone::CancerTarget)
            .map(|channel| channel.channel_id.clone())
            .collect::<Vec<_>>();
        ids.extend(
            self.venturi_placements
                .iter()
                .map(|placement| placement.target_channel_id.clone()),
        );
        ids.sort();
        ids.dedup();
        ids
    }

    /// Derives a channel ID from a split stage and branch label.
    #[must_use]
    pub fn branch_channel_id(stage_id: &str, branch_label: &str) -> String {
        format!("{stage_id}_{branch_label}")
    }
}
