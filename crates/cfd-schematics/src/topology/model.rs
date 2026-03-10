use crate::domain::therapy_metadata::TherapyZone;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SplitKind {
    NFurcation(usize),
}

impl SplitKind {
    #[must_use]
    pub const fn branch_count(self) -> usize {
        match self {
            Self::NFurcation(n) => n,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BranchRole {
    Treatment,
    WbcCollection,
    RbcBypass,
    Neutral,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TreatmentActuationMode {
    UltrasoundOnly,
    VenturiCavitation,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VenturiPlacementMode {
    StraightSegment,
    CurvaturePeakDeanNumber,
    DiffuserShoulder,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SerpentineSpec {
    pub segments: usize,
    pub bend_radius_m: f64,
    pub segment_length_m: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChannelRouteSpec {
    pub length_m: f64,
    pub width_m: f64,
    pub height_m: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub serpentine: Option<SerpentineSpec>,
    #[serde(default = "default_therapy_zone")]
    pub therapy_zone: TherapyZone,
}

fn default_therapy_zone() -> TherapyZone {
    TherapyZone::MixedFlow
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BranchSpec {
    pub label: String,
    pub role: BranchRole,
    pub treatment_path: bool,
    pub route: ChannelRouteSpec,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SplitStageSpec {
    pub stage_id: String,
    pub split_kind: SplitKind,
    pub branches: Vec<BranchSpec>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyChannelSpec {
    pub channel_id: String,
    pub route: ChannelRouteSpec,
}

pub type SeriesChannelSpec = TopologyChannelSpec;
pub type ParallelChannelSpec = TopologyChannelSpec;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThroatGeometrySpec {
    pub throat_width_m: f64,
    pub throat_height_m: f64,
    pub throat_length_m: f64,
    pub inlet_width_m: f64,
    pub outlet_width_m: f64,
    pub convergent_half_angle_deg: f64,
    pub divergent_half_angle_deg: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiConfig {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub target_channel_ids: Vec<String>,
    pub serial_throat_count: u8,
    pub throat_geometry: ThroatGeometrySpec,
    pub placement_mode: VenturiPlacementMode,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VenturiPlacementSpec {
    pub placement_id: String,
    pub target_channel_id: String,
    pub serial_throat_count: u8,
    pub throat_geometry: ThroatGeometrySpec,
    pub placement_mode: VenturiPlacementMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub struct DeanSiteEstimate {
    pub dean_number: f64,
    pub curvature_radius_m: f64,
    pub arc_length_m: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BlueprintTopologySpec {
    pub topology_id: String,
    pub design_name: String,
    pub box_dims_mm: (f64, f64),
    pub inlet_width_m: f64,
    pub outlet_width_m: f64,
    pub trunk_length_m: f64,
    pub outlet_tail_length_m: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub series_channels: Vec<SeriesChannelSpec>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parallel_channels: Vec<ParallelChannelSpec>,
    pub split_stages: Vec<SplitStageSpec>,
    #[serde(default)]
    pub venturi_placements: Vec<VenturiPlacementSpec>,
    pub treatment_mode: TreatmentActuationMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TopologyOptimizationStage {
    AsymmetricSplitResidenceSeparation,
    AsymmetricSplitVenturiCavitationSelectivity,
    InPlaceDeanSerpentineRefinement,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyLineageEvent {
    pub stage: TopologyOptimizationStage,
    pub mutation: String,
    pub source_blueprint: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TopologyLineageMetadata {
    pub root_blueprint_name: String,
    pub current_stage: TopologyOptimizationStage,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub option1_source_blueprint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub option2_source_blueprint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ga_seed_blueprint: Option<String>,
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

    #[must_use]
    pub fn branch_channel_id(stage_id: &str, branch_label: &str) -> String {
        format!("{stage_id}_{branch_label}")
    }
}
