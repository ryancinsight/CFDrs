use super::{ChannelSpec, NodeSpec};
use crate::geometry::metadata::{BlueprintRenderHints, MetadataContainer};
use crate::topology::{BlueprintTopologySpec, TopologyLineageMetadata};
use serde::{Deserialize, Serialize};
use std::fmt;

mod analysis_impl;
mod metadata_impl;
mod venturi_impl;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NetworkBlueprint {
    pub name: String,
    pub box_dims: (f64, f64),
    pub box_outline: Vec<((f64, f64), (f64, f64))>,
    pub nodes: Vec<NodeSpec>,
    pub channels: Vec<ChannelSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub render_hints: Option<BlueprintRenderHints>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub topology: Option<BlueprintTopologySpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lineage: Option<TopologyLineageMetadata>,
    #[serde(skip)]
    pub metadata: Option<MetadataContainer>,
    #[serde(default)]
    pub geometry_authored: bool,
}

/// Summary of physical channel footprint overlap in a blueprint.
#[derive(Debug, Clone, Copy, Default)]
pub struct ChannelOverlapAnalysis {
    /// Maximum overlap fraction across all channel pairs [0, 1].
    pub max_overlap_fraction: f64,
    /// Width ratio (wider / narrower) at the most-overlapping channel pair.
    pub width_ratio_at_worst: f64,
    /// Number of channel pairs with any physical footprint overlap.
    pub overlap_pair_count: usize,
}

impl fmt::Display for NetworkBlueprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.describe())
    }
}

#[cfg(test)]
mod tests {
    use crate::domain::therapy_metadata::TherapyZone;
    use crate::geometry::metadata::ChannelVenturiSpec;
    use crate::topology::presets::parallel_path_spec;
    use crate::topology::{
        ChannelRouteSpec, ParallelChannelSpec, VenturiConfig, VenturiPlacementMode,
    };
    use crate::{BlueprintTopologyFactory, TreatmentActuationMode};

    #[test]
    fn add_venturi_attaches_metadata_to_existing_parallel_channel() {
        let topology = parallel_path_spec(
            "venturi-blueprint",
            2.0e-3,
            2.0e-3,
            12.0e-3,
            12.0e-3,
            vec![ParallelChannelSpec {
                channel_id: "treatment_lane".to_string(),
                route: ChannelRouteSpec {
                    length_m: 10.0e-3,
                    width_m: 1.6e-3,
                    height_m: 1.0e-3,
                    serpentine: None,
                    therapy_zone: TherapyZone::CancerTarget,
                },
            }],
            TreatmentActuationMode::UltrasoundOnly,
        );
        let mut blueprint =
            BlueprintTopologyFactory::build(&topology).expect("parallel topology should build");

        blueprint
            .add_venturi(&VenturiConfig {
                target_channel_ids: vec!["treatment_lane".to_string()],
                serial_throat_count: 2,
                throat_geometry: crate::topology::ThroatGeometrySpec {
                    throat_width_m: 80.0e-6,
                    throat_height_m: 1.0e-3,
                    throat_length_m: 300.0e-6,
                    inlet_width_m: 0.0,
                    outlet_width_m: 0.0,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                },
                placement_mode: VenturiPlacementMode::StraightSegment,
            })
            .expect("venturi should attach");

        let channel = blueprint
            .channels
            .iter()
            .find(|channel| channel.id.as_str() == "treatment_lane")
            .expect("target channel must exist");
        assert!(channel.venturi_geometry.is_some());
        assert_eq!(
            channel
                .metadata
                .as_ref()
                .and_then(|metadata| metadata.get::<ChannelVenturiSpec>())
                .expect("ChannelVenturiSpec must be inserted")
                .n_throats,
            2
        );
        assert_eq!(
            blueprint
                .topology
                .as_ref()
                .expect("topology metadata must be preserved")
                .venturi_placements
                .len(),
            1
        );
    }
}
