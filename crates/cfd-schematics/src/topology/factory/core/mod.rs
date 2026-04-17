//! Core orchestrator for `BlueprintTopologyFactory`.
//!
//! This module provides the canonical entry point for converting declarative
//! [`BlueprintTopologySpec`]s into [`NetworkBlueprint`] graphs.  Instead of
//! maintaining its own geometry builders (SSOT violation), it delegates to the
//! canonical [`create_geometry`](crate::geometry::generator::create_geometry)
//! pipeline via [`GeometryGeneratorBuilder`].

use crate::domain::model::NetworkBlueprint;
use crate::topology::model::{
    BlueprintTopologySpec, SerpentineSpec, SplitKind, TreatmentActuationMode,
    ThroatGeometrySpec, VenturiPlacementMode, VenturiPlacementSpec,
};


/// High-level mutations available to the GA optimization engine.
#[derive(Debug, Clone, PartialEq)]
pub enum BlueprintTopologyMutation {
    UpdateBranchWidth {
        stage_id: String,
        branch_label: String,
        new_width_m: f64,
    },
    ReplaceSplitKind {
        stage_id: String,
        split_kind: SplitKind,
    },
    InsertStage {
        stage_index: usize,
        split_kind: SplitKind,
    },
    UpdateVenturiConfiguration {
        placements: Vec<VenturiPlacementSpec>,
        treatment_mode: TreatmentActuationMode,
    },
    SetTreatmentSerpentine {
        stage_id: String,
        branch_label: String,
        serpentine: Option<SerpentineSpec>,
    },
    SetTreatmentChannelVenturi {
        target_channel_id: String,
        serial_throat_count: u8,
        throat_geometry: ThroatGeometrySpec,
        placement_mode: VenturiPlacementMode,
    },
    SetTreatmentChannelSerpentine {
        target_channel_id: String,
        serpentine: Option<SerpentineSpec>,
    },
    InsertTreatmentSplitMerge {
        target_channel_id: String,
        split_kind: SplitKind,
        treatment_serpentine: Option<SerpentineSpec>,
        venturi_serial_throat_count: Option<u8>,
        venturi_throat_geometry: Option<ThroatGeometrySpec>,
        venturi_placement_mode: VenturiPlacementMode,
    },
}

/// Core interface for turning declarative [`BlueprintTopologySpec`]s into
/// [`NetworkBlueprint`] graphs.
///
/// ## SSOT Architecture
///
/// This factory is a **thin facade** that delegates all geometry generation
/// to the canonical [`GeometryGeneratorBuilder`] pipeline.  No ad-hoc
/// node/channel construction is performed here — that logic lives exclusively
/// in [`GeometryGenerator`](crate::geometry::generator::GeometryGenerator).
pub struct BlueprintTopologyFactory;


mod build_impl;
mod mutation_impl;
mod spec_analysis_impl;

impl BlueprintTopologyFactory {
    /// Entrypoint: builds a fully detailed blueprint graph from a declarative
    /// topology spec by delegating to the canonical `create_geometry` pipeline.
    ///
    /// # Pipeline
    ///
    /// 1. Validate the spec (`validation::validate_spec`)
    /// 2. Convert `BlueprintTopologySpec` → `SplitType[]` + `GeometryConfig`
    /// 3. Delegate to `GeometryGeneratorBuilder` (canonical pipeline)
    /// 4. Apply venturi placements post-hoc
    /// 5. Attach topology + lineage metadata
    ///
    /// # Errors
    ///
    /// Returns a descriptive error string if the spec violates any geometric
    /// or structural constraint.
    pub fn build(spec: &BlueprintTopologySpec) -> Result<NetworkBlueprint, String> {
        super::validation::validate_spec(spec)?;

        let lineage = Self::lineage_for_spec(spec);

        let mut blueprint = if spec.has_series_path() && spec.split_stages.is_empty() {
            Self::build_series_path(spec, lineage)
        } else if spec.has_parallel_paths() && spec.split_stages.is_empty() {
            Self::build_parallel_path(spec, lineage)
        } else {
            Self::build_split_tree(spec, lineage)?
        };
        let resolved_spec = Self::resolve_materialized_venturi_targets(&blueprint, spec);
        blueprint.topology = Some(resolved_spec.clone());

        // Post-process: apply venturi placements
        if resolved_spec.has_venturi()
            && !Self::has_materialized_venturi_geometry(&blueprint, &resolved_spec)
        {
            super::modifiers::venturi::apply_venturi_placements(&mut blueprint, &resolved_spec)?;
        }

        Ok(blueprint)
    }


    // build_impl.rs: build_series_path, build_parallel_path, build_split_tree,
    // reconcile_channel_ids, venturi helpers, leading_merge_side_treatment_channels

    // mutation_impl.rs: validate_spec, mutate

    // spec_analysis_impl.rs: estimate_dean_site, lineage_for_spec, spec queries
}

#[cfg(test)]
mod tests {
    use super::{BlueprintTopologyFactory, BlueprintTopologyMutation};
    use crate::topology::presets::enumerate_milestone12_topologies;
    use crate::{
        BlueprintTopologySpec, SerpentineSpec, SplitKind, TopologyOptimizationStage,
        TreatmentActuationMode, VenturiPlacementMode,
    };
    use crate::domain::therapy_metadata::TherapyZone;

    fn base_blueprint() -> crate::NetworkBlueprint {
        let request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Tri-BASE")
            .expect("tri base request");
        crate::build_milestone12_blueprint(&request).expect("base blueprint")
    }

    fn treatment_venturi_channel_count(blueprint: &crate::NetworkBlueprint) -> usize {
        blueprint
            .channels
            .iter()
            .filter(|channel| {
                channel.venturi_geometry.is_some()
                    && channel.therapy_zone == Some(TherapyZone::CancerTarget)
            })
            .count()
    }

    fn inserted_treatment_split_merge_blueprint(split_kind: SplitKind) -> crate::NetworkBlueprint {
        let blueprint = base_blueprint();
        let target_channel_id = blueprint
            .treatment_channel_ids()
            .into_iter()
            .next()
            .expect("treatment channel");
        let topology = blueprint.topology_spec().expect("topology");
        let route = topology
            .channel_route(&target_channel_id)
            .expect("treatment route");

        BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                target_channel_id,
                split_kind,
                treatment_serpentine: None,
                venturi_serial_throat_count: Some(1),
                venturi_throat_geometry: Some(crate::ThroatGeometrySpec {
                    throat_width_m: 65.0e-6,
                    throat_height_m: route.height_m,
                    throat_length_m: 240.0e-6,
                    inlet_width_m: route.width_m,
                    outlet_width_m: route.width_m,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                }),
                venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect("split-merge mutation")
    }

    #[test]
    fn treatment_channel_mutations_target_cancer_path_only() {
        let blueprint = base_blueprint();
        let target_channel_id = blueprint
            .treatment_channel_ids()
            .into_iter()
            .next()
            .expect("treatment channel");
        let topology = blueprint.topology_spec().expect("topology");
        let route = topology
            .channel_route(&target_channel_id)
            .expect("route");

        let mutated = BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                target_channel_id: target_channel_id.clone(),
                serial_throat_count: 2,
                throat_geometry: crate::ThroatGeometrySpec {
                    throat_width_m: 80.0e-6,
                    throat_height_m: route.height_m,
                    throat_length_m: 300.0e-6,
                    inlet_width_m: route.width_m,
                    outlet_width_m: route.width_m,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                },
                placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect("venturi mutation");

        assert!(mutated
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi));
        assert!(mutated
            .channels
            .iter()
            .any(|channel| channel.therapy_zone == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)));
    }

    #[test]
    fn split_merge_insertion_preserves_geometry_authored_validation() {
        let mut request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Quad-Y")
            .expect("quad mirror request");
        request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
        request.venturi_throat_count = 1;
        let blueprint = crate::build_milestone12_blueprint(&request).expect("quad blueprint");
        let target_channel_id = blueprint
            .treatment_channel_ids()
            .into_iter()
            .next()
            .expect("treatment channel");

        let mutated = BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                target_channel_id,
                split_kind: SplitKind::NFurcation(3),
                treatment_serpentine: Some(SerpentineSpec {
                    wave_type: crate::topology::SerpentineWaveType::Sine,
                    segments: 4,
                    bend_radius_m: 1.2e-3,
                    segment_length_m: 4.0e-3,
                }),
                venturi_serial_throat_count: Some(2),
                venturi_throat_geometry: Some(crate::ThroatGeometrySpec {
                    throat_width_m: 65.0e-6,
                    throat_height_m: 1.0e-3,
                    throat_length_m: 240.0e-6,
                    inlet_width_m: 1.4e-3,
                    outlet_width_m: 1.4e-3,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                }),
                venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect("split-merge mutation");

        assert!(mutated.is_geometry_authored());
        assert!(mutated.validate().is_ok());
        let topology = mutated.topology_spec().expect("mutated topology");
        assert!(topology.split_stages.len() >= 2);
        assert!(
            topology.venturi_placements.iter().all(|placement| mutated.channels.iter().any(
                |channel| {
                    (channel.id.as_str() == placement.target_channel_id
                        || channel.id.as_str().starts_with(&placement.target_channel_id))
                        && channel.venturi_geometry.is_some()
                }
            )),
            "every declared venturi placement must materialize venturi geometry on a matching channel"
        );
    }

    #[test]
    fn inserted_center_bifurcation_expands_venturis_across_both_child_channels() {
        let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(2));
        let topology = blueprint.topology_spec().expect("resolved topology");

        assert_eq!(treatment_venturi_channel_count(&blueprint), 2);
        assert_eq!(topology.venturi_placements.len(), 2);
        assert!(topology.venturi_placements.iter().all(|placement| {
            blueprint.channels.iter().any(|channel| {
                channel.id.as_str() == placement.target_channel_id
                    && channel.therapy_zone == Some(TherapyZone::CancerTarget)
                    && channel.venturi_geometry.is_some()
            })
        }));
    }

    #[test]
    fn inserted_center_trifurcation_expands_venturis_across_all_child_channels() {
        let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(3));
        let topology = blueprint.topology_spec().expect("resolved topology");

        assert_eq!(treatment_venturi_channel_count(&blueprint), 3);
        assert_eq!(topology.venturi_placements.len(), 3);
        assert!(topology.venturi_placements.iter().all(|placement| {
            blueprint.channels.iter().any(|channel| {
                channel.id.as_str() == placement.target_channel_id
                    && channel.therapy_zone == Some(TherapyZone::CancerTarget)
                    && channel.venturi_geometry.is_some()
            })
        }));
    }

    #[test]
    fn venturis_never_materialize_on_healthy_bypass_channels() {
        let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(3));

        assert!(blueprint.channels.iter().all(|channel| {
            !(channel.therapy_zone == Some(TherapyZone::HealthyBypass)
                && channel.venturi_geometry.is_some())
        }));
    }
}
