//! Core orchestrator for `BlueprintTopologyFactory`.
//!
//! This module provides the canonical entry point for converting declarative
//! [`BlueprintTopologySpec`]s into [`NetworkBlueprint`] graphs.  Instead of
//! maintaining its own geometry builders (SSOT violation), it delegates to the
//! canonical [`create_geometry`](crate::geometry::generator::create_geometry)
//! pipeline via [`GeometryGeneratorBuilder`].

use crate::domain::model::NetworkBlueprint;
use crate::topology::model::{
    BlueprintTopologySpec, SerpentineSpec, SplitKind, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
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

/// Core interface for turning declarative [`BlueprintTopologySpec`]s into
/// [`NetworkBlueprint`] graphs.
///
/// ## SSOT Architecture
///
/// This factory is a **thin facade** that delegates all geometry generation
/// to the canonical
/// [`GeometryGeneratorBuilder`](crate::geometry::generator::GeometryGeneratorBuilder)
/// pipeline. No ad-hoc
/// node/channel construction is performed here — that logic lives exclusively
/// in the private `GeometryGenerator` implementation.
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
    use crate::domain::therapy_metadata::TherapyZone;
    use crate::topology::presets::enumerate_milestone12_topologies;
    use crate::{
        BlueprintTopologySpec, SerpentineSpec, SplitKind, TopologyOptimizationStage,
        TreatmentActuationMode, VenturiPlacementMode,
    };
    use aequitas::systems::si::quantities::{Angle, Length};

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
                    throat_width_m: Length::from_base(65.0e-6),
                    throat_height_m: route.height_m,
                    throat_length_m: Length::from_base(240.0e-6),
                    inlet_width_m: route.width_m,
                    outlet_width_m: route.width_m,
                    convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                    divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
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
        let route = topology.channel_route(&target_channel_id).expect("route");

        let mutated = BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                target_channel_id: target_channel_id.clone(),
                serial_throat_count: 2,
                throat_geometry: crate::ThroatGeometrySpec {
                    throat_width_m: Length::from_base(80.0e-6),
                    throat_height_m: route.height_m,
                    throat_length_m: Length::from_base(300.0e-6),
                    inlet_width_m: route.width_m,
                    outlet_width_m: route.width_m,
                    convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                    divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                },
                placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect("venturi mutation");

        assert!(mutated
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi));
        assert!(mutated.channels.iter().any(|channel| channel.therapy_zone
            == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)));
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
                    bend_radius_m: Length::from_base(1.2e-3),
                    segment_length_m: Length::from_base(4.0e-3),
                }),
                venturi_serial_throat_count: Some(2),
                venturi_throat_geometry: Some(crate::ThroatGeometrySpec {
                    throat_width_m: Length::from_base(65.0e-6),
                    throat_height_m: Length::from_base(1.0e-3),
                    throat_length_m: Length::from_base(240.0e-6),
                    inlet_width_m: Length::from_base(1.4e-3),
                    outlet_width_m: Length::from_base(1.4e-3),
                    convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                    divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
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

    #[test]
    fn dean_site_prefers_target_route_serpentine_without_changing_values() {
        let mut blueprint = base_blueprint();
        let target_channel_id = blueprint
            .treatment_channel_ids()
            .into_iter()
            .next()
            .expect("treatment channel");
        blueprint
            .channels
            .iter_mut()
            .find(|channel| channel.therapy_zone == Some(TherapyZone::CancerTarget))
            .expect("materialized treatment channel")
            .id
            .0
            .clone_from(&target_channel_id);
        let mut topology = blueprint.topology_spec().expect("topology").clone();
        let expected_serpentine = SerpentineSpec {
            wave_type: crate::topology::SerpentineWaveType::Sine,
            segments: 3,
            bend_radius_m: Length::from_base(1.2e-3),
            segment_length_m: Length::from_base(4.0e-3),
        };
        for stage in &mut topology.split_stages {
            for branch in &mut stage.branches {
                if format!("{}_{}", stage.stage_id, branch.label) == target_channel_id {
                    branch.route.serpentine = Some(expected_serpentine.clone());
                }
            }
        }
        let route = topology
            .channel_route(&target_channel_id)
            .expect("target route");
        let channel = blueprint
            .channels
            .iter()
            .find(|channel| channel.id.as_str() == target_channel_id)
            .expect("target channel");
        let area = channel.cross_section.area().into_base();
        let hydraulic_diameter = channel.cross_section.hydraulic_diameter().into_base();
        let average_velocity = 2.0e-9 / area;
        let reynolds = average_velocity * hydraulic_diameter / 4.0e-6;
        let serpentine = route.serpentine.as_ref().expect("target serpentine");
        let expected_radius = serpentine.bend_radius_m.into_base();
        let expected_arc_length = (serpentine.segments.max(1) as f64
            * (serpentine.segment_length_m.into_base().max(0.0)
                + std::f64::consts::PI * serpentine.bend_radius_m.into_base().max(0.0)))
        .max(channel.length_m.into_base());
        let route_height = route.height_m;
        let route_width = route.width_m;
        blueprint.topology = Some(topology);

        let estimate = BlueprintTopologyFactory::estimate_dean_site(
            &blueprint,
            &crate::VenturiPlacementSpec {
                placement_id: "target".to_string(),
                target_channel_id,
                serial_throat_count: 1,
                throat_geometry: crate::ThroatGeometrySpec {
                    throat_width_m: Length::from_base(50.0e-6),
                    throat_height_m: route_height,
                    throat_length_m: Length::from_base(100.0e-6),
                    inlet_width_m: route_width,
                    outlet_width_m: route_width,
                    convergent_half_angle: Angle::from_base(0.1),
                    divergent_half_angle: Angle::from_base(0.1),
                },
                placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            2.0e-9,
            4.0e-6,
        )
        .expect("target route estimate");

        let expected_dean = reynolds * (hydraulic_diameter / (2.0 * expected_radius)).sqrt();
        assert_eq!(estimate.curvature_radius_m.into_base(), expected_radius);
        assert_eq!(estimate.arc_length_m.into_base(), expected_arc_length);
        assert_eq!(estimate.dean_number.into_base(), expected_dean);
    }
}
