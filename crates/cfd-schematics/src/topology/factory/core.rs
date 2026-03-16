//! Core orchestrator for `BlueprintTopologyFactory`.
//!
//! This module provides the canonical entry point for converting declarative
//! [`BlueprintTopologySpec`]s into [`NetworkBlueprint`] graphs.  Instead of
//! maintaining its own geometry builders (SSOT violation), it delegates to the
//! canonical [`create_geometry`](crate::geometry::generator::create_geometry)
//! pipeline via [`GeometryGeneratorBuilder`].

use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::therapy_metadata::TherapyZone;
use crate::domain::model::{EdgeId, NetworkBlueprint};
use crate::geometry::generator::{
    create_parallel_geometry_from_spec, create_primitive_selective_tree_geometry_from_spec,
    create_series_geometry_from_spec, GeometryGeneratorBuilder,
};
use crate::geometry::metadata::{BlueprintRenderHints, ChannelVisualRole};
use crate::geometry::types::SplitType;
use crate::topology::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec, DeanSiteEstimate,
    SerpentineSpec, SplitKind, SplitStageSpec, ThroatGeometrySpec, TopologyLineageMetadata,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiPlacementMode,
    VenturiPlacementSpec,
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
        serpentine: Option<crate::topology::model::SerpentineSpec>,
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
            Self::build_series_path(spec, lineage)?
        } else if spec.has_parallel_paths() && spec.split_stages.is_empty() {
            Self::build_parallel_path(spec, lineage)?
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

    fn build_series_path(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> Result<NetworkBlueprint, String> {
        let mut blueprint = create_series_geometry_from_spec(spec);
        blueprint.lineage = Some(lineage);
        blueprint.render_hints = Some(Self::render_hints_for_spec(spec));
        Ok(blueprint)
    }

    fn build_parallel_path(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> Result<NetworkBlueprint, String> {
        let mut blueprint = create_parallel_geometry_from_spec(spec);
        blueprint.lineage = Some(lineage);
        blueprint.render_hints = Some(Self::render_hints_for_spec(spec));
        Ok(blueprint)
    }

    /// Construct a split-tree blueprint via the canonical GeometryGeneratorBuilder.
    fn build_split_tree(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> Result<NetworkBlueprint, String> {
        if let Some(mut blueprint) = create_primitive_selective_tree_geometry_from_spec(spec)? {
            blueprint.lineage = Some(lineage);
            return Ok(blueprint);
        }

        let splits = Self::spec_to_split_types(spec)?;

        let representative_width_m = Self::representative_channel_width_m(spec);
        let representative_height_m = Self::representative_channel_height_m(spec);

        let geo_config = GeometryConfig {
            channel_width: representative_width_m * 1e3,   // m → mm
            channel_height: representative_height_m * 1e3, // m → mm
            ..GeometryConfig::default()
        };

        let mut blueprint = GeometryGeneratorBuilder::new(
            spec.box_dims_mm,
            &splits,
            &geo_config,
            &ChannelTypeConfig::AllStraight,
        )
        .with_topology_spec(spec.clone())
        .with_lineage(lineage)
        .with_render_hints(Self::render_hints_for_spec(spec))
        .build();

        blueprint.name = spec.design_name.clone();

        // Reconcile auto-generated channel IDs with the spec's semantic names
        Self::reconcile_channel_ids(&mut blueprint, spec);

        Ok(blueprint)
    }

    /// Reconcile auto-generated channel IDs (`ch_0`, `ch_1`, ...) with the
    /// topology spec's semantic channel names.
    ///
    /// For series-path specs, channels are renamed 1:1 in positional order.
    /// For split-tree specs, only the trunk and outlet channels can be
    /// reconciled positionally; branch channels keep their auto-generated IDs.
    fn reconcile_channel_ids(blueprint: &mut NetworkBlueprint, spec: &BlueprintTopologySpec) {
        if !spec.series_channels.is_empty() {
            // Series path: channels map 1:1 in order
            let named_ids: Vec<String> = spec
                .series_channels
                .iter()
                .map(|ch| ch.channel_id.clone())
                .collect();
            for (i, name) in named_ids.into_iter().enumerate() {
                if let Some(channel) = blueprint.channels.get_mut(i) {
                    channel.id = EdgeId(name);
                }
            }
        } else if !spec.parallel_channels.is_empty() {
            let named_ids: Vec<String> = spec
                .parallel_channels
                .iter()
                .map(|ch| ch.channel_id.clone())
                .collect();
            for (i, name) in named_ids.into_iter().enumerate() {
                if let Some(channel) = blueprint.channels.get_mut(i) {
                    channel.id = EdgeId(name);
                }
            }
        }
    }

    fn has_materialized_venturi_geometry(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> bool {
        !spec.venturi_placements.is_empty()
            && spec.venturi_placements.iter().all(|placement| {
                blueprint.channels.iter().any(|channel| {
                    (channel.id.as_str() == placement.target_channel_id
                        || channel.id.as_str().starts_with(&placement.target_channel_id))
                        && channel.venturi_geometry.is_some()
                })
            })
    }

    fn resolve_materialized_venturi_targets(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> BlueprintTopologySpec {
        let expanded = Self::selective_materialized_venturi_placements(blueprint, spec);
        if expanded.is_empty() {
            return spec.clone();
        }

        let mut resolved = spec.clone();
        resolved.venturi_placements = expanded;
        resolved
    }

    fn selective_materialized_venturi_placements(
        blueprint: &NetworkBlueprint,
        spec: &BlueprintTopologySpec,
    ) -> Vec<VenturiPlacementSpec> {
        if spec.split_stages.is_empty() || spec.venturi_placements.is_empty() {
            return Vec::new();
        }

        let preferred = Self::leading_merge_side_treatment_channels(blueprint, true);
        let fallback = Self::leading_merge_side_treatment_channels(blueprint, false);
        let actual_ids = if fallback.len() > preferred.len() {
            fallback
        } else if preferred.is_empty() {
            fallback
        } else {
            preferred
        };

        if actual_ids.is_empty() {
            return Vec::new();
        }

        let template = spec
            .venturi_placements
            .iter()
            .max_by_key(|placement| placement.serial_throat_count)
            .or_else(|| spec.venturi_placements.first())
            .expect("non-empty venturi placements required for selective expansion");

        actual_ids
            .into_iter()
            .enumerate()
            .map(|(idx, materialized_id)| VenturiPlacementSpec {
                placement_id: format!("venturi_{idx}"),
                target_channel_id: materialized_id,
                serial_throat_count: template.serial_throat_count,
                throat_geometry: template.throat_geometry.clone(),
                placement_mode: template.placement_mode,
            })
            .collect()
    }

    fn leading_merge_side_treatment_channels(
        blueprint: &NetworkBlueprint,
        require_venturi: bool,
    ) -> Vec<String> {
        let mut channels = blueprint
            .channels
            .iter()
            .filter(|channel| {
                channel.therapy_zone == Some(TherapyZone::CancerTarget)
                    && channel.visual_role != Some(ChannelVisualRole::Trunk)
                    && (!require_venturi || channel.venturi_geometry.is_some())
            })
            .filter_map(|channel| {
                let max_x = channel
                    .path
                    .iter()
                    .map(|(x, _)| *x)
                    .fold(f64::NEG_INFINITY, f64::max);
                (max_x > blueprint.box_dims.0 * 0.5 + 1.0e-6).then(|| {
                    let centroid_x =
                        channel.path.iter().map(|(x, _)| *x).sum::<f64>() / channel.path.len() as f64;
                    let centroid_y =
                        channel.path.iter().map(|(_, y)| *y).sum::<f64>() / channel.path.len() as f64;
                    (channel.id.0.clone(), centroid_x, centroid_y)
                })
            })
            .collect::<Vec<_>>();

        channels.sort_by(|left, right| {
            left.1
                .total_cmp(&right.1)
                .then(left.2.total_cmp(&right.2))
                .then(left.0.cmp(&right.0))
        });
        let Some((_, lead_x, _)) = channels.first() else {
            return Vec::new();
        };
        let lead_x = *lead_x;
        channels
            .into_iter()
            .filter(|(_, centroid_x, _)| (*centroid_x - lead_x).abs() <= 1.0e-6)
            .map(|(channel_id, _, _)| channel_id)
            .collect()
    }

    /// Validate a spec without building (called by `NetworkBlueprint::validate`).
    pub fn validate_spec(spec: &BlueprintTopologySpec) -> Result<(), String> {
        super::validation::validate_spec(spec)
    }

    /// Derives a new blueprint variant based on a valid GA mutation command.
    pub fn mutate(
        blueprint: &NetworkBlueprint,
        mutation: BlueprintTopologyMutation,
        next_stage: TopologyOptimizationStage,
    ) -> Result<NetworkBlueprint, String> {
        let spec = blueprint
            .topology
            .as_ref()
            .ok_or("Blueprint has no topology spec to mutate")?;
        if spec.is_selective_routing() && !blueprint.is_geometry_authored() {
            return Err(format!(
                "selective-routing mutation requires create_geometry-authored provenance on blueprint '{}'",
                blueprint.name
            ));
        }

        let mut new_spec = spec.clone();
        match mutation {
            BlueprintTopologyMutation::UpdateBranchWidth {
                stage_id,
                branch_label,
                new_width_m,
            } => {
                let branch = new_spec
                    .split_stages
                    .iter_mut()
                    .find(|s| s.stage_id == stage_id)
                    .and_then(|s| s.branches.iter_mut().find(|b| b.label == branch_label))
                    .ok_or_else(|| format!("Branch {stage_id}/{branch_label} not found"))?;
                branch.route.width_m = new_width_m;
            }
            BlueprintTopologyMutation::ReplaceSplitKind {
                stage_id,
                split_kind,
            } => {
                let stage = new_spec
                    .split_stages
                    .iter_mut()
                    .find(|s| s.stage_id == stage_id)
                    .ok_or_else(|| format!("Stage {stage_id} not found"))?;
                stage.split_kind = split_kind;
            }
            BlueprintTopologyMutation::InsertStage {
                stage_index,
                split_kind,
            } => {
                let new_stage = crate::topology::model::SplitStageSpec {
                    stage_id: format!("stage_{}", new_spec.split_stages.len()),
                    split_kind,
                    branches: (0..split_kind.branch_count())
                        .map(|i| crate::topology::model::BranchSpec {
                            label: format!("branch_{i}"),
                            role: crate::topology::model::BranchRole::Neutral,
                            treatment_path: i == 0,
                            route: crate::topology::model::ChannelRouteSpec {
                                length_m: 10.0e-3,
                                width_m: 1.0e-3,
                                height_m: 0.5e-3,
                                serpentine: None,
                                therapy_zone:
                                    crate::domain::therapy_metadata::TherapyZone::MixedFlow,
                            },
                            recovery_sub_split: None,
                        })
                        .collect(),
                };
                let idx = stage_index.min(new_spec.split_stages.len());
                new_spec.split_stages.insert(idx, new_stage);
            }
            BlueprintTopologyMutation::UpdateVenturiConfiguration {
                placements,
                treatment_mode,
            } => {
                new_spec.venturi_placements = placements;
                new_spec.treatment_mode = treatment_mode;
            }
            BlueprintTopologyMutation::SetTreatmentSerpentine {
                stage_id,
                branch_label,
                serpentine,
            } => {
                let branch = new_spec
                    .split_stages
                    .iter_mut()
                    .find(|s| s.stage_id == stage_id)
                    .and_then(|s| s.branches.iter_mut().find(|b| b.label == branch_label))
                    .ok_or_else(|| {
                        format!("Branch {stage_id}/{branch_label} not found for serpentine update")
                    })?;
                branch.route.serpentine = serpentine;
            }
            BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                target_channel_id,
                serial_throat_count,
                throat_geometry,
                placement_mode,
            } => {
                let route = new_spec
                    .channel_route(&target_channel_id)
                    .ok_or_else(|| format!("Treatment channel '{}' not found", target_channel_id))?;
                if route.therapy_zone != TherapyZone::CancerTarget {
                    return Err(format!(
                        "venturi mutation requires a CancerTarget channel, but '{}' is {:?}",
                        target_channel_id, route.therapy_zone
                    ));
                }
                let placement_id = format!("mut_vt_{}", new_spec.venturi_placements.len());
                let serial_throat_count = serial_throat_count.max(1);
                let mut replaced = false;
                for placement in &mut new_spec.venturi_placements {
                    if placement.target_channel_id == target_channel_id {
                        placement.serial_throat_count = serial_throat_count;
                        placement.throat_geometry = throat_geometry.clone();
                        placement.placement_mode = placement_mode;
                        replaced = true;
                    }
                }
                if !replaced {
                    new_spec.venturi_placements.push(VenturiPlacementSpec {
                        placement_id,
                        target_channel_id: target_channel_id.clone(),
                        serial_throat_count,
                        throat_geometry,
                        placement_mode,
                    });
                }
                new_spec.treatment_mode = TreatmentActuationMode::VenturiCavitation;
            }
            BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
                target_channel_id,
                serpentine,
            } => {
                let (stage_id, branch_label) = resolve_stage_branch_for_channel(&new_spec, &target_channel_id)
                    .ok_or_else(|| {
                        format!(
                            "treatment-channel serpentine update requires split-stage channel '{}'",
                            target_channel_id
                        )
                    })?;
                let branch = new_spec
                    .split_stages
                    .iter_mut()
                    .find(|stage| stage.stage_id == stage_id)
                    .and_then(|stage| stage.branches.iter_mut().find(|branch| branch.label == branch_label))
                    .ok_or_else(|| {
                        format!(
                            "branch '{}' in stage '{}' not found for channel '{}'",
                            branch_label, stage_id, target_channel_id
                        )
                    })?;
                branch.route.serpentine = serpentine;
            }
            BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                target_channel_id,
                split_kind,
                treatment_serpentine,
                venturi_serial_throat_count,
                venturi_throat_geometry,
                venturi_placement_mode,
            } => {
                let (stage_id, branch_label) = resolve_stage_branch_for_channel(&new_spec, &target_channel_id)
                    .ok_or_else(|| {
                        format!(
                            "split-merge insertion requires split-stage treatment channel '{}'",
                            target_channel_id
                        )
                    })?;
                let stage_index = new_spec
                    .split_stages
                    .iter()
                    .position(|stage| stage.stage_id == stage_id)
                    .ok_or_else(|| format!("stage '{}' not found", stage_id))?;
                let parent_branch = new_spec.split_stages[stage_index]
                    .branches
                    .iter()
                    .find(|branch| branch.label == branch_label)
                    .ok_or_else(|| format!("branch '{}' not found in '{}'", branch_label, stage_id))?
                    .clone();
                if !parent_branch.treatment_path || parent_branch.route.therapy_zone != TherapyZone::CancerTarget {
                    return Err(format!(
                        "split-merge insertion requires CancerTarget treatment channel '{}'",
                        target_channel_id
                    ));
                }

                let inserted_stage_id = format!("stage_{}", stage_index + 1);
                let inserted_stage = build_inserted_treatment_stage(
                    &inserted_stage_id,
                    split_kind,
                    &parent_branch.route,
                    treatment_serpentine,
                )?;
                if let Some(serial_throat_count) = venturi_serial_throat_count {
                    let throat_geometry = venturi_throat_geometry.clone().ok_or_else(|| {
                        "venturi split-merge insertion requires throat geometry".to_string()
                    })?;
                    let treatment_branch = inserted_stage
                        .branches
                        .iter()
                        .find(|branch| branch.treatment_path)
                        .ok_or_else(|| "inserted split stage lost treatment branch".to_string())?;
                    new_spec.venturi_placements.push(VenturiPlacementSpec {
                        placement_id: format!("mut_vt_{}", new_spec.venturi_placements.len()),
                        target_channel_id: BlueprintTopologySpec::branch_channel_id(
                            &inserted_stage_id,
                            &treatment_branch.label,
                        ),
                        serial_throat_count: serial_throat_count.max(1),
                        throat_geometry,
                        placement_mode: venturi_placement_mode,
                    });
                    new_spec.treatment_mode = TreatmentActuationMode::VenturiCavitation;
                }
                new_spec
                    .split_stages
                    .insert(stage_index + 1, inserted_stage);
                for (index, stage) in new_spec.split_stages.iter_mut().enumerate() {
                    stage.stage_id = format!("stage_{index}");
                }
            }
        }

        let mut result = Self::build(&new_spec)?;
        if let Some(source_hints) = blueprint.render_hints() {
            if source_hints.mirror_x || source_hints.mirror_y {
                Self::mirror_blueprint_geometry(
                    &mut result,
                    new_spec.box_dims_mm,
                    source_hints.mirror_x,
                    source_hints.mirror_y,
                );
                if let Some(render_hints) = result.render_hints.as_mut() {
                    render_hints.mirror_x = source_hints.mirror_x;
                    render_hints.mirror_y = source_hints.mirror_y;
                }
            }
        }
        if let Some(ref mut lineage) = result.lineage {
            lineage.current_stage = next_stage;
            lineage
                .mutations
                .push(crate::topology::model::TopologyLineageEvent {
                    stage: next_stage,
                    mutation: format!("{:?}", "mutation applied"),
                    source_blueprint: Some(blueprint.name.clone()),
                });
        }
        Ok(result)
    }

    /// Computes Dean Number estimate for a venturi placement site.
    pub fn estimate_dean_site(
        blueprint: &NetworkBlueprint,
        placement: &VenturiPlacementSpec,
        flow_m3_s: f64,
        kinematic_viscosity_m2_s: f64,
    ) -> Option<DeanSiteEstimate> {
        let channel = blueprint
            .channels
            .iter()
            .find(|ch| ch.id.as_str() == placement.target_channel_id)?;

        let d_h = channel.cross_section.hydraulic_diameter();
        let area = channel.cross_section.area();
        if area <= 0.0 || d_h <= 0.0 {
            return None;
        }

        let v_avg = flow_m3_s / area;
        let reynolds = v_avg * d_h / kinematic_viscosity_m2_s;

        // Estimate curvature radius from polyline path (3-point circumradius)
        let curve_radius_m = channel
            .path
            .windows(3)
            .filter_map(|pts| {
                let (ax, ay) = pts[0];
                let (bx, by) = pts[1];
                let (cx, cy) = pts[2];
                let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
                if d.abs() < 1e-15 {
                    return None;
                }
                let ux = ((ax * ax + ay * ay) * (by - cy)
                    + (bx * bx + by * by) * (cy - ay)
                    + (cx * cx + cy * cy) * (ay - by))
                    / d;
                let uy = ((ax * ax + ay * ay) * (cx - bx)
                    + (bx * bx + by * by) * (ax - cx)
                    + (cx * cx + cy * cy) * (bx - ax))
                    / d;
                let r = ((ax - ux).powi(2) + (ay - uy).powi(2)).sqrt();
                Some(r * 1e-3) // mm → m
            })
            .next()
            .unwrap_or(5.0e-3);

        // De = Re √(D_h / (2 R_c))
        let dean_num = if curve_radius_m > 0.0 {
            reynolds * (d_h / (2.0 * curve_radius_m)).sqrt()
        } else {
            0.0
        };

        Some(DeanSiteEstimate {
            dean_number: dean_num,
            curvature_radius_m: curve_radius_m,
            arc_length_m: channel.length_m,
        })
    }

    /// Extract standard lineage tracing metadata.
    ///
    /// The optimization stage is derived from the spec's treatment mode:
    /// - `VenturiCavitation` → `AsymmetricSplitVenturiCavitationSelectivity`
    /// - `UltrasoundOnly` → `AsymmetricSplitResidenceSeparation`
    pub fn lineage_for_spec(spec: &BlueprintTopologySpec) -> TopologyLineageMetadata {
        let stage = match spec.treatment_mode {
            TreatmentActuationMode::VenturiCavitation => {
                TopologyOptimizationStage::AsymmetricSplitVenturiCavitationSelectivity
            }
            TreatmentActuationMode::UltrasoundOnly => {
                TopologyOptimizationStage::AsymmetricSplitResidenceSeparation
            }
        };
        TopologyLineageMetadata {
            root_blueprint_name: spec.design_name.clone(),
            current_stage: stage,
            option1_source_blueprint: None,
            option2_source_blueprint: None,
            ga_seed_blueprint: None,
            mutations: Vec::new(),
        }
    }

    // ── Private helpers ─────────────────────────────────────────────────

    /// Convert `BlueprintTopologySpec.split_stages` → `Vec<SplitType>`.
    ///
    /// For specs with no split stages (series/parallel), returns an empty
    /// split array so that `create_geometry` generates a linear channel.
    fn spec_to_split_types(spec: &BlueprintTopologySpec) -> Result<Vec<SplitType>, String> {
        spec.split_stages
            .iter()
            .map(|stage| match stage.split_kind {
                SplitKind::NFurcation(2) => Ok(SplitType::Bifurcation),
                SplitKind::NFurcation(3) => Ok(SplitType::Trifurcation),
                SplitKind::NFurcation(4) => Ok(SplitType::Quadfurcation),
                SplitKind::NFurcation(5) => Ok(SplitType::Pentafurcation),
                SplitKind::NFurcation(other) => Err(format!(
                    "split-tree geometry generation supports only N=2,3,4,5; stage '{}' requested N={other}",
                    stage.stage_id
                )),
            })
            .collect()
    }

    /// Derive a representative channel width from the spec for `GeometryConfig`.
    fn representative_channel_width_m(spec: &BlueprintTopologySpec) -> f64 {
        if let Some(first_series) = spec.series_channels.first() {
            return first_series.route.width_m;
        }
        if let Some(first_parallel) = spec.parallel_channels.first() {
            return first_parallel.route.width_m;
        }
        if let Some(first_stage) = spec.split_stages.first() {
            if let Some(first_branch) = first_stage.branches.first() {
                return first_branch.route.width_m;
            }
        }
        1.0e-3 // Default 1mm
    }

    /// Derive a representative channel height from the spec for `GeometryConfig`.
    fn representative_channel_height_m(spec: &BlueprintTopologySpec) -> f64 {
        if let Some(first_series) = spec.series_channels.first() {
            return first_series.route.height_m;
        }
        if let Some(first_parallel) = spec.parallel_channels.first() {
            return first_parallel.route.height_m;
        }
        if let Some(first_stage) = spec.split_stages.first() {
            if let Some(first_branch) = first_stage.branches.first() {
                return first_branch.route.height_m;
            }
        }
        0.5e-3 // Default 0.5mm
    }

    fn render_hints_for_spec(spec: &BlueprintTopologySpec) -> BlueprintRenderHints {
        BlueprintRenderHints {
            stage_sequence: spec.stage_sequence_label(),
            split_layers: spec.visible_split_layers(),
            throat_count_hint: spec.venturi_count(),
            treatment_label: if spec.treatment_mode == TreatmentActuationMode::VenturiCavitation {
                "venturi".to_string()
            } else {
                "ultrasound".to_string()
            },
            mirror_x: false,
            mirror_y: false,
        }
    }

    fn mirror_blueprint_geometry(
        blueprint: &mut NetworkBlueprint,
        box_dims_mm: (f64, f64),
        mirror_x: bool,
        mirror_y: bool,
    ) {
        let reflect = |(x, y): (f64, f64)| -> (f64, f64) {
            let mapped_x = if mirror_x { box_dims_mm.0 - x } else { x };
            let mapped_y = if mirror_y { box_dims_mm.1 - y } else { y };
            (mapped_x, mapped_y)
        };

        for node in &mut blueprint.nodes {
            node.point = reflect(node.point);
        }
        for channel in &mut blueprint.channels {
            channel.path = channel.path.iter().copied().map(reflect).collect();
        }
        blueprint.box_outline = blueprint
            .box_outline
            .iter()
            .map(|(start, end)| (reflect(*start), reflect(*end)))
            .collect();
    }
}

fn resolve_stage_branch_for_channel(
    spec: &BlueprintTopologySpec,
    target_channel_id: &str,
) -> Option<(String, String)> {
    spec.split_stages.iter().find_map(|stage| {
        stage.branches.iter().find_map(|branch| {
            let channel_id = BlueprintTopologySpec::branch_channel_id(&stage.stage_id, &branch.label);
            (channel_id == target_channel_id).then(|| (stage.stage_id.clone(), branch.label.clone()))
        })
    })
}

fn build_inserted_treatment_stage(
    stage_id: &str,
    split_kind: SplitKind,
    parent_route: &ChannelRouteSpec,
    treatment_serpentine: Option<SerpentineSpec>,
) -> Result<SplitStageSpec, String> {
    let parent_width_m = parent_route.width_m;
    let height_m = parent_route.height_m;
    let length_m = parent_route.length_m.max(1.0e-3);
    let branches = match split_kind {
        SplitKind::NFurcation(2) => {
            let treatment_width = parent_width_m * 0.68;
            vec![
                BranchSpec {
                    label: "upper".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    route: ChannelRouteSpec {
                        length_m,
                        width_m: treatment_width,
                        height_m,
                        serpentine: treatment_serpentine.clone(),
                        therapy_zone: TherapyZone::CancerTarget,
                    },
                    recovery_sub_split: None,
                },
                BranchSpec {
                    label: "lower".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    route: ChannelRouteSpec {
                        length_m,
                        width_m: (parent_width_m - treatment_width).max(f64::EPSILON),
                        height_m,
                        serpentine: None,
                        therapy_zone: TherapyZone::HealthyBypass,
                    },
                    recovery_sub_split: None,
                },
            ]
        }
        SplitKind::NFurcation(3) => {
            let center_width = parent_width_m * 0.45;
            let side_width = (parent_width_m - center_width) * 0.5;
            vec![
                bypass_branch("left", side_width, height_m, length_m, BranchRole::WbcCollection),
                treatment_branch("center", center_width, height_m, length_m, treatment_serpentine.clone()),
                bypass_branch("right", side_width, height_m, length_m, BranchRole::RbcBypass),
            ]
        }
        SplitKind::NFurcation(4) => {
            let treatment_width = parent_width_m * 0.40;
            let side_width = (parent_width_m - treatment_width) / 3.0;
            vec![
                bypass_branch("arm_0", side_width, height_m, length_m, BranchRole::WbcCollection),
                treatment_branch("arm_1", treatment_width, height_m, length_m, treatment_serpentine.clone()),
                bypass_branch("arm_2", side_width, height_m, length_m, BranchRole::RbcBypass),
                bypass_branch("arm_3", side_width, height_m, length_m, BranchRole::RbcBypass),
            ]
        }
        SplitKind::NFurcation(5) => {
            let treatment_width = parent_width_m * 0.36;
            let side_width = (parent_width_m - treatment_width) / 4.0;
            vec![
                bypass_branch("arm_0", side_width, height_m, length_m, BranchRole::WbcCollection),
                bypass_branch("arm_1", side_width, height_m, length_m, BranchRole::Neutral),
                treatment_branch("center", treatment_width, height_m, length_m, treatment_serpentine),
                bypass_branch("arm_3", side_width, height_m, length_m, BranchRole::RbcBypass),
                bypass_branch("arm_4", side_width, height_m, length_m, BranchRole::RbcBypass),
            ]
        }
        SplitKind::NFurcation(other) => {
            return Err(format!("unsupported split-merge insertion N={other}"));
        }
    };

    Ok(SplitStageSpec {
        stage_id: stage_id.to_string(),
        split_kind,
        branches,
    })
}

fn treatment_branch(
    label: &str,
    width_m: f64,
    height_m: f64,
    length_m: f64,
    serpentine: Option<SerpentineSpec>,
) -> BranchSpec {
    BranchSpec {
        label: label.to_string(),
        role: BranchRole::Treatment,
        treatment_path: true,
        route: ChannelRouteSpec {
            length_m,
            width_m,
            height_m,
            serpentine,
            therapy_zone: TherapyZone::CancerTarget,
        },
        recovery_sub_split: None,
    }
}

fn bypass_branch(
    label: &str,
    width_m: f64,
    height_m: f64,
    length_m: f64,
    role: BranchRole,
) -> BranchSpec {
    BranchSpec {
        label: label.to_string(),
        role,
        treatment_path: false,
        route: ChannelRouteSpec {
            length_m,
            width_m,
            height_m,
            serpentine: None,
            therapy_zone: TherapyZone::HealthyBypass,
        },
        recovery_sub_split: None,
    }
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
                    segments: 5,
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
