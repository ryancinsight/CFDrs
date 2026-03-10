//! Core orchestrator for `BlueprintTopologyFactory`.
//!
//! This module provides the canonical entry point for converting declarative
//! [`BlueprintTopologySpec`]s into [`NetworkBlueprint`] graphs.  Instead of
//! maintaining its own geometry builders (SSOT violation), it delegates to the
//! canonical [`create_geometry`](crate::geometry::generator::create_geometry)
//! pipeline via [`GeometryGeneratorBuilder`].

use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{EdgeId, NetworkBlueprint};
use crate::geometry::generator::{
    create_primitive_selective_tree_geometry_from_spec, GeometryGeneratorBuilder,
};
use crate::geometry::metadata::BlueprintRenderHints;
use crate::geometry::types::SplitType;
use crate::topology::model::{
    BlueprintTopologySpec, DeanSiteEstimate, SplitKind, TopologyLineageMetadata,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiPlacementSpec,
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

        // Dispatch to the appropriate construction path:
        // - Series: direct linear chain construction (create_geometry cannot
        //   express multi-segment series paths)
        // - Split-tree: delegate to the canonical GeometryGeneratorBuilder
        let mut blueprint = if spec.has_series_path() && spec.split_stages.is_empty() {
            Self::build_series_path(spec, lineage)?
        } else {
            Self::build_split_tree(spec, lineage)?
        };

        // Post-process: apply venturi placements
        if spec.has_venturi() && !Self::has_materialized_venturi_geometry(&blueprint) {
            super::modifiers::venturi::apply_venturi_placements(&mut blueprint, spec)?;
        }

        Ok(blueprint)
    }

    /// Construct a linear series-path blueprint directly from the spec.
    ///
    /// Series topologies (inlet → section₁ → section₂ → … → outlet) cannot
    /// be expressed as `SplitType[]`, so we build nodes and channels directly
    /// from the topology spec's `series_channels` array.
    fn build_series_path(
        spec: &BlueprintTopologySpec,
        lineage: TopologyLineageMetadata,
    ) -> Result<NetworkBlueprint, String> {
        use crate::domain::model::{ChannelSpec, NodeKind, NodeSpec};
        use crate::domain::therapy_metadata::TherapyZoneMetadata;

        let mut bp = NetworkBlueprint {
            name: spec.design_name.clone(),
            box_dims: (127.76, 85.47),
            box_outline: Vec::new(),
            nodes: Vec::new(),
            channels: Vec::new(),
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
        };

        // Layout: place nodes linearly across the box width
        let n_nodes = spec.series_channels.len() + 1;
        let (box_w, box_h) = spec.box_dims_mm;
        let y_mid = box_h * 0.5;
        let x_step = box_w / n_nodes as f64;

        // Create inlet node
        let inlet_name = "inlet";
        bp.add_node(NodeSpec::new_at(
            inlet_name,
            NodeKind::Inlet,
            (x_step * 0.5, y_mid),
        ));

        // Create intermediate junction nodes and outlet
        let mut node_names: Vec<String> = vec![inlet_name.to_string()];
        for i in 0..spec.series_channels.len() {
            let name = if i == spec.series_channels.len() - 1 {
                "outlet".to_string()
            } else {
                format!("junction_{}", i)
            };
            let kind = if i == spec.series_channels.len() - 1 {
                NodeKind::Outlet
            } else {
                NodeKind::Junction
            };
            bp.add_node(NodeSpec::new_at(
                &name,
                kind,
                (x_step * (i as f64 + 1.5), y_mid),
            ));
            node_names.push(name);
        }

        // Create channels from spec
        for (i, ch) in spec.series_channels.iter().enumerate() {
            let from = &node_names[i];
            let to = &node_names[i + 1];
            bp.add_channel(
                ChannelSpec::new_pipe_rect(
                    &ch.channel_id,
                    from,
                    to,
                    ch.route.length_m,
                    ch.route.width_m,
                    ch.route.height_m,
                    0.0,
                    0.0,
                )
                .with_metadata(TherapyZoneMetadata::new(ch.route.therapy_zone)),
            );
        }

        // Attach topology and lineage
        bp.topology = Some(spec.clone());
        bp.lineage = Some(lineage);
        bp.box_dims = spec.box_dims_mm;

        Ok(bp)
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

        let splits = Self::spec_to_split_types(spec);

        let representative_width_m = Self::representative_channel_width_m(spec);
        let representative_height_m = Self::representative_channel_height_m(spec);

        let geo_config = GeometryConfig {
            channel_width: representative_width_m * 1e3,  // m → mm
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

    fn has_materialized_venturi_geometry(blueprint: &NetworkBlueprint) -> bool {
        blueprint
            .channels
            .iter()
            .any(|channel| channel.venturi_geometry.is_some())
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
        }

        let mut result = Self::build(&new_spec)?;
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
    /// - `VenturiCavitation` → `SelectiveVenturiCavitation`
    /// - `UltrasoundOnly` → `SelectiveAcousticResidenceSeparation`
    pub fn lineage_for_spec(spec: &BlueprintTopologySpec) -> TopologyLineageMetadata {
        let stage = match spec.treatment_mode {
            TreatmentActuationMode::VenturiCavitation => {
                TopologyOptimizationStage::SelectiveVenturiCavitation
            }
            TreatmentActuationMode::UltrasoundOnly => {
                TopologyOptimizationStage::SelectiveAcousticResidenceSeparation
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
    fn spec_to_split_types(spec: &BlueprintTopologySpec) -> Vec<SplitType> {
        spec.split_stages
            .iter()
            .map(|stage| match stage.split_kind {
                SplitKind::Bifurcation => SplitType::Bifurcation,
                SplitKind::Trifurcation => SplitType::Trifurcation,
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
        }
    }
}
