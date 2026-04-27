//! Spec analysis, Dean number estimation, and query methods for BlueprintTopologyFactory.
use super::BlueprintTopologyFactory;
use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::BlueprintRenderHints;
use crate::geometry::types::SplitType;
use crate::topology::model::{
    BlueprintTopologySpec, DeanSiteEstimate, SplitKind, TopologyLineageMetadata,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiPlacementSpec,
};

impl BlueprintTopologyFactory {
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

        let spec_serpentine = blueprint.topology_spec().and_then(|topology| {
            topology
                .channel_route(&placement.target_channel_id)
                .and_then(|route| route.serpentine.clone())
                .or_else(|| {
                    topology
                        .treatment_channel_ids()
                        .into_iter()
                        .filter_map(|channel_id| topology.channel_route(&channel_id))
                        .find_map(|route| route.serpentine.clone())
                })
        });

        // Estimate curvature radius from polyline path (3-point circumradius)
        let path_curve_radius_m = channel.path.windows(3).find_map(|pts| {
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
        });
        let spec_curve_radius_m = spec_serpentine
            .as_ref()
            .map(|serpentine| serpentine.bend_radius_m);
        let curve_radius_m = spec_curve_radius_m
            .or(path_curve_radius_m)
            .unwrap_or(5.0e-3);

        let spec_arc_length_m = spec_serpentine.as_ref().map(|serpentine| {
            let segments = serpentine.segments.max(1) as f64;
            let straight_length = segments * serpentine.segment_length_m.max(0.0);
            let bend_length = segments * std::f64::consts::PI * serpentine.bend_radius_m.max(0.0);
            (straight_length + bend_length).max(channel.length_m)
        });
        let arc_length_m = spec_arc_length_m.unwrap_or(channel.length_m);

        // De = Re √(D_h / (2 R_c))
        let dean_num = if curve_radius_m > 0.0 {
            reynolds * (d_h / (2.0 * curve_radius_m)).sqrt()
        } else {
            0.0
        };

        Some(DeanSiteEstimate {
            dean_number: dean_num,
            curvature_radius_m: curve_radius_m,
            arc_length_m,
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
    pub(super) fn spec_to_split_types(
        spec: &BlueprintTopologySpec,
    ) -> Result<Vec<SplitType>, String> {
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
    pub(super) fn representative_channel_width_m(spec: &BlueprintTopologySpec) -> f64 {
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
    pub(super) fn representative_channel_height_m(spec: &BlueprintTopologySpec) -> f64 {
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

    pub(super) fn render_hints_for_spec(spec: &BlueprintTopologySpec) -> BlueprintRenderHints {
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

    pub(super) fn mirror_blueprint_geometry(
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
