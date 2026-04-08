//! Primitive Selective Tree (PST) geometry generation.
//!
//! Provides entry-points for generating PST network blueprints from either a
//! `PrimitiveSelectiveTreeRequest` or a declarative `BlueprintTopologySpec`.
use super::super::super::types::{Point2D, SplitType};
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::BlueprintRenderHints;
use crate::topology::{
    BlueprintTopologyFactory, BlueprintTopologySpec, SplitKind, TreatmentActuationMode,
};

mod annotation;

use annotation::annotate_primitive_tree;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimitiveSelectiveSplitKind {
    Bi,
    Tri,
    Quad,
    Penta,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrimitiveSelectiveTreeRequest {
    pub name: String,
    pub box_dims_mm: (f64, f64),
    pub split_sequence: Vec<PrimitiveSelectiveSplitKind>,
    pub main_width_m: f64,
    pub throat_width_m: f64,
    pub throat_length_m: f64,
    pub channel_height_m: f64,
    pub first_trifurcation_center_frac: f64,
    pub later_trifurcation_center_frac: f64,
    pub bifurcation_treatment_frac: f64,
    pub treatment_branch_venturi_enabled: bool,
    pub treatment_branch_throat_count: u8,
    pub center_serpentine: Option<super::super::CenterSerpentinePathSpec>,
}

pub fn create_primitive_selective_tree_geometry(
    request: &PrimitiveSelectiveTreeRequest,
) -> NetworkBlueprint {
    let splits: Vec<SplitType> = request
        .split_sequence
        .iter()
        .enumerate()
        .map(|(idx, split)| match split {
            PrimitiveSelectiveSplitKind::Bi => SplitType::Bifurcation,
            PrimitiveSelectiveSplitKind::Tri => SplitType::SymmetricTrifurcation {
                center_ratio: if idx == 0 {
                    request.first_trifurcation_center_frac
                } else {
                    request.later_trifurcation_center_frac
                },
            },
            PrimitiveSelectiveSplitKind::Quad => SplitType::Quadfurcation,
            PrimitiveSelectiveSplitKind::Penta => SplitType::Pentafurcation,
        })
        .collect();

    let geometry_config = GeometryConfig {
        wall_clearance: 1.0,
        channel_width: (request.main_width_m * 1.0e3).clamp(0.2, 12.0),
        channel_height: (request.channel_height_m * 1.0e3).clamp(0.2, 5.0),
        ..Default::default()
    };

    let split_depth = splits.len().max(1);
    let total_branches = splits
        .iter()
        .fold(1usize, |product, split| product.saturating_mul(split.branch_count()));
    let retry_scales = [1.0, 1.2, 1.45, 1.75, 2.1, 2.5];
    let mut last_blueprint = None;

    for retry_scale in retry_scales {
        let internal_dims = expanded_generation_box_dims(
            request.box_dims_mm,
            split_depth,
            total_branches,
            retry_scale,
        );
        let mut blueprint = super::super::create_geometry(
            internal_dims,
            &splits,
            &geometry_config,
            &ChannelTypeConfig::AllStraight,
        );
        if internal_dims != request.box_dims_mm {
            scale_blueprint_geometry(&mut blueprint, request.box_dims_mm);
        }
        annotate_primitive_tree(&mut blueprint, request);
        if blueprint.unresolved_channel_overlap_count() == 0 {
            return blueprint;
        }
        last_blueprint = Some(blueprint);
    }

    last_blueprint.expect("primitive selective tree generation should produce a blueprint")
}

fn expanded_generation_box_dims(
    target_dims: (f64, f64),
    split_depth: usize,
    total_branches: usize,
    retry_scale: f64,
) -> (f64, f64) {
    let depth_scale = (1.0 + 0.18 * split_depth.saturating_sub(1) as f64) * retry_scale;
    let branch_scale = ((total_branches as f64) / 9.0).sqrt().clamp(1.0, 3.2) * retry_scale;
    let width = target_dims.0 * depth_scale.max(1.0);
    let height = target_dims.1 * branch_scale.max(1.0);
    (width, height)
}

fn scale_blueprint_geometry(blueprint: &mut NetworkBlueprint, target_dims: (f64, f64)) {
    let source_dims = blueprint.box_dims;
    if source_dims.0 <= 0.0 || source_dims.1 <= 0.0 {
        return;
    }
    let scale_x = target_dims.0 / source_dims.0;
    let scale_y = target_dims.1 / source_dims.1;

    let scale_point = |(x, y): Point2D| (x * scale_x, y * scale_y);

    for node in &mut blueprint.nodes {
        node.point = scale_point(node.point);
    }
    for channel in &mut blueprint.channels {
        channel.path = channel.path.iter().copied().map(scale_point).collect();
    }
    blueprint.box_outline = blueprint
        .box_outline
        .iter()
        .map(|(start, end)| (scale_point(*start), scale_point(*end)))
        .collect();
    blueprint.box_dims = target_dims;
}

/// Build a primitive selective tree directly from a declarative topology spec
/// when the spec matches the canonical selective-routing contract.
pub fn create_primitive_selective_tree_geometry_from_spec(
    spec: &BlueprintTopologySpec,
) -> Result<Option<NetworkBlueprint>, String> {
    if spec.split_stages.is_empty() || spec.has_series_path() || spec.has_parallel_paths() {
        return Ok(None);
    }

    let mut split_sequence = Vec::with_capacity(spec.split_stages.len());
    let mut parent_width_m = spec.inlet_width_m;
    let mut first_trifurcation_center_frac = 0.5;
    let mut later_trifurcation_center_frac = 0.5;
    let mut bifurcation_treatment_frac = 0.5;
    let mut saw_later_trifurcation = false;
    let mut center_serpentine = None;

    for (stage_index, stage) in spec.split_stages.iter().enumerate() {
        let treatment_branches: Vec<_> = stage
            .branches
            .iter()
            .filter(|branch| branch.treatment_path)
            .collect();
        if treatment_branches.len() != 1 {
            return Ok(None);
        }
        let treatment_branch = treatment_branches[0];
        if parent_width_m <= 0.0 || treatment_branch.route.width_m <= 0.0 {
            return Ok(None);
        }

        let treatment_fraction = (treatment_branch.route.width_m / parent_width_m).clamp(0.0, 1.0);
        center_serpentine = center_serpentine.or(
            treatment_branch
                .route
                .serpentine
                .as_ref()
                .map(|serpentine| super::super::CenterSerpentinePathSpec {
                    segments: serpentine.segments,
                    bend_radius_m: serpentine.bend_radius_m,
                    wave_type: serpentine.wave_type,
                }),
        );

        match stage.split_kind {
            SplitKind::NFurcation(2) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Bi);
                bifurcation_treatment_frac = treatment_fraction;
            }
            SplitKind::NFurcation(3) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Tri);
                if stage_index == 0 {
                    first_trifurcation_center_frac = treatment_fraction;
                } else {
                    later_trifurcation_center_frac = treatment_fraction;
                    saw_later_trifurcation = true;
                }
            }
            SplitKind::NFurcation(4) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Quad);
            }
            SplitKind::NFurcation(5) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Penta);
            }
            SplitKind::NFurcation(_) => return Ok(None),
        }

        parent_width_m = treatment_branch.route.width_m;
    }

    if !saw_later_trifurcation {
        later_trifurcation_center_frac = first_trifurcation_center_frac;
    }

    let representative_height_m = spec
        .split_stages
        .first()
        .and_then(|stage| stage.branches.first())
        .map_or(1.0e-3, |branch| branch.route.height_m);
    let strongest_venturi = spec
        .venturi_placements
        .iter()
        .max_by_key(|placement| placement.serial_throat_count);

    let request = PrimitiveSelectiveTreeRequest {
        name: spec.design_name.clone(),
        box_dims_mm: spec.box_dims_mm,
        split_sequence,
        main_width_m: spec.inlet_width_m,
        throat_width_m: strongest_venturi
            .map_or(parent_width_m, |placement| placement.throat_geometry.throat_width_m),
        throat_length_m: strongest_venturi.map_or(spec.trunk_length_m / 8.0, |placement| {
            placement.throat_geometry.throat_length_m
        }),
        channel_height_m: representative_height_m,
        first_trifurcation_center_frac,
        later_trifurcation_center_frac,
        bifurcation_treatment_frac,
        treatment_branch_venturi_enabled: spec.treatment_mode
            == TreatmentActuationMode::VenturiCavitation
            && !spec.venturi_placements.is_empty(),
        treatment_branch_throat_count: strongest_venturi
            .map_or(0, |placement| placement.serial_throat_count),
        center_serpentine,
    };

    let mut blueprint = create_primitive_selective_tree_geometry(&request);
    blueprint.name.clone_from(&spec.design_name);
    blueprint.topology = Some(spec.clone());
    blueprint.lineage = Some(BlueprintTopologyFactory::lineage_for_spec(spec));
    blueprint.render_hints = Some(BlueprintRenderHints {
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
    });
    Ok(Some(blueprint))
}