use crate::domain::model::NetworkBlueprint;
use crate::topology::model::{
    BlueprintTopologySpec, SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiConfig, VenturiPlacementMode,
};
use crate::topology::BlueprintTopologyFactory;

use super::super::modifiers::with_venturi;
use super::support::{
    apply_request_mirror, bypass_branch, milestone12_default_stage_layouts,
    mirror_blueprint_geometry, treatment_branch,
};
use super::{Milestone12PrimitiveSelectiveSpec, Milestone12TopologyRequest};

#[must_use]
pub fn milestone12_primitive_selective_tree_spec(
    request: &Milestone12PrimitiveSelectiveSpec,
) -> BlueprintTopologySpec {
    let stage_layouts = if request.stage_layouts.is_empty() {
        milestone12_default_stage_layouts(request)
    } else {
        request.stage_layouts.clone()
    };
    assert_eq!(
        stage_layouts.len(),
        request.split_kinds.len(),
        "Milestone 12 stage_layouts must match split_kinds length"
    );
    let mut split_stages = Vec::with_capacity(request.split_kinds.len());
    let mut outlet_treatment_width_m = request.inlet_width_m;

    for (stage_index, (split_kind, stage_layout)) in request
        .split_kinds
        .iter()
        .copied()
        .zip(stage_layouts.iter())
        .enumerate()
    {
        let is_last = stage_index + 1 == request.split_kinds.len();
        assert_eq!(
            split_kind.branch_count(),
            stage_layout.branches.len(),
            "Milestone 12 stage {stage_index} branch count must match {split_kind:?}"
        );
        assert_eq!(
            stage_layout.split_kind, split_kind,
            "Milestone 12 stage {stage_index} split_kind must match the sequence SSOT"
        );
        let serpentine = is_last.then_some(request.center_serpentine.clone()).flatten();
        let branches = stage_layout
            .branches
            .iter()
            .map(|branch| {
                if branch.treatment_path {
                    treatment_branch(
                        &branch.label,
                        branch.role,
                        branch.width_m,
                        request.channel_height_m,
                        request.branch_length_m,
                        serpentine.clone(),
                    )
                } else {
                    bypass_branch(
                        &branch.label,
                        branch.width_m,
                        request.channel_height_m,
                        request.branch_length_m,
                        branch.role,
                    )
                }
            })
            .collect::<Vec<_>>();

        outlet_treatment_width_m = branches
            .iter()
            .filter(|branch| branch.treatment_path)
            .map(|branch| branch.route.width_m)
            .sum::<f64>()
            .max(outlet_treatment_width_m.min(request.inlet_width_m));

        split_stages.push(SplitStageSpec {
            stage_id: format!("stage_{stage_index}"),
            split_kind,
            branches,
        });
    }

    let spec = BlueprintTopologySpec {
        topology_id: request.topology_id.clone(),
        design_name: request.design_name.clone(),
        box_dims_mm: request.box_dims_mm,
        inlet_width_m: request.inlet_width_m,
        outlet_width_m: outlet_treatment_width_m,
        trunk_length_m: request.branch_length_m,
        outlet_tail_length_m: request.outlet_tail_length_m,
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        split_stages,
        venturi_placements: Vec::new(),
        treatment_mode: request.treatment_mode,
    };

    if request.treatment_mode == TreatmentActuationMode::VenturiCavitation
        && request.venturi_throat_count > 0
    {
        return with_venturi(
            spec,
            VenturiConfig {
                target_channel_ids: request.venturi_target_channel_ids.clone(),
                serial_throat_count: request.venturi_throat_count,
                throat_geometry: ThroatGeometrySpec {
                    throat_width_m: request.venturi_throat_width_m,
                    throat_height_m: request.channel_height_m,
                    throat_length_m: request.venturi_throat_length_m,
                    inlet_width_m: outlet_treatment_width_m,
                    outlet_width_m: outlet_treatment_width_m,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                },
                placement_mode: request.venturi_placement_mode,
            },
        );
    }

    spec
}

/// Build the canonical Milestone 12 topology spec from a declarative request.
#[must_use]
pub fn build_milestone12_topology_spec(
    request: &Milestone12TopologyRequest,
) -> BlueprintTopologySpec {
    milestone12_primitive_selective_tree_spec(request)
}

/// Build a canonical Milestone 12 blueprint through the topology factory.
///
/// The returned blueprint is authored through the canonical `create_geometry`
/// pipeline and therefore preserves split-tree autolayout and geometry
/// provenance.
///
/// # Errors
///
/// Returns an error if the request yields an invalid topology or geometry.
pub fn build_milestone12_blueprint(
    request: &Milestone12TopologyRequest,
) -> Result<NetworkBlueprint, String> {
    let mut blueprint = BlueprintTopologyFactory::build(&build_milestone12_topology_spec(request))?;
    apply_request_mirror(&mut blueprint, request);
    Ok(blueprint)
}

/// Promote a canonical Milestone 12 Option 1 selective-routing blueprint into
/// a canonical Option 2 venturi-capable blueprint.
///
/// The returned blueprint is rebuilt through the topology factory so geometry
/// provenance and auto-layout remain authoritative.
///
/// # Errors
///
/// Returns an error when the source blueprint is not a geometry-authored
/// selective Milestone 12 design or when the treatment route cannot supply a
/// physically valid venturi throat geometry.
pub fn promote_milestone12_option1_to_option2(
    blueprint: &NetworkBlueprint,
    serial_throat_count: u8,
    placement_mode: VenturiPlacementMode,
) -> Result<NetworkBlueprint, String> {
    let topology = blueprint
        .topology_spec()
        .ok_or_else(|| {
            format!(
                "Milestone 12 promotion requires topology metadata on blueprint '{}'",
                blueprint.name
            )
        })?
        .clone();
    if !topology.is_selective_routing() {
        return Err(format!(
            "Milestone 12 promotion requires selective-routing topology, but '{}' is '{}'",
            blueprint.name,
            topology.stage_sequence_label()
        ));
    }
    if !blueprint.is_geometry_authored() {
        return Err(format!(
            "Milestone 12 promotion requires create_geometry provenance; '{}' is not geometry-authored",
            blueprint.name
        ));
    }

    let treatment_channel_ids = topology.treatment_channel_ids();
    let representative_id = treatment_channel_ids
        .first()
        .ok_or_else(|| {
            format!(
                "Milestone 12 promotion requires at least one treatment channel in '{}'",
                blueprint.name
            )
        })?
        .clone();
    let representative_route = topology.channel_route(&representative_id).ok_or_else(|| {
        format!(
            "Milestone 12 promotion could not resolve treatment channel '{}' in '{}'",
            representative_id, blueprint.name
        )
    })?;

    let throat_width_m =
        (representative_route.width_m * 0.4).clamp(60.0e-6, representative_route.width_m * 0.85);
    let throat_length_m =
        (representative_route.length_m / 8.0).clamp(300.0e-6, representative_route.length_m * 0.5);

    let throat_height_m = representative_route.height_m;
    let inlet_width_m = representative_route.width_m;

    let spec = with_venturi(
        topology.clone(),
        VenturiConfig {
            target_channel_ids: treatment_channel_ids,
            serial_throat_count: serial_throat_count.max(1),
            throat_geometry: ThroatGeometrySpec {
                throat_width_m,
                throat_height_m,
                throat_length_m,
                inlet_width_m,
                outlet_width_m: inlet_width_m,
                convergent_half_angle_deg: 7.0,
                divergent_half_angle_deg: 7.0,
            },
            placement_mode,
        },
    );

    let mut promoted = BlueprintTopologyFactory::build(&spec)?;
    if let Some(render_hints) = blueprint.render_hints() {
        mirror_blueprint_geometry(
            &mut promoted,
            topology.box_dims_mm,
            render_hints.mirror_x,
            render_hints.mirror_y,
        );
        if let Some(promoted_hints) = promoted.render_hints.as_mut() {
            promoted_hints.mirror_x = render_hints.mirror_x;
            promoted_hints.mirror_y = render_hints.mirror_y;
        }
    }
    if let Some(lineage) = promoted.lineage.as_mut() {
        lineage.option1_source_blueprint = Some(blueprint.name.clone());
        lineage.option2_source_blueprint = Some(promoted.name.clone());
    }
    Ok(promoted)
}