//! Mutation and validation methods for BlueprintTopologyFactory.
use crate::domain::therapy_metadata::TherapyZone;
use crate::domain::model::NetworkBlueprint;
use crate::topology::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec,
    SerpentineSpec, SplitKind, SplitStageSpec, TopologyOptimizationStage,
    TreatmentActuationMode, VenturiPlacementSpec,
};
use super::{BlueprintTopologyFactory, BlueprintTopologyMutation};

impl BlueprintTopologyFactory {
    /// Validate a spec without building (called by `NetworkBlueprint::validate`).
    pub fn validate_spec(spec: &BlueprintTopologySpec) -> Result<(), String> {
        super::super::validation::validate_spec(spec)
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
                    .ok_or_else(|| format!("Treatment channel '{target_channel_id}' not found"))?;
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
                            "treatment-channel serpentine update requires split-stage channel '{target_channel_id}'"
                        )
                    })?;
                let branch = new_spec
                    .split_stages
                    .iter_mut()
                    .find(|stage| stage.stage_id == stage_id)
                    .and_then(|stage| stage.branches.iter_mut().find(|branch| branch.label == branch_label))
                    .ok_or_else(|| {
                        format!(
                            "branch '{branch_label}' in stage '{stage_id}' not found for channel '{target_channel_id}'"
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
                            "split-merge insertion requires split-stage treatment channel '{target_channel_id}'"
                        )
                    })?;
                let stage_index = new_spec
                    .split_stages
                    .iter()
                    .position(|stage| stage.stage_id == stage_id)
                    .ok_or_else(|| format!("stage '{stage_id}' not found"))?;
                let parent_branch = new_spec.split_stages[stage_index]
                    .branches
                    .iter()
                    .find(|branch| branch.label == branch_label)
                    .ok_or_else(|| format!("branch '{branch_label}' not found in '{stage_id}'"))?
                    .clone();
                if !parent_branch.treatment_path || parent_branch.route.therapy_zone != TherapyZone::CancerTarget {
                    return Err(format!(
                        "split-merge insertion requires CancerTarget treatment channel '{target_channel_id}'"
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
        if let Some(existing_lineage) = blueprint.lineage.clone() {
            result.lineage = Some(existing_lineage);
        }
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

