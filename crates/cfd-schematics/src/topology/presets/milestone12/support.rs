use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::BlueprintRenderHints;
use crate::topology::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec, SerpentineSpec, SplitKind,
};

use super::{
    Milestone12PrimitiveSelectiveSpec, Milestone12StageBranchSpec, Milestone12StageLayout,
    Milestone12TopologyRequest,
};

pub(super) fn mirror_blueprint_geometry(
    blueprint: &mut NetworkBlueprint,
    box_dims_mm: (f64, f64),
    mirror_x: bool,
    mirror_y: bool,
) {
    if !mirror_x && !mirror_y {
        return;
    }

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

pub(super) fn apply_request_mirror(
    blueprint: &mut NetworkBlueprint,
    request: &Milestone12TopologyRequest,
) {
    mirror_blueprint_geometry(
        blueprint,
        request.box_dims_mm,
        request.mirror_x,
        request.mirror_y,
    );
    if let Some(hints) = blueprint.render_hints.as_mut() {
        hints.mirror_x = request.mirror_x;
        hints.mirror_y = request.mirror_y;
    } else {
        let topology = blueprint.topology_spec();
        blueprint.render_hints = Some(BlueprintRenderHints {
            stage_sequence: topology
                .map_or_else(String::new, BlueprintTopologySpec::stage_sequence_label),
            split_layers: topology.map_or(0, BlueprintTopologySpec::visible_split_layers),
            throat_count_hint: topology.map_or(0, BlueprintTopologySpec::venturi_count),
            treatment_label: if topology.is_some_and(BlueprintTopologySpec::has_venturi) {
                "venturi".to_string()
            } else {
                "ultrasound".to_string()
            },
            mirror_x: request.mirror_x,
            mirror_y: request.mirror_y,
        });
    }
}

pub(super) fn treatment_branch(
    label: &str,
    role: BranchRole,
    width_m: f64,
    height_m: f64,
    length_m: f64,
    serpentine: Option<SerpentineSpec>,
) -> BranchSpec {
    BranchSpec {
        label: label.to_string(),
        role,
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

pub(super) fn bypass_branch(
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

fn default_stage_branch_specs(
    split_kind: SplitKind,
    stage_index: usize,
    parent_width_m: f64,
    request: &Milestone12PrimitiveSelectiveSpec,
) -> Vec<Milestone12StageBranchSpec> {
    match split_kind {
        SplitKind::NFurcation(2) => {
            let treatment_frac = request.bifurcation_treatment_frac.clamp(0.15, 0.85);
            let treatment_width = parent_width_m * treatment_frac;
            let bypass_width = parent_width_m - treatment_width;
            vec![
                Milestone12StageBranchSpec {
                    label: "upper".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: treatment_width,
                },
                Milestone12StageBranchSpec {
                    label: "lower".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: bypass_width,
                },
            ]
        }
        SplitKind::NFurcation(3) => {
            let center_frac = if stage_index == 0 {
                request.first_trifurcation_center_frac
            } else {
                request.later_trifurcation_center_frac
            }
            .clamp(0.15, 0.70);
            let center_width = parent_width_m * center_frac;
            let side_width = (parent_width_m - center_width) * 0.5;
            vec![
                Milestone12StageBranchSpec {
                    label: "left".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: side_width,
                },
                Milestone12StageBranchSpec {
                    label: "center".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: center_width,
                },
                Milestone12StageBranchSpec {
                    label: "right".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: side_width,
                },
            ]
        }
        SplitKind::NFurcation(4) => {
            let center_frac = if stage_index == 0 {
                request.first_trifurcation_center_frac
            } else {
                request.later_trifurcation_center_frac
            }
            .clamp(0.15, 0.55);
            let treatment_width = parent_width_m * center_frac;
            let bypass_width = (parent_width_m - treatment_width) / 3.0;
            vec![
                Milestone12StageBranchSpec {
                    label: "arm_0".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: bypass_width,
                },
                Milestone12StageBranchSpec {
                    label: "arm_1".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: treatment_width,
                },
                Milestone12StageBranchSpec {
                    label: "arm_2".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: bypass_width,
                },
                Milestone12StageBranchSpec {
                    label: "arm_3".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: bypass_width,
                },
            ]
        }
        SplitKind::NFurcation(5) => {
            let center_frac = if stage_index == 0 {
                request.first_trifurcation_center_frac
            } else {
                request.later_trifurcation_center_frac
            }
            .clamp(0.15, 0.55);
            let center_width = parent_width_m * center_frac;
            let bypass_width = (parent_width_m - center_width) / 4.0;
            (0..5)
                .map(|index| {
                    if index == 2 {
                        Milestone12StageBranchSpec {
                            label: "center".to_string(),
                            role: BranchRole::Treatment,
                            treatment_path: true,
                            width_m: center_width,
                        }
                    } else {
                        let role = if index < 2 {
                            BranchRole::Neutral
                        } else {
                            BranchRole::RbcBypass
                        };
                        Milestone12StageBranchSpec {
                            label: format!("arm_{index}"),
                            role,
                            treatment_path: false,
                            width_m: bypass_width,
                        }
                    }
                })
                .collect()
        }
        SplitKind::NFurcation(other) => {
            panic!("Milestone 12 primitive selective tree only supports Bi/Tri/Quad/Penta, got N={other}");
        }
    }
}

#[must_use]
pub fn milestone12_default_stage_layouts(
    request: &Milestone12PrimitiveSelectiveSpec,
) -> Vec<Milestone12StageLayout> {
    let mut parent_width_m = request.inlet_width_m;
    request
        .split_kinds
        .iter()
        .copied()
        .enumerate()
        .map(|(stage_index, split_kind)| {
            let branches =
                default_stage_branch_specs(split_kind, stage_index, parent_width_m, request);
            let treatment_width_m: f64 = branches
                .iter()
                .filter(|branch| branch.treatment_path)
                .map(|branch| branch.width_m)
                .sum();
            parent_width_m = treatment_width_m.max(f64::EPSILON);
            Milestone12StageLayout {
                split_kind,
                branches,
            }
        })
        .collect()
}
