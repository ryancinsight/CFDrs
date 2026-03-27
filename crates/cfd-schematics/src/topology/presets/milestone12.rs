use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::BlueprintRenderHints;
use crate::topology::BlueprintTopologyFactory;

use super::super::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec, SerpentineSpec, SplitKind,
    SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode, VenturiConfig,
    VenturiPlacementMode,
};
use super::helpers::{PLATE_HEIGHT_MM, PLATE_WIDTH_MM};
use super::modifiers::with_venturi;
use super::sequence::MILESTONE12_SWEEP_SEQUENCES;

#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageBranchSpec {
    pub label: String,
    pub role: BranchRole,
    pub treatment_path: bool,
    pub width_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageLayout {
    pub split_kind: SplitKind,
    pub branches: Vec<Milestone12StageBranchSpec>,
}

#[derive(Debug, Clone)]
pub struct Milestone12PrimitiveSelectiveSpec {
    pub topology_id: String,
    pub design_name: String,
    pub mirror_x: bool,
    pub mirror_y: bool,
    pub box_dims_mm: (f64, f64),
    pub split_kinds: Vec<SplitKind>,
    pub inlet_width_m: f64,
    pub channel_height_m: f64,
    pub branch_length_m: f64,
    pub outlet_tail_length_m: f64,
    pub stage_layouts: Vec<Milestone12StageLayout>,
    pub first_trifurcation_center_frac: f64,
    pub later_trifurcation_center_frac: f64,
    pub bifurcation_treatment_frac: f64,
    pub treatment_mode: TreatmentActuationMode,
    pub venturi_throat_count: u8,
    pub venturi_throat_width_m: f64,
    pub venturi_throat_length_m: f64,
    pub center_serpentine: Option<SerpentineSpec>,
    pub venturi_placement_mode: VenturiPlacementMode,
    pub venturi_target_channel_ids: Vec<String>,
}

impl Milestone12PrimitiveSelectiveSpec {
    #[must_use]
    pub fn new(
        topology_id: impl Into<String>,
        design_name: impl Into<String>,
        split_kinds: Vec<SplitKind>,
        inlet_width_m: f64,
        channel_height_m: f64,
        branch_length_m: f64,
        outlet_tail_length_m: f64,
    ) -> Self {
        Self {
            topology_id: topology_id.into(),
            design_name: design_name.into(),
            mirror_x: false,
            mirror_y: false,
            box_dims_mm: (PLATE_WIDTH_MM, PLATE_HEIGHT_MM),
            split_kinds,
            inlet_width_m,
            channel_height_m,
            branch_length_m,
            outlet_tail_length_m,
            stage_layouts: Vec::new(),
            first_trifurcation_center_frac: 0.45,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
            venturi_throat_count: 0,
            venturi_throat_width_m: inlet_width_m,
            venturi_throat_length_m: branch_length_m / 8.0,
            center_serpentine: None,
            venturi_placement_mode: VenturiPlacementMode::StraightSegment,
            venturi_target_channel_ids: Vec::new(),
        }
    }
}

/// Canonical declarative Milestone 12 topology request.
///
/// This is the single-source request shape consumed by
/// [`build_milestone12_blueprint`] and
/// [`build_milestone12_topology_spec`].
pub type Milestone12TopologyRequest = Milestone12PrimitiveSelectiveSpec;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MirrorVariant {
    suffix: &'static str,
    mirror_x: bool,
    mirror_y: bool,
}

const MILESTONE12_MIRROR_VARIANTS: [MirrorVariant; 4] = [
    MirrorVariant {
        suffix: "base",
        mirror_x: false,
        mirror_y: false,
    },
    MirrorVariant {
        suffix: "x",
        mirror_x: true,
        mirror_y: false,
    },
    MirrorVariant {
        suffix: "y",
        mirror_x: false,
        mirror_y: true,
    },
    MirrorVariant {
        suffix: "xy",
        mirror_x: true,
        mirror_y: true,
    },
];

fn mirror_blueprint_geometry(
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

fn apply_request_mirror(
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
            stage_sequence: topology.map_or_else(String::new, super::super::model::BlueprintTopologySpec::stage_sequence_label),
            split_layers: topology.map_or(0, super::super::model::BlueprintTopologySpec::visible_split_layers),
            throat_count_hint: topology.map_or(0, super::super::model::BlueprintTopologySpec::venturi_count),
            treatment_label: if topology.is_some_and(super::super::model::BlueprintTopologySpec::has_venturi) {
                "venturi".to_string()
            } else {
                "ultrasound".to_string()
            },
            mirror_x: request.mirror_x,
            mirror_y: request.mirror_y,
        });
    }
}

fn treatment_branch(
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
        let serpentine = is_last
            .then_some(request.center_serpentine.clone())
            .flatten();
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

/// Enumerate the canonical Milestone 12 split-tree scaffolds.
///
/// The catalog fixes only the split-stage lineage; callers remain responsible
/// for sweeping operating conditions and Milestone-12-specific geometric
/// parameters such as asymmetric width fractions, venturi throat geometry, and
/// terminal serpentine settings.
#[must_use]
pub fn enumerate_milestone12_topologies() -> Vec<Milestone12TopologyRequest> {
    let mut requests = Vec::with_capacity(
        MILESTONE12_SWEEP_SEQUENCES.len() * MILESTONE12_MIRROR_VARIANTS.len(),
    );
    for seq in &MILESTONE12_SWEEP_SEQUENCES {
        let split_kinds = seq.to_split_kinds();
        let family = seq.label().replace('\u{2192}', "");
        for variant in &MILESTONE12_MIRROR_VARIANTS {
            let mut request = Milestone12TopologyRequest::new(
                format!("pst-{}-{}", family.to_ascii_lowercase(), variant.suffix),
                format!("{family}-{}", variant.suffix.to_ascii_uppercase()),
                split_kinds.clone(),
                6.0e-3,
                1.0e-3,
                8.0e-3,
                8.0e-3,
            );
            request.mirror_x = variant.mirror_x;
            request.mirror_y = variant.mirror_y;
            requests.push(request);
        }
    }
    requests
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{BlueprintTopologyFactory, BlueprintTopologyMutation, TopologyOptimizationStage};

    #[test]
    fn default_stage_layouts_cover_bi_tri_quad_penta_roots() {
        let request = Milestone12PrimitiveSelectiveSpec::new(
            "m12",
            "m12",
            vec![
                SplitKind::NFurcation(2),
                SplitKind::NFurcation(3),
                SplitKind::NFurcation(4),
                SplitKind::NFurcation(5),
            ],
            6.0e-3,
            1.0e-3,
            8.0e-3,
            8.0e-3,
        );
        let layouts = milestone12_default_stage_layouts(&request);
        assert_eq!(layouts.len(), 4);
        assert_eq!(layouts[0].branches.len(), 2);
        assert_eq!(layouts[1].branches.len(), 3);
        assert_eq!(layouts[2].branches.len(), 4);
        assert_eq!(layouts[3].branches.len(), 5);
    }

    #[test]
    fn explicit_even_way_treatment_lane_sets_are_preserved() {
        let mut request = Milestone12PrimitiveSelectiveSpec::new(
            "quad",
            "quad",
            vec![SplitKind::NFurcation(4)],
            6.0e-3,
            1.0e-3,
            8.0e-3,
            8.0e-3,
        );
        request.stage_layouts = vec![Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(4),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "arm_0".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: 1.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "arm_1".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 2.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "arm_2".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 2.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "arm_3".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: 1.0e-3,
                },
            ],
        }];
        request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
        request.venturi_throat_count = 1;
        request.venturi_target_channel_ids = vec![
            BlueprintTopologySpec::branch_channel_id("stage_0", "arm_1"),
            BlueprintTopologySpec::branch_channel_id("stage_0", "arm_2"),
        ];

        let spec = milestone12_primitive_selective_tree_spec(&request);
        assert_eq!(spec.treatment_channel_ids().len(), 2);
        assert_eq!(spec.venturi_placements.len(), 2);
    }

    #[test]
    fn asymmetric_branch_widths_are_conserved() {
        let mut request = Milestone12PrimitiveSelectiveSpec::new(
            "tri",
            "tri",
            vec![SplitKind::NFurcation(3)],
            6.0e-3,
            1.0e-3,
            8.0e-3,
            8.0e-3,
        );
        request.stage_layouts = vec![Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "left".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: 1.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "center".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 2.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "right".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: 3.0e-3,
                },
            ],
        }];

        let spec = milestone12_primitive_selective_tree_spec(&request);
        let total: f64 = spec.split_stages[0]
            .branches
            .iter()
            .map(|branch| branch.route.width_m)
            .sum();
        assert!((total - request.inlet_width_m).abs() < 1.0e-12);
    }

    #[test]
    fn enumerate_milestone12_topologies_covers_mirrored_root_catalog() {
        let catalog = enumerate_milestone12_topologies();
        assert_eq!(
            catalog.len(),
            MILESTONE12_SWEEP_SEQUENCES.len() * MILESTONE12_MIRROR_VARIANTS.len(),
        );
        // Single-stage roots
        assert!(catalog.iter().any(|request| request.design_name == "Bi-BASE"));
        assert!(catalog.iter().any(|request| request.design_name == "Quad-Y"));
        assert!(catalog.iter().any(|request| request.design_name == "Penta-XY"));
        // Multi-layer Quad/Penta cascades
        assert!(catalog.iter().any(|request| request.design_name == "QuadTriBi-BASE"));
        assert!(catalog.iter().any(|request| request.design_name == "PentaQuadBi-BASE"));
        assert!(catalog.iter().any(|request| request.design_name == "PentaQuadTri-XY"));
        assert!(catalog.iter().any(|request| request.design_name == "PentaTriBi-XY"));
        assert!(catalog.iter().any(|request| request.design_name == "QuadBi-X"));
        assert!(catalog.iter().any(|request| request.design_name == "PentaTri-Y"));
    }

    #[test]
    fn build_milestone12_blueprint_uses_canonical_geometry_authoring() {
        let request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Tri-XY")
            .expect("mirrored Tri scaffold should exist");

        let blueprint =
            build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");

        assert!(blueprint.is_geometry_authored());
        assert!(blueprint.render_hints().is_some_and(|hints| hints.mirror_x && hints.mirror_y));
        blueprint
            .validate()
            .expect("Milestone 12 blueprint should validate");
    }

    #[test]
    fn promote_option1_to_option2_rebuilds_geometry_authored_blueprint() {
        let mut request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Quad-Y")
            .expect("mirrored Quad scaffold should exist");
        request.treatment_mode = TreatmentActuationMode::UltrasoundOnly;
        request.venturi_throat_count = 0;
        let option1 =
            build_milestone12_blueprint(&request).expect("Milestone 12 Option 1 should build");

        let promoted = promote_milestone12_option1_to_option2(
            &option1,
            1,
            VenturiPlacementMode::StraightSegment,
        )
        .expect("promotion should succeed");

        assert!(promoted.is_geometry_authored());
        assert!(promoted
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi));
        assert!(promoted.render_hints().is_some_and(|hints| !hints.mirror_x && hints.mirror_y));
        assert_eq!(
            promoted
                .lineage()
                .and_then(|lineage| lineage.option1_source_blueprint.as_deref()),
            Some(option1.name.as_str())
        );
    }

    #[test]
    fn selective_mutation_rejects_non_geometry_authored_blueprint() {
        let request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Tri-BASE")
            .expect("Tri scaffold should exist");
        let mut blueprint =
            build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");
        blueprint.metadata = None;
        blueprint.geometry_authored = false;

        let error = BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::UpdateBranchWidth {
                stage_id: "stage_0".to_string(),
                branch_label: "center".to_string(),
                new_width_m: 1.8e-3,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect_err("non-canonical selective mutation must fail");
        assert!(error.contains("create_geometry-authored provenance"));
    }
}
