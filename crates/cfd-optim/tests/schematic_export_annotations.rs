use cfd_optim::save_blueprint_schematic_svg;
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::{
    BlueprintTopologyFactory, BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec,
    NetworkBlueprint, SplitKind, SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_svg_path(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}.svg"))
}

fn selective_blueprint(
    id: &str,
    split_kinds: &[SplitKind],
    treatment_mode: TreatmentActuationMode,
) -> NetworkBlueprint {
    let split_stages = split_kinds
        .iter()
        .enumerate()
        .map(|(index, split_kind)| SplitStageSpec {
            stage_id: format!("stage{index}"),
            split_kind: *split_kind,
            branches: match split_kind {
                SplitKind::Bifurcation => vec![
                    BranchSpec {
                        label: "ctc".to_string(),
                        role: BranchRole::Treatment,
                        treatment_path: true,
                        route: ChannelRouteSpec {
                            length_m: 22.0e-3,
                            width_m: if index == 0 { 1.4e-3 } else { 1.0e-3 },
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::CancerTarget,
                        },
                    },
                    BranchSpec {
                        label: "waste".to_string(),
                        role: BranchRole::Neutral,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: 18.0e-3,
                            width_m: if index == 0 { 0.6e-3 } else { 0.7e-3 },
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                    },
                ],
                SplitKind::Trifurcation => vec![
                    BranchSpec {
                        label: "wbc".to_string(),
                        role: BranchRole::WbcCollection,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: 24.0e-3,
                            width_m: if index == 0 { 1.0e-3 } else { 0.45e-3 },
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                    },
                    BranchSpec {
                        label: "ctc".to_string(),
                        role: BranchRole::Treatment,
                        treatment_path: true,
                        route: ChannelRouteSpec {
                            length_m: 28.0e-3,
                            width_m: if index == 0 { 2.0e-3 } else { 1.0e-3 },
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::CancerTarget,
                        },
                    },
                    BranchSpec {
                        label: "rbc".to_string(),
                        role: BranchRole::RbcBypass,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: 24.0e-3,
                            width_m: if index == 0 { 1.0e-3 } else { 0.45e-3 },
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                    },
                ],
            },
        })
        .collect::<Vec<_>>();
    let last_stage_id = format!("stage{}", split_kinds.len().saturating_sub(1));
    let venturi_placements = if treatment_mode == TreatmentActuationMode::VenturiCavitation {
        vec![VenturiPlacementSpec {
            placement_id: format!("{id}-venturi"),
            target_channel_id: BlueprintTopologySpec::branch_channel_id(&last_stage_id, "ctc"),
            serial_throat_count: 2,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m: 45e-6,
                throat_height_m: 1.0e-3,
                throat_length_m: 250e-6,
                inlet_width_m: 1.4e-3,
                outlet_width_m: 1.4e-3,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 9.0,
            },
            placement_mode: VenturiPlacementMode::StraightSegment,
        }]
    } else {
        Vec::new()
    };
    let spec = BlueprintTopologySpec {
        topology_id: id.to_string(),
        design_name: id.to_string(),
        box_dims_mm: (127.76, 85.47),
        inlet_width_m: 5.0e-3,
        outlet_width_m: 4.0e-3,
        trunk_length_m: 20.0e-3,
        outlet_tail_length_m: 14.0e-3,
        split_stages,
        venturi_placements,
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        treatment_mode,
    };
    BlueprintTopologyFactory::build(&spec).expect("spec should build")
}

#[test]
fn save_schematic_svg_adds_throat_markers_for_venturi_blueprint() {
    let blueprint = selective_blueprint(
        "venturi_test",
        &[SplitKind::Trifurcation],
        TreatmentActuationMode::VenturiCavitation,
    );
    let out = unique_svg_path("cfd_optim_venturi_blueprint");

    save_blueprint_schematic_svg(&blueprint, &out)
        .expect("venturi blueprint schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("IN"));
    assert!(svg.contains("OUT"));
    assert!(svg.contains("TH1"));
}

#[test]
fn save_schematic_svg_keeps_node_markers_without_throats_for_nonventuri_blueprint() {
    let blueprint = selective_blueprint(
        "nonventuri_test",
        &[SplitKind::Trifurcation],
        TreatmentActuationMode::UltrasoundOnly,
    );
    let out = unique_svg_path("cfd_optim_nonventuri_blueprint");

    save_blueprint_schematic_svg(&blueprint, &out)
        .expect("non-venturi blueprint schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("IN"));
    assert!(svg.contains("OUT"));
    assert!(!svg.contains("TH1"));
}

#[test]
fn save_schematic_svg_renders_selective_blueprint_topology_from_blueprint_ssot() {
    let blueprint = selective_blueprint(
        "selective_cif",
        &[SplitKind::Trifurcation, SplitKind::Bifurcation],
        TreatmentActuationMode::VenturiCavitation,
    );
    let out = unique_svg_path("cfd_optim_selective_cif");

    save_blueprint_schematic_svg(&blueprint, &out)
        .expect("selective blueprint schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("TH1"));
    assert!(svg.contains("S1"));
    assert!(svg.contains("S2"));
}

#[test]
fn save_schematic_svg_renders_full_tree_split_markers_for_dtcv() {
    let blueprint = selective_blueprint(
        "selective_dtcv",
        &[SplitKind::Trifurcation, SplitKind::Trifurcation],
        TreatmentActuationMode::VenturiCavitation,
    );
    let out = unique_svg_path("cfd_optim_selective_dtcv");

    assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
    blueprint
        .validate()
        .expect("report-grade selective venturi blueprint must be planar before export");

    save_blueprint_schematic_svg(&blueprint, &out).expect("dtcv blueprint export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("TH1"));
    assert!(svg.contains("S2"));
}

#[test]
fn save_schematic_svg_renders_acoustic_dtcv_without_throats() {
    let blueprint = selective_blueprint(
        "selective_dtcv_acoustic",
        &[SplitKind::Trifurcation, SplitKind::Trifurcation],
        TreatmentActuationMode::UltrasoundOnly,
    );
    let out = unique_svg_path("cfd_optim_selective_dtcv_acoustic");

    save_blueprint_schematic_svg(&blueprint, &out)
        .expect("acoustic dtcv blueprint export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("S2"));
    assert!(!svg.contains("TH1"));
}

#[test]
fn save_schematic_svg_rejects_non_geometry_authored_blueprint() {
    let mut blueprint = selective_blueprint(
        "selective_noncanonical",
        &[SplitKind::Trifurcation, SplitKind::Trifurcation],
        TreatmentActuationMode::VenturiCavitation,
    );
    let out = unique_svg_path("cfd_optim_noncanonical");
    blueprint.metadata = None;

    let error = save_blueprint_schematic_svg(&blueprint, &out)
        .expect_err("canonical export must reject non-geometry-authored blueprints");
    assert!(
        error
            .to_string()
            .contains("geometry-authored blueprint provenance"),
        "unexpected error: {error}"
    );
}
