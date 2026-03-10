use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::{
    BlueprintTopologyFactory, BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec,
    SplitKind, SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
};

fn selective_spec() -> BlueprintTopologySpec {
    BlueprintTopologySpec {
        topology_id: "tri_tri_canonical".to_string(),
        design_name: "tri_tri_canonical".to_string(),
        box_dims_mm: (127.76, 85.47),
        inlet_width_m: 5.5e-3,
        outlet_width_m: 4.0e-3,
        trunk_length_m: 20.0e-3,
        outlet_tail_length_m: 14.0e-3,
        split_stages: vec![
            SplitStageSpec {
                stage_id: "stage0".to_string(),
                split_kind: SplitKind::Trifurcation,
                branches: vec![
                    BranchSpec {
                        label: "wbc".to_string(),
                        role: BranchRole::WbcCollection,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: 24.0e-3,
                            width_m: 1.0e-3,
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
                            width_m: 2.0e-3,
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
                            width_m: 1.0e-3,
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                    },
                ],
            },
            SplitStageSpec {
                stage_id: "stage1".to_string(),
                split_kind: SplitKind::Trifurcation,
                branches: vec![
                    BranchSpec {
                        label: "wbc".to_string(),
                        role: BranchRole::WbcCollection,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: 20.0e-3,
                            width_m: 0.45e-3,
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
                            length_m: 24.0e-3,
                            width_m: 1.3e-3,
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
                            length_m: 20.0e-3,
                            width_m: 0.25e-3,
                            height_m: 1.0e-3,
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                    },
                ],
            },
        ],
        venturi_placements: vec![VenturiPlacementSpec {
            placement_id: "tri_tri_venturi".to_string(),
            target_channel_id: BlueprintTopologySpec::branch_channel_id("stage1", "ctc"),
            serial_throat_count: 2,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m: 45e-6,
                throat_height_m: 1.0e-3,
                throat_length_m: 250e-6,
                inlet_width_m: 1.3e-3,
                outlet_width_m: 1.3e-3,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 9.0,
            },
            placement_mode: VenturiPlacementMode::StraightSegment,
        }],
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        treatment_mode: TreatmentActuationMode::VenturiCavitation,
    }
}

#[test]
fn selective_factory_build_uses_canonical_geometry_authoring() {
    let blueprint = BlueprintTopologyFactory::build(&selective_spec())
        .expect("selective spec should build through canonical geometry path");

    assert!(
        blueprint.is_geometry_authored(),
        "selective factory builds must carry create_geometry provenance"
    );
    assert!(
        blueprint.channels.iter().all(|channel| channel.path.len() >= 2),
        "report-grade channels must have explicit routed paths"
    );
    assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
    blueprint
        .validate()
        .expect("canonical selective blueprint should remain structurally valid");
}
