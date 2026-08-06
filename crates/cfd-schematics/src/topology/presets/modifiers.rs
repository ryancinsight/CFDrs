//! Builder modifiers for topology specs (venturi, serpentine, Dean placement).

use super::super::model::{
    BlueprintTopologySpec, SerpentineSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiConfig, VenturiPlacementMode, VenturiPlacementSpec,
};
use super::plate_presets::default_venturi_half_angle;
use aequitas::systems::si::quantities::Length;

fn resolved_venturi_geometry(
    spec: &BlueprintTopologySpec,
    channel_id: &str,
    geometry: &ThroatGeometrySpec,
) -> ThroatGeometrySpec {
    let resolved_width_m = spec
        .channel_route(channel_id)
        .map_or(geometry.throat_width_m * 3.0, |route| route.width_m);

    ThroatGeometrySpec {
        throat_width_m: geometry.throat_width_m,
        throat_height_m: geometry.throat_height_m,
        throat_length_m: geometry.throat_length_m,
        inlet_width_m: if geometry.inlet_width_m.into_base() > 0.0 {
            geometry.inlet_width_m
        } else {
            resolved_width_m
        },
        outlet_width_m: if geometry.outlet_width_m.into_base() > 0.0 {
            geometry.outlet_width_m
        } else {
            resolved_width_m
        },
        convergent_half_angle: if geometry.convergent_half_angle.into_base() > 0.0 {
            geometry.convergent_half_angle
        } else {
            default_venturi_half_angle()
        },
        divergent_half_angle: if geometry.divergent_half_angle.into_base() > 0.0 {
            geometry.divergent_half_angle
        } else {
            default_venturi_half_angle()
        },
    }
}

/// Builder: attach a generic venturi augmentation to a topology spec.
///
/// When `config.target_channel_ids` is empty, the modifier attaches the
/// venturi to every treatment-path channel in the spec.
#[must_use]
pub fn with_venturi(
    mut spec: BlueprintTopologySpec,
    config: VenturiConfig,
) -> BlueprintTopologySpec {
    let target_channel_ids = if config.target_channel_ids.is_empty() {
        spec.treatment_channel_ids()
    } else {
        config.target_channel_ids
    };

    spec.venturi_placements = target_channel_ids
        .into_iter()
        .enumerate()
        .map(|(idx, channel_id)| VenturiPlacementSpec {
            placement_id: format!("venturi_{idx}"),
            target_channel_id: channel_id.clone(),
            serial_throat_count: config.serial_throat_count,
            throat_geometry: resolved_venturi_geometry(&spec, &channel_id, &config.throat_geometry),
            placement_mode: config.placement_mode,
        })
        .collect();
    spec.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    spec
}

/// Builder: attach venturi throat placements to a spec (Option 1 → Option 2).
///
/// Targets all treatment-path channels with identical throat geometry.
/// Returns a new spec with `treatment_mode = VenturiCavitation`.
#[must_use]
pub fn with_venturi_placements(
    spec: BlueprintTopologySpec,
    throat_width_m: f64,
    throat_height_m: f64,
    throat_length_m: f64,
    serial_throat_count: u8,
    placement_mode: VenturiPlacementMode,
) -> BlueprintTopologySpec {
    with_venturi(
        spec,
        VenturiConfig {
            target_channel_ids: Vec::new(),
            serial_throat_count,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m: Length::from_base(throat_width_m),
                throat_height_m: Length::from_base(throat_height_m),
                throat_length_m: Length::from_base(throat_length_m),
                inlet_width_m: Length::from_base(0.0),
                outlet_width_m: Length::from_base(0.0),
                convergent_half_angle: default_venturi_half_angle(),
                divergent_half_angle: default_venturi_half_angle(),
            },
            placement_mode,
        },
    )
}

/// Builder: add serpentine path to a treatment-path branch (Option 3 mutation).
///
/// Modifies the specified branch's route to include a serpentine specification.
/// Returns `None` if the stage/branch combination is not found.
#[must_use]
pub fn with_branch_serpentine(
    mut spec: BlueprintTopologySpec,
    stage_id: &str,
    branch_label: &str,
    segments: usize,
    bend_radius_m: f64,
    segment_length_m: f64,
) -> Option<BlueprintTopologySpec> {
    let stage = spec
        .split_stages
        .iter_mut()
        .find(|s| s.stage_id == stage_id)?;
    let branch = stage
        .branches
        .iter_mut()
        .find(|b| b.label == branch_label)?;
    branch.route.serpentine = Some(SerpentineSpec {
        segments,
        bend_radius_m: Length::from_base(bend_radius_m),
        segment_length_m: Length::from_base(segment_length_m),
        wave_type: crate::topology::SerpentineWaveType::Sine,
    });
    Some(spec)
}

/// Builder: relocate venturi placements to Dean vortex peak sites.
///
/// Changes all venturi placement modes to
/// [`VenturiPlacementMode::CurvaturePeakDeanNumber`], placing throats at
/// the apex of serpentine bends where secondary Dean circulation maximises
/// curvature-induced inertial focusing of stiff CTCs.
#[must_use]
pub fn with_dean_venturi_placement(mut spec: BlueprintTopologySpec) -> BlueprintTopologySpec {
    for placement in &mut spec.venturi_placements {
        placement.placement_mode = VenturiPlacementMode::CurvaturePeakDeanNumber;
    }
    spec
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::therapy_metadata::TherapyZone;
    use crate::topology::model::{
        BranchRole, BranchSpec, ChannelRouteSpec, ParallelChannelSpec, SplitKind, SplitStageSpec,
    };
    use aequitas::systems::si::quantities::Angle;

    #[test]
    fn venturi_modifier_targets_explicit_parallel_channel_ids() {
        let spec = BlueprintTopologySpec {
            topology_id: "parallel".to_string(),
            design_name: "parallel".to_string(),
            box_dims_mm: (127.76, 85.47),
            inlet_width_m: 2.0e-3,
            outlet_width_m: 2.0e-3,
            trunk_length_m: 12.0e-3,
            outlet_tail_length_m: 12.0e-3,
            series_channels: Vec::new(),
            parallel_channels: vec![
                ParallelChannelSpec {
                    channel_id: "straight_lane".to_string(),
                    route: ChannelRouteSpec {
                        length_m: Length::from_base(10.0e-3),
                        width_m: Length::from_base(2.0e-3),
                        height_m: Length::from_base(1.0e-3),
                        serpentine: None,
                        therapy_zone: TherapyZone::CancerTarget,
                    },
                },
                ParallelChannelSpec {
                    channel_id: "serpentine_lane".to_string(),
                    route: ChannelRouteSpec {
                        length_m: Length::from_base(18.0e-3),
                        width_m: Length::from_base(1.5e-3),
                        height_m: Length::from_base(1.0e-3),
                        serpentine: Some(SerpentineSpec {
                            wave_type: crate::topology::SerpentineWaveType::Sine,
                            segments: 4,
                            bend_radius_m: Length::from_base(1.2e-3),
                            segment_length_m: Length::from_base(4.5e-3),
                        }),
                        therapy_zone: TherapyZone::CancerTarget,
                    },
                },
            ],
            split_stages: Vec::new(),
            venturi_placements: Vec::new(),
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
        };

        let updated = with_venturi(
            spec,
            VenturiConfig {
                target_channel_ids: vec![
                    "straight_lane".to_string(),
                    "serpentine_lane".to_string(),
                ],
                serial_throat_count: 2,
                throat_geometry: ThroatGeometrySpec {
                    throat_width_m: Length::from_base(80.0e-6),
                    throat_height_m: Length::from_base(1.0e-3),
                    throat_length_m: Length::from_base(300.0e-6),
                    inlet_width_m: Length::from_base(0.0),
                    outlet_width_m: Length::from_base(0.0),
                    convergent_half_angle: Angle::from_base(0.0),
                    divergent_half_angle: Angle::from_base(0.0),
                },
                placement_mode: VenturiPlacementMode::StraightSegment,
            },
        );

        assert_eq!(updated.venturi_placements.len(), 2);
        assert_eq!(
            updated.venturi_placements[0].target_channel_id,
            "straight_lane"
        );
        assert!(
            (updated.venturi_placements[1]
                .throat_geometry
                .inlet_width_m
                .into_base()
                - 1.5e-3)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn typed_venturi_geometry_round_trips_si_base_units() {
        let geometry = ThroatGeometrySpec {
            throat_width_m: Length::from_base(80.0e-6),
            throat_height_m: Length::from_base(1.0e-3),
            throat_length_m: Length::from_base(300.0e-6),
            inlet_width_m: Length::from_base(1.5e-3),
            outlet_width_m: Length::from_base(1.5e-3),
            convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
            divergent_half_angle: Angle::from_base(9.0_f64.to_radians()),
        };

        let encoded = serde_json::to_string(&geometry).expect("serialize typed geometry");
        let decoded: ThroatGeometrySpec =
            serde_json::from_str(&encoded).expect("deserialize typed geometry");

        assert!((decoded.throat_width_m.into_base() - 80.0e-6).abs() <= f64::EPSILON);
        assert!((decoded.throat_height_m.into_base() - 1.0e-3).abs() <= f64::EPSILON);
        assert!((decoded.throat_length_m.into_base() - 300.0e-6).abs() <= f64::EPSILON);
        assert!((decoded.inlet_width_m.into_base() - 1.5e-3).abs() <= f64::EPSILON);
        assert!((decoded.outlet_width_m.into_base() - 1.5e-3).abs() <= f64::EPSILON);
        // JSON performs one binary64 decimal round-trip; one machine epsilon
        // at the value's scale bounds the resulting representation drift.
        let angle_tolerance = f64::EPSILON * 7.0_f64.to_radians().abs().max(1.0);
        assert!(
            (decoded.convergent_half_angle.into_base() - 7.0_f64.to_radians()).abs()
                <= angle_tolerance
        );
        assert!(
            (decoded.divergent_half_angle.into_base() - 9.0_f64.to_radians()).abs() <= f64::EPSILON
        );
    }

    #[test]
    fn typed_channel_route_round_trips_si_base_units() {
        let route = ChannelRouteSpec {
            length_m: Length::from_base(12.0e-3),
            width_m: Length::from_base(1.5e-3),
            height_m: Length::from_base(0.75e-3),
            serpentine: None,
            therapy_zone: TherapyZone::CancerTarget,
        };

        let encoded = serde_json::to_string(&route).expect("serialize typed channel route");
        let decoded: ChannelRouteSpec =
            serde_json::from_str(&encoded).expect("deserialize typed channel route");

        // JSON performs one binary64 decimal round-trip; one machine epsilon
        // at the value scale bounds the resulting representation drift.
        let tolerance = |value: f64| f64::EPSILON * value.abs().max(1.0);
        assert!((decoded.length_m.into_base() - 12.0e-3).abs() <= tolerance(12.0e-3));
        assert!((decoded.width_m.into_base() - 1.5e-3).abs() <= tolerance(1.5e-3));
        assert!((decoded.height_m.into_base() - 0.75e-3).abs() <= tolerance(0.75e-3));
        assert_eq!(decoded.serpentine, route.serpentine);
        assert_eq!(decoded.therapy_zone, route.therapy_zone);
    }

    #[test]
    fn branch_serpentine_modifier_updates_named_branch() {
        let spec = BlueprintTopologySpec {
            topology_id: "split".to_string(),
            design_name: "split".to_string(),
            box_dims_mm: (127.76, 85.47),
            inlet_width_m: 6.0e-3,
            outlet_width_m: 2.0e-3,
            trunk_length_m: 8.0e-3,
            outlet_tail_length_m: 8.0e-3,
            series_channels: Vec::new(),
            parallel_channels: Vec::new(),
            split_stages: vec![SplitStageSpec {
                stage_id: "stage_0".to_string(),
                split_kind: SplitKind::NFurcation(3),
                branches: vec![
                    BranchSpec {
                        label: "left".to_string(),
                        role: BranchRole::WbcCollection,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: Length::from_base(8.0e-3),
                            width_m: Length::from_base(2.0e-3),
                            height_m: Length::from_base(1.0e-3),
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                        recovery_sub_split: None,
                    },
                    BranchSpec {
                        label: "center".to_string(),
                        role: BranchRole::Treatment,
                        treatment_path: true,
                        route: ChannelRouteSpec {
                            length_m: Length::from_base(8.0e-3),
                            width_m: Length::from_base(2.0e-3),
                            height_m: Length::from_base(1.0e-3),
                            serpentine: None,
                            therapy_zone: TherapyZone::CancerTarget,
                        },
                        recovery_sub_split: None,
                    },
                    BranchSpec {
                        label: "right".to_string(),
                        role: BranchRole::RbcBypass,
                        treatment_path: false,
                        route: ChannelRouteSpec {
                            length_m: Length::from_base(8.0e-3),
                            width_m: Length::from_base(2.0e-3),
                            height_m: Length::from_base(1.0e-3),
                            serpentine: None,
                            therapy_zone: TherapyZone::HealthyBypass,
                        },
                        recovery_sub_split: None,
                    },
                ],
            }],
            venturi_placements: Vec::new(),
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
        };

        let updated = with_branch_serpentine(spec, "stage_0", "center", 5, 1.8e-3, 4.0e-3)
            .expect("treatment branch must exist");
        assert_eq!(
            updated.split_stages[0].branches[1].route.serpentine,
            Some(SerpentineSpec {
                wave_type: crate::topology::SerpentineWaveType::Sine,
                segments: 5,
                bend_radius_m: Length::from_base(1.8e-3),
                segment_length_m: Length::from_base(4.0e-3),
            })
        );
    }
}
