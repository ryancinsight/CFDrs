//! Primitive selective blueprint-candidate construction.
use crate::constraints::{
    CHANNEL_HEIGHT_M, SERPENTINE_BEND_RADIUS_M, TREATMENT_WIDTH_MM, VENTURI_INLET_DIAM_M,
};
use crate::design::sequence_catalog::PrimitiveSplitSequence;
use crate::design::space::params::TreatmentZoneMode;
use crate::domain::{BlueprintCandidate, OperatingPoint};
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::topology::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec, SerpentineSpec, SplitKind,
    SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode, VenturiPlacementMode,
    VenturiPlacementSpec,
};

#[must_use]
pub(crate) fn primitive_selective_candidate(
    id: String,
    sequence: PrimitiveSplitSequence,
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    throat_diameter_m: f64,
    throat_length_m: f64,
    channel_width_m: f64,
    serpentine_segments: usize,
    cif_pretri_center_frac: f64,
    cif_terminal_tri_center_frac: f64,
    cif_terminal_bi_treat_frac: f64,
    treatment_zone_mode: TreatmentZoneMode,
    centerline_venturi_throat_count: u8,
) -> BlueprintCandidate {
    let split_kinds = primitive_sequence_to_split_kinds(sequence);
    let mut spec = BlueprintTopologySpec {
        topology_id: format!(
            "pst-{}",
            sequence.label().replace('→', "").to_ascii_lowercase()
        ),
        design_name: sequence.label().to_string(),
        box_dims_mm: (127.76, 85.47),
        inlet_width_m: channel_width_m,
        outlet_width_m: channel_width_m,
        trunk_length_m: TREATMENT_WIDTH_MM * 1e-3,
        outlet_tail_length_m: TREATMENT_WIDTH_MM * 1e-3,
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        split_stages: Vec::new(),
        venturi_placements: Vec::new(),
        treatment_mode: match treatment_zone_mode {
            TreatmentZoneMode::UltrasoundOnly => TreatmentActuationMode::UltrasoundOnly,
            TreatmentZoneMode::VenturiThroats => TreatmentActuationMode::VenturiCavitation,
        },
    };

    let mut parent_width = channel_width_m;
    for (index, split_kind) in split_kinds.iter().copied().enumerate() {
        let is_last = index + 1 == split_kinds.len();
        let branches = match split_kind {
            SplitKind::Trifurcation => {
                let center_frac = if index == 0 {
                    cif_pretri_center_frac
                } else {
                    cif_terminal_tri_center_frac
                };
                let center_width = parent_width * center_frac;
                let bypass_width = parent_width * (1.0 - center_frac) / 2.0;
                vec![
                    treatment_branch("center", center_width, serpentine_segments, is_last),
                    bypass_branch("left", bypass_width, BranchRole::WbcCollection),
                    bypass_branch("right", bypass_width, BranchRole::RbcBypass),
                ]
            }
            SplitKind::Bifurcation => {
                let treatment_width = parent_width * cif_terminal_bi_treat_frac;
                let bypass_width = parent_width * (1.0 - cif_terminal_bi_treat_frac);
                vec![
                    treatment_branch("upper", treatment_width, serpentine_segments, is_last),
                    bypass_branch("lower", bypass_width, BranchRole::RbcBypass),
                ]
            }
        };

        parent_width = branches
            .iter()
            .find(|branch| branch.treatment_path)
            .map_or(parent_width, |branch| branch.route.width_m);
        spec.split_stages.push(SplitStageSpec {
            stage_id: format!("stage_{index}"),
            split_kind,
            branches,
        });
    }

    if matches!(treatment_zone_mode, TreatmentZoneMode::VenturiThroats) {
        spec.venturi_placements = spec
            .treatment_channel_ids()
            .into_iter()
            .enumerate()
            .map(|(index, channel_id)| VenturiPlacementSpec {
                placement_id: format!("v{index}"),
                target_channel_id: channel_id,
                serial_throat_count: centerline_venturi_throat_count,
                throat_geometry: ThroatGeometrySpec {
                    throat_width_m: throat_diameter_m,
                    throat_height_m: CHANNEL_HEIGHT_M,
                    throat_length_m,
                    inlet_width_m: VENTURI_INLET_DIAM_M,
                    outlet_width_m: VENTURI_INLET_DIAM_M,
                    convergent_half_angle_deg: 7.0,
                    divergent_half_angle_deg: 7.0,
                },
                placement_mode: VenturiPlacementMode::StraightSegment,
            })
            .collect();
    }

    BlueprintCandidate::from_topology_spec(
        id,
        &spec,
        OperatingPoint {
            flow_rate_m3_s,
            inlet_gauge_pa,
            feed_hematocrit: 0.45,
            patient_context: None,
        },
    )
    .expect("primitive selective blueprint candidate should be constructible")
}

fn primitive_sequence_to_split_kinds(sequence: PrimitiveSplitSequence) -> Vec<SplitKind> {
    sequence
        .label()
        .split('→')
        .map(|segment| match segment {
            "Tri" => SplitKind::Trifurcation,
            _ => SplitKind::Bifurcation,
        })
        .collect()
}

fn treatment_branch(
    label: &str,
    width_m: f64,
    serpentine_segments: usize,
    is_last: bool,
) -> BranchSpec {
    BranchSpec {
        label: label.to_string(),
        role: BranchRole::Treatment,
        treatment_path: true,
        route: ChannelRouteSpec {
            length_m: TREATMENT_WIDTH_MM * 1e-3,
            width_m,
            height_m: CHANNEL_HEIGHT_M,
            serpentine: (is_last && serpentine_segments >= 2).then(|| SerpentineSpec {
                segments: serpentine_segments,
                bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
            }),
            therapy_zone: TherapyZone::CancerTarget,
        },
    }
}

fn bypass_branch(label: &str, width_m: f64, role: BranchRole) -> BranchSpec {
    BranchSpec {
        label: label.to_string(),
        role,
        treatment_path: false,
        route: ChannelRouteSpec {
            length_m: TREATMENT_WIDTH_MM * 1e-3,
            width_m,
            height_m: CHANNEL_HEIGHT_M,
            serpentine: None,
            therapy_zone: TherapyZone::HealthyBypass,
        },
    }
}
