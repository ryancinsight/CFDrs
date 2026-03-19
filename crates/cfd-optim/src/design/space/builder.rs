//! Primitive selective blueprint-candidate construction.
use crate::constraints::{SERPENTINE_BEND_RADIUS_M, TREATMENT_WIDTH_MM};
use crate::domain::{BlueprintCandidate, OperatingPoint};
use cfd_schematics::topology::{
    presets::{build_milestone12_blueprint, Milestone12TopologyRequest},
    SerpentineSpec, TreatmentActuationMode, VenturiPlacementMode,
};

#[must_use]
pub(crate) fn primitive_selective_candidate(
    id: String,
    request: &Milestone12TopologyRequest,
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    throat_diameter_m: f64,
    throat_length_m: f64,
    channel_width_m: f64,
    serpentine_segments: usize,
    cif_pretri_center_frac: f64,
    cif_terminal_tri_center_frac: f64,
    cif_terminal_bi_treat_frac: f64,
    treatment_mode: TreatmentActuationMode,
    centerline_venturi_throat_count: u8,
) -> BlueprintCandidate {
    let request = Milestone12TopologyRequest {
        topology_id: request.topology_id.clone(),
        design_name: request.design_name.clone(),
        mirror_x: request.mirror_x,
        mirror_y: request.mirror_y,
        box_dims_mm: request.box_dims_mm,
        split_kinds: request.split_kinds.clone(),
        inlet_width_m: channel_width_m,
        channel_height_m: request.channel_height_m,
        branch_length_m: TREATMENT_WIDTH_MM * 1e-3,
        outlet_tail_length_m: TREATMENT_WIDTH_MM * 1e-3,
        stage_layouts: request.stage_layouts.clone(),
        first_trifurcation_center_frac: cif_pretri_center_frac,
        later_trifurcation_center_frac: cif_terminal_tri_center_frac,
        bifurcation_treatment_frac: cif_terminal_bi_treat_frac,
        treatment_mode,
        venturi_throat_count: centerline_venturi_throat_count,
        venturi_throat_width_m: throat_diameter_m,
        venturi_throat_length_m: throat_length_m,
        center_serpentine: (serpentine_segments >= 2).then_some(SerpentineSpec {
            segments: serpentine_segments,
            bend_radius_m: SERPENTINE_BEND_RADIUS_M,
            segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
            wave_type: cfd_schematics::SerpentineWaveType::Sine,
        }),
        venturi_placement_mode: VenturiPlacementMode::StraightSegment,
        venturi_target_channel_ids: Vec::new(),
    };
    let blueprint = build_milestone12_blueprint(&request)
        .expect("primitive selective blueprint candidate should be constructible");

    BlueprintCandidate::new(
        id,
        blueprint,
        OperatingPoint {
            flow_rate_m3_s,
            inlet_gauge_pa,
            feed_hematocrit: 0.45,
            patient_context: None,
        },
    )
}
