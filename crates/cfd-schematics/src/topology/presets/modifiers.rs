//! Builder modifiers for topology specs (venturi, serpentine, Dean placement).

use super::super::model::{
    BlueprintTopologySpec, SerpentineSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
};
use super::helpers::VENTURI_HALF_ANGLE_DEG;

/// Builder: attach venturi throat placements to a spec (Option 1 → Option 2).
///
/// Targets all treatment-path channels with identical throat geometry.
/// Returns a new spec with `treatment_mode = VenturiCavitation`.
#[must_use]
pub fn with_venturi_placements(
    mut spec: BlueprintTopologySpec,
    throat_width_m: f64,
    throat_height_m: f64,
    throat_length_m: f64,
    serial_throat_count: u8,
    placement_mode: VenturiPlacementMode,
) -> BlueprintTopologySpec {
    let treatment_ids = spec.treatment_channel_ids();
    spec.venturi_placements = treatment_ids
        .into_iter()
        .enumerate()
        .map(|(idx, channel_id)| {
            let inlet_width = spec
                .split_stages
                .iter()
                .flat_map(|s| s.branches.iter())
                .find(|b| {
                    spec.split_stages.iter().any(|stage| {
                        BlueprintTopologySpec::branch_channel_id(&stage.stage_id, &b.label)
                            == channel_id
                    })
                })
                .map_or(throat_width_m * 3.0, |b| b.route.width_m);

            VenturiPlacementSpec {
                placement_id: format!("venturi_{idx}"),
                target_channel_id: channel_id,
                serial_throat_count,
                throat_geometry: ThroatGeometrySpec {
                    throat_width_m,
                    throat_height_m,
                    throat_length_m,
                    inlet_width_m: inlet_width,
                    outlet_width_m: inlet_width,
                    convergent_half_angle_deg: VENTURI_HALF_ANGLE_DEG,
                    divergent_half_angle_deg: VENTURI_HALF_ANGLE_DEG,
                },
                placement_mode,
            }
        })
        .collect();
    spec.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    spec
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
        bend_radius_m,
        segment_length_m,
    });
    Some(spec)
}

/// Builder: relocate venturi placements to Dean vortex peak sites.
///
/// Changes all venturi placement modes to
/// [`VenturiPlacementMode::CurvaturePeakDeanNumber`], placing throats at
/// the apex of serpentine bends where secondary Dean circulation maximises
/// centrifugal focusing of stiff CTCs.
#[must_use]
pub fn with_dean_venturi_placement(mut spec: BlueprintTopologySpec) -> BlueprintTopologySpec {
    for placement in &mut spec.venturi_placements {
        placement.placement_mode = VenturiPlacementMode::CurvaturePeakDeanNumber;
    }
    spec
}
