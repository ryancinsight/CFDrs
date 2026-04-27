//! Internal parallel-lane geometry helpers.
use super::super::super::finalize_preset_blueprint;
use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::generator::CenterSerpentinePathSpec;
use crate::topology::presets::{parallel_path_spec, with_venturi};
use crate::topology::{
    ChannelRouteSpec, ParallelChannelSpec, SerpentineSpec, ThroatGeometrySpec,
    TreatmentActuationMode, VenturiConfig, VenturiPlacementMode,
};
use crate::BlueprintTopologyFactory;

/// Optional serpentine geometry applied only to center treatment lanes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenterSerpentineSpec {
    pub segments: usize,
    pub bend_radius_m: f64,
    pub segment_length_m: f64,
}

pub(super) fn normalized_center_serpentine(
    center_serpentine: Option<CenterSerpentineSpec>,
) -> Option<CenterSerpentineSpec> {
    center_serpentine.and_then(|spec| {
        if spec.segments < 2 || spec.bend_radius_m <= 0.0 || spec.segment_length_m <= 0.0 {
            None
        } else {
            Some(spec)
        }
    })
}

pub(super) fn generator_center_serpentine(
    center_serpentine: Option<CenterSerpentineSpec>,
) -> Option<CenterSerpentinePathSpec> {
    normalized_center_serpentine(center_serpentine).map(|spec| CenterSerpentinePathSpec {
        segments: spec.segments,
        bend_radius_m: spec.bend_radius_m,
        wave_type: crate::SerpentineWaveType::default(),
    })
}

pub(super) fn parallel_lane(
    channel_id: impl Into<String>,
    length_m: f64,
    width_m: f64,
    height_m: f64,
    therapy_zone: TherapyZone,
    serpentine: Option<SerpentineSpec>,
) -> ParallelChannelSpec {
    ParallelChannelSpec {
        channel_id: channel_id.into(),
        route: ChannelRouteSpec {
            length_m,
            width_m,
            height_m,
            serpentine,
            therapy_zone,
        },
    }
}

pub(super) fn canonical_parallel_blueprint(
    name: String,
    inlet_width_m: f64,
    trunk_length_m: f64,
    outlet_tail_length_m: f64,
    parallel_channels: Vec<ParallelChannelSpec>,
    treatment_mode: TreatmentActuationMode,
) -> NetworkBlueprint {
    let topology = parallel_path_spec(
        &name,
        inlet_width_m,
        inlet_width_m,
        trunk_length_m,
        outlet_tail_length_m,
        parallel_channels,
        treatment_mode,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid canonical parallel topology"),
    )
}

pub(super) fn canonical_parallel_venturi_blueprint(
    name: String,
    inlet_width_m: f64,
    trunk_length_m: f64,
    outlet_tail_length_m: f64,
    throat_width_m: f64,
    throat_height_m: f64,
    throat_length_m: f64,
    target_channel_ids: Vec<String>,
    parallel_channels: Vec<ParallelChannelSpec>,
) -> NetworkBlueprint {
    let topology = with_venturi(
        parallel_path_spec(
            &name,
            inlet_width_m,
            inlet_width_m,
            trunk_length_m,
            outlet_tail_length_m,
            parallel_channels,
            TreatmentActuationMode::VenturiCavitation,
        ),
        VenturiConfig {
            target_channel_ids,
            serial_throat_count: 1,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m,
                throat_height_m,
                throat_length_m,
                inlet_width_m: 0.0,
                outlet_width_m: 0.0,
                convergent_half_angle_deg: 7.0,
                divergent_half_angle_deg: 7.0,
            },
            placement_mode: VenturiPlacementMode::StraightSegment,
        },
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid canonical venturi topology"),
    )
}
