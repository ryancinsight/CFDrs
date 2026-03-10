use super::finalize_preset_blueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::topology::presets::{parallel_path_spec, with_venturi_placements};
use crate::topology::{
    ChannelRouteSpec, ParallelChannelSpec, SerpentineSpec, TreatmentActuationMode,
    VenturiPlacementMode,
};
use crate::BlueprintTopologyFactory;

fn parallel_channels_for_n_furcation(
    splits: usize,
    length_m: f64,
    width_m: f64,
    height_m: f64,
    serpentine: Option<SerpentineSpec>,
) -> Vec<ParallelChannelSpec> {
    (0..splits.max(1))
        .map(|index| ParallelChannelSpec {
            channel_id: if index == 0 {
                "throat_section".to_string()
            } else {
                format!("arm_{}", index + 1)
            },
            route: ChannelRouteSpec {
                length_m,
                width_m,
                height_m,
                serpentine: serpentine.clone(),
                therapy_zone: TherapyZone::CancerTarget,
            },
        })
        .collect()
}

#[must_use]
pub fn n_furcation_venturi_rect(
    name: impl Into<String>,
    splits: usize,
    trunk_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let base = parallel_path_spec(
        &name,
        main_width_m,
        main_width_m,
        trunk_length_m,
        trunk_length_m,
        parallel_channels_for_n_furcation(splits, trunk_length_m, main_width_m, height_m, None),
        TreatmentActuationMode::VenturiCavitation,
    );
    let spec = with_venturi_placements(
        base,
        throat_width_m,
        height_m,
        throat_length_m,
        1,
        VenturiPlacementMode::StraightSegment,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&spec).expect("valid n_furcation venturi topology spec"),
    )
}

#[must_use]
pub fn n_furcation_serpentine_rect(
    name: impl Into<String>,
    splits: usize,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let spec = parallel_path_spec(
        &name,
        main_width_m,
        main_width_m,
        trunk_length_m,
        trunk_length_m,
        parallel_channels_for_n_furcation(
            splits,
            segment_length_m,
            main_width_m,
            height_m,
            Some(SerpentineSpec {
                segments: segments.max(2),
                bend_radius_m: main_width_m * 0.5,
                segment_length_m,
            }),
        ),
        TreatmentActuationMode::UltrasoundOnly,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&spec).expect("valid n_furcation serpentine topology spec"),
    )
}
