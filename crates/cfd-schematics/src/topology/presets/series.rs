//! Series-path topology presets.

use crate::domain::therapy_metadata::TherapyZone;

use super::super::model::{
    BlueprintTopologySpec, SeriesChannelSpec, SerpentineSpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
};
use super::helpers::{series_channel, throat_geometry, PLATE_HEIGHT_MM, PLATE_WIDTH_MM};

/// Create a canonical linear series-path topology spec.
#[must_use]
pub fn series_path_spec(
    name: &str,
    inlet_width_m: f64,
    outlet_width_m: f64,
    trunk_length_m: f64,
    outlet_tail_length_m: f64,
    series_channels: Vec<SeriesChannelSpec>,
    treatment_mode: TreatmentActuationMode,
    venturi_placements: Vec<VenturiPlacementSpec>,
) -> BlueprintTopologySpec {
    BlueprintTopologySpec {
        topology_id: format!("{name}_topology"),
        design_name: name.to_string(),
        box_dims_mm: (PLATE_WIDTH_MM, PLATE_HEIGHT_MM),
        inlet_width_m,
        outlet_width_m,
        trunk_length_m,
        outlet_tail_length_m,
        series_channels,
        parallel_channels: Vec::new(),
        split_stages: Vec::new(),
        venturi_placements,
        treatment_mode,
    }
}

/// Canonical single-venturi series path.
#[must_use]
pub fn single_venturi_series_spec(
    name: &str,
    inlet_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> BlueprintTopologySpec {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (inlet_width_m + throat_width_m) * 0.5;
    series_path_spec(
        name,
        inlet_width_m,
        inlet_width_m,
        l_taper,
        l_taper,
        vec![
            series_channel(
                "inlet_section",
                l_taper,
                inlet_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
            series_channel(
                "throat_section",
                l_throat,
                throat_width_m,
                height_m,
                TherapyZone::CancerTarget,
                None,
            ),
            series_channel(
                "diffuser_section",
                l_taper,
                inlet_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
        ],
        TreatmentActuationMode::VenturiCavitation,
        vec![VenturiPlacementSpec {
            placement_id: "venturi_0".to_string(),
            target_channel_id: "throat_section".to_string(),
            serial_throat_count: 1,
            throat_geometry: throat_geometry(inlet_width_m, throat_width_m, height_m, l_throat),
            placement_mode: VenturiPlacementMode::StraightSegment,
        }],
    )
}

/// Canonical two-throat serial venturi path.
#[must_use]
pub fn serial_double_venturi_series_spec(
    name: &str,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    inter_length_m: f64,
) -> BlueprintTopologySpec {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    series_path_spec(
        name,
        main_width_m,
        main_width_m,
        l_taper,
        l_taper,
        vec![
            series_channel(
                "inlet_section",
                l_taper,
                main_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
            series_channel(
                "throat_section",
                l_throat,
                throat_width_m,
                height_m,
                TherapyZone::CancerTarget,
                None,
            ),
            series_channel(
                "diffuser_1",
                l_taper,
                main_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
            series_channel(
                "mid_section",
                inter_length_m,
                main_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
            series_channel(
                "throat_section_2",
                l_throat,
                throat_width_m,
                height_m,
                TherapyZone::CancerTarget,
                None,
            ),
            series_channel(
                "diffuser_section",
                l_taper,
                main_width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ),
        ],
        TreatmentActuationMode::VenturiCavitation,
        vec![
            VenturiPlacementSpec {
                placement_id: "venturi_0".to_string(),
                target_channel_id: "throat_section".to_string(),
                serial_throat_count: 1,
                throat_geometry: throat_geometry(main_width_m, throat_width_m, height_m, l_throat),
                placement_mode: VenturiPlacementMode::StraightSegment,
            },
            VenturiPlacementSpec {
                placement_id: "venturi_1".to_string(),
                target_channel_id: "throat_section_2".to_string(),
                serial_throat_count: 1,
                throat_geometry: throat_geometry(main_width_m, throat_width_m, height_m, l_throat),
                placement_mode: VenturiPlacementMode::StraightSegment,
            },
        ],
    )
}

/// Canonical serpentine treatment path.
#[must_use]
pub fn serpentine_series_spec(
    name: &str,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    height_m: f64,
) -> BlueprintTopologySpec {
    let segment_count = segments.max(1);
    let serpentine = (segment_count > 1).then_some(SerpentineSpec {
        segments: segment_count,
        bend_radius_m: width_m * 0.5,
        segment_length_m,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    });
    let series_channels = (0..segment_count)
        .map(|index| {
            series_channel(
                format!("segment_{}", index + 1),
                segment_length_m,
                width_m,
                height_m,
                TherapyZone::CancerTarget,
                serpentine.clone(),
            )
        })
        .collect();
    series_path_spec(
        name,
        width_m,
        width_m,
        segment_length_m,
        segment_length_m,
        series_channels,
        TreatmentActuationMode::UltrasoundOnly,
        Vec::new(),
    )
}

/// Canonical venturi followed by treatment-zone serpentine path.
#[must_use]
pub fn venturi_serpentine_series_spec(
    name: &str,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    segments: usize,
    segment_length_m: f64,
) -> BlueprintTopologySpec {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let segment_count = segments.max(1);
    let serpentine = (segment_count > 1).then_some(SerpentineSpec {
        segments: segment_count,
        bend_radius_m: main_width_m * 0.5,
        segment_length_m,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    });
    let mut series_channels = vec![
        series_channel(
            "inlet_section",
            l_taper,
            main_width_m,
            height_m,
            TherapyZone::MixedFlow,
            None,
        ),
        series_channel(
            "throat_section",
            l_throat,
            throat_width_m,
            height_m,
            TherapyZone::CancerTarget,
            None,
        ),
        series_channel(
            "diffuser_section",
            l_taper,
            main_width_m,
            height_m,
            TherapyZone::MixedFlow,
            None,
        ),
    ];
    series_channels.extend((0..segment_count).map(|index| {
        series_channel(
            format!("segment_{}", index + 1),
            segment_length_m,
            main_width_m,
            height_m,
            TherapyZone::CancerTarget,
            serpentine.clone(),
        )
    }));
    series_path_spec(
        name,
        main_width_m,
        main_width_m,
        l_taper,
        segment_length_m,
        series_channels,
        TreatmentActuationMode::VenturiCavitation,
        vec![VenturiPlacementSpec {
            placement_id: "venturi_0".to_string(),
            target_channel_id: "throat_section".to_string(),
            serial_throat_count: 1,
            throat_geometry: throat_geometry(main_width_m, throat_width_m, height_m, l_throat),
            placement_mode: VenturiPlacementMode::StraightSegment,
        }],
    )
}

/// Canonical bend-apex venturi placement along a serpentine path.
#[must_use]
pub fn serpentine_bend_venturi_series_spec(
    name: &str,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    bend_radius_m: f64,
) -> BlueprintTopologySpec {
    let segment_count = segments.max(2);
    let n_venturis = segment_count - 1;
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (width_m + throat_width_m) * 0.5;
    let serpentine = Some(SerpentineSpec {
        segments: segment_count,
        bend_radius_m,
        segment_length_m,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    });
    let mut series_channels = Vec::with_capacity(segment_count + n_venturis * 3);
    let mut venturi_placements = Vec::with_capacity(n_venturis);
    for index in 0..segment_count {
        series_channels.push(series_channel(
            format!("segment_{}", index + 1),
            segment_length_m,
            width_m,
            height_m,
            TherapyZone::MixedFlow,
            serpentine.clone(),
        ));
        if index < n_venturis {
            let venturi_index = index + 1;
            series_channels.push(series_channel(
                format!("approach_{venturi_index}"),
                l_taper,
                width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ));
            let throat_id = format!("throat_{venturi_index}");
            series_channels.push(series_channel(
                throat_id.clone(),
                l_throat,
                throat_width_m,
                height_m,
                TherapyZone::CancerTarget,
                None,
            ));
            series_channels.push(series_channel(
                format!("recovery_{venturi_index}"),
                l_taper,
                width_m,
                height_m,
                TherapyZone::MixedFlow,
                None,
            ));
            venturi_placements.push(VenturiPlacementSpec {
                placement_id: format!("venturi_{index}"),
                target_channel_id: throat_id,
                serial_throat_count: 1,
                throat_geometry: throat_geometry(width_m, throat_width_m, height_m, l_throat),
                placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            });
        }
    }

    series_path_spec(
        name,
        width_m,
        width_m,
        segment_length_m,
        segment_length_m,
        series_channels,
        TreatmentActuationMode::VenturiCavitation,
        venturi_placements,
    )
}

/// Canonical constriction-expansion leukapheresis array.
#[must_use]
pub fn constriction_expansion_series_spec(
    name: &str,
    n_cycles: usize,
    wide_length_m: f64,
    narrow_length_m: f64,
    wide_width_m: f64,
    narrow_width_m: f64,
    height_m: f64,
) -> BlueprintTopologySpec {
    let cycle_count = n_cycles.max(1);
    let mut series_channels = Vec::with_capacity(cycle_count * 2);
    for index in 0..cycle_count {
        series_channels.push(series_channel(
            format!("wide_{index}"),
            wide_length_m,
            wide_width_m,
            height_m,
            TherapyZone::MixedFlow,
            None,
        ));
        series_channels.push(series_channel(
            format!("narrow_{index}"),
            narrow_length_m,
            narrow_width_m,
            height_m,
            TherapyZone::CancerTarget,
            None,
        ));
    }
    series_path_spec(
        name,
        wide_width_m,
        wide_width_m,
        wide_length_m,
        wide_length_m,
        series_channels,
        TreatmentActuationMode::UltrasoundOnly,
        Vec::new(),
    )
}

/// Canonical spiral / Dean-dominant treatment path.
#[must_use]
pub fn spiral_serpentine_series_spec(
    name: &str,
    n_turns: usize,
    turn_length_m: f64,
    width_m: f64,
    height_m: f64,
) -> BlueprintTopologySpec {
    let turn_count = n_turns.max(1);
    let bend_radius_m = (turn_length_m / (2.0 * std::f64::consts::PI)).max(width_m * 0.5);
    let series_channels = (0..turn_count)
        .map(|index| {
            series_channel(
                format!("spiral_{index}"),
                turn_length_m,
                width_m,
                height_m,
                TherapyZone::CancerTarget,
                Some(SerpentineSpec {
                    segments: 4,
                    bend_radius_m,
                    segment_length_m: turn_length_m * 0.25,
                    wave_type: crate::topology::SerpentineWaveType::Sine,
                }),
            )
        })
        .collect();
    series_path_spec(
        name,
        width_m,
        width_m,
        turn_length_m,
        turn_length_m,
        series_channels,
        TreatmentActuationMode::UltrasoundOnly,
        Vec::new(),
    )
}
