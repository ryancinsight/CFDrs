//! Series-topology composite presets: single flow path with multiple features.
#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use super::{shah_london, BLOOD_MU};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use crate::topology::presets::{serial_double_venturi_series_spec, venturi_serpentine_series_spec};
use crate::BlueprintTopologyFactory;

/// Rectangular venturi followed immediately by a serpentine — all in series.
///
/// **Topology**:
/// `inlet → inlet_section → [venturi] → segment_1 → … → segment_n → outlet`
///
/// Provides a single cavitation site (venturi throat) upstream of a full
/// serpentine sweep across the treatment zone.
///
/// # Channel names
/// - `"inlet_section"` — venturi inlet approach (`main_width_m × height_m`)
/// - `"throat_section"` — venturi throat (`throat_width_m × height_m`, [`TherapyZone::CancerTarget`])
/// - `"diffuser_section"` — venturi recovery (`main_width_m × height_m`)
/// - `"segment_1"` … `"segment_n"` — serpentine straights
#[must_use]
pub fn venturi_serpentine_rect(
    name: impl Into<String>,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    segments: usize,
    segment_length_m: f64,
) -> NetworkBlueprint {
    let name = name.into();
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let topology = venturi_serpentine_series_spec(
        &name,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        segments,
        segment_length_m,
    );
    let lineage = BlueprintTopologyFactory::lineage_for_spec(&topology);
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_jn", NodeKind::Junction));
    for i in 0..segments {
        bp.add_node(NodeSpec::new(format!("serp_{i}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "contraction",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "throat_section",
            "contraction",
            "throat_jn",
            l_throat,
            throat_width_m,
            height_m,
            shah_london(throat_width_m, height_m, l_throat, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );
    let diffuser_to = if segments == 0 {
        "outlet".to_string()
    } else {
        "serp_0".to_string()
    };
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "diffuser_section",
            "throat_jn",
            diffuser_to,
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    for i in 0..segments {
        let from = format!("serp_{i}");
        let to = if i + 1 == segments {
            "outlet".to_string()
        } else {
            format!("serp_{}", i + 1)
        };
        let mut spec = ChannelSpec::new_pipe_rect(
            format!("segment_{}", i + 1),
            from,
            to,
            segment_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, segment_length_m, BLOOD_MU),
            0.0,
        );
        spec.channel_shape = ChannelShape::Serpentine {
            segments,
            bend_radius_m: main_width_m * 0.5,
        };
        bp.add_channel(spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)));
    }

    bp.with_topology_spec(topology).with_lineage(lineage)
}

/// Two venturi throats in series on the same flow path — closed loop.
///
/// **Topology**:
/// `inlet → [venturi-1] → mid_jn → [venturi-2] → outlet`
///
/// Provides two sequential cavitation events on the same fluid parcel.
///
/// # Channel names
/// - `"inlet_section"` — first venturi approach
/// - `"throat_section"` — first venturi throat ([`TherapyZone::CancerTarget`])
/// - `"mid_section"` — inter-venturi recovery / approach segment
/// - `"throat_section_2"` — second venturi throat ([`TherapyZone::CancerTarget`])
/// - `"diffuser_section"` — final recovery
#[must_use]
pub fn serial_double_venturi_rect(
    name: impl Into<String>,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    inter_length_m: f64,
) -> NetworkBlueprint {
    let name = name.into();
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let topology = serial_double_venturi_series_spec(
        &name,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        inter_length_m,
    );
    let lineage = BlueprintTopologyFactory::lineage_for_spec(&topology);
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("mid_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("contraction_2", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_2", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "contraction_1",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "throat_section",
            "contraction_1",
            "throat_1",
            l_throat,
            throat_width_m,
            height_m,
            shah_london(throat_width_m, height_m, l_throat, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "diffuser_1",
            "throat_1",
            "mid_jn",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "mid_section",
            "mid_jn",
            "contraction_2",
            inter_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, inter_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "throat_section_2",
            "contraction_2",
            "throat_2",
            l_throat,
            throat_width_m,
            height_m,
            shah_london(throat_width_m, height_m, l_throat, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "diffuser_section",
            "throat_2",
            "outlet",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp.with_topology_spec(topology).with_lineage(lineage)
}

#[cfg(test)]
mod tests {
    use super::{serial_double_venturi_rect, venturi_serpentine_rect};

    #[test]
    fn venturi_serpentine_rect_embeds_canonical_topology_metadata() {
        let blueprint =
            venturi_serpentine_rect("vs-meta", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3, 4, 8.0e-3);
        let topology = blueprint
            .topology_spec()
            .expect("venturi_serpentine_rect should attach topology metadata");

        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 1);
        assert!(topology.has_serpentine());
        assert!(topology
            .treatment_channel_ids()
            .iter()
            .any(|channel_id| channel_id == "throat_section"));
    }

    #[test]
    fn serial_double_venturi_rect_embeds_serial_path_metadata() {
        let blueprint =
            serial_double_venturi_rect("sdv-meta", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3, 4.0e-3);
        let topology = blueprint
            .topology_spec()
            .expect("serial_double_venturi_rect should attach topology metadata");

        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 2);
        assert_eq!(topology.treatment_channel_ids().len(), 2);
    }
}
