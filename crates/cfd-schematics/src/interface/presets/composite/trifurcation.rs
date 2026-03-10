//! Trifurcation-based composite presets: symmetric 3-branch topologies.
#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use super::{shah_london, BLOOD_MU};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use crate::geometry::metadata::VenturiGeometryMetadata;

/// Symmetric trifurcation with a venturi throat in each of three branches.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [3 × venturi branch] → merge_jn → trunk_out → outlet`
///
/// Each branch carries `Q/3`; venturi throat in each provides cavitation.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"throat_section"`, `"throat_section_2"`, `"throat_section_3"` — three throats
/// - `"trunk_out"` — downstream trunk
#[must_use]
pub fn trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for b in 1..=3 {
        bp.add_node(NodeSpec::new(
            format!("contraction_{b}"),
            NodeKind::Junction,
        ));
        bp.add_node(NodeSpec::new(format!("throat_{b}"), NodeKind::Junction));
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_jn",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    for b in 1..=3 {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("branch_{b}_in"),
                "split_jn",
                format!("contraction_{b}"),
                l_taper,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        let throat_name = if b == 1 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{b}")
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                throat_name,
                format!("contraction_{b}"),
                format!("throat_{b}"),
                l_throat,
                throat_width_m,
                height_m,
                shah_london(throat_width_m, height_m, l_throat, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget))
            .with_venturi_geometry(VenturiGeometryMetadata {
                throat_width_m,
                throat_height_m: height_m,
                throat_length_m: l_throat,
                inlet_width_m: main_width_m,
                outlet_width_m: main_width_m,
                convergent_half_angle_deg: ((main_width_m - throat_width_m).abs() * 0.5)
                    .atan2(l_taper)
                    .to_degrees(),
                divergent_half_angle_deg: ((main_width_m - throat_width_m).abs() * 0.5)
                    .atan2(l_taper)
                    .to_degrees(),
                throat_position: 0.5,
            })
            .with_metadata(VenturiGeometryMetadata {
                throat_width_m,
                throat_height_m: height_m,
                throat_length_m: l_throat,
                inlet_width_m: main_width_m,
                outlet_width_m: main_width_m,
                convergent_half_angle_deg: ((main_width_m - throat_width_m).abs() * 0.5)
                    .atan2(l_taper)
                    .to_degrees(),
                divergent_half_angle_deg: ((main_width_m - throat_width_m).abs() * 0.5)
                    .atan2(l_taper)
                    .to_degrees(),
                throat_position: 0.5,
            }),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("branch_{b}_out"),
                format!("throat_{b}"),
                "merge_jn",
                l_taper,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_jn",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}

/// Symmetric trifurcation with a full serpentine in each of three arms.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [3 × serpentine arm] → merge_jn → trunk_out → outlet`
///
/// Three parallel serpentine arms provide maximum coverage with `Q/3` in each.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"segment_1"` … `"segment_n"` — arm-1 serpentine
/// - `"arm2_seg_1"` … `"arm2_seg_n"` — arm-2
/// - `"arm3_seg_1"` … `"arm3_seg_n"` — arm-3
/// - `"trunk_out"` — downstream trunk
#[must_use]
pub fn trifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for arm in 1..=3 {
        for i in 1..segments {
            bp.add_node(NodeSpec::new(format!("ar{arm}_{i}"), NodeKind::Junction));
        }
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_jn",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    for arm in 1..=3u32 {
        for i in 0..segments {
            let actual_from = if i == 0 {
                "split_jn".to_string()
            } else {
                format!("ar{arm}_{i}")
            };
            let to = if i + 1 == segments {
                "merge_jn".to_string()
            } else {
                format!("ar{arm}_{}", i + 1)
            };
            let seg_name = match arm {
                1 => format!("segment_{}", i + 1),
                2 => format!("arm2_seg_{}", i + 1),
                _ => format!("arm3_seg_{}", i + 1),
            };
            bp.add_channel({
                let mut spec = ChannelSpec::new_pipe_rect(
                    seg_name,
                    actual_from,
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
                spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow))
            });
        }
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_jn",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}
