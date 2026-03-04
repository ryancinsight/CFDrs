//! Bifurcation-based composite presets: symmetric 2-branch split topologies.

use super::{shah_london, BLOOD_MU};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

/// Symmetric bifurcation with a venturi throat in each branch — closed loop.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [branch_1 venturi, branch_2 venturi] → merge_jn → trunk_out → outlet`
///
/// Both branches carry `Q/2`; venturi throat in each provides cavitation.
/// All paths reconverge at a single outlet.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk (`main_width_m × height_m`)
/// - `"throat_section"` — branch-1 venturi throat (primary, [`TherapyZone::CancerTarget`])
/// - `"throat_section_2"` — branch-2 venturi throat
/// - `"trunk_out"` — downstream trunk
#[must_use]
pub fn bifurcation_venturi_rect(
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
    bp.add_node(NodeSpec::new("contraction_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("contraction_2", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_2", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // Upstream trunk
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

    // Branch 1 venturi
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "branch_1_in",
            "split_jn",
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
            "branch_1_out",
            "throat_1",
            "merge_jn",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Branch 2 venturi
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "branch_2_in",
            "split_jn",
            "contraction_2",
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
            "branch_2_out",
            "throat_2",
            "merge_jn",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Downstream trunk
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

/// Symmetric bifurcation with a full serpentine in each arm — closed loop.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [arm_1 serpentine, arm_2 serpentine] → merge_jn → trunk_out → outlet`
///
/// Two parallel serpentine arms provide uniform exposure with full coverage.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"segment_1"` … `"segment_n"` — arm-1 serpentine segments
/// - `"arm2_seg_1"` … `"arm2_seg_n"` — arm-2 serpentine segments
#[must_use]
pub fn bifurcation_serpentine_rect(
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
    // Intermediate junction nodes between segments (arm*_0 is split_jn itself; skip it)
    for i in 1..segments {
        bp.add_node(NodeSpec::new(format!("arm1_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("arm2_{i}"), NodeKind::Junction));
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

    // Arm 1 serpentine
    for i in 0..segments {
        let from = format!("arm1_{i}");
        let to = if i + 1 == segments {
            "merge_jn".to_string()
        } else {
            format!("arm1_{}", i + 1)
        };
        // Connect split_jn to arm1_0 via the first channel
        let actual_from = if i == 0 { "split_jn".to_string() } else { from };
        let mut spec = ChannelSpec::new_pipe_rect(
            format!("segment_{}", i + 1),
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
        bp.add_channel(spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)));
    }

    // Arm 2 serpentine
    for i in 0..segments {
        let from = format!("arm2_{i}");
        let to = if i + 1 == segments {
            "merge_jn".to_string()
        } else {
            format!("arm2_{}", i + 1)
        };
        let actual_from = if i == 0 { "split_jn".to_string() } else { from };
        let mut spec = ChannelSpec::new_pipe_rect(
            format!("arm2_seg_{}", i + 1),
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
        bp.add_channel(spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)));
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

/// Double-level symmetric bifurcation with 4 parallel serpentine arms — closed loop.
///
/// **Topology**:
/// `inlet → inlet_section → split_1 → [split_2a → {arm1, arm2} → merge_2a,`
/// `                                    split_2b → {arm3, arm4} → merge_2b]`
/// `                               → merge_1 → trunk_out → outlet`
///
/// Four parallel serpentine arms (Q/4 each) provide superior flow uniformity
/// and well-plate coverage compared to a single-level bifurcation.  Power-of-2
/// branching gives exact equal flow splitting at every junction.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"branch_a_in"`, `"branch_b_in"` — level-1 to level-2 connectors
/// - `"arm1_seg_1"` … `"arm1_seg_n"` — arm-1 serpentine segments
/// - `"arm2_seg_1"` … `"arm2_seg_n"` — arm-2 serpentine segments
/// - `"arm3_seg_1"` … `"arm3_seg_n"` — arm-3 serpentine segments
/// - `"arm4_seg_1"` … `"arm4_seg_n"` — arm-4 serpentine segments
/// - `"branch_a_out"`, `"branch_b_out"` — level-2 to level-1 merge connectors
/// - `"trunk_out"` — downstream trunk
#[must_use]
pub fn double_bifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let sub_trunk = trunk_length_m * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    // ── Nodes ──
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("split_2a", NodeKind::Junction));
    bp.add_node(NodeSpec::new("split_2b", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_2a", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_2b", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_1", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for i in 1..segments {
        bp.add_node(NodeSpec::new(format!("arm1_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("arm2_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("arm3_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("arm4_{i}"), NodeKind::Junction));
    }

    // ── Inlet trunk ──
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_1",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // ── Level-1 → level-2 connectors ──
    for (name_ch, to) in [("branch_a_in", "split_2a"), ("branch_b_in", "split_2b")] {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                name_ch,
                "split_1",
                to,
                sub_trunk,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, sub_trunk, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // ── Four serpentine arms ──
    let arm_defs: [(&str, &str, &str, &str); 4] = [
        ("split_2a", "merge_2a", "arm1", "arm1_seg"),
        ("split_2a", "merge_2a", "arm2", "arm2_seg"),
        ("split_2b", "merge_2b", "arm3", "arm3_seg"),
        ("split_2b", "merge_2b", "arm4", "arm4_seg"),
    ];
    for (split_node, merge_node, node_prefix, seg_prefix) in arm_defs {
        for i in 0..segments {
            let actual_from = if i == 0 {
                split_node.to_string()
            } else {
                format!("{node_prefix}_{i}")
            };
            let to = if i + 1 == segments {
                merge_node.to_string()
            } else {
                format!("{node_prefix}_{}", i + 1)
            };
            let mut spec = ChannelSpec::new_pipe_rect(
                format!("{seg_prefix}_{}", i + 1),
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
            bp.add_channel(spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)));
        }
    }

    // ── Level-2 merge → level-1 merge connectors ──
    for (from, name_ch) in [("merge_2a", "branch_a_out"), ("merge_2b", "branch_b_out")] {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                name_ch,
                from,
                "merge_1",
                sub_trunk,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, sub_trunk, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // ── Outlet trunk ──
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_1",
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
