//! Mixed-tree composite presets: topologies combining bifurcation and trifurcation at different levels.

use super::{shah_london, BLOOD_MU};
use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

/// Two-level symmetric trifurcation with a venturi in each of the 9 leaves — closed loop.
///
/// **Topology**: 2-level ternary tree → 9 venturis → converging ternary tree → outlet.
///
/// Nine cavitation sites in a 3×3 grid covering the full treatment zone.
#[must_use]
pub fn double_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    let l1_labels = ["A", "B", "C"];
    // Level-2 labels: AA, AB, AC, BA, BB, BC, CA, CB, CC
    let l2_labels: Vec<String> = l1_labels
        .iter()
        .flat_map(|a| l1_labels.iter().map(move |b| format!("{a}{b}")))
        .collect();

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for l1 in &l1_labels {
        bp.add_node(NodeSpec::new(format!("split_{l1}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l1}"), NodeKind::Junction));
    }
    for l2 in &l2_labels {
        bp.add_node(NodeSpec::new(format!("ctr_{l2}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{l2}"), NodeKind::Junction));
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_0",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_0",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    for l1 in &l1_labels {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_in"),
                "split_0",
                format!("split_{l1}"),
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_out"),
                format!("merge_{l1}"),
                "merge_0",
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    for (idx, l2) in l2_labels.iter().enumerate() {
        let parent = &l2[0..1];
        let throat_name = if idx == 0 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{l2}")
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("leaf_{l2}_in"),
                format!("split_{parent}"),
                format!("ctr_{l2}"),
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
                throat_name,
                format!("ctr_{l2}"),
                format!("throat_{l2}"),
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
                format!("leaf_{l2}_out"),
                format!("throat_{l2}"),
                format!("merge_{parent}"),
                l_taper,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp
}

/// Mixed split: bifurcation at level-1, trifurcation at level-2 → 6 venturis — closed loop.
///
/// **Topology**: bi split → [tri, tri] → 6 venturis → [tri-merge, tri-merge] → bi-merge → outlet.
///
/// Six cavitation sites in a 2×3 grid.
#[must_use]
pub fn bifurcation_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    // Level-1: binary (L, R)
    let l1 = ["L", "R"];
    // Level-2: ternary per branch (A, B, C)
    let l2_sub = ["A", "B", "C"];

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for s in &l1 {
        bp.add_node(NodeSpec::new(format!("split_{s}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{s}"), NodeKind::Junction));
        for t in &l2_sub {
            bp.add_node(NodeSpec::new(format!("ctr_{s}{t}"), NodeKind::Junction));
            bp.add_node(NodeSpec::new(format!("throat_{s}{t}"), NodeKind::Junction));
        }
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_0",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_0",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    let mut throat_idx = 0usize;
    for s in &l1 {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{s}_in"),
                "split_0",
                format!("split_{s}"),
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{s}_out"),
                format!("merge_{s}"),
                "merge_0",
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        for t in &l2_sub {
            let label = format!("{s}{t}");
            let throat_name = if throat_idx == 0 {
                "throat_section".to_string()
            } else {
                format!("throat_section_{label}")
            };
            throat_idx += 1;
            bp.add_channel(
                ChannelSpec::new_pipe_rect(
                    format!("leaf_{label}_in"),
                    format!("split_{s}"),
                    format!("ctr_{label}"),
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
                    throat_name,
                    format!("ctr_{label}"),
                    format!("throat_{label}"),
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
                    format!("leaf_{label}_out"),
                    format!("throat_{label}"),
                    format!("merge_{s}"),
                    l_taper,
                    main_width_m,
                    height_m,
                    shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
                    0.0,
                )
                .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
            );
        }
    }

    bp
}

/// Mixed split: trifurcation at level-1, bifurcation at level-2 → 6 venturis — closed loop.
///
/// **Topology**: tri split → [bi, bi, bi] → 6 venturis → [bi-merge ×3] → tri-merge → outlet.
///
/// Six cavitation sites in a 3×2 grid.
///
/// # Channel names
/// - `"inlet_section"` — level-0 trunk
/// - `"throat_section"` — leaf-AL venturi throat (primary)
/// - `"throat_section_AR"`, `"throat_section_BL"`, … — remaining throats
#[must_use]
pub fn trifurcation_bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    // Level-1: ternary (A, B, C); level-2: binary (L, R)
    let l1 = ["A", "B", "C"];
    let l2_sub = ["L", "R"];

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for s in &l1 {
        bp.add_node(NodeSpec::new(format!("split_{s}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{s}"), NodeKind::Junction));
        for t in &l2_sub {
            bp.add_node(NodeSpec::new(format!("ctr_{s}{t}"), NodeKind::Junction));
            bp.add_node(NodeSpec::new(format!("throat_{s}{t}"), NodeKind::Junction));
        }
    }

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_0",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_0",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    let mut throat_idx = 0usize;
    for s in &l1 {
        // Level-1 trifurcation distribution branches
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{s}_in"),
                "split_0",
                format!("split_{s}"),
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{s}_out"),
                format!("merge_{s}"),
                "merge_0",
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );

        // Level-2 bifurcation + venturi at each leaf
        for t in &l2_sub {
            let label = format!("{s}{t}");
            let throat_name = if throat_idx == 0 {
                "throat_section".to_string()
            } else {
                format!("throat_section_{label}")
            };
            throat_idx += 1;
            bp.add_channel(
                ChannelSpec::new_pipe_rect(
                    format!("leaf_{label}_in"),
                    format!("split_{s}"),
                    format!("ctr_{label}"),
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
                    throat_name,
                    format!("ctr_{label}"),
                    format!("throat_{label}"),
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
                    format!("leaf_{label}_out"),
                    format!("throat_{label}"),
                    format!("merge_{s}"),
                    l_taper,
                    main_width_m,
                    height_m,
                    shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
                    0.0,
                )
                .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
            );
        }
    }

    bp
}

/// Trifurcation → Bifurcation → Bifurcation → 12 parallel venturi throats.
///
/// **Topology**:
/// `inlet → split_0 → [A/B/C] → [AL/AR, BL/BR, CL/CR] → [ALL/ALR, … CRR] → 12 venturis → merging tree → outlet`
///
/// Mixed 3-level tree: one ternary split at level 1, two binary splits at levels 2–3.
/// Width-scaled at the trifurcation split using `center_frac`; bifurcation arms
/// halve width at each subsequent level.
///
/// # Channel names
/// - `"inlet_section"` — trunk
/// - `"throat_section"` — first leaf (label `"ALL"`)
/// - `"throat_section_{LABEL}"` — remaining 11 throats
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn trifurcation_bifurcation_bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    branch2_length_m: f64,
    main_width_m: f64,
    center_frac: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let center_frac = center_frac.clamp(0.20, 0.70);
    let periph_frac = (1.0 - center_frac) * 0.5;

    // Level-1 (trifurcation): A, B (center), C
    let tri_labels = ["A", "B", "C"];
    // Level-2 (bifurcation of each tri arm): {A,B,C}{L,R}
    let bi_labels: Vec<String> = tri_labels
        .iter()
        .flat_map(|t| ["L", "R"].iter().map(move |b| format!("{t}{b}")))
        .collect();
    // Level-3 (bifurcation of each level-2 arm): {A,B,C}{L,R}{L,R}
    let leaf_labels: Vec<String> = bi_labels
        .iter()
        .flat_map(|t| ["L", "R"].iter().map(move |b| format!("{t}{b}")))
        .collect();

    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for l1 in &tri_labels {
        bp.add_node(NodeSpec::new(format!("split_{l1}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l1}"), NodeKind::Junction));
    }
    for l2 in &bi_labels {
        bp.add_node(NodeSpec::new(format!("split_{l2}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l2}"), NodeKind::Junction));
    }
    for leaf in &leaf_labels {
        bp.add_node(NodeSpec::new(format!("ctr_{leaf}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{leaf}"), NodeKind::Junction));
    }

    // Width at each level:
    //   Level-1: trifurcation width scaling
    //   Level-2: halved from level-1 (bifurcation — equal split)
    //   Level-3: halved from level-2
    let w_l1 =
        |l1: &str| -> f64 { main_width_m * if l1 == "B" { center_frac } else { periph_frac } };
    let w_l2 = |l2: &str| -> f64 { w_l1(&l2[0..1]) * 0.5 };
    let w_leaf = |leaf: &str| -> f64 { w_l2(&leaf[0..2]) * 0.5 };

    // Trunk in/out
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_0",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_0",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Level-1 trifurcation branches
    for l1 in &tri_labels {
        let w = w_l1(l1);
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_in"),
                "split_0",
                format!("split_{l1}"),
                branch1_length_m,
                w,
                height_m,
                shah_london(w, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_out"),
                format!("merge_{l1}"),
                "merge_0",
                branch1_length_m,
                w,
                height_m,
                shah_london(w, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // Level-2 bifurcation branches
    for l2 in &bi_labels {
        let parent = &l2[0..1];
        let w = w_l2(l2);
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b2_{l2}_in"),
                format!("split_{parent}"),
                format!("split_{l2}"),
                branch2_length_m,
                w,
                height_m,
                shah_london(w, height_m, branch2_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b2_{l2}_out"),
                format!("merge_{l2}"),
                format!("merge_{parent}"),
                branch2_length_m,
                w,
                height_m,
                shah_london(w, height_m, branch2_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // Level-3 venturis (12 leaves)
    for (idx, leaf) in leaf_labels.iter().enumerate() {
        let parent = &leaf[0..2];
        let throat_name = if idx == 0 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{leaf}")
        };
        let w = w_leaf(leaf);
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("leaf_{leaf}_in"),
                format!("split_{parent}"),
                format!("ctr_{leaf}"),
                l_taper,
                w,
                height_m,
                shah_london(w, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                throat_name,
                format!("ctr_{leaf}"),
                format!("throat_{leaf}"),
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
                format!("leaf_{leaf}_out"),
                format!("throat_{leaf}"),
                format!("merge_{parent}"),
                l_taper,
                w,
                height_m,
                shah_london(w, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp
}
