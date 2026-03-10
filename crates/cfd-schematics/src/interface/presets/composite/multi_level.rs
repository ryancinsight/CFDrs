//! Multi-level bifurcation-tree composite presets: 2- and 3-level binary branching.
#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use super::{shah_london, BLOOD_MU};
use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

/// Two-level symmetric bifurcation with a venturi in each of the 4 leaves â€” closed loop.
///
/// **Topology**:
/// `inlet â†’ split_0 â†’ [split_L, split_R] â†’ 4Ã— venturis â†’ [merge_L, merge_R] â†’ merge_0 â†’ outlet`
///
/// Four cavitation sites distributed over a 2Ã—2 grid.
///
/// # Channel names
/// - `"inlet_section"` â€” level-0 trunk
/// - `"throat_section"` â€” leaf-LL venturi throat (primary)
/// - `"throat_section_LR"`, `"throat_section_RL"`, `"throat_section_RR"` â€” remaining throats
#[must_use]
pub fn double_bifurcation_venturi_rect(
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

    // Level-0 junctions
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    // Level-1 junctions
    for side in &["L", "R"] {
        bp.add_node(NodeSpec::new(format!("split_{side}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{side}"), NodeKind::Junction));
        // Level-2 venturi junctions
        for leaf in &["L", "R"] {
            bp.add_node(NodeSpec::new(
                format!("ctr_{side}{leaf}"),
                NodeKind::Junction,
            ));
            bp.add_node(NodeSpec::new(
                format!("throat_{side}{leaf}"),
                NodeKind::Junction,
            ));
        }
    }

    // Level-0 trunk
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

    for side in &["L", "R"] {
        // Level-1 distribution branches
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("branch_{side}_in"),
                "split_0",
                format!("split_{side}"),
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
                format!("branch_{side}_out"),
                format!("merge_{side}"),
                "merge_0",
                branch1_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch1_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        // Level-2 venturi in each leaf (LL, LR or RL, RR)
        for leaf in &["L", "R"] {
            let label = format!("{side}{leaf}");
            let throat_name = if label == "LL" {
                "throat_section".to_string()
            } else {
                format!("throat_section_{label}")
            };
            bp.add_channel(
                ChannelSpec::new_pipe_rect(
                    format!("leaf_{label}_in"),
                    format!("split_{side}"),
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
                    format!("merge_{side}"),
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

/// Three-level symmetric bifurcation with a venturi in each of the 8 leaves â€” closed loop.
///
/// **Topology**: 3-level binary tree â†’ 8 venturis â†’ converging binary tree â†’ outlet.
///
/// Eight cavitation sites covering all octant zones of the treatment area.
#[must_use]
pub fn triple_bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    branch2_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    let labels_l1 = ["L", "R"];
    let labels_l2 = ["LL", "LR", "RL", "RR"];
    let labels_l3 = ["LLL", "LLR", "LRL", "LRR", "RLL", "RLR", "RRL", "RRR"];

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for l in &labels_l1 {
        bp.add_node(NodeSpec::new(format!("split_{l}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l}"), NodeKind::Junction));
    }
    for l in &labels_l2 {
        bp.add_node(NodeSpec::new(format!("split_{l}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l}"), NodeKind::Junction));
    }
    for l in &labels_l3 {
        bp.add_node(NodeSpec::new(format!("ctr_{l}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{l}"), NodeKind::Junction));
    }

    // Trunk
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

    // Level-1 branches (L, R)
    for l1 in &labels_l1 {
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

    // Level-2 branches (LL, LR, RL, RR)
    for l2 in &labels_l2 {
        let parent = &l2[0..1]; // "L" or "R"
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b2_{l2}_in"),
                format!("split_{parent}"),
                format!("split_{l2}"),
                branch2_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch2_length_m, BLOOD_MU),
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
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, branch2_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // Level-3 venturis (LLL â€¦ RRR)
    for (idx, l3) in labels_l3.iter().enumerate() {
        let parent = &l3[0..2]; // "LL", "LR", "RL", "RR"
        let throat_name = if idx == 0 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{l3}")
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("leaf_{l3}_in"),
                format!("split_{parent}"),
                format!("ctr_{l3}"),
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
                format!("ctr_{l3}"),
                format!("throat_{l3}"),
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
                format!("leaf_{l3}_out"),
                format!("throat_{l3}"),
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

/// Three-level symmetric trifurcation with a venturi in each of the 27 leaves.
///
/// **Topology**: `inlet → [Tri] → [Tri×3] → [Tri×9] → 27 venturis → converging tree → outlet`
///
/// Channel widths are **scaled at each split** using `center_frac`:
/// - Center arm (`"B"`) = parent_width × center_frac
/// - Each peripheral arm (`"A"`, `"C"`) = parent_width × (1 − center_frac) / 2
///
/// With `center_frac = 1/3` (symmetric), all arms are equal width.
/// With `center_frac > 1/3`, the center arm is wider — enabling Zweifach-Fung
/// junction routing of large/stiff cells (cancer, WBC) toward the center arm.
///
/// # Channel names
/// - `"inlet_section"` — trunk
/// - `"throat_section"` — first leaf venturi (label `"AAA"`)
/// - `"throat_section_{LABEL}"` — remaining 26 throats (AAB … CCC)
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn triple_trifurcation_venturi_rect(
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

    let labels_l1 = ["A", "B", "C"]; // center = B
    let labels_l2: Vec<String> = labels_l1
        .iter()
        .flat_map(|a| labels_l1.iter().map(move |b| format!("{a}{b}")))
        .collect();
    let labels_l3: Vec<String> = labels_l2
        .iter()
        .flat_map(|a| labels_l1.iter().map(move |b| format!("{a}{b}")))
        .collect();

    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for l1 in &labels_l1 {
        bp.add_node(NodeSpec::new(format!("split_{l1}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l1}"), NodeKind::Junction));
    }
    for l2 in &labels_l2 {
        bp.add_node(NodeSpec::new(format!("split_{l2}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l2}"), NodeKind::Junction));
    }
    for l3 in &labels_l3 {
        bp.add_node(NodeSpec::new(format!("ctr_{l3}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{l3}"), NodeKind::Junction));
    }

    // Arm width: 'B' = center_frac × parent, 'A'/'C' = periph_frac × parent
    let w_for_label = |label: &str| -> f64 {
        label.chars().fold(main_width_m, |w, ch| {
            w * if ch == 'B' { center_frac } else { periph_frac }
        })
    };

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

    // Level-1 branches
    for l1 in &labels_l1 {
        let w = w_for_label(l1);
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

    // Level-2 branches
    for l2 in &labels_l2 {
        let w = w_for_label(l2);
        let parent = &l2[0..1];
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

    // Level-3 venturis (27 leaves)
    for (idx, l3) in labels_l3.iter().enumerate() {
        let parent = &l3[0..2];
        let throat_name = if idx == 0 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{l3}")
        };
        let w_leaf = w_for_label(l3);
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("leaf_{l3}_in"),
                format!("split_{parent}"),
                format!("ctr_{l3}"),
                l_taper,
                w_leaf,
                height_m,
                shah_london(w_leaf, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                throat_name,
                format!("ctr_{l3}"),
                format!("throat_{l3}"),
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
                format!("leaf_{l3}_out"),
                format!("throat_{l3}"),
                format!("merge_{parent}"),
                l_taper,
                w_leaf,
                height_m,
                shah_london(w_leaf, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp
}

/// Four-level symmetric trifurcation with a venturi in each of the 81 leaves.
///
/// **Topology**: `inlet → [Tri^4] → 81 venturis → converging Tri^4 → outlet`
///
/// Channel widths scale at each level using `center_frac`.
/// With `center_frac ≈ 1/3`, symmetric split → terminal channels reach
/// `w ≈ main_width × (1/3)^4 ≈ 0.012 × main_width` (possibly below fab limit;
/// optimizer will discard via `plate_fits` check).
///
/// # Channel names
/// - `"inlet_section"` — trunk
/// - `"throat_section"` — first leaf (label `"AAAA"`)
/// - `"throat_section_{LABEL}"` — remaining 80 throats
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn quad_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    branch2_length_m: f64,
    branch3_length_m: f64,
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

    let labels_l1 = ["A", "B", "C"];
    let labels_l2: Vec<String> = labels_l1
        .iter()
        .flat_map(|a| labels_l1.iter().map(move |b| format!("{a}{b}")))
        .collect();
    let labels_l3: Vec<String> = labels_l2
        .iter()
        .flat_map(|a| labels_l1.iter().map(move |b| format!("{a}{b}")))
        .collect();
    let labels_l4: Vec<String> = labels_l3
        .iter()
        .flat_map(|a| labels_l1.iter().map(move |b| format!("{a}{b}")))
        .collect();

    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_0", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    for l1 in &labels_l1 {
        bp.add_node(NodeSpec::new(format!("split_{l1}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l1}"), NodeKind::Junction));
    }
    for l2 in &labels_l2 {
        bp.add_node(NodeSpec::new(format!("split_{l2}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l2}"), NodeKind::Junction));
    }
    for l3 in &labels_l3 {
        bp.add_node(NodeSpec::new(format!("split_{l3}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("merge_{l3}"), NodeKind::Junction));
    }
    for l4 in &labels_l4 {
        bp.add_node(NodeSpec::new(format!("ctr_{l4}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{l4}"), NodeKind::Junction));
    }

    // Arm width: 'B' = center_frac × parent, 'A'/'C' = periph_frac × parent
    let w_for_label = |label: &str| -> f64 {
        label.chars().fold(main_width_m, |w, ch| {
            w * if ch == 'B' { center_frac } else { periph_frac }
        })
    };

    let branch_lens = [branch1_length_m, branch2_length_m, branch3_length_m];

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

    // Levels 1–3 internal branches
    for l1 in &labels_l1 {
        let w = w_for_label(l1);
        let len = branch_lens[0];
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_in"),
                "split_0",
                format!("split_{l1}"),
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b1_{l1}_out"),
                format!("merge_{l1}"),
                "merge_0",
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }
    for l2 in &labels_l2 {
        let w = w_for_label(l2);
        let len = branch_lens[1];
        let parent = &l2[0..1];
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b2_{l2}_in"),
                format!("split_{parent}"),
                format!("split_{l2}"),
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b2_{l2}_out"),
                format!("merge_{l2}"),
                format!("merge_{parent}"),
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }
    for l3 in &labels_l3 {
        let w = w_for_label(l3);
        let len = branch_lens[2];
        let parent = &l3[0..2];
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b3_{l3}_in"),
                format!("split_{parent}"),
                format!("split_{l3}"),
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("b3_{l3}_out"),
                format!("merge_{l3}"),
                format!("merge_{parent}"),
                len,
                w,
                height_m,
                shah_london(w, height_m, len, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // Level-4 venturis (81 leaves)
    for (idx, l4) in labels_l4.iter().enumerate() {
        let parent = &l4[0..3];
        let throat_name = if idx == 0 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{l4}")
        };
        let w_leaf = w_for_label(l4);
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("leaf_{l4}_in"),
                format!("split_{parent}"),
                format!("ctr_{l4}"),
                l_taper,
                w_leaf,
                height_m,
                shah_london(w_leaf, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                throat_name,
                format!("ctr_{l4}"),
                format!("throat_{l4}"),
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
                format!("leaf_{l4}_out"),
                format!("throat_{l4}"),
                format!("merge_{parent}"),
                l_taper,
                w_leaf,
                height_m,
                shah_london(w_leaf, height_m, l_taper, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp
}
