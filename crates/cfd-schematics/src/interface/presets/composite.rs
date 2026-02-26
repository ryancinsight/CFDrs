//! Composite closed-loop preset networks for millifluidic SDT device designs.
//!
//! Every function in this module produces a [`NetworkBlueprint`] with exactly
//! **one `NodeKind::Inlet`** and **one `NodeKind::Outlet`**.  Every split
//! junction is paired with a corresponding converging junction so the full
//! 127.76 × 85.47 mm 96-well-plate cuboid is always traversed completely.
//!
//! # Channel-name conventions (relied on by downstream consumers)
//! - `"inlet_section"` — first inlet / trunk approach channel
//! - `"throat_section"` — primary venturi throat (may also be `"throat_section_N"`)
//! - `"segment_1"`, `"segment_2"`, … — serpentine straight segments
//! - `"trunk_out"` — final trunk channel before the outlet

use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

const BLOOD_MU: f64 = 3.5e-3; // Pa·s (whole blood, high-shear Newtonian approx.)

// ── Resistance helpers ────────────────────────────────────────────────────────

/// Shah-London hydraulic resistance for a rectangular duct [Pa·s/m³].
fn shah_london(w: f64, h: f64, l: f64, mu: f64) -> f64 {
    let alpha = h.min(w) / h.max(w);
    let po = 96.0
        * (1.0
            - 1.3553 * alpha
            + 1.9467 * alpha.powi(2)
            - 1.7012 * alpha.powi(3)
            + 0.9564 * alpha.powi(4)
            - 0.2537 * alpha.powi(5));
    let d_h = 2.0 * w * h / (w + h);
    let area = w * h;
    po * mu * l / (d_h * d_h * area)
}

// ── 1. Venturi + Serpentine (series) ─────────────────────────────────────────

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
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    // ── Nodes ───────────────────────────────────────────────────────────────
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_jn", NodeKind::Junction));
    // n junction nodes between the venturi diffuser and the final outlet
    // serp_0 .. serp_{n-1}  (each is a 180° bend in the serpentine)
    for i in 0..segments {
        bp.add_node(NodeSpec::new(format!("serp_{i}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // ── Venturi channels ────────────────────────────────────────────────────
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

    // ── Serpentine channels ─────────────────────────────────────────────────
    for i in 0..segments {
        let from = format!("serp_{i}");
        let to = if i + 1 == segments {
            "outlet".to_string()
        } else {
            format!("serp_{}", i + 1)
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("segment_{}", i + 1),
                from,
                to,
                segment_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    bp
}

// ── 2. Bifurcation + Venturi in each branch ───────────────────────────────────

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

// ── 3. Trifurcation + Venturi in each branch ─────────────────────────────────

/// Symmetric trifurcation with a venturi throat in each branch — closed loop.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [3 × branch venturi] → merge_jn → trunk_out → outlet`
///
/// Three cavitation sites distributed across the three treatment-zone columns.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"throat_section"` — branch-1 throat (primary)
/// - `"throat_section_2"`, `"throat_section_3"` — branch 2–3 throats
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
    for i in 1..=3 {
        bp.add_node(NodeSpec::new(format!("contraction_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("throat_{i}"), NodeKind::Junction));
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

    for i in 1..=3 {
        let throat_name = if i == 1 {
            "throat_section".to_string()
        } else {
            format!("throat_section_{i}")
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("branch_{i}_in"),
                "split_jn",
                format!("contraction_{i}"),
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
                format!("contraction_{i}"),
                format!("throat_{i}"),
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
                format!("branch_{i}_out"),
                format!("throat_{i}"),
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

// ── 4. Bifurcation + Serpentine in each arm ───────────────────────────────────

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
        let actual_from = if i == 0 {
            "split_jn".to_string()
        } else {
            from
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("segment_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
    }

    // Arm 2 serpentine
    for i in 0..segments {
        let from = format!("arm2_{i}");
        let to = if i + 1 == segments {
            "merge_jn".to_string()
        } else {
            format!("arm2_{}", i + 1)
        };
        let actual_from = if i == 0 {
            "split_jn".to_string()
        } else {
            from
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("arm2_seg_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                main_width_m,
                height_m,
                shah_london(main_width_m, height_m, segment_length_m, BLOOD_MU),
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

// ── 5. Trifurcation + Serpentine in each arm ─────────────────────────────────

/// Symmetric trifurcation with a full serpentine in each arm — closed loop.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [3 × serpentine arm] → merge_jn → trunk_out → outlet`
///
/// Three parallel serpentine arms provide full 6×6 coverage.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"segment_1"` … `"segment_n"` — arm-1 serpentine segments
/// - `"arm2_seg_1"` … `"arm3_seg_n"` — arm-2, arm-3 segments
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
    // Intermediate junction nodes between segments (arm*_0 is split_jn itself; skip it)
    for arm in 1..=3 {
        for i in 1..segments {
            bp.add_node(NodeSpec::new(format!("arm{arm}_{i}"), NodeKind::Junction));
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

    for arm in 1..=3 {
        for i in 0..segments {
            let to = if i + 1 == segments {
                "merge_jn".to_string()
            } else {
                format!("arm{arm}_{}", i + 1)
            };
            let actual_from = if i == 0 {
                "split_jn".to_string()
            } else {
                format!("arm{arm}_{i}")
            };
            let ch_name = if arm == 1 {
                format!("segment_{}", i + 1)
            } else {
                format!("arm{arm}_seg_{}", i + 1)
            };
            bp.add_channel(
                ChannelSpec::new_pipe_rect(
                    ch_name,
                    actual_from,
                    to,
                    segment_length_m,
                    main_width_m,
                    height_m,
                    shah_london(main_width_m, height_m, segment_length_m, BLOOD_MU),
                    0.0,
                )
                .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
            );
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

// ── 6. Cell-separation (center venturi + peripheral bypass) ───────────────────

/// Inertial-focusing cell-separation network with closed loop.
///
/// **Topology**:
/// `inlet → approach → split_jn → [center: venturi (CancerTarget), peripheral_1 + peripheral_2 (HealthyBypass)] → merge_jn → outlet`
///
/// The curved approach section drives Dean-flow margination so that large, stiff
/// cells (cancer cells, WBCs) focus to the channel center while small, deformable
/// RBCs migrate to the walls.  At the split junction:
/// - Centre stream → venturi throat (sonodynamic treatment).
/// - Two peripheral streams → low-shear bypass (RBC protection).
///
/// # Channel names
/// - `"inlet_section"` — upstream approach (`main_width_m × height_m`)
/// - `"throat_section"` — centre venturi throat ([`TherapyZone::CancerTarget`])
/// - `"peripheral_1"`, `"peripheral_2"` — wall-side bypass ([`TherapyZone::HealthyBypass`])
#[must_use]
pub fn cell_separation_rect(
    name: impl Into<String>,
    approach_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    bypass_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // Approach
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_jn",
            approach_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, approach_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Centre venturi (cancer/WBC target)
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "center_approach",
            "split_jn",
            "contraction",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
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
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "center_diffuser",
            "throat_jn",
            "merge_jn",
            l_taper,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Peripheral bypass streams (×2, for RBC protection)
    for i in 1..=2 {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("peripheral_{i}"),
                "split_jn",
                "merge_jn",
                bypass_length_m,
                main_width_m * 0.5, // peripheral channels are narrower
                height_m,
                shah_london(main_width_m * 0.5, height_m, bypass_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
        );
    }

    bp
}

// ── 7. Serial double venturi ─────────────────────────────────────────────────

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
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
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

    bp
}

// ── 8. Double bifurcation + venturi (2-level Bi tree, 4 venturis) ─────────────

/// Two-level symmetric bifurcation with a venturi in each of the 4 leaves — closed loop.
///
/// **Topology**:
/// `inlet → split_0 → [split_L, split_R] → 4× venturis → [merge_L, merge_R] → merge_0 → outlet`
///
/// Four cavitation sites distributed over a 2×2 grid.
///
/// # Channel names
/// - `"inlet_section"` — level-0 trunk
/// - `"throat_section"` — leaf-LL venturi throat (primary)
/// - `"throat_section_LR"`, `"throat_section_RL"`, `"throat_section_RR"` — remaining throats
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
            bp.add_node(NodeSpec::new(format!("ctr_{side}{leaf}"), NodeKind::Junction));
            bp.add_node(NodeSpec::new(format!("throat_{side}{leaf}"), NodeKind::Junction));
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

// ── 9. Triple bifurcation + venturi (3-level Bi tree, 8 venturis) ─────────────

/// Three-level symmetric bifurcation with a venturi in each of the 8 leaves — closed loop.
///
/// **Topology**: 3-level binary tree → 8 venturis → converging binary tree → outlet.
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

    // Level-3 venturis (LLL … RRR)
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

// ── 10. Double trifurcation + venturi (2-level Tri tree, 9 venturis) ──────────

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

// ── 11. Bifurcation then trifurcation + venturi (2-level mixed, 6 venturis) ───

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

// ── 12. Trifurcation then bifurcation + venturi (2-level mixed, 6 venturis) ──

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

// ── 13. Asymmetric bifurcation + serpentine arms (Zweifach-Fung cell routing) ─

/// Asymmetric bifurcation — wide serpentine arm (cancer/WBC) + narrow bypass (RBC).
///
/// Exploits the **Zweifach-Fung bifurcation law**: at a T-junction, large stiff
/// cells (cancer cells, WBCs) preferentially enter the arm carrying the greater
/// volumetric flow (the wider, lower-resistance arm).  Small deformable RBCs
/// distribute more uniformly and partly enter the narrow bypass arm.
///
/// **Topology**:
/// `inlet → inlet_section → split_jn → [wide serpentine (CancerTarget), narrow bypass (HealthyBypass)] → merge_jn → trunk_out → outlet`
///
/// - Wide arm: `wide_width_m`; handles ≈ 89 % of flow (8× lower resistance).
/// - Narrow arm: `wide_width_m × 0.5`; handles ≈ 11 % of flow.
/// - Both arms have the same segment count and segment length.
///
/// # Channel names
/// - `"inlet_section"` — upstream trunk
/// - `"wide_seg_1"` … `"wide_seg_n"` — wide arm serpentine (cancer/WBC route)
/// - `"narrow_seg_1"` … `"narrow_seg_n"` — narrow arm bypass (RBC route)
/// - `"trunk_out"` — downstream trunk
#[must_use]
pub fn asymmetric_bifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    wide_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let narrow_width_m = wide_width_m * 0.5;
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("merge_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));
    // Intermediate nodes between segments (split_jn is "wide_0"/"narrow_0"; skip index 0)
    for i in 1..segments {
        bp.add_node(NodeSpec::new(format!("wide_{i}"), NodeKind::Junction));
        bp.add_node(NodeSpec::new(format!("narrow_{i}"), NodeKind::Junction));
    }

    // Upstream trunk (wide cross-section)
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "split_jn",
            trunk_length_m,
            wide_width_m,
            height_m,
            shah_london(wide_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Wide arm — carries ~89 % of total flow; cancer cells / WBCs route here.
    for i in 0..segments {
        let actual_from = if i == 0 {
            "split_jn".to_string()
        } else {
            format!("wide_{i}")
        };
        let to = if i + 1 == segments {
            "merge_jn".to_string()
        } else {
            format!("wide_{}", i + 1)
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("wide_seg_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                wide_width_m,
                height_m,
                shah_london(wide_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
        );
    }

    // Narrow arm — carries ~11 % of total flow; RBCs route here.
    for i in 0..segments {
        let actual_from = if i == 0 {
            "split_jn".to_string()
        } else {
            format!("narrow_{i}")
        };
        let to = if i + 1 == segments {
            "merge_jn".to_string()
        } else {
            format!("narrow_{}", i + 1)
        };
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("narrow_seg_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                narrow_width_m,
                height_m,
                shah_london(narrow_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
        );
    }

    // Downstream trunk (wide cross-section)
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "merge_jn",
            "outlet",
            trunk_length_m,
            wide_width_m,
            height_m,
            shah_london(wide_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}

// ── 14. Constriction-expansion array (N cycles, rectangular) ─────────────────

/// N alternating wide→narrow constriction-expansion cycles for cumulative WBC margination.
///
/// **Topology**:
/// `inlet → wide_0 → narrow_0 → wide_1 → narrow_1 → … → outlet`
///
/// Each narrow constriction forces cells to re-equilibrate laterally, enhancing
/// WBC margination through repeated inertial focusing events.  The cumulative
/// margination enhancement after N cycles follows:
///
/// ```text
/// E_N = 1 − (1 − E_1)^N
/// ```
///
/// where E_1 is the single-cycle WBC enrichment fraction.
///
/// # Arguments
/// - `n_cycles` — number of wide→narrow cycles (2–20 recommended)
/// - `wide_length_m` — length of each wide section [m]
/// - `narrow_length_m` — length of each narrow constriction [m]
/// - `wide_width_m` — width of wide sections [m]
/// - `narrow_width_m` — width of constrictions [m] (typically 0.3–0.7 × wide_width_m)
/// - `height_m` — channel height [m] (same for all sections)
///
/// # Channel names
/// - `"wide_0"` … `"wide_{n-1}"` — wide inlet sections ([`TherapyZone::MixedFlow`])
/// - `"narrow_0"` … `"narrow_{n-1}"` — constriction sections ([`TherapyZone::CancerTarget`])
///
/// # References
/// Wu, Z. et al. (2019). Integrated multifunctional microfluidics for simultaneous
/// isolation and detection of CTC. *Sci. Rep.* 9, 7356.
#[must_use]
pub fn constriction_expansion_array_rect(
    name: impl Into<String>,
    n_cycles: usize,
    wide_length_m: f64,
    narrow_length_m: f64,
    wide_width_m: f64,
    narrow_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    assert!(n_cycles >= 1, "n_cycles must be ≥ 1");
    let mut bp = NetworkBlueprint::new(name);

    // ── Nodes ────────────────────────────────────────────────────────────────
    // inlet, j_0 … j_{2n-2}, outlet  (2n-1 intermediate junctions)
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    for k in 0..2 * n_cycles - 1 {
        bp.add_node(NodeSpec::new(format!("j_{k}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // Helper: channel endpoint at position p in the sequence
    // p=0 → "inlet", p=1 … 2n-1 → "j_{p-1}", p=2n → "outlet"
    let node_at = |p: usize| -> String {
        if p == 0 {
            "inlet".to_string()
        } else if p == 2 * n_cycles {
            "outlet".to_string()
        } else {
            format!("j_{}", p - 1)
        }
    };

    // ── Channels ─────────────────────────────────────────────────────────────
    for i in 0..n_cycles {
        // Wide section: position 2i → 2i+1
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("wide_{i}"),
                node_at(2 * i),
                node_at(2 * i + 1),
                wide_length_m,
                wide_width_m,
                height_m,
                shah_london(wide_width_m, height_m, wide_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
        );
        // Narrow constriction: position 2i+1 → 2i+2
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("narrow_{i}"),
                node_at(2 * i + 1),
                node_at(2 * i + 2),
                narrow_length_m,
                narrow_width_m,
                height_m,
                shah_london(narrow_width_m, height_m, narrow_length_m, BLOOD_MU),
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
        );
    }

    bp
}

// ── 15. Spiral channel (N turns, rectangular) ────────────────────────────────

/// Tight spiral channel with N turns for Dean-flow dominant WBC/RBC separation.
///
/// **Topology**:
/// `inlet → spiral_0 → spiral_1 → … → spiral_{n-1} → outlet`
///
/// N identical arc segments, each representing one 360° spiral turn.  The
/// topology is a linear chain; the geometry layer renders the tight spiral.
/// All segments carry the same flow rate Q, giving a uniform Dean number:
///
/// ```text
/// De = Re · √(D_h / (2 · R_mean))
/// ```
///
/// # Arguments
/// - `n_turns` — number of complete 360° turns (2–20 recommended)
/// - `turn_length_m` — arc length per turn [m] (≈ 2π × R_mean)
/// - `width_m` — channel width [m]
/// - `height_m` — channel height [m]
///
/// # Channel names
/// - `"spiral_0"` … `"spiral_{n-1}"` — one turn per channel ([`TherapyZone::CancerTarget`])
///
/// # References
/// Nivedita, N. & Papautsky, I. (2013). Continuous separation of blood cells in spiral
/// microfluidic devices. *Biomicrofluidics*, 7, 054101.
#[must_use]
pub fn spiral_channel_rect(
    name: impl Into<String>,
    n_turns: usize,
    turn_length_m: f64,
    width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    assert!(n_turns >= 1, "n_turns must be ≥ 1");
    let mut bp = NetworkBlueprint::new(name);

    // ── Nodes ────────────────────────────────────────────────────────────────
    // inlet, turn_0 … turn_{n-2}, outlet  (n-1 intermediate junctions)
    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    for k in 0..n_turns.saturating_sub(1) {
        bp.add_node(NodeSpec::new(format!("turn_{k}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // Helper: node name at position p (0 = inlet, n = outlet)
    let node_at = |p: usize| -> String {
        if p == 0 {
            "inlet".to_string()
        } else if p == n_turns {
            "outlet".to_string()
        } else {
            format!("turn_{}", p - 1)
        }
    };

    // ── Channels ─────────────────────────────────────────────────────────────
    let r_seg = shah_london(width_m, height_m, turn_length_m, BLOOD_MU);
    for i in 0..n_turns {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("spiral_{i}"),
                node_at(i),
                node_at(i + 1),
                turn_length_m,
                width_m,
                height_m,
                r_seg,
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
        );
    }

    bp
}

// ── 16. Parallel microchannel array (N channels, rectangular) ────────────────

/// N identical parallel microchannels all connecting a single Inlet to a single Outlet.
///
/// **Topology**: `inlet ──[ch_0, ch_1, …, ch_{n-1}]──▶ outlet` (star topology)
///
/// All N channels have identical cross-section, length, and hydraulic resistance.
/// The 1D solver distributes flow equally: Q_per_channel = Q_total / N.
///
/// Enables micro-scale inertial focusing (D_h < 150 µm, κ_WBC > 0.15) at
/// clinical throughput (≥ 10 mL/min) by replicating N identical units within
/// a single 96-well-plate millifluidic chip.
///
/// # Theorem — Parallel Resistance
/// For N identical channels in parallel, each with resistance R_ch:
/// ```text
/// R_total = R_ch / N
/// Q_i = Q_total / N   (uniform distribution by symmetry)
/// ```
///
/// # Arguments
/// - `n_channels` — number of parallel microchannels (10–500 recommended)
/// - `channel_length_m` — length of each channel [m]
/// - `channel_width_m` — width of each channel [m] (micro-scale, < 200 µm for focusing)
/// - `channel_height_m` — height of each channel [m]
///
/// # Channel names
/// - `"ch_0"` … `"ch_{n-1}"` — N identical parallel channels ([`TherapyZone::CancerTarget`])
///
/// # References
/// Bhagat, A. A. S. et al. (2010). Dean flow fractionation of Cryptosporidium oocysts.
/// *Lab Chip*, 10, 2605–2614.
#[must_use]
pub fn parallel_microchannel_array_rect(
    name: impl Into<String>,
    n_channels: usize,
    channel_length_m: f64,
    channel_width_m: f64,
    channel_height_m: f64,
) -> NetworkBlueprint {
    assert!(n_channels >= 1, "n_channels must be ≥ 1");
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    let r_ch = shah_london(channel_width_m, channel_height_m, channel_length_m, BLOOD_MU);
    for i in 0..n_channels {
        bp.add_channel(
            ChannelSpec::new_pipe_rect(
                format!("ch_{i}"),
                "inlet",
                "outlet",
                channel_length_m,
                channel_width_m,
                channel_height_m,
                r_ch,
                0.0,
            )
            .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
        );
    }

    bp
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::NodeKind;

    fn count_inlets_outlets(bp: &NetworkBlueprint) -> (usize, usize) {
        let inlets = bp.nodes.iter().filter(|n| n.kind == NodeKind::Inlet).count();
        let outlets = bp.nodes.iter().filter(|n| n.kind == NodeKind::Outlet).count();
        (inlets, outlets)
    }

    #[test]
    fn venturi_serpentine_has_single_inlet_outlet() {
        let bp = venturi_serpentine_rect("t", 2e-3, 0.5e-3, 0.5e-3, 1e-3, 6, 7.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
    }

    #[test]
    fn bifurcation_venturi_has_single_inlet_outlet() {
        let bp = bifurcation_venturi_rect("t", 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn trifurcation_venturi_has_single_inlet_outlet() {
        let bp = trifurcation_venturi_rect("t", 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn bifurcation_serpentine_has_single_inlet_outlet() {
        let bp = bifurcation_serpentine_rect("t", 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn trifurcation_serpentine_has_single_inlet_outlet() {
        let bp = trifurcation_serpentine_rect("t", 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn cell_separation_has_single_inlet_outlet() {
        let bp = cell_separation_rect("t", 22.5e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3, 22.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn serial_double_venturi_has_single_inlet_outlet() {
        let bp = serial_double_venturi_rect("t", 2e-3, 0.5e-3, 0.5e-3, 1e-3, 20e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn double_bifurcation_venturi_has_single_inlet_outlet() {
        let bp = double_bifurcation_venturi_rect("t", 20e-3, 10e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn triple_bifurcation_venturi_has_single_inlet_outlet() {
        let bp =
            triple_bifurcation_venturi_rect("t", 15e-3, 8e-3, 5e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn double_trifurcation_venturi_has_single_inlet_outlet() {
        let bp = double_trifurcation_venturi_rect("t", 20e-3, 10e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn bifurcation_trifurcation_venturi_has_single_inlet_outlet() {
        let bp =
            bifurcation_trifurcation_venturi_rect("t", 20e-3, 10e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn trifurcation_bifurcation_venturi_has_single_inlet_outlet() {
        let bp =
            trifurcation_bifurcation_venturi_rect("t", 20e-3, 10e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
    }

    #[test]
    fn asymmetric_bifurcation_serpentine_has_single_inlet_outlet() {
        let bp = asymmetric_bifurcation_serpentine_rect("t", 20e-3, 6, 7.5e-3, 2e-3, 0.5e-3);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
    }

    #[test]
    fn constriction_expansion_array_has_single_inlet_outlet() {
        let bp = constriction_expansion_array_rect("t", 10, 500e-6, 250e-6, 300e-6, 150e-6, 60e-6);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
        // 10 cycles → 10 wide + 10 narrow = 20 channels
        assert_eq!(bp.channels.len(), 20, "10 cycles → 20 channels");
        // 2*10 - 1 = 19 intermediate junctions + 1 inlet + 1 outlet = 21 nodes
        assert_eq!(bp.nodes.len(), 21, "10 cycles → 21 nodes");
    }

    #[test]
    fn spiral_channel_has_single_inlet_outlet() {
        let bp = spiral_channel_rect("t", 8, 11e-3, 400e-6, 80e-6);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
        // 8 turns → 8 channels
        assert_eq!(bp.channels.len(), 8, "8 turns → 8 channels");
        // 8-1 = 7 intermediate junctions + inlet + outlet = 9 nodes
        assert_eq!(bp.nodes.len(), 9, "8 turns → 9 nodes");
    }

    #[test]
    fn parallel_microchannel_array_has_single_inlet_outlet() {
        let bp = parallel_microchannel_array_rect("t", 50, 30e-3, 100e-6, 60e-6);
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1, "must have exactly 1 inlet");
        assert_eq!(o, 1, "must have exactly 1 outlet");
        assert_eq!(bp.channels.len(), 50, "N=50 → 50 channels");
        assert_eq!(bp.nodes.len(), 2, "parallel array → only inlet + outlet");
    }

    #[test]
    fn constriction_array_single_cycle() {
        // Minimal n_cycles=1: inlet → j_0 → outlet (1 wide + 1 narrow)
        let bp = constriction_expansion_array_rect("t", 1, 500e-6, 250e-6, 300e-6, 150e-6, 60e-6);
        assert_eq!(bp.channels.len(), 2, "1 cycle → 2 channels");
        assert_eq!(bp.nodes.len(), 3, "1 cycle → inlet + j_0 + outlet");
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn spiral_single_turn() {
        // n_turns=1: inlet → outlet directly
        let bp = spiral_channel_rect("t", 1, 11e-3, 400e-6, 80e-6);
        assert_eq!(bp.channels.len(), 1, "1 turn → 1 channel");
        assert_eq!(bp.nodes.len(), 2, "1 turn → inlet + outlet");
        let (i, o) = count_inlets_outlets(&bp);
        assert_eq!(i, 1);
        assert_eq!(o, 1);
    }

    #[test]
    fn all_nodes_referenced_by_channels() {
        // Every node ID that appears in channels should exist as a node
        let bp = bifurcation_venturi_rect("t", 20e-3, 2e-3, 0.5e-3, 0.5e-3, 1e-3);
        let node_ids: std::collections::HashSet<&str> =
            bp.nodes.iter().map(|n| n.id.as_str()).collect();
        for ch in &bp.channels {
            assert!(
                node_ids.contains(ch.from.as_str()),
                "channel '{}' references missing from-node '{}'",
                ch.id.as_str(),
                ch.from.as_str()
            );
            assert!(
                node_ids.contains(ch.to.as_str()),
                "channel '{}' references missing to-node '{}'",
                ch.id.as_str(),
                ch.to.as_str()
            );
        }
    }
}
