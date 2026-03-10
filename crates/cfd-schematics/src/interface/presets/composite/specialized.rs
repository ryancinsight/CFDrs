//! Specialized composite presets: cell separation, asymmetric, constriction,
//! spiral, and parallel microchannel array topologies.
#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use super::{shah_london, specialized_layout, BLOOD_MU};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use crate::geometry::generator::{
    create_primitive_selective_tree_geometry, create_selective_tree_geometry,
    CenterSerpentinePathSpec, PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
    SelectiveTreeRequest, SelectiveTreeTopology,
};
use crate::geometry::metadata::{AsymmetricTrifurcationParams, VenturiGeometryMetadata};
use crate::topology::presets::{
    constriction_expansion_series_spec, parallel_microchannel_array_spec,
    spiral_serpentine_series_spec,
};
use crate::BlueprintTopologyFactory;

/// Optional serpentine geometry applied only to center treatment lanes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenterSerpentineSpec {
    pub segments: usize,
    pub bend_radius_m: f64,
    pub segment_length_m: f64,
}

fn normalized_center_serpentine(
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

fn generator_center_serpentine(
    center_serpentine: Option<CenterSerpentineSpec>,
) -> Option<CenterSerpentinePathSpec> {
    normalized_center_serpentine(center_serpentine).map(|spec| CenterSerpentinePathSpec {
        segments: spec.segments,
        bend_radius_m: spec.bend_radius_m,
    })
}

fn venturi_metadata(
    throat_width_m: f64,
    throat_height_m: f64,
    throat_length_m: f64,
    inlet_width_m: f64,
    outlet_width_m: f64,
    taper_length_m: f64,
) -> VenturiGeometryMetadata {
    let half_delta_in = ((inlet_width_m - throat_width_m) * 0.5).abs();
    let half_delta_out = ((outlet_width_m - throat_width_m) * 0.5).abs();
    let safe_taper = taper_length_m.max(1e-9);
    let conv_angle = half_delta_in.atan2(safe_taper).to_degrees().max(1.0);
    let diff_angle = half_delta_out.atan2(safe_taper).to_degrees().max(1.0);
    VenturiGeometryMetadata {
        throat_width_m,
        throat_height_m,
        throat_length_m,
        inlet_width_m,
        outlet_width_m,
        convergent_half_angle_deg: conv_angle,
        divergent_half_angle_deg: diff_angle,
        throat_position: 0.5,
    }
}

fn with_venturi_geometry_metadata(
    channel: ChannelSpec,
    metadata: VenturiGeometryMetadata,
) -> ChannelSpec {
    channel
        .with_venturi_geometry(metadata.clone())
        .with_metadata(metadata)
}

#[must_use]
pub fn primitive_selective_split_tree_rect(
    name: impl Into<String>,
    box_dims_mm: (f64, f64),
    split_sequence: &[PrimitiveSelectiveSplitKind],
    main_width_m: f64,
    first_trifurcation_center_frac: f64,
    later_trifurcation_center_frac: f64,
    bifurcation_treatment_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
    venturi_treatment_enabled: bool,
    treatment_branch_throat_count: u8,
    center_serpentine: Option<CenterSerpentineSpec>,
) -> NetworkBlueprint {
    let request = PrimitiveSelectiveTreeRequest {
        name: name.into(),
        box_dims_mm,
        split_sequence: split_sequence.to_vec(),
        main_width_m,
        throat_width_m,
        throat_length_m,
        channel_height_m: height_m,
        first_trifurcation_center_frac,
        later_trifurcation_center_frac,
        bifurcation_treatment_frac,
        treatment_branch_venturi_enabled: venturi_treatment_enabled,
        treatment_branch_throat_count,
        center_serpentine: generator_center_serpentine(center_serpentine),
    };
    create_primitive_selective_tree_geometry(&request)
}

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
    if let Some(channel) = bp
        .channels
        .iter_mut()
        .find(|channel| channel.id.as_str() == "throat_section")
    {
        let metadata = venturi_metadata(
            throat_width_m,
            height_m,
            l_throat,
            main_width_m,
            main_width_m,
            l_taper,
        );
        channel.venturi_geometry = Some(metadata.clone());
        channel
            .metadata
            .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
            .insert(metadata);
    }
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
    narrow_frac: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let narrow_width_m = wide_width_m * narrow_frac.clamp(0.10, 0.90);
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
        bp.add_channel({
            let mut spec = ChannelSpec::new_pipe_rect(
                format!("wide_seg_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                wide_width_m,
                height_m,
                shah_london(wide_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            );
            spec.channel_shape = ChannelShape::Serpentine {
                segments,
                bend_radius_m: wide_width_m * 0.5,
            };
            spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget))
        });
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
        bp.add_channel({
            let mut spec = ChannelSpec::new_pipe_rect(
                format!("narrow_seg_{}", i + 1),
                actual_from,
                to,
                segment_length_m,
                narrow_width_m,
                height_m,
                shah_london(narrow_width_m, height_m, segment_length_m, BLOOD_MU),
                0.0,
            );
            spec.channel_shape = ChannelShape::Serpentine {
                segments,
                bend_radius_m: narrow_width_m * 0.5,
            };
            spec.with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass))
        });
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
    let name = name.into();
    assert!(n_cycles >= 1, "n_cycles must be ≥ 1");
    let topology = constriction_expansion_series_spec(
        &name,
        n_cycles,
        wide_length_m,
        narrow_length_m,
        wide_width_m,
        narrow_width_m,
        height_m,
    );
    let lineage = BlueprintTopologyFactory::lineage_for_spec(&topology);
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

    bp.with_topology_spec(topology).with_lineage(lineage)
}

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
    let name = name.into();
    assert!(n_turns >= 1, "n_turns must be ≥ 1");
    let topology = spiral_serpentine_series_spec(&name, n_turns, turn_length_m, width_m, height_m);
    let lineage = BlueprintTopologyFactory::lineage_for_spec(&topology);
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

    bp.with_topology_spec(topology).with_lineage(lineage)
}

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
    let name = name.into();
    assert!(n_channels >= 1, "n_channels must be ≥ 1");
    let topology = parallel_microchannel_array_spec(
        &name,
        n_channels,
        channel_length_m,
        channel_width_m,
        channel_height_m,
    );
    let lineage = BlueprintTopologyFactory::lineage_for_spec(&topology);
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    let r_ch = shah_london(
        channel_width_m,
        channel_height_m,
        channel_length_m,
        BLOOD_MU,
    );
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

    bp.with_topology_spec(topology).with_lineage(lineage)
}

/// Selective-routing cascade center-trifurcation separator with a single outlet.
///
/// **Topology** (n_levels = 2 example):
/// ```text
/// inlet → split_lv0 → [center → split_lv1 → [center → throat_in → throat_out → merge_pool]
///                                           → L1/R1 → merge_pool]
///                    → L0/R0 → merge_pool
///                    merge_pool → outlet
/// ```
///
/// At each cascade level, ONLY the center arm is re-trifurcated.
/// All peripheral side arms (L_i / R_i) flow directly into a shared
/// `periph_merge` node, then join the post-venturi center stream at `outlet_merge`
/// before reaching the single `outlet`.
///
/// The deepest center arm passes through a venturi throat (`TherapyZone::CancerTarget`).
/// All side arms carry `TherapyZone::HealthyBypass` — wide, low-shear channels
/// that protect RBCs while cancer cells are concentrated toward the center by
/// Zweifach-Fung junction routing.
///
/// **Single outlet**: all arms converge before the outlet node.
///
/// # Parameters
/// - `n_levels` — cascade depth ∈ \[1, 5\]; 1 = single trifurcation, no re-splitting
/// - `center_frac` — center-arm width fraction (0.25–0.65)
///
/// # Channel names
/// - `"inlet_section"` — inlet trunk
/// - `"center_lv{N}"` — center-arm straight at level N
/// - `"L_lv{N}"`, `"R_lv{N}"` — peripheral side arms at level N
/// - `"throat_section"` — venturi throat on the deepest center arm
/// - `"periph_collect_{N}"` — branches feeding peripheral arms into the merge node
/// - `"center_to_merge"` — post-venturi channel joining outlet_merge
/// - `"trunk_out"` — outlet merge → outlet
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn cascade_center_trifurcation_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch_length_m: f64,
    n_levels: u8,
    main_width_m: f64,
    center_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
    venturi_treatment_enabled: bool,
    center_serpentine: Option<CenterSerpentineSpec>,
) -> NetworkBlueprint {
    let request = SelectiveTreeRequest {
        name: name.into(),
        box_dims_mm: (127.76, 85.47),
        trunk_length_m,
        branch_length_m,
        hybrid_branch_length_m: branch_length_m,
        main_width_m,
        throat_width_m,
        throat_length_m,
        channel_height_m: height_m,
        topology: SelectiveTreeTopology::CascadeCenterTrifurcation {
            n_levels: n_levels.clamp(1, 5) as usize,
            center_frac: center_frac.clamp(0.20, 0.70),
            venturi_treatment_enabled,
            center_serpentine: generator_center_serpentine(center_serpentine),
        },
    };
    create_selective_tree_geometry(&request)
}

/// Selective routing: pre-trifurcation cascade followed
/// by terminal trifurcation + bifurcation skimming, with a single venturi
/// treatment arm and a single external outlet.
///
/// Backward-compatible wrapper over
/// [`incremental_filtration_tri_bi_rect_staged`], using the same center fraction
/// for both pre-trifurcation and terminal-trifurcation stages.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn incremental_filtration_tri_bi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    pretri_branch_length_m: f64,
    hybrid_branch_length_m: f64,
    n_pretri: u8,
    main_width_m: f64,
    center_frac: f64,
    bi_treat_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    incremental_filtration_tri_bi_rect_staged(
        name,
        trunk_length_m,
        pretri_branch_length_m,
        hybrid_branch_length_m,
        n_pretri,
        main_width_m,
        center_frac,
        center_frac,
        bi_treat_frac,
        throat_width_m,
        throat_length_m,
        height_m,
        true,
        None,
    )
}

/// Selective routing with staged trifurcation control:
/// separate pre-trifurcation and terminal-trifurcation center fractions.
///
/// **Topology**:
/// `inlet → pre-tri cascade (center arm only) → hybrid_tri → hybrid_bi`
/// where:
/// - pre-tri side arms route to peripheral bypass merge (`periph_merge`),
/// - hybrid-tri side arms also route to `periph_merge`,
/// - hybrid-bi low-flow arm routes to `periph_merge`,
/// - hybrid-bi treatment arm routes through `"throat_section"` and then merges.
///
/// This encodes "first trifurcations, then trifurcation/bifurcation" with
/// progressive RBC skimming to periphery.
///
/// # Parameters
/// - `n_pretri` — number of initial center-cascade trifurcation levels (1–3)
/// - `pretri_center_frac` — center-arm width fraction for pre-trifurcation levels
/// - `terminal_tri_center_frac` — center-arm width fraction for terminal trifurcation
/// - `bi_treat_frac` — treatment-arm width fraction at terminal bifurcation
/// - `outlet_tail_length_m` — distance from remerge node to outlet. Use short
///   values for "remerge near outlet" layouts.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn incremental_filtration_tri_bi_rect_staged(
    name: impl Into<String>,
    trunk_length_m: f64,
    pretri_branch_length_m: f64,
    hybrid_branch_length_m: f64,
    n_pretri: u8,
    main_width_m: f64,
    pretri_center_frac: f64,
    terminal_tri_center_frac: f64,
    bi_treat_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
    venturi_treatment_enabled: bool,
    center_serpentine: Option<CenterSerpentineSpec>,
) -> NetworkBlueprint {
    incremental_filtration_tri_bi_rect_staged_remerge(
        name,
        trunk_length_m,
        pretri_branch_length_m,
        hybrid_branch_length_m,
        n_pretri,
        main_width_m,
        pretri_center_frac,
        terminal_tri_center_frac,
        bi_treat_frac,
        throat_width_m,
        throat_length_m,
        trunk_length_m,
        height_m,
        venturi_treatment_enabled,
        center_serpentine,
    )
}

/// Selective-routing staged preset with explicit
/// post-remerge outlet-tail control.
///
/// This additive variant allows "remerge right before outlet" geometry by
/// shortening `outlet_tail_length_m` while keeping all upstream staged
/// trifurcation/bifurcation controls unchanged.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn incremental_filtration_tri_bi_rect_staged_remerge(
    name: impl Into<String>,
    trunk_length_m: f64,
    pretri_branch_length_m: f64,
    hybrid_branch_length_m: f64,
    n_pretri: u8,
    main_width_m: f64,
    pretri_center_frac: f64,
    terminal_tri_center_frac: f64,
    bi_treat_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    outlet_tail_length_m: f64,
    height_m: f64,
    venturi_treatment_enabled: bool,
    center_serpentine: Option<CenterSerpentineSpec>,
) -> NetworkBlueprint {
    let request = SelectiveTreeRequest {
        name: name.into(),
        box_dims_mm: (127.76, 85.47),
        trunk_length_m,
        branch_length_m: pretri_branch_length_m,
        hybrid_branch_length_m,
        main_width_m,
        throat_width_m,
        throat_length_m,
        channel_height_m: height_m,
        topology: SelectiveTreeTopology::IncrementalFiltrationTriBi {
            n_pretri: n_pretri.clamp(1, 3) as usize,
            pretri_center_frac: pretri_center_frac.clamp(0.20, 0.70),
            terminal_tri_center_frac: terminal_tri_center_frac.clamp(0.20, 0.70),
            bi_treat_frac: bi_treat_frac.clamp(0.50, 0.85),
            venturi_treatment_enabled,
            center_serpentine: generator_center_serpentine(center_serpentine),
            outlet_tail_length_m: outlet_tail_length_m.max(throat_width_m),
        },
    };
    create_selective_tree_geometry(&request)
}

/// Asymmetric 3-stream trifurcation with center-only venturi for selective SDT.
///
/// **Topology**:
/// `inlet → split_jn → [center→throat_in→venturi→outlet_merge,
///                       left→periph_merge (WBC collection),
///                       right→periph_merge (RBC waste)] → outlet`
///
/// The center arm (cancer-enriched by inertial focusing) receives venturi
/// cavitation treatment while the left (wide, WBC collection) and right
/// (narrow, RBC waste via CFL layer margination) arms bypass treatment.
///
/// # Parameters
/// - `center_frac` — center arm width fraction ∈ [0.20, 0.60]
/// - `left_frac` — left arm (WBC collection) width fraction ∈ [0.15, 0.50]
/// - `right_frac` is derived: `1.0 - center_frac - left_frac`, clamped ≥ 0.05
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn asymmetric_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch_length_m: f64,
    main_width_m: f64,
    center_frac: f64,
    left_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let center_frac = center_frac.clamp(0.20, 0.60);
    let left_frac = left_frac.clamp(0.15, 0.50);
    let right_frac = (1.0 - center_frac - left_frac).clamp(0.05, 0.50);
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (main_width_m + throat_width_m) * 0.5;
    let w_center = main_width_m * center_frac;
    let w_left = main_width_m * left_frac;
    let w_right = main_width_m * right_frac;

    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("split_jn", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_in", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat_out", NodeKind::Junction));
    bp.add_node(NodeSpec::new("periph_merge", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet_merge", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    // Inlet trunk
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
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow))
        .with_metadata(AsymmetricTrifurcationParams {
            center_frac,
            left_frac,
            right_frac,
        }),
    );

    // Center arm → venturi
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "center_arm",
            "split_jn",
            "throat_in",
            branch_length_m,
            w_center,
            height_m,
            shah_london(w_center, height_m, branch_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );

    // Left arm (WBC collection) → periph_merge
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "left_arm",
            "split_jn",
            "periph_merge",
            branch_length_m,
            w_left,
            height_m,
            shah_london(w_left, height_m, branch_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
    );

    // Right arm (RBC waste via CFL margination) → periph_merge
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "right_arm",
            "split_jn",
            "periph_merge",
            branch_length_m,
            w_right,
            height_m,
            shah_london(w_right, height_m, branch_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
    );

    // Venturi throat
    bp.add_channel(with_venturi_geometry_metadata(
        ChannelSpec::new_pipe_rect(
            "throat_section",
            "throat_in",
            "throat_out",
            l_throat,
            throat_width_m,
            height_m,
            shah_london(throat_width_m, height_m, l_throat, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
        venturi_metadata(
            throat_width_m,
            height_m,
            l_throat,
            w_center,
            w_center,
            l_taper,
        ),
    ));

    // Diffuser: throat_out → outlet_merge
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "center_to_merge",
            "throat_out",
            "outlet_merge",
            l_taper,
            w_center,
            height_m,
            shah_london(w_center, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    // Peripheral merge → outlet_merge
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "periph_to_merge",
            "periph_merge",
            "outlet_merge",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
    );

    // Outlet trunk
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "trunk_out",
            "outlet_merge",
            "outlet",
            trunk_length_m,
            main_width_m,
            height_m,
            shah_london(main_width_m, height_m, trunk_length_m, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    specialized_layout::attach_asymmetric_trifurcation_layout(&mut bp);
    bp
}

/// Cascade Tri→Bi→Tri with center-only venturi (selective SDT treatment).
///
/// **Topology**:
/// `inlet → tri_stage1 → [L/R→periph_merge, center→bi_stage2]`
/// `→ [bypass→periph_merge, treat→tri_stage3]`
/// `→ [L/R→periph_merge, center→throat_in→venturi] → outlet`
///
/// Three progressive focusing stages ensure only the most cancer-enriched
/// center stream reaches the venturi cavitation zone, while RBCs and WBCs
/// are progressively skimmed to peripheral bypass channels.
///
/// # Parameters
/// - `tri1_center_frac` — stage-1 trifurcation center fraction ∈ [0.25, 0.65]
/// - `bi_treat_frac` — stage-2 bifurcation treatment arm fraction ∈ [0.50, 0.85]
/// - `tri3_center_frac` — stage-3 trifurcation center fraction ∈ [0.25, 0.65]
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn cascade_tri_bi_tri_selective_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch_length_m: f64,
    main_width_m: f64,
    tri1_center_frac: f64,
    bi_treat_frac: f64,
    tri3_center_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let request = SelectiveTreeRequest {
        name: name.into(),
        box_dims_mm: (127.76, 85.47),
        trunk_length_m,
        branch_length_m,
        hybrid_branch_length_m: branch_length_m,
        main_width_m,
        throat_width_m,
        throat_length_m,
        channel_height_m: height_m,
        topology: SelectiveTreeTopology::TriBiTriSelective {
            first_center_frac: tri1_center_frac.clamp(0.25, 0.65),
            bi_treat_frac: bi_treat_frac.clamp(0.50, 0.85),
            second_center_frac: tri3_center_frac.clamp(0.25, 0.65),
        },
    };
    create_selective_tree_geometry(&request)
}

/// Double-trifurcation selective-routing topology with differential venturi throat counts.
///
/// Implements the SDT Millifluidic Device layout from the
/// `treatment_zone_plate_trifurcation` figure:
///
/// ```text
/// Inlet → 1st split (trifurcation)
///   ├── outer_upper → 2nd split → ch 1–3 (RBC/WBC bypass, 0 throats)
///   ├── center      → 2nd split → ch 4–6 (CTC treatment lanes)
///   │                                  ├── center_L  (center_throat_count throats or acoustic)
///   │                                  ├── center_C  (center_throat_count throats or acoustic)
///   │                                  └── center_R  (center_throat_count throats or acoustic)
///   └── outer_lower → 2nd split → ch 7–9 (RBC/WBC bypass, 0 throats)
///         → all streams → merge → outlet
/// ```
///
/// **Hemolysis partitioning invariant**: bypass channels accumulate shear-only
/// HI (no cavitation contribution) because they bypass the sonication zone.
/// Center channels accumulate both shear HI and cavitation HI from
/// `center_throat_count` serial venturi throats.
///
/// **FDA MI compliance**: each throat on the center channels is tagged with
/// [`ChannelVenturiSpec`] so that `compute_metrics` can enforce
/// `MI_equiv < 1.9` independently per throat.
///
/// # Parameters
///
/// - `split1_center_frac`   — center-arm width fraction at the first trifurcation ∈ [0.25, 0.65].
/// - `split2_center_frac`   — center-arm width fraction at the second trifurcation ∈ [0.25, 0.65].
/// - `center_throat_count`  — serial venturi throats per center channel (0–4).
/// - `throat_width_m`       — throat constriction width [m].
/// - `throat_length_m`      — length of each throat segment [m].
/// - `inter_throat_spacing_m` — re-development length between throats; must be
///                              `> 10 × D_h` (Shah & London 1978, §2-3).
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn double_trifurcation_cif_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch_length_m: f64,
    main_width_m: f64,
    split1_center_frac: f64,
    split2_center_frac: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    inter_throat_spacing_m: f64,
    center_throat_count: u8,
    height_m: f64,
) -> NetworkBlueprint {
    let request = SelectiveTreeRequest {
        name: name.into(),
        box_dims_mm: (127.76, 85.47),
        trunk_length_m,
        branch_length_m,
        hybrid_branch_length_m: branch_length_m,
        main_width_m,
        throat_width_m,
        throat_length_m,
        channel_height_m: height_m,
        topology: SelectiveTreeTopology::DoubleTrifurcationCif {
            split1_center_frac: split1_center_frac.clamp(0.25, 0.65),
            split2_center_frac: split2_center_frac.clamp(0.25, 0.65),
            center_throat_count: center_throat_count.min(4),
            inter_throat_spacing_m,
        },
    };
    create_selective_tree_geometry(&request)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::CrossSectionSpec;
    use crate::geometry::metadata::IncrementalFiltrationParams;

    fn rect_width(bp: &NetworkBlueprint, id: &str) -> f64 {
        let ch = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == id)
            .expect("channel must exist");
        match ch.cross_section {
            CrossSectionSpec::Rectangular { width_m, .. } => width_m,
            CrossSectionSpec::Circular { .. } => panic!("expected rectangular cross-section"),
        }
    }

    fn channel_length(bp: &NetworkBlueprint, id: &str) -> f64 {
        bp.channels
            .iter()
            .find(|c| c.id.as_str() == id)
            .map_or_else(|| panic!("channel {id} must exist"), |c| c.length_m)
    }

    #[test]
    fn staged_cif_terminal_fraction_controls_hybrid_center_width() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-staged",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.33,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            true,
            None,
        );

        let center_lv1 = rect_width(&bp, "center_lv1");
        let hy_tri_center = rect_width(&bp, "hy_tri_center");
        let expected = center_lv1 * 0.55;

        assert!(
            (hy_tri_center - expected).abs() < 1e-12,
            "terminal staged fraction must set hy_tri_center width"
        );
    }

    #[test]
    fn staged_cif_metadata_encodes_all_stage_parameters() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-meta", 12e-3, 8e-3, 6e-3, 3, 2.0e-3, 0.45, 0.55, 0.76, 100e-6, 300e-6, 1.0e-3,
            true, None,
        );
        let inlet = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "inlet_section")
            .expect("inlet_section must exist");
        let params = inlet
            .metadata
            .as_ref()
            .and_then(
                crate::geometry::metadata::MetadataContainer::get::<IncrementalFiltrationParams>,
            )
            .expect("IncrementalFiltrationParams metadata must exist");

        assert_eq!(params.n_pretri, 3);
        assert!((params.pretri_center_frac - 0.45).abs() < 1e-12);
        assert!((params.terminal_tri_center_frac - 0.55).abs() < 1e-12);
        assert!((params.bi_treat_frac - 0.76).abs() < 1e-12);
        assert!((params.outlet_tail_length_m - 12e-3).abs() < 1e-12);
    }

    #[test]
    fn staged_cif_remerge_tail_controls_trunk_out_length() {
        let bp = incremental_filtration_tri_bi_rect_staged_remerge(
            "cif-remerge",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.45,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.5e-3,
            1.0e-3,
            true,
            None,
        );

        assert!((channel_length(&bp, "periph_to_merge") - 12e-3).abs() < 1e-12);
        assert!((channel_length(&bp, "trunk_out") - 1.5e-3).abs() < 1e-12);
    }

    #[test]
    fn staged_cif_pretri_levels_ramp_center_bias() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-stage-ramp",
            12e-3,
            8e-3,
            6e-3,
            3,
            2.0e-3,
            0.45,
            0.60,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            true,
            None,
        );

        let c0 = rect_width(&bp, "center_lv0");
        let c1 = rect_width(&bp, "center_lv1");
        let c2 = rect_width(&bp, "center_lv2");
        let l0 = rect_width(&bp, "L_lv0");
        let l1 = rect_width(&bp, "L_lv1");
        let l2 = rect_width(&bp, "L_lv2");

        let f0 = c0 / (c0 + 2.0 * l0);
        let f1 = c1 / (c1 + 2.0 * l1);
        let f2 = c2 / (c2 + 2.0 * l2);

        assert!(f1 >= f0 - 1e-12);
        assert!(f2 >= f1 - 1e-12);
        assert!(f2 > f0 + 1e-6);
    }

    #[test]
    fn staged_cif_acoustic_serpentine_omits_throat_and_keeps_bypass_straight() {
        let serpentine = CenterSerpentineSpec {
            segments: 5,
            bend_radius_m: 3.5e-3,
            segment_length_m: 8.0e-3,
        };
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-acoustic-serp",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.45,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            false,
            Some(serpentine),
        );

        assert!(
            bp.channels
                .iter()
                .all(|c| c.id.as_str() != "throat_section"),
            "acoustic-only selective-routing layout must not emit a venturi throat section"
        );
        assert!(
            bp.channels
                .iter()
                .any(|c| c.id.as_str() == "treatment_section"),
            "acoustic-only selective-routing layout must emit a center treatment section"
        );

        let center_channel = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "center_lv0")
            .expect("center_lv0 must exist");
        assert!(matches!(
            center_channel.channel_shape,
            ChannelShape::Serpentine { .. }
        ));

        let bypass_channel = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "L_lv0")
            .expect("L_lv0 must exist");
        assert_eq!(bypass_channel.channel_shape, ChannelShape::Straight);
    }

    #[test]
    fn constriction_expansion_array_embeds_series_topology_metadata() {
        let bp =
            constriction_expansion_array_rect("ce-meta", 4, 3.0e-3, 1.5e-3, 200e-6, 80e-6, 60e-6);
        let topology = bp
            .topology_spec()
            .expect("constriction_expansion_array_rect should attach topology metadata");

        assert_eq!(topology.treatment_channel_ids().len(), 4);
        assert!(topology.is_leukapheresis_topology());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }

    #[test]
    fn spiral_channel_embeds_series_topology_metadata() {
        let bp = spiral_channel_rect("spiral-meta", 5, 6.0e-3, 200e-6, 60e-6);
        let topology = bp
            .topology_spec()
            .expect("spiral_channel_rect should attach topology metadata");

        assert_eq!(topology.treatment_channel_ids().len(), 5);
        assert!(topology.has_serpentine());
        assert!(topology.is_leukapheresis_topology());
    }

    #[test]
    fn parallel_microchannel_array_embeds_parallel_topology_metadata() {
        let bp = parallel_microchannel_array_rect("parallel-meta", 8, 12.0e-3, 140e-6, 55e-6);
        let topology = bp
            .topology_spec()
            .expect("parallel_microchannel_array_rect should attach topology metadata");

        assert!(topology.has_parallel_paths());
        assert_eq!(topology.parallel_channels.len(), 8);
        assert_eq!(topology.treatment_channel_ids().len(), 8);
        assert!(topology.is_leukapheresis_topology());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }
}
