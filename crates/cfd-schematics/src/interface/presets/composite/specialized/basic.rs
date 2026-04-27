//! Basic composite preset factory functions.
use super::super::super::finalize_preset_blueprint;
use super::parallel_lane::{
    canonical_parallel_blueprint, canonical_parallel_venturi_blueprint,
    generator_center_serpentine, parallel_lane, CenterSerpentineSpec,
};
use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::generator::{
    create_primitive_selective_tree_geometry, PrimitiveSelectiveSplitKind,
    PrimitiveSelectiveTreeRequest,
};
use crate::topology::presets::{
    constriction_expansion_series_spec, parallel_microchannel_array_spec,
    spiral_serpentine_series_spec,
};
use crate::topology::{SerpentineSpec, TreatmentActuationMode};
use crate::BlueprintTopologyFactory;

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
    let l_taper = 2.5 * (main_width_m + throat_width_m);
    let name = name.into();
    let channels = vec![
        parallel_lane(
            "throat_section",
            2.0 * l_taper + l_throat,
            main_width_m,
            height_m,
            TherapyZone::CancerTarget,
            None,
        ),
        parallel_lane(
            "peripheral_1",
            bypass_length_m,
            main_width_m * 0.5,
            height_m,
            TherapyZone::HealthyBypass,
            None,
        ),
        parallel_lane(
            "peripheral_2",
            bypass_length_m,
            main_width_m * 0.5,
            height_m,
            TherapyZone::HealthyBypass,
            None,
        ),
    ];
    canonical_parallel_venturi_blueprint(
        name,
        main_width_m,
        approach_length_m,
        bypass_length_m,
        throat_width_m,
        height_m,
        l_throat,
        vec!["throat_section".to_string()],
        channels,
    )
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
    let name = name.into();
    let serpentine = Some(SerpentineSpec {
        segments: segments.max(2),
        bend_radius_m: wide_width_m * 0.5,
        segment_length_m,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    });
    canonical_parallel_blueprint(
        name,
        wide_width_m,
        trunk_length_m,
        trunk_length_m,
        vec![
            parallel_lane(
                "wide_arm",
                segment_length_m * segments.max(1) as f64,
                wide_width_m,
                height_m,
                TherapyZone::CancerTarget,
                serpentine.clone(),
            ),
            parallel_lane(
                "narrow_arm",
                segment_length_m * segments.max(1) as f64,
                narrow_width_m,
                height_m,
                TherapyZone::HealthyBypass,
                Some(SerpentineSpec {
                    segments: segments.max(2),
                    bend_radius_m: narrow_width_m * 0.5,
                    segment_length_m,
                    wave_type: crate::topology::SerpentineWaveType::Sine,
                }),
            ),
        ],
        TreatmentActuationMode::UltrasoundOnly,
    )
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
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology)
            .expect("valid constriction_expansion_array topology spec"),
    )
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
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid spiral_channel topology spec"),
    )
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
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology)
            .expect("valid parallel_microchannel_array topology spec"),
    )
}
