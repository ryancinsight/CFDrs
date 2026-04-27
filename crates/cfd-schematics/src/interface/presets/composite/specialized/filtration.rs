//! Cascade and incremental filtration composite preset functions.
use super::parallel_lane::{generator_center_serpentine, CenterSerpentineSpec};
use crate::domain::model::NetworkBlueprint;
use crate::geometry::generator::{
    create_selective_tree_geometry, SelectiveTreeRequest, SelectiveTreeTopology,
};

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
