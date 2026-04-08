//! Asymmetric trifurcation and selective tree composite preset functions.
use super::parallel_lane::{canonical_parallel_venturi_blueprint, parallel_lane};
use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::generator::{create_selective_tree_geometry, SelectiveTreeRequest, SelectiveTreeTopology};

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
    let w_center = main_width_m * center_frac;
    let w_left = main_width_m * left_frac;
    let w_right = main_width_m * right_frac;
    let name = name.into();
    canonical_parallel_venturi_blueprint(
        name,
        main_width_m,
        trunk_length_m,
        trunk_length_m,
        throat_width_m,
        height_m,
        l_throat,
        vec!["center_arm".to_string()],
        vec![
            parallel_lane(
                "left_arm",
                branch_length_m,
                w_left,
                height_m,
                TherapyZone::HealthyBypass,
                None,
            ),
            parallel_lane(
                "center_arm",
                branch_length_m,
                w_center,
                height_m,
                TherapyZone::CancerTarget,
                None,
            ),
            parallel_lane(
                "right_arm",
                branch_length_m,
                w_right,
                height_m,
                TherapyZone::HealthyBypass,
                None,
            ),
        ],
    )
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
///   `> 10 * D_h` (Shah & London 1978, §2-3).
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

