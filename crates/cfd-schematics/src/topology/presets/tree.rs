//! Split-tree topology presets (asymmetric / symmetric bifurcation / trifurcation).

use crate::domain::therapy_metadata::TherapyZone;

use super::super::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec, SplitKind, SplitStageSpec,
    TreatmentActuationMode,
};
use super::helpers::{PLATE_HEIGHT_MM, PLATE_WIDTH_MM};

/// Create an N-level asymmetric split tree spec (Option 1 seed).
///
/// Each stage independently specifies its split kind (Bi or Tri) and
/// per-branch width allocation.  The center branch carries the treatment
/// path; peripheral branches carry bypass flow for RBC/WBC separation.
///
/// # Parameters
///
/// - `name`: unique design name
/// - `split_kinds`: sequence of [`SplitKind`] for each level (1–4 stages)
/// - `branch_widths_m`: per-stage Vec of `(width_m, role, is_treatment)` per branch
/// - `channel_height_m`: uniform channel depth
/// - `trunk_length_m`: inlet trunk channel length
///
/// # Panics
///
/// Panics if `split_kinds.len() != branch_widths_m.len()` or any branch count
/// mismatches its `SplitKind`.
#[must_use]
pub fn asymmetric_split_tree_spec(
    name: &str,
    split_kinds: &[SplitKind],
    branch_widths_m: &[Vec<(f64, BranchRole, bool)>],
    channel_height_m: f64,
    trunk_length_m: f64,
) -> BlueprintTopologySpec {
    assert_eq!(
        split_kinds.len(),
        branch_widths_m.len(),
        "split_kinds and branch_widths_m must match in length"
    );

    let stages: Vec<SplitStageSpec> = split_kinds
        .iter()
        .zip(branch_widths_m.iter())
        .enumerate()
        .map(|(idx, (&kind, widths))| {
            assert_eq!(
                widths.len(),
                kind.branch_count(),
                "stage {idx} expects {} branches, got {}",
                kind.branch_count(),
                widths.len()
            );

            let branches = widths
                .iter()
                .enumerate()
                .map(|(bi, &(w, role, treatment))| {
                    let label = match (kind, bi) {
                        (SplitKind::Trifurcation, 0) => "center".to_string(),
                        (SplitKind::Trifurcation, 1) => "left".to_string(),
                        (SplitKind::Trifurcation, 2) => "right".to_string(),
                        (SplitKind::Bifurcation, 0) => "upper".to_string(),
                        (SplitKind::Bifurcation, _) => "lower".to_string(),
                        _ => format!("arm_{bi}"),
                    };

                    let therapy_zone = if treatment {
                        TherapyZone::CancerTarget
                    } else {
                        TherapyZone::HealthyBypass
                    };

                    BranchSpec {
                        label,
                        role,
                        treatment_path: treatment,
                        route: ChannelRouteSpec {
                            length_m: trunk_length_m * 0.8,
                            width_m: w,
                            height_m: channel_height_m,
                            serpentine: None,
                            therapy_zone,
                        },
                    }
                })
                .collect();

            SplitStageSpec {
                stage_id: format!("stage_{idx}"),
                split_kind: kind,
                branches,
            }
        })
        .collect();

    BlueprintTopologySpec {
        topology_id: format!("{name}_topology"),
        design_name: name.to_string(),
        box_dims_mm: (PLATE_WIDTH_MM, PLATE_HEIGHT_MM),
        inlet_width_m: branch_widths_m[0].iter().map(|(w, _, _)| *w).sum::<f64>(),
        outlet_width_m: branch_widths_m[0].iter().map(|(w, _, _)| *w).sum::<f64>(),
        trunk_length_m,
        outlet_tail_length_m: trunk_length_m,
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        split_stages: stages,
        venturi_placements: Vec::new(),
        treatment_mode: TreatmentActuationMode::UltrasoundOnly,
    }
}

/// Symmetric bifurcation spec: equal-width branches at every level.
///
/// Convenience wrapper over [`asymmetric_split_tree_spec`] with
/// `center_frac = 0.5` at every bifurcation (equal halves).
#[must_use]
pub fn symmetric_bifurcation_spec(
    name: &str,
    n_levels: usize,
    inlet_width_m: f64,
    channel_height_m: f64,
    trunk_length_m: f64,
) -> BlueprintTopologySpec {
    let kinds = vec![SplitKind::Bifurcation; n_levels];
    let mut widths = Vec::with_capacity(n_levels);
    let mut parent_width = inlet_width_m;
    for _ in 0..n_levels {
        let half = parent_width * 0.5;
        widths.push(vec![
            (half, BranchRole::Treatment, true),
            (half, BranchRole::Neutral, false),
        ]);
        parent_width = half;
    }
    asymmetric_split_tree_spec(name, &kinds, &widths, channel_height_m, trunk_length_m)
}

/// Symmetric trifurcation spec: per-arm fractions from center_frac.
///
/// Center arm gets `center_frac × parent_width`, each peripheral arm
/// gets `(1 − center_frac) / 2 × parent_width`.
#[must_use]
pub fn asymmetric_trifurcation_spec(
    name: &str,
    n_levels: usize,
    inlet_width_m: f64,
    center_frac: f64,
    channel_height_m: f64,
    trunk_length_m: f64,
) -> BlueprintTopologySpec {
    let kinds = vec![SplitKind::Trifurcation; n_levels];
    let periph_frac = (1.0 - center_frac) / 2.0;
    let mut widths = Vec::with_capacity(n_levels);
    let mut parent_width = inlet_width_m;
    for _ in 0..n_levels {
        widths.push(vec![
            (parent_width * center_frac, BranchRole::Treatment, true),
            (parent_width * periph_frac, BranchRole::WbcCollection, false),
            (parent_width * periph_frac, BranchRole::RbcBypass, false),
        ]);
        parent_width *= center_frac;
    }
    asymmetric_split_tree_spec(name, &kinds, &widths, channel_height_m, trunk_length_m)
}
