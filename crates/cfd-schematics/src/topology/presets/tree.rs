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
                        (SplitKind::NFurcation(2), 0) => "upper".to_string(),
                        (SplitKind::NFurcation(2), 1) => "lower".to_string(),
                        (SplitKind::NFurcation(3), 0) => "left".to_string(),
                        (SplitKind::NFurcation(3), 1) => "center".to_string(),
                        (SplitKind::NFurcation(3), 2) => "right".to_string(),
                        (SplitKind::NFurcation(n), idx) if n % 2 != 0 && idx == n / 2 => {
                            "center".to_string()
                        }
                        (_, _) => format!("arm_{bi}"),
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
                        recovery_sub_split: None,
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

/// Symmetric N-furcation spec: split `parent_width` based on whether N is even or odd,
/// creating `splits` branches. Center branch (if N is odd) or first branch (if N even) is the treatment path.
/// All branches divide the total parent flow equally.
#[must_use]
pub fn symmetric_n_furcation_spec(
    name: &str,
    n_levels: usize,
    splits: usize,
    inlet_width_m: f64,
    channel_height_m: f64,
    trunk_length_m: f64,
) -> BlueprintTopologySpec {
    let kinds = vec![SplitKind::NFurcation(splits); n_levels];
    let frac = 1.0 / (splits as f64);
    let mut widths = Vec::with_capacity(n_levels);
    let mut parent_width = inlet_width_m;
    for _ in 0..n_levels {
        let mut level_branches = Vec::with_capacity(splits);
        let branch_width = parent_width * frac;

        for i in 0..splits {
            let (role, treatment) =
                if (splits % 2 != 0 && i == splits / 2) || (splits % 2 == 0 && i == 0) {
                    (BranchRole::Treatment, true)
                } else {
                    (BranchRole::Neutral, false)
                };
            level_branches.push((branch_width, role, treatment));
        }

        widths.push(level_branches);
        parent_width *= frac;
    }
    asymmetric_split_tree_spec(name, &kinds, &widths, channel_height_m, trunk_length_m)
}
