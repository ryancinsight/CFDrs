//! Incremental filtration (CIF) staged selective-routing functions.
//!
//! Implements the staged filtration model where cells pass through:
//! 1. Pre-trifurcation center-cascade skimming stages,
//! 2. One terminal trifurcation skimming stage,
//! 3. One terminal asymmetric bifurcation selecting the treatment arm.

use super::cascade_routing::{tri_center_q_frac, tri_center_q_frac_cross_junction};
use super::routing_probability::{p_center, p_treat_bifurcation, SE_CANCER, SE_RBC, SE_WBC};
use super::IncrementalFiltrationResult;

/// Internal incremental filtration computation from pre-solved flow fractions.
pub(crate) fn incremental_from_q_fractions(
    pretri_q_center_fracs: &[f64],
    terminal_tri_q_center_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    let mut f_cancer = 1.0_f64;
    let mut f_wbc = 1.0_f64;
    let mut f_rbc = 1.0_f64;
    let mut q_total_treat = 1.0_f64;

    for &q_center in pretri_q_center_fracs {
        let q = q_center.clamp(1.0e-9, 1.0 - 1.0e-9);
        f_cancer *= p_center(q, SE_CANCER);
        f_wbc *= p_center(q, SE_WBC);
        f_rbc *= p_center(q, SE_RBC);
        q_total_treat *= q;
    }

    let q_tri = terminal_tri_q_center_frac.clamp(1.0e-9, 1.0 - 1.0e-9);
    f_cancer *= p_center(q_tri, SE_CANCER);
    f_wbc *= p_center(q_tri, SE_WBC);
    f_rbc *= p_center(q_tri, SE_RBC);
    q_total_treat *= q_tri;

    let q_bi = terminal_bi_treat_frac.clamp(0.50, 0.85);
    f_cancer *= p_treat_bifurcation(q_bi, SE_CANCER);
    f_wbc *= p_treat_bifurcation(q_bi, SE_WBC);
    f_rbc *= p_treat_bifurcation(q_bi, SE_RBC);
    q_total_treat *= q_bi;

    let rbc_periph = (1.0 - f_rbc).clamp(0.0, 1.0);
    let sep_eff = (f_cancer - f_rbc).abs().clamp(0.0, 1.0);
    let center_hematocrit_ratio = if q_total_treat > 1.0e-12 {
        (f_rbc / q_total_treat).clamp(0.0, 2.0)
    } else {
        1.0
    };

    IncrementalFiltrationResult {
        cancer_center_fraction: f_cancer.clamp(0.0, 1.0),
        wbc_center_fraction: f_wbc.clamp(0.0, 1.0),
        rbc_center_fraction: f_rbc.clamp(0.0, 1.0),
        rbc_peripheral_fraction: rbc_periph,
        separation_efficiency: sep_eff,
        center_hematocrit_ratio,
    }
}

/// Stage-wise pre-trifurcation center-arm width fractions for selective routing.
///
/// This models progressive asymmetric focusing across pre-trifurcation levels:
/// deeper levels are allowed to bias further toward the center arm so that
/// large, stiff cells remain in the treatment path while deformable RBCs are
/// progressively skimmed to peripheral bypass branches.
///
/// The returned fractions are:
/// - clamped to `[0.20, 0.70]`,
/// - monotone non-decreasing over stage index,
/// - equal to `pretri_center_frac` when `n_pretri == 1`.
///   Constructed via a monotone interpolation from `start` to
///   `max(start, terminal_tri_frac)` with final clamping.
pub fn cif_pretri_stage_center_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Vec<f64> {
    let n = n_pretri.clamp(1, 3) as usize;
    let start = pretri_center_frac.clamp(0.20, 0.70);
    if n == 1 {
        return vec![start];
    }

    let terminal = terminal_tri_frac.clamp(0.20, 0.70);
    let target = start.max(terminal);

    (0..n)
        .map(|i| {
            let t = ((i + 1) as f64 / n as f64).powi(2);
            (start + (target - start) * t).clamp(0.20, 0.70)
        })
        .collect()
}

/// Stage-wise pre-trifurcation center-flow fractions for selective routing.
///
/// Each value is `tri_center_q_frac(stage_center_frac)` for the corresponding
/// stage returned by [`cif_pretri_stage_center_fracs`].
pub fn cif_pretri_stage_q_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Vec<f64> {
    cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)
        .into_iter()
        .map(tri_center_q_frac)
        .collect()
}

/// Compute staged selective-routing separation.
///
/// This is the canonical staged selective-routing model for:
/// - pre-trifurcation skimming (`pretri_center_frac`),
/// - terminal-trifurcation skimming (`terminal_tri_frac`),
/// - terminal treatment bifurcation (`terminal_bi_treat_frac`).
pub fn incremental_filtration_separation_staged(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    let pretri_q_fracs = cif_pretri_stage_q_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
    let q_tri = tri_center_q_frac(terminal_tri_frac);
    incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
}

/// Compute staged selective-routing separation from solved per-stage flow fractions.
///
/// This additive API is intended for coupling with solved-network extraction
/// where each pre-trifurcation stage can carry a different center-flow fraction.
pub fn incremental_filtration_separation_from_qfracs(
    pretri_q_center_fracs: &[f64],
    terminal_tri_q_center_frac: f64,
    terminal_bi_treat_frac: f64,
) -> IncrementalFiltrationResult {
    incremental_from_q_fractions(
        pretri_q_center_fracs,
        terminal_tri_q_center_frac,
        terminal_bi_treat_frac,
    )
}

/// Compute staged selective-routing separation with cross-junction diameter effects.
///
/// Identical to [`incremental_filtration_separation_staged`] but uses
/// [`tri_center_q_frac_cross_junction`] to account for cross-junction
/// minor losses at each trifurcation.  When channels of different widths
/// intersect, the junction K-factor correction steers additional flow to
/// the wider center arm, enhancing RBC peripheral skimming.
///
/// # Arguments
/// * `n_pretri`              — number of initial center-cascade trifurcations (1–3).
/// * `pretri_center_frac`    — center-arm width fraction for the pre-cascade trifurcations.
/// * `terminal_tri_frac`     — center-arm width fraction for the terminal trifurcation.
/// * `terminal_bi_treat_frac`— flow fraction into the treatment arm of the terminal bifurcation.
/// * `parent_width_m`        — parent trunk channel width [m].
/// * `channel_height_m`      — channel height [m].
pub fn incremental_filtration_separation_cross_junction(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> IncrementalFiltrationResult {
    let stage_fracs =
        cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);

    // Convert width fractions to flow fractions using cross-junction model.
    // Each stage re-splits from the current center width, which narrows by
    // center_frac at each level.
    let mut current_parent_w = parent_width_m;
    let pretri_q_fracs: Vec<f64> = stage_fracs
        .iter()
        .map(|&frac| {
            let q = tri_center_q_frac_cross_junction(frac, current_parent_w, channel_height_m);
            current_parent_w *= frac; // center arm becomes parent for next level
            q
        })
        .collect();

    let q_tri =
        tri_center_q_frac_cross_junction(terminal_tri_frac, current_parent_w, channel_height_m);
    incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
}
