//! Incremental filtration (CIF) staged selective-routing functions.
//!
//! Implements the staged filtration model where cells pass through:
//! 1. Pre-trifurcation center-cascade skimming stages,
//! 2. One terminal trifurcation skimming stage,
//! 3. One terminal asymmetric bifurcation selecting the treatment arm.

use super::cascade_routing::{tri_center_q_frac, tri_center_q_frac_cross_junction};
use super::routing_probability::{p_center, p_treat_bifurcation, SE_CANCER, SE_RBC, SE_WBC};
use super::IncrementalFiltrationResult;
use cfd_core::error::{Error, Result};

const PRETRI_CENTER_FRAC_MIN: f64 = 0.20;
const PRETRI_CENTER_FRAC_MAX: f64 = 0.70;
const TERMINAL_BI_TREAT_FRAC_MIN: f64 = 0.50;
const TERMINAL_BI_TREAT_FRAC_MAX: f64 = 0.85;

fn validate_pretri_count(n_pretri: u8) -> Result<usize> {
    if !(1..=3).contains(&n_pretri) {
        return Err(Error::InvalidConfiguration(
            "Incremental filtration pre-trifurcation stage count must lie in [1, 3]".to_string(),
        ));
    }

    Ok(n_pretri as usize)
}

fn validate_center_fraction(name: &str, value: f64) -> Result<f64> {
    if !value.is_finite() || !(PRETRI_CENTER_FRAC_MIN..=PRETRI_CENTER_FRAC_MAX).contains(&value) {
        return Err(Error::InvalidConfiguration(format!(
            "Incremental filtration {name} must lie in [{PRETRI_CENTER_FRAC_MIN}, {PRETRI_CENTER_FRAC_MAX}]"
        )));
    }

    Ok(value)
}

fn validate_terminal_bifurcation_fraction(value: f64) -> Result<f64> {
    if !value.is_finite()
        || !(TERMINAL_BI_TREAT_FRAC_MIN..=TERMINAL_BI_TREAT_FRAC_MAX).contains(&value)
    {
        return Err(Error::InvalidConfiguration(format!(
            "Incremental filtration terminal bifurcation treatment fraction must lie in [{TERMINAL_BI_TREAT_FRAC_MIN}, {TERMINAL_BI_TREAT_FRAC_MAX}]"
        )));
    }

    Ok(value)
}

fn validate_unit_open_fraction(name: &str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 || value >= 1.0 {
        return Err(Error::InvalidConfiguration(format!(
            "Incremental filtration {name} must lie in the open interval (0, 1)"
        )));
    }

    Ok(value)
}

fn validate_positive_geometry(name: &str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "Incremental filtration {name} must be finite and positive"
        )));
    }

    Ok(value)
}

fn cif_pretri_stage_center_fracs_impl(n: usize, start: f64, terminal: f64) -> Vec<f64> {
    if n == 1 {
        return vec![start];
    }

    let target = start.max(terminal);

    (0..n)
        .map(|i| {
            let t = ((i + 1) as f64 / n as f64).powi(2);
            start + (target - start) * t
        })
        .collect()
}

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
    checked_cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)
        .unwrap_or_else(|_| {
            cif_pretri_stage_center_fracs_impl(
                n_pretri.clamp(1, 3) as usize,
                pretri_center_frac.clamp(PRETRI_CENTER_FRAC_MIN, PRETRI_CENTER_FRAC_MAX),
                terminal_tri_frac.clamp(PRETRI_CENTER_FRAC_MIN, PRETRI_CENTER_FRAC_MAX),
            )
        })
}

/// Checked stage-wise pre-trifurcation center-arm width fractions.
pub fn checked_cif_pretri_stage_center_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Result<Vec<f64>> {
    let n = validate_pretri_count(n_pretri)?;
    let start = validate_center_fraction("pre-trifurcation center fraction", pretri_center_frac)?;
    let terminal =
        validate_center_fraction("terminal trifurcation center fraction", terminal_tri_frac)?;

    Ok(cif_pretri_stage_center_fracs_impl(n, start, terminal))
}

/// Stage-wise pre-trifurcation center-flow fractions for selective routing.
///
/// Each value is `tri_center_q_frac(stage_center_frac)` for the corresponding
/// stage returned by [`cif_pretri_stage_center_fracs`].
///
/// This is a legacy width-fraction-only surrogate. When the parent width and
/// channel height are known, prefer [`cif_pretri_stage_q_fracs_cross_junction`]
/// so the staged split follows the geometry-aware rectangular conductance model.
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

/// Checked stage-wise pre-trifurcation center-flow fractions.
pub fn checked_cif_pretri_stage_q_fracs(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
) -> Result<Vec<f64>> {
    Ok(
        checked_cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)?
            .into_iter()
            .map(tri_center_q_frac)
            .collect(),
    )
}

/// Geometry-aware stage-wise pre-trifurcation center-flow fractions.
///
/// Each pre-trifurcation stage uses [`tri_center_q_frac_cross_junction`] with
/// the current parent width and shared channel height, then narrows the parent
/// width by the authored center fraction before evaluating the next stage.
pub fn cif_pretri_stage_q_fracs_cross_junction(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> Vec<f64> {
    checked_cif_pretri_stage_q_fracs_cross_junction(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        parent_width_m,
        channel_height_m,
    )
    .unwrap_or_else(|_| {
        let stage_fracs =
            cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
        let mut current_parent_w = parent_width_m.max(1.0e-9);
        stage_fracs
            .into_iter()
            .map(|stage_center_frac| {
                let q_frac = tri_center_q_frac_cross_junction(
                    stage_center_frac,
                    current_parent_w,
                    channel_height_m.max(1.0e-9),
                );
                current_parent_w *= stage_center_frac;
                q_frac
            })
            .collect()
    })
}

/// Checked geometry-aware stage-wise pre-trifurcation center-flow fractions.
pub fn checked_cif_pretri_stage_q_fracs_cross_junction(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> Result<Vec<f64>> {
    let stage_fracs =
        checked_cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)?;
    let mut current_parent_w = validate_positive_geometry("parent width", parent_width_m)?;
    let channel_height_m = validate_positive_geometry("channel height", channel_height_m)?;

    Ok(stage_fracs
        .into_iter()
        .map(|stage_center_frac| {
            let q_frac = tri_center_q_frac_cross_junction(
                stage_center_frac,
                current_parent_w,
                channel_height_m,
            );
            current_parent_w *= stage_center_frac;
            q_frac
        })
        .collect())
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
    checked_incremental_filtration_separation_staged(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        terminal_bi_treat_frac,
    )
    .unwrap_or_else(|_| {
        let pretri_q_fracs =
            cif_pretri_stage_q_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
        let q_tri = tri_center_q_frac(terminal_tri_frac);
        incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
    })
}

/// Checked staged selective-routing separation.
pub fn checked_incremental_filtration_separation_staged(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
) -> Result<IncrementalFiltrationResult> {
    validate_terminal_bifurcation_fraction(terminal_bi_treat_frac)?;
    let pretri_q_fracs =
        checked_cif_pretri_stage_q_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)?;
    let q_tri = tri_center_q_frac(validate_center_fraction(
        "terminal trifurcation center fraction",
        terminal_tri_frac,
    )?);
    Ok(incremental_from_q_fractions(
        &pretri_q_fracs,
        q_tri,
        terminal_bi_treat_frac,
    ))
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
    checked_incremental_filtration_separation_from_qfracs(
        pretri_q_center_fracs,
        terminal_tri_q_center_frac,
        terminal_bi_treat_frac,
    )
    .unwrap_or_else(|_| {
        incremental_from_q_fractions(
            pretri_q_center_fracs,
            terminal_tri_q_center_frac,
            terminal_bi_treat_frac,
        )
    })
}

/// Checked staged selective-routing separation from solved per-stage flow fractions.
pub fn checked_incremental_filtration_separation_from_qfracs(
    pretri_q_center_fracs: &[f64],
    terminal_tri_q_center_frac: f64,
    terminal_bi_treat_frac: f64,
) -> Result<IncrementalFiltrationResult> {
    for &q_frac in pretri_q_center_fracs {
        validate_unit_open_fraction("pre-trifurcation center flow fraction", q_frac)?;
    }
    validate_unit_open_fraction(
        "terminal trifurcation center flow fraction",
        terminal_tri_q_center_frac,
    )?;
    validate_terminal_bifurcation_fraction(terminal_bi_treat_frac)?;

    Ok(incremental_from_q_fractions(
        pretri_q_center_fracs,
        terminal_tri_q_center_frac,
        terminal_bi_treat_frac,
    ))
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
    checked_incremental_filtration_separation_cross_junction(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        terminal_bi_treat_frac,
        parent_width_m,
        channel_height_m,
    )
    .unwrap_or_else(|_| {
        let pretri_q_fracs = cif_pretri_stage_q_fracs_cross_junction(
            n_pretri,
            pretri_center_frac,
            terminal_tri_frac,
            parent_width_m,
            channel_height_m,
        );

        let stage_fracs =
            cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac);
        let current_parent_w = stage_fracs
            .iter()
            .fold(parent_width_m, |width_m, stage_center_frac| {
                width_m * stage_center_frac
            });

        let q_tri =
            tri_center_q_frac_cross_junction(terminal_tri_frac, current_parent_w, channel_height_m);
        incremental_from_q_fractions(&pretri_q_fracs, q_tri, terminal_bi_treat_frac)
    })
}

/// Checked staged selective-routing separation with cross-junction diameter effects.
pub fn checked_incremental_filtration_separation_cross_junction(
    n_pretri: u8,
    pretri_center_frac: f64,
    terminal_tri_frac: f64,
    terminal_bi_treat_frac: f64,
    parent_width_m: f64,
    channel_height_m: f64,
) -> Result<IncrementalFiltrationResult> {
    validate_positive_geometry("parent width", parent_width_m)?;
    validate_positive_geometry("channel height", channel_height_m)?;
    validate_terminal_bifurcation_fraction(terminal_bi_treat_frac)?;
    let stage_fracs =
        checked_cif_pretri_stage_center_fracs(n_pretri, pretri_center_frac, terminal_tri_frac)?;
    let pretri_q_fracs = checked_cif_pretri_stage_q_fracs_cross_junction(
        n_pretri,
        pretri_center_frac,
        terminal_tri_frac,
        parent_width_m,
        channel_height_m,
    )?;

    let current_parent_w = stage_fracs
        .iter()
        .fold(parent_width_m, |width_m, stage_center_frac| {
            width_m * stage_center_frac
        });

    let q_tri =
        tri_center_q_frac_cross_junction(terminal_tri_frac, current_parent_w, channel_height_m);
    Ok(incremental_from_q_fractions(
        &pretri_q_fracs,
        q_tri,
        terminal_bi_treat_frac,
    ))
}
