//! Pure math helpers for report-metrics computation.

use crate::constraints::ACOUSTIC_HALF_WAVELENGTH_M;
use cfd_1d::cascade_treatment_flow_fractions;
use cfd_schematics::topology::SplitStageSpec;

pub(super) fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

pub(super) fn coefficient_of_variation(values: &[f64], mean: f64) -> f64 {
    if values.len() <= 1 || mean <= 1.0e-18 {
        return 0.0;
    }
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / values.len() as f64;
    variance.sqrt() / mean
}

/// Compute the p-th percentile by sorting the mutable slice in-place,
/// avoiding a `.to_vec()` heap clone.  The caller's slice is left sorted.
pub(super) fn percentile_mut(values: &mut [f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_unstable_by(|a, b| a.total_cmp(b));
    let index = ((values.len() - 1) as f64 * p.clamp(0.0, 1.0)).round() as usize;
    values[index]
}

pub(super) fn direct_linear_risk(value: f64, low_risk: f64, high_risk: f64) -> f64 {
    if value <= low_risk {
        0.0
    } else if value >= high_risk {
        1.0
    } else {
        (value - low_risk) / (high_risk - low_risk)
    }
}

pub(super) fn inverse_linear_risk(value: f64, high_risk: f64, low_risk: f64) -> f64 {
    if value <= high_risk {
        1.0
    } else if value >= low_risk {
        0.0
    } else {
        (low_risk - value) / (low_risk - high_risk)
    }
}

pub(super) fn log_risk(value: f64, low_risk: f64, high_risk: f64) -> f64 {
    if value <= low_risk {
        0.0
    } else if value >= high_risk {
        1.0
    } else {
        let log_value = value.ln();
        let log_low = low_risk.ln();
        let log_high = high_risk.ln();
        (log_value - log_low) / (log_high - log_low)
    }
}

pub(super) fn cumulative_pass_damage(pass_damage: f64, passes: f64) -> f64 {
    if pass_damage <= 0.0 || passes <= 0.0 {
        0.0
    } else {
        1.0 - (1.0 - pass_damage.clamp(0.0, 1.0)).powf(passes)
    }
}

pub(super) fn resonance_match(hydraulic_diameter_m: f64) -> f64 {
    if hydraulic_diameter_m <= 0.0 {
        0.0
    } else {
        (1.0 - ((hydraulic_diameter_m - ACOUSTIC_HALF_WAVELENGTH_M).abs()
            / ACOUSTIC_HALF_WAVELENGTH_M.max(1.0e-18)))
        .clamp(0.0, 1.0)
    }
}

/// Per-stage treatment-path flow fractions via Hagen–Poiseuille conductance.
///
/// Thin wrapper: maps [`SplitStageSpec`] branches to `(width, is_treatment)`
/// pairs and delegates to [`cfd_1d::cascade_treatment_flow_fractions`] for the
/// actual Q ∝ w³ conductance computation.
///
/// Uses a fixed-capacity stack buffer for the slice-of-slices reference array
/// (topologies have ≤ 8 stages), eliminating the second heap `Vec<&[...]>`.
pub(super) fn split_stage_flow_fractions(stages: &[SplitStageSpec]) -> (Vec<f64>, f64) {
    let stage_data: Vec<Vec<(f64, bool)>> = stages
        .iter()
        .map(|s| {
            s.branches
                .iter()
                .map(|b| (b.route.width_m, b.treatment_path))
                .collect()
        })
        .collect();
    // Stack buffer for the slice-of-slices (topologies have ≤ 8 stages).
    let mut refs_buf: [&[(f64, bool)]; 8] = [&[]; 8];
    let n = stage_data.len().min(refs_buf.len());
    for (slot, data) in refs_buf[..n].iter_mut().zip(stage_data.iter()) {
        *slot = data.as_slice();
    }
    cascade_treatment_flow_fractions(&refs_buf[..n])
}
