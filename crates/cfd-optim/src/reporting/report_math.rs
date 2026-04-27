//! Pure math helpers for report-metrics computation.

use crate::constraints::ACOUSTIC_HALF_WAVELENGTH_M;
use cfd_1d::cascade_treatment_flow_fractions;
use cfd_schematics::topology::SplitStageSpec;

const GA_CONVERGENCE_TAIL_EPSILON: f64 = 1.0e-3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum GaConvergenceTrend {
    Improving {
        tail_len: usize,
        tail_delta: f64,
    },
    Regressing {
        tail_len: usize,
        tail_delta: f64,
    },
    NearPlateau {
        tail_len: usize,
        tail_delta_abs: f64,
    },
}

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

pub(super) fn percentile(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let index = ((sorted.len() - 1) as f64 * p.clamp(0.0, 1.0)).round() as usize;
    sorted[index]
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

/// Per-stage treatment-path flow fractions via rectangular laminar conductance.
///
/// Thin wrapper: maps [`SplitStageSpec`] branches to
/// `(width, height, is_treatment)` triples and delegates to
/// [`cfd_1d::cascade_treatment_flow_fractions`] for the actual rectangular
/// branch-conductance weighting.
pub(super) fn split_stage_flow_fractions(stages: &[SplitStageSpec]) -> (Vec<f64>, f64) {
    let stage_data: Vec<Vec<(f64, f64, bool)>> = stages
        .iter()
        .map(|s| {
            s.branches
                .iter()
                .map(|b| (b.route.width_m, b.route.height_m, b.treatment_path))
                .collect()
        })
        .collect();
    let refs: Vec<&[(f64, f64, bool)]> = stage_data.iter().map(|v| v.as_slice()).collect();
    cascade_treatment_flow_fractions(&refs)
}

/// Classify the trailing GA fitness window used by the Milestone 12 report.
///
/// The helper centralizes the trailing-delta logic so the convergence figure
/// and narrative prose describe the same window and threshold.
pub(super) fn ga_convergence_trend(best_per_gen: &[f64]) -> GaConvergenceTrend {
    let tail_len = best_per_gen.len().min(5);
    let tail_start = best_per_gen.len().saturating_sub(tail_len);
    let tail_delta = if tail_len >= 2 {
        best_per_gen[best_per_gen.len() - 1] - best_per_gen[tail_start]
    } else {
        0.0
    };

    if tail_delta > GA_CONVERGENCE_TAIL_EPSILON {
        GaConvergenceTrend::Improving {
            tail_len,
            tail_delta,
        }
    } else if tail_delta < -GA_CONVERGENCE_TAIL_EPSILON {
        GaConvergenceTrend::Regressing {
            tail_len,
            tail_delta,
        }
    } else {
        GaConvergenceTrend::NearPlateau {
            tail_len,
            tail_delta_abs: tail_delta.abs(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ga_convergence_trend_reports_trailing_improvement() {
        match ga_convergence_trend(&[0.70, 0.71, 0.72, 0.73, 0.74, 0.79]) {
            GaConvergenceTrend::Improving {
                tail_len,
                tail_delta,
            } => {
                assert_eq!(tail_len, 5);
                assert!((tail_delta - 0.08).abs() < 1.0e-12);
            }
            other => panic!("unexpected trend: {other:?}"),
        }
    }

    #[test]
    fn ga_convergence_trend_reports_near_plateau() {
        match ga_convergence_trend(&[0.70, 0.7002, 0.7001]) {
            GaConvergenceTrend::NearPlateau {
                tail_len,
                tail_delta_abs,
            } => {
                assert_eq!(tail_len, 3);
                assert!((tail_delta_abs - 0.0001).abs() < 1.0e-12);
            }
            other => panic!("unexpected trend: {other:?}"),
        }
    }

    #[test]
    fn ga_convergence_trend_reports_regression() {
        match ga_convergence_trend(&[0.79, 0.76, 0.75, 0.74]) {
            GaConvergenceTrend::Regressing {
                tail_len,
                tail_delta,
            } => {
                assert_eq!(tail_len, 4);
                assert!((tail_delta + 0.05).abs() < 1.0e-12);
            }
            other => panic!("unexpected trend: {other:?}"),
        }
    }
}
