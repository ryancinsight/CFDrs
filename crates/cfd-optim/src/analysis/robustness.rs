//! Design robustness analysis via parametric perturbation sweeps.
//!
//! Evaluates how sensitive a blueprint's evaluation is to ±10 % / ±20 % changes in
//! key operating parameters (low_speed, inlet_gauge, etc.). The coefficient
//! of variation (CV = σ/μ) across all perturbed evaluations quantifies overall robustness.
//!
//! A design is considered *robust* when CV < 0.10 (< 10 % relative variation).

use serde::{Deserialize, Serialize};

use crate::{
    application::objectives::evaluate_goal,
    domain::{BlueprintCandidate, OptimizationGoal},
};

// ── Constants ─────────────────────────────────────────────────────────────────

/// Default perturbation fractions: ±10 % and ±20 %.
pub const STANDARD_PERTURBATIONS: [f64; 4] = [-0.20, -0.10, 0.10, 0.20];

/// CV threshold below which a design is classified as "robust".
pub const ROBUSTNESS_CV_THRESHOLD: f64 = 0.10;

// ── Output types ──────────────────────────────────────────────────────────────

/// Robustness summary for a single blueprint candidate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobustnessReport {
    pub candidate_id: String,
    pub score_nominal: f64,
    pub score_min: f64,
    pub score_max: f64,
    pub score_cv: f64,
    pub worst_case_param: String,
    pub is_robust: bool,
}

// ── Core function ─────────────────────────────────────────────────────────────

/// Run parametric perturbation sweeps and produce a [RobustnessReport].
#[must_use]
pub fn robustness_sweep_blueprint(
    candidate: &BlueprintCandidate,
    goal: OptimizationGoal,
    perturbation_fracs: &[f64],
) -> RobustnessReport {
    let score_nominal = evaluate_goal(candidate, goal).map_or(0.0, |e| e.score_or_zero());

    let mut scored: Vec<(f64, String)> = vec![(score_nominal, "nominal".into())];

    // ── Inlet pressure perturbations ─────────────────────────────────────────
    for &frac in perturbation_fracs {
        let mut c = candidate.clone();
        c.operating_point.inlet_gauge_pa *= 1.0 + frac;
        let s = evaluate_goal(&c, goal).map_or(0.0, |e| e.score_or_zero());
        scored.push((s, format!("inlet_pressure{:+.0}%", frac * 100.0)));
    }

    // ── Flow rate perturbations ───────────────────────────────────────────────
    // The blueprint assumes inlet flow rate as derived or explicitly set.
    // We adjust flow rate by perturbing the target flow speed limit constraint in operating_point.
    for &frac in perturbation_fracs {
        let mut c = candidate.clone();
        if true {
            c.operating_point.flow_rate_m3_s *= 1.0 + frac;
            let s = evaluate_goal(&c, goal).map_or(0.0, |e| e.score_or_zero());
            scored.push((s, format!("flow_rate{:+.0}%", frac * 100.0)));
        }
    }

    // ── Statistics ────────────────────────────────────────────────────────────
    let vals: Vec<f64> = scored.iter().map(|(s, _)| *s).collect();
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let variance = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let std_dev = variance.sqrt();
    let score_cv = if mean.abs() > 1e-9 {
        std_dev / mean
    } else {
        0.0
    };

    let (score_min, worst_label) = scored
        .iter()
        .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .cloned()
        .unwrap_or((0.0, "none".into()));

    let score_max = vals.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    RobustnessReport {
        candidate_id: candidate.id.clone(),
        score_nominal,
        score_min,
        score_max,
        score_cv,
        worst_case_param: worst_label,
        is_robust: score_cv < ROBUSTNESS_CV_THRESHOLD,
    }
}
