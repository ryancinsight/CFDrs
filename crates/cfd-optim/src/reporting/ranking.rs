//! Deterministic ranking helpers for Milestone 12 report generation.
//!
//! Provides a 5-tier deterministic tie-break sort, the oncology priority score
//! composite, and the shortlisting function used by all three report tracks
//! (Option 1, Option 2, RBC-protected).

use crate::metrics::SdtMetrics;
use crate::RankedDesign;

/// Absolute percentage difference between two values.
///
/// Returns `NAN` when `a` is below numerical zero (`|a| < 1e-12`).
pub fn pct_diff(a: f64, b: f64) -> f64 {
    let denom = a.abs();
    if denom < 1e-12 {
        f64::NAN
    } else {
        (a - b).abs() / denom * 100.0
    }
}

/// Composite oncology priority score used in the 5-tier deterministic tie-break.
///
/// Cancer-directed focus:
/// - maximize cancer-targeted cavitation and oncology selectivity,
/// - prefer lower RBC venturi exposure,
/// - penalize elevated pediatric cumulative hemolysis.
pub fn oncology_priority_score(metrics: &SdtMetrics) -> f64 {
    let cancer = metrics.cancer_targeted_cavitation.clamp(0.0, 1.0);
    let selectivity = metrics.oncology_selectivity_index.clamp(0.0, 1.0);
    let rbc_shield = (1.0 - metrics.rbc_venturi_exposure_fraction).clamp(0.0, 1.0);
    let hi_gate = (1.0 - metrics.projected_hemolysis_15min_pediatric_3kg / 0.01).clamp(0.0, 1.0);
    0.45 * cancer + 0.30 * selectivity + 0.15 * rbc_shield + 0.10 * hi_gate
}

/// Sort a ranked-design slice by the 5-tier deterministic tie-break:
///
/// 1. Score descending
/// 2. Oncology-priority score descending
/// 3. RBC venturi exposure ascending
/// 4. Clot-risk index ascending
/// 5. Candidate ID lexical ascending
pub fn sort_by_report_priority(ranked: &mut [RankedDesign]) {
    ranked.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then_with(|| {
                oncology_priority_score(&b.metrics).total_cmp(&oncology_priority_score(&a.metrics))
            })
            .then_with(|| {
                a.metrics
                    .rbc_venturi_exposure_fraction
                    .total_cmp(&b.metrics.rbc_venturi_exposure_fraction)
            })
            .then_with(|| {
                a.metrics
                    .clotting_risk_index
                    .total_cmp(&b.metrics.clotting_risk_index)
            })
            .then_with(|| a.candidate.id.cmp(&b.candidate.id))
    });
}

/// Apply the deterministic tie-break sort and return the top-`n` designs.
///
/// Reassigns `rank` fields to 1-based positions after sorting.
///
/// # Errors
/// Returns an error if `ranked.len() < n`.
pub fn shortlist_report(
    mut ranked: Vec<RankedDesign>,
    n: usize,
    label: &str,
) -> Result<Vec<RankedDesign>, Box<dyn std::error::Error>> {
    sort_by_report_priority(&mut ranked);
    if ranked.len() < n {
        return Err(format!("{label}: required {n} candidates, found {}", ranked.len()).into());
    }
    ranked.truncate(n);
    for (i, d) in ranked.iter_mut().enumerate() {
        d.rank = i + 1;
    }
    Ok(ranked)
}
