//! Deterministic ranking helpers for Milestone 12 report generation.

use crate::metrics::SdtMetrics;
use crate::reporting::Milestone12ReportDesign;

/// Absolute percentage difference between two values.
#[must_use]
pub fn pct_diff(a: f64, b: f64) -> f64 {
    let denom = a.abs();
    if denom < 1.0e-12 {
        f64::NAN
    } else {
        (a - b).abs() / denom * 100.0
    }
}

/// Composite oncology priority score used in deterministic report tie-breaks.
#[must_use]
pub fn oncology_priority_score(metrics: &SdtMetrics) -> f64 {
    let cancer = metrics.cancer_targeted_cavitation.clamp(0.0, 1.0);
    let selectivity = metrics.oncology_selectivity_index.clamp(0.0, 1.0);
    let rbc_shield = (1.0 - metrics.rbc_venturi_exposure_fraction).clamp(0.0, 1.0);
    let hi_gate = (1.0 - metrics.projected_hemolysis_15min_pediatric_3kg / 0.01).clamp(0.0, 1.0);
    0.45 * cancer + 0.30 * selectivity + 0.15 * rbc_shield + 0.10 * hi_gate
}

/// Sort blueprint-native report records by score and deterministic tie-breaks.
pub fn sort_report_designs(ranked: &mut [Milestone12ReportDesign]) {
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

/// Truncate blueprint-native Milestone 12 report designs after ordering.
///
/// # Errors
/// Returns an error if `ranked.len() < n`.
pub fn shortlist_report_designs(
    mut ranked: Vec<Milestone12ReportDesign>,
    n: usize,
    label: &str,
) -> Result<Vec<Milestone12ReportDesign>, Box<dyn std::error::Error>> {
    sort_report_designs(&mut ranked);
    if ranked.len() < n {
        return Err(format!("{label}: required {n} candidates, found {}", ranked.len()).into());
    }
    ranked.truncate(n);
    for (index, design) in ranked.iter_mut().enumerate() {
        design.rank = index + 1;
    }
    Ok(ranked)
}
