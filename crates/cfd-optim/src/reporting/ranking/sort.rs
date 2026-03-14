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
///
/// Pre-computes `oncology_priority_score` for each element once (O(n)) instead
/// of recomputing it O(n log n) times inside the comparator.
pub fn sort_report_designs(ranked: &mut [Milestone12ReportDesign]) {
    // Pre-compute oncology scores to avoid redundant recalculation per comparison.
    let oncology_scores: Vec<f64> = ranked
        .iter()
        .map(|d| oncology_priority_score(&d.metrics))
        .collect();
    // Build index array for indirect sort, then permute in-place.
    let mut indices: Vec<usize> = (0..ranked.len()).collect();
    indices.sort_unstable_by(|&i, &j| {
        ranked[j]
            .score
            .total_cmp(&ranked[i].score)
            .then_with(|| oncology_scores[j].total_cmp(&oncology_scores[i]))
            .then_with(|| {
                ranked[i]
                    .metrics
                    .rbc_venturi_exposure_fraction
                    .total_cmp(&ranked[j].metrics.rbc_venturi_exposure_fraction)
            })
            .then_with(|| {
                ranked[i]
                    .metrics
                    .clotting_risk_index
                    .total_cmp(&ranked[j].metrics.clotting_risk_index)
            })
            .then_with(|| ranked[i].candidate.id.cmp(&ranked[j].candidate.id))
    });
    // Permute ranked in-place following the sorted index order.
    permute_in_place(ranked, &mut indices);
}

/// Permute `data` in-place according to `indices` (cycle-leader algorithm).
fn permute_in_place<T>(data: &mut [T], indices: &mut [usize]) {
    for i in 0..indices.len() {
        let mut current = i;
        while indices[current] != i {
            let target = indices[current];
            data.swap(current, target);
            indices[current] = current;
            current = target;
        }
        indices[current] = current;
    }
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
