//! Multi-objective Pareto front computation for SDT design candidates.
//!
//! Uses the NSGA-II non-dominated sort from `cfd_math::statistics::pareto` to
//! identify designs that simultaneously optimise cancer-targeted cavitation,
//! RBC safety (lysis risk), and three-population separation efficiency.
//!
//! # Objectives
//!
//! | # | Metric | Direction |
//! |---|--------|-----------|
//! | 1 | `cancer_targeted_cavitation` | maximise |
//! | 2 | `lysis_risk_index` | **minimise** (inverted: 1 − lysis) |
//! | 3 | `three_pop_sep_efficiency` | maximise |
//!
//! # References
//! - Deb, K. et al. (2002) "A fast and elitist multiobjective genetic algorithm: NSGA-II"

use cfd_math::statistics::{crowding_distances, pareto_front_nd};

use crate::optimizer::RankedDesign;

// ── Output types ──────────────────────────────────────────────────────────────

/// Pareto-optimal set for the three primary SDT objectives.
///
/// Members are non-dominated across (cancer_cav, lysis_safety, sep3_efficiency).
/// Crowding distances quantify how isolated each member is in objective space —
/// higher distance = more diverse = preferred under NSGA-II selection pressure.
#[derive(Debug, Clone)]
pub struct SdtParetoFront {
    /// Non-dominated candidates (Pareto front members), in input order.
    pub members: Vec<RankedDesign>,
    /// NSGA-II crowding distance per member (same index as `members`).
    /// Boundary members on each objective axis receive [`f64::INFINITY`].
    pub crowding_distances: Vec<f64>,
}

impl SdtParetoFront {
    /// Number of designs on the Pareto front.
    #[must_use]
    pub fn len(&self) -> usize {
        self.members.len()
    }

    /// Returns `true` if the Pareto front is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// Candidate with the highest `cancer_targeted_cavitation` on the front.
    #[must_use]
    pub fn best_cancer_cav(&self) -> Option<&RankedDesign> {
        self.members.iter().max_by(|a, b| {
            a.metrics
                .cancer_targeted_cavitation
                .partial_cmp(&b.metrics.cancer_targeted_cavitation)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Candidate with the lowest `lysis_risk_index` (safest for RBCs).
    #[must_use]
    pub fn safest_rbc(&self) -> Option<&RankedDesign> {
        self.members.iter().min_by(|a, b| {
            a.metrics
                .lysis_risk_index
                .partial_cmp(&b.metrics.lysis_risk_index)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Candidate with the highest `three_pop_sep_efficiency`.
    #[must_use]
    pub fn best_sep3(&self) -> Option<&RankedDesign> {
        self.members.iter().max_by(|a, b| {
            a.metrics
                .three_pop_sep_efficiency
                .partial_cmp(&b.metrics.three_pop_sep_efficiency)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Top-k members by crowding distance (most isolated / diverse first).
    ///
    /// Returns up to `k` members sorted by descending crowding distance.
    /// Members with `f64::INFINITY` distance sort to the front.
    #[must_use]
    pub fn top_k_by_crowding(&self, k: usize) -> Vec<&RankedDesign> {
        let mut indexed: Vec<(usize, f64)> = self
            .crowding_distances
            .iter()
            .enumerate()
            .map(|(i, &d)| (i, d))
            .collect();

        // Sort descending: INFINITY > finite, NaN last
        indexed.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Less));

        indexed
            .into_iter()
            .take(k)
            .map(|(i, _)| &self.members[i])
            .collect()
    }
}

// ── Main computation ──────────────────────────────────────────────────────────

/// Compute the 3-objective Pareto front from a slice of feasible candidates.
///
/// **Feasibility filter should be applied by the caller** — pass only candidates
/// where `metrics.pressure_feasible && metrics.fda_main_compliant` are both
/// `true` (or whichever fitness gate is appropriate for your analysis).
///
/// The three objectives are:
/// 1. Maximise `cancer_targeted_cavitation`
/// 2. Maximise `1 − lysis_risk_index.clamp(0, 1)` (= minimise lysis risk)
/// 3. Maximise `three_pop_sep_efficiency`
///
/// Returns an [`SdtParetoFront`] with members in the same relative order as
/// `feasible` and corresponding NSGA-II crowding distances.
#[must_use]
pub fn compute_sdt_pareto_front(feasible: &[RankedDesign]) -> SdtParetoFront {
    if feasible.is_empty() {
        return SdtParetoFront {
            members: Vec::new(),
            crowding_distances: Vec::new(),
        };
    }

    // Build objective vectors for non-dominated sort.
    // All three are cast to "maximise" by inverting lysis_risk.
    let objectives: Vec<Vec<f64>> = feasible
        .iter()
        .map(|d| {
            vec![
                d.metrics.cancer_targeted_cavitation,
                1.0 - d.metrics.lysis_risk_index.clamp(0.0, 1.0),
                d.metrics.three_pop_sep_efficiency,
            ]
        })
        .collect();

    let is_maximized = [true, true, true];
    let front_indices = pareto_front_nd(&objectives, &is_maximized);

    // Collect front objectives in the same index order for crowding distance.
    let front_objs: Vec<Vec<f64>> = front_indices
        .iter()
        .map(|&i| objectives[i].clone())
        .collect();

    let dists = crowding_distances(&front_objs);

    let members: Vec<RankedDesign> = front_indices.iter().map(|&i| feasible[i].clone()).collect();

    SdtParetoFront {
        members,
        crowding_distances: dists,
    }
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        design::{CrossSectionShape, DesignCandidate, DesignTopology},
        metrics::SdtMetrics,
    };

    /// Construct a minimal `RankedDesign` with only the three Pareto metrics set.
    fn make_design(cancer_cav: f64, lysis_risk: f64, sep3: f64) -> RankedDesign {
        let mut metrics = SdtMetrics::default();
        metrics.cancer_targeted_cavitation = cancer_cav;
        metrics.lysis_risk_index = lysis_risk;
        metrics.three_pop_sep_efficiency = sep3;
        metrics.pressure_feasible = true;
        metrics.fda_main_compliant = true;

        RankedDesign {
            rank: 1,
            candidate: DesignCandidate {
                id: "test".into(),
                topology: DesignTopology::SingleVenturi,
                flow_rate_m3_s: 5e-6,
                inlet_gauge_pa: 300_000.0,
                throat_diameter_m: 100e-6,
                inlet_diameter_m: 4e-3,
                throat_length_m: 200e-6,
                channel_width_m: 4e-3,
                channel_height_m: 1e-3,
                serpentine_segments: 6,
                segment_length_m: 0.045,
                bend_radius_m: 4.5e-3,
                feed_hematocrit: 0.45,
                trifurcation_center_frac: 1.0 / 3.0,
                cif_pretri_center_frac: 1.0 / 3.0,
                cif_terminal_tri_center_frac: 1.0 / 3.0,
                cif_terminal_bi_treat_frac: 0.68,
                asymmetric_narrow_frac: 0.5,
                trifurcation_left_frac: 1.0 / 3.0,
                cross_section_shape: CrossSectionShape::Rectangular,
            },
            metrics,
            score: 0.5,
        }
    }

    #[test]
    fn empty_input_returns_empty_front() {
        let front = compute_sdt_pareto_front(&[]);
        assert!(front.is_empty());
        assert_eq!(front.len(), 0);
    }

    #[test]
    fn single_candidate_is_pareto_optimal() {
        let d = make_design(0.5, 0.001, 0.6);
        let front = compute_sdt_pareto_front(std::slice::from_ref(&d));
        assert_eq!(front.len(), 1);
    }

    #[test]
    fn dominated_design_excluded() {
        // d0 dominates d1 on all three objectives
        let d0 = make_design(0.8, 0.001, 0.7); // high cav, low lysis, high sep
        let d1 = make_design(0.5, 0.005, 0.4); // dominated on all axes
        let front = compute_sdt_pareto_front(&[d0, d1]);
        assert_eq!(front.len(), 1);
        assert!((front.members[0].metrics.cancer_targeted_cavitation - 0.8).abs() < 1e-10);
    }

    #[test]
    fn two_incomparable_designs_both_on_front() {
        // d0 better on cancer_cav; d1 better on lysis safety — trade-off
        let d0 = make_design(0.9, 0.010, 0.5); // high cav, higher lysis
        let d1 = make_design(0.3, 0.001, 0.5); // low cav, safe lysis
        let front = compute_sdt_pareto_front(&[d0, d1]);
        assert_eq!(front.len(), 2, "both designs should be non-dominated");
    }

    #[test]
    fn best_helpers_return_correct_extremes() {
        let d0 = make_design(0.9, 0.010, 0.4); // best cancer_cav
        let d1 = make_design(0.3, 0.001, 0.4); // best lysis safety (safest_rbc)
        let d2 = make_design(0.5, 0.005, 0.9); // best sep3
        let front = compute_sdt_pareto_front(&[d0, d1, d2]);
        assert_eq!(front.len(), 3);

        let best_cav = front.best_cancer_cav().unwrap();
        assert!((best_cav.metrics.cancer_targeted_cavitation - 0.9).abs() < 1e-10);

        let safest = front.safest_rbc().unwrap();
        assert!((safest.metrics.lysis_risk_index - 0.001).abs() < 1e-10);

        let best_sep = front.best_sep3().unwrap();
        assert!((best_sep.metrics.three_pop_sep_efficiency - 0.9).abs() < 1e-10);
    }

    #[test]
    fn top_k_by_crowding_returns_at_most_k() {
        let designs: Vec<RankedDesign> = (0..5)
            .map(|i| {
                let t = i as f64 / 4.0;
                make_design(t, 0.01 - t * 0.009, 0.5)
            })
            .collect();
        let front = compute_sdt_pareto_front(&designs);
        let top3 = front.top_k_by_crowding(3);
        assert!(top3.len() <= 3);
    }
}
