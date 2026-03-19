use std::sync::OnceLock;

use crate::{BlueprintCandidate, OptimError, SdtMetrics};
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// ParetoPoint — lightweight struct for Pareto-front visualization
// ---------------------------------------------------------------------------

/// Lightweight data point for Pareto-front scatter plots.
///
/// ~32 bytes per entry vs ~16-20 KB for a full [`Milestone12ReportDesign`].
/// Saved to disk as `option2_pool_all.json` / `ga_pool_all.json` so the
/// report phase can render the background cloud without loading full
/// candidates + metrics into memory.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ParetoPoint {
    pub cancer_targeted_cavitation: f64,
    pub rbc_venturi_protection: f64,
    pub score: f64,
    pub tag: ParetoTag,
}

/// Which optimization track produced this Pareto point.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ParetoTag {
    Option2,
    Ga,
}

const GA_NON_HYDROSDT_PENALTY: f64 = 0.25;
const GA_NO_CAVITATION_GAIN_PENALTY: f64 = 0.10;
const GA_CAVITATION_GAIN_MARGIN_DEFAULT: f64 = 0.01;

fn hydrosdt_cavitation_gain_margin() -> f64 {
    static MARGIN: OnceLock<f64> = OnceLock::new();
    *MARGIN.get_or_init(|| {
        std::env::var("M12_GA_CAVITATION_GAIN_MARGIN")
            .ok()
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| value.is_finite() && *value >= 0.0)
            .unwrap_or(GA_CAVITATION_GAIN_MARGIN_DEFAULT)
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Milestone12ReportDesign {
    pub rank: usize,
    pub candidate: BlueprintCandidate,
    pub metrics: SdtMetrics,
    pub score: f64,
}

impl Milestone12ReportDesign {
    #[must_use]
    pub fn new(
        rank: usize,
        candidate: BlueprintCandidate,
        metrics: SdtMetrics,
        score: f64,
    ) -> Self {
        Self {
            rank,
            candidate,
            metrics,
            score,
        }
    }

    pub fn from_blueprint_candidate(
        rank: usize,
        candidate: BlueprintCandidate,
        score: f64,
    ) -> Result<Self, OptimError> {
        let metrics = compute_blueprint_report_metrics(&candidate)?;
        Ok(Self::new(rank, candidate, metrics, score))
    }

    #[must_use]
    pub fn topology_display_name(&self) -> String {
        self.candidate.blueprint().topology_spec().map_or_else(
            || self.candidate.blueprint().name.clone(),
            |spec| spec.display_name(),
        )
    }

    #[must_use]
    pub fn topology_short_code(&self) -> String {
        self.candidate.blueprint().topology_spec().map_or_else(
            || self.candidate.blueprint().name.clone(),
            |spec| spec.short_code(),
        )
    }

    #[must_use]
    pub fn stage_sequence_label(&self) -> String {
        self.candidate
            .blueprint()
            .topology_spec()
            .map_or_else(String::new, |spec| spec.stage_sequence_label())
    }

    #[must_use]
    pub fn visible_split_layers(&self) -> usize {
        self.candidate
            .blueprint()
            .topology_spec()
            .map_or(0, |spec| spec.visible_split_layers())
    }

    #[must_use]
    pub fn flow_rate_ml_min(&self) -> f64 {
        self.candidate.operating_point.flow_rate_m3_s * 6.0e7
    }

    #[must_use]
    pub fn inlet_gauge_kpa(&self) -> f64 {
        self.candidate.operating_point.inlet_gauge_pa * 1.0e-3
    }

    #[must_use]
    pub fn throat_width_um(&self) -> Option<f64> {
        self.candidate
            .blueprint()
            .topology_spec()
            .and_then(|spec| spec.venturi_placements.first())
            .map(|placement| placement.throat_geometry.throat_width_m * 1.0e6)
    }
    #[must_use]
    pub fn throat_length_um(&self) -> Option<f64> {
        self.candidate
            .blueprint()
            .topology_spec()
            .and_then(|spec| spec.venturi_placements.first())
            .map(|placement| placement.throat_geometry.throat_length_m * 1.0e6)
    }
}

#[must_use]
pub fn pareto_point_from_report_design(
    design: &Milestone12ReportDesign,
    tag: ParetoTag,
) -> ParetoPoint {
    ParetoPoint {
        cancer_targeted_cavitation: design.metrics.cancer_targeted_cavitation,
        rbc_venturi_protection: design.metrics.rbc_venturi_protection,
        score: design.score,
        tag,
    }
}

#[must_use]
pub fn pareto_pool_from_report_designs(
    designs: &[Milestone12ReportDesign],
    tag: ParetoTag,
    limit: usize,
) -> Vec<ParetoPoint> {
    designs
        .iter()
        .filter(|design| {
            design.metrics.cancer_targeted_cavitation.is_finite()
                && design.metrics.rbc_venturi_protection.is_finite()
        })
        .take(limit)
        .map(|design| pareto_point_from_report_design(design, tag))
        .collect()
}

#[must_use]
pub fn is_hydrosdt_venturi_report_candidate(design: &Milestone12ReportDesign) -> bool {
    design.metrics.active_venturi_throat_count > 0
        && design.metrics.serial_venturi_stages_per_path > 0
        && design.metrics.cavitation_number.is_finite()
        && design.metrics.cavitation_number < 1.0
        && design.metrics.cancer_targeted_cavitation > 0.0
}

#[must_use]
pub fn cavitation_gain_over_baseline(
    design: &Milestone12ReportDesign,
    baseline_option2: &Milestone12ReportDesign,
) -> f64 {
    design.metrics.serial_cavitation_dose_fraction
        - baseline_option2.metrics.serial_cavitation_dose_fraction
}

#[must_use]
pub fn improves_hydrosdt_cavitation(
    design: &Milestone12ReportDesign,
    baseline_option2: &Milestone12ReportDesign,
) -> bool {
    is_hydrosdt_venturi_report_candidate(design)
        && cavitation_gain_over_baseline(design, baseline_option2)
            > hydrosdt_cavitation_gain_margin()
}

#[must_use]
pub fn ga_report_selection_score(
    design: &Milestone12ReportDesign,
    baseline_option2: &Milestone12ReportDesign,
) -> f64 {
    if !is_hydrosdt_venturi_report_candidate(design) {
        return design.score - GA_NON_HYDROSDT_PENALTY;
    }

    let cavitation_gain = cavitation_gain_over_baseline(design, baseline_option2);
    let required_margin = hydrosdt_cavitation_gain_margin();
    if cavitation_gain > required_margin {
        design.score
    } else {
        design.score
            - GA_NO_CAVITATION_GAIN_PENALTY
            - (required_margin - cavitation_gain).max(0.0)
    }
}

#[must_use]
pub fn rank_ga_hydrosdt_report_designs(
    raw_ranked: &[Milestone12ReportDesign],
    baseline_option2: &Milestone12ReportDesign,
) -> Vec<Milestone12ReportDesign> {
    let mut hydrosdt_ranked: Vec<Milestone12ReportDesign> = raw_ranked
        .iter()
        .filter(|design| is_hydrosdt_venturi_report_candidate(design))
        .cloned()
        .collect();

    hydrosdt_ranked.sort_by(|left, right| {
        ga_report_selection_score(right, baseline_option2)
            .total_cmp(&ga_report_selection_score(left, baseline_option2))
            .then_with(|| right.score.total_cmp(&left.score))
            .then_with(|| {
                right
                    .metrics
                    .serial_cavitation_dose_fraction
                    .total_cmp(&left.metrics.serial_cavitation_dose_fraction)
            })
            .then_with(|| left.candidate.id.cmp(&right.candidate.id))
    });

    let improving_ranked: Vec<Milestone12ReportDesign> = hydrosdt_ranked
        .iter()
        .filter(|design| improves_hydrosdt_cavitation(design, baseline_option2))
        .cloned()
        .collect();

    let mut ranked = if improving_ranked.is_empty() {
        hydrosdt_ranked
    } else {
        improving_ranked
    };
    for (index, design) in ranked.iter_mut().enumerate() {
        design.rank = index + 1;
    }
    ranked
}

pub fn compute_blueprint_report_metrics(
    candidate: &BlueprintCandidate,
) -> Result<SdtMetrics, OptimError> {
    crate::reporting::report_metrics::compute_blueprint_report_metrics(candidate)
}
