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

pub fn compute_blueprint_report_metrics(
    candidate: &BlueprintCandidate,
) -> Result<SdtMetrics, OptimError> {
    crate::reporting::report_metrics::compute_blueprint_report_metrics(candidate)
}
