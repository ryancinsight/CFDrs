use serde::{Deserialize, Serialize};

use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

use super::blueprint_graph::solve_blueprint_candidate;
use super::blueprint_separation::{
    BlueprintSeparationMetrics, compute_blueprint_separation_metrics,
};
use super::residence::{ResidenceMetrics, compute_residence_metrics};
use super::safety::{BlueprintSafetyMetrics, compute_blueprint_safety_metrics};
use super::venturi::{BlueprintVenturiMetrics, compute_blueprint_venturi_metrics};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintEvaluation {
    pub residence: ResidenceMetrics,
    pub separation: BlueprintSeparationMetrics,
    pub venturi: BlueprintVenturiMetrics,
    pub safety: BlueprintSafetyMetrics,
}

pub fn evaluate_blueprint_candidate(
    candidate: &BlueprintCandidate,
) -> Result<BlueprintEvaluation, OptimError> {
    let solve = solve_blueprint_candidate(candidate)?;
    let residence = compute_residence_metrics(candidate, &solve);
    let separation = compute_blueprint_separation_metrics(candidate)?;
    let venturi = compute_blueprint_venturi_metrics(candidate, &solve, &separation)?;
    let safety = compute_blueprint_safety_metrics(candidate, &solve);
    Ok(BlueprintEvaluation {
        residence,
        separation,
        venturi,
        safety,
    })
}
