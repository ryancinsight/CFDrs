//! Physics-based metric computation for SDT millifluidic design candidates.

mod blueprint_eval;
mod blueprint_graph;
mod blueprint_separation;
mod healthy_cell_protection;
mod residence;
mod safety;
mod sdt_metrics;
mod venturi;

pub use blueprint_eval::{evaluate_blueprint_candidate, BlueprintEvaluation};
pub use blueprint_graph::{solve_blueprint_candidate, BlueprintSolveSample, BlueprintSolveSummary};
pub use blueprint_separation::{
    compute_blueprint_separation_metrics, BlueprintSeparationMetrics,
    StageBlueprintSeparationSummary,
};
pub use cfd_1d::physics::hemolysis::giersiepen_hi;
pub use healthy_cell_protection::healthy_cell_protection_index;
pub use residence::{compute_residence_metrics, ResidenceMetrics};
pub use safety::{compute_blueprint_safety_metrics, BlueprintSafetyMetrics};
pub use sdt_metrics::{ChannelHemolysis, SdtMetrics};
pub use venturi::{
    compute_blueprint_venturi_metrics, BlueprintVenturiMetrics, VenturiPlacementMetrics,
};
