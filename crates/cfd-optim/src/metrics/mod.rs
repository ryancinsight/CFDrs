//! Physics-based metric computation for SDT millifluidic design candidates.

mod blueprint_eval;
mod blueprint_graph;
mod blueprint_separation;
mod healthy_cell_protection;
mod residence;
mod safety;
mod sdt_metrics;
mod venturi;

pub use blueprint_eval::{BlueprintEvaluation, evaluate_blueprint_candidate};
pub use blueprint_graph::{BlueprintSolveSample, BlueprintSolveSummary, solve_blueprint_candidate};
pub use blueprint_separation::{
    BlueprintSeparationMetrics, StageBlueprintSeparationSummary,
    compute_blueprint_separation_metrics,
};
pub use cfd_1d::physics::hemolysis::giersiepen_hi;
pub use healthy_cell_protection::healthy_cell_protection_index;
pub(crate) use residence::compute_typed_residence_metrics;
pub use residence::{ResidenceMetrics, compute_residence_metrics};
pub(crate) use safety::compute_typed_blueprint_safety_metrics;
pub use safety::{BlueprintSafetyMetrics, compute_blueprint_safety_metrics};
pub use sdt_metrics::{ChannelHemolysis, SdtMetrics};
pub use venturi::{
    BlueprintVenturiMetrics, VenturiPlacementMetrics, compute_blueprint_venturi_metrics,
};
