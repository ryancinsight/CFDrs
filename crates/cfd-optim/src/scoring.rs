//! Multi-objective scoring functions for SDT design candidates.
//!
//! Two optimisation modes are provided:
//!
//! - **[`OptimMode::SdtCavitation`]** — maximise cavitation intensity at the
//!   venturi throat while maintaining FDA haemolysis compliance in the main
//!   channels.
//!
//! - **[`OptimMode::UniformExposure`]** — maximise spatial uniformity of flow
//!   across all 36 treatment wells and maximise residence time in the exposure
//!   zone.
//!
//! - **[`OptimMode::Combined`]** — weighted combination of both objectives.
//!
//! Any candidate that violates a hard constraint — non-feasible pressure drop,
//! FDA main-channel shear exceedance, or plate-boundary overflow — receives a
//! hard score of `0.0`.  The smooth penalty mode provides the alternative
//! search signal when gradient-preserving feasibility shaping is required.

mod constraints;
mod description;
mod modes;
mod score;
mod types;

pub use constraints::sigmoid_penalty;
pub use description::score_description;
pub use score::score_candidate;
pub use types::{OptimMode, ScoreMode, SdtWeights};

/// Score assigned to hard-constraint violations.
const INFEASIBILITY_SCORE: f64 = 0.0;

#[cfg(test)]
mod tests;
