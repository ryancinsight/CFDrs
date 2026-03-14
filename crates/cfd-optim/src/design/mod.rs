//! Milestone 12 design-space helpers.
//!
//! ## Sub-modules
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`space`] | `build_milestone12_blueprint_candidate_space()` |

pub(crate) mod space;

pub use space::{build_milestone12_blueprint_candidate_space, build_milestone12_candidate_params, CandidateParams};
