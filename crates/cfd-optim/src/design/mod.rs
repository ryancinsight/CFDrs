//! Milestone 12 design-space helpers.
//!
//! ## Sub-modules
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`space`] | `build_milestone12_blueprint_candidate_space()` |

mod sequence_catalog;
pub(crate) mod space;
pub(crate) use sequence_catalog::{
    primitive_sequence_metadata, TRI_FIRST_PRIMITIVE_SELECTIVE_SEQUENCES,
};

pub use space::build_milestone12_blueprint_candidate_space;
