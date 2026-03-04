//! Design topology and candidate generation for SDT millifluidic devices.
//!
//! This module defines the [`DesignTopology`] enum, the
//! [`DesignCandidate`] parameter struct, blueprint conversion methods,
//! and the parametric sweep / random-sampling functions that populate
//! the initial candidate space.
//!
//! ## Sub-modules
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`topology`] | `DesignTopology` enum + query methods |
//! | [`candidate`] | `DesignCandidate` struct + accessor helpers |
//! | [`blueprint`] | `to_blueprint()`, `to_channel_system()`, `total_path_length_mm()` |
//! | [`space`] | `build_candidate_space()`, `sample_random_candidates()` |
//! | [`cif_schematic`] | Design Option 1 ESDT CIF SVG schematic generator |

mod blueprint;
mod candidate;
pub mod cif_schematic;
mod space;
mod topology;

pub use candidate::{CrossSectionShape, DesignCandidate, TreatmentZoneMode};
pub use space::{build_candidate_space, sample_random_candidates};
pub use topology::DesignTopology;
