//! Re-exports for the `resistance::traits` facade.
//!
//! ## Canonical Source
//!
//! All trait definitions are the **canonical implementations** in
//! [`crate::resistance::models::traits`]. This facade makes them accessible
//! via the shorter path `cfd_1d::resistance::traits::ResistanceModel`.
//!
//! ## Re-exported Traits & Types
//!
//! - [`ResistanceModel`]: Core trait defining the hydraulic resistance interface.
//! - [`FlowConditions`]: Value object describing flow state for resistance calculation.
//!
//! ## Design Rationale (DIP + SSOT)
//!
//! Resistance models depend on the abstract `ResistanceModel` trait, not on
//! concrete model types. The single definition in `models::traits` is the
//! authoritative source; this facade provides a shorter import path without
//! duplication.
pub use crate::resistance::models::traits::{FlowConditions, ResistanceModel};
