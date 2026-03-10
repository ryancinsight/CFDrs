//! Shared orchestration helpers for Milestone 12 example pipelines.
//!
//! Consolidates tracing setup, workspace resolution, fast-mode parsing,
//! progress tracking, figure export, and topology filtering that were
//! previously duplicated across `milestone12_option1`, `milestone12_option2`,
//! `milestone12_ga`, and `milestone12_report`.

mod helpers;
mod progress;
mod setup;

pub use helpers::{
    blueprint_lineage_key, ga_matches_lineage_sequence, is_selective_report_topology,
    report_eligible_venturi_oncology, save_figure,
};
pub use progress::ScanProgress;
pub use setup::{fast_env, fast_mode, init_tracing, option2_mode, resolve_output_directories};
