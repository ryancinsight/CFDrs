//! Shared orchestration helpers for Milestone 12 example pipelines.
//!
//! Consolidates tracing setup, workspace resolution, fast-mode parsing,
//! progress tracking, figure export, and topology filtering that were
//! previously duplicated across `milestone12_option1`, `milestone12_option2`,
//! `milestone12_ga`, and `milestone12_report`.

mod helpers;
mod milestone12;
mod progress;
mod setup;

pub use helpers::{
    blueprint_lineage_key, ga_matches_lineage_sequence, is_selective_report_topology,
    report_eligible_venturi_oncology, save_figure,
};
pub use milestone12::{
    refresh_milestone12_reports, run_milestone12_ga, run_milestone12_option1,
    run_milestone12_option2, run_milestone12_report, run_milestone12_validation,
    Milestone12GaRun, Milestone12Option1Run, Milestone12Option2Run, Milestone12RequestedStage,
    Milestone12StageArtifact, Milestone12ValidationRun,
};
pub use progress::ScanProgress;
pub use setup::{
    ensure_release_reports, fast_env, fast_mode, init_tracing, milestone12_ranked_pool_size,
    option2_mode, resolve_output_directories,
};
