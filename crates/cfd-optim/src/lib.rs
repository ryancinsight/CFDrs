//! # cfd-optim — Millifluidic SDT Design Optimiser
//!
//! Generates and ranks the **top-5 millifluidic designs** for Sonodynamic
//! Therapy (SDT) using either:
//!
//! - **Hydrodynamic cavitation** through venturi throats — collapses microbubbles
//!   to activate sonosensitisers in targeted tissue wells.
//! - **Uniform light / ultrasound exposure** across a 6 × 6 well-plate centre
//!   zone — ensures every well receives equal treatment dose.
//!
//! All designs are constrained to fit within a standard **96-well plate footprint**
//! (ANSI/SLAS 1-2004, 127.76 × 85.47 mm) and comply with the FDA haemolysis
//! guideline for blood-contacting devices (≤ 150 Pa sustained wall shear).
//!
//! ## Quick start
//!
//! ```rust,ignore
//! use cfd_optim::build_milestone12_blueprint_candidate_space;
//!
//! let candidates = build_milestone12_blueprint_candidate_space()
//!     .expect("Milestone 12 candidate space should build");
//! tracing::info!("{} blueprint candidates", candidates.len());
//! ```
//!
//! ## Module structure
//!
//! | Module | Purpose |
//! |--------|---------|
//! [`domain`]        | Blueprint-native candidate and optimisation goal types |
//! [`application`]   | Goal evaluation, deterministic search, GA refinement, reporting |
//! [`metrics`]       | Blueprint-native and legacy physics evaluation layers |
//! [`constraints`]   | 96-well plate geometry, blood properties, FDA limits |
//! [`delivery`]      | JSON/SVG export |
//! [`design`]        | Milestone 12 blueprint candidate-space helpers |
//! [`error`]         | Error types |

#![warn(clippy::all)]
#![warn(clippy::pedantic)]
// Numerical cast patterns common in optimisation scoring.
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
// Physics variable names (x, y, k, n, …) are conventional.
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::manual_let_else)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::useless_conversion)]
#![allow(clippy::inline_always)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::format_push_string)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::new_without_default)]
#![allow(clippy::trivially_copy_pass_by_ref)]
#![allow(clippy::float_cmp)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::unused_self)]
#![allow(clippy::struct_excessive_bools)]

// ── SDT-specific modules ─────────────────────────────────────────────────────

mod analysis;
mod application;
pub mod constraints;
mod delivery;
mod design;
mod domain;
mod error;
mod metrics;
mod reporting;
mod scoring;

// ── Top-level re-exports ─────────────────────────────────────────────────────

pub use analysis::{robustness_sweep_blueprint, RobustnessReport, STANDARD_PERTURBATIONS};
pub use application::objectives::{
    evaluate_blueprint_genetic_refinement, evaluate_goal,
    evaluate_selective_acoustic_residence_separation, evaluate_selective_venturi_cavitation,
    BlueprintEvaluationStatus, BlueprintObjectiveEvaluation,
};
pub use application::orchestration::{
    blueprint_lineage_key as orchestration_lineage_key, fast_env, fast_mode,
    ga_matches_lineage_sequence, init_tracing, is_selective_report_topology, option2_mode,
    refresh_milestone12_reports, report_eligible_venturi_oncology, resolve_output_directories,
    run_milestone12_ga, run_milestone12_option1, run_milestone12_option2,
    run_milestone12_report, run_milestone12_validation, save_figure, Milestone12GaRun,
    Milestone12Option1Run, Milestone12Option2Run, Milestone12RequestedStage,
    Milestone12StageArtifact, Milestone12ValidationRun, ScanProgress,
};
pub use application::reporting::evidence::{
    build_goal_evidence, render_canonical_report, render_exploratory_report,
    validate_canonical_manifest, write_canonical_report, write_exploratory_report,
    EvidenceRunManifest, GoalEvidence, ValidationEvidence,
};
pub use application::reporting::figures::{required_figure_ids, FigureManifestEntry};
pub use application::reporting::narrative::{
    build_component_audit, render_goal_narrative, ComponentAuditEntry,
};
pub use application::search::deterministic::{
    build_blueprint_candidates_from_specs, rank_blueprint_candidates,
};
pub use application::search::genetic::{
    BlueprintEvolutionResult, BlueprintGeneticOptimizer, BlueprintRankedCandidate,
};
pub use application::search::mutations::{
    build_milestone12_ga_seed_pair, generate_ga_mutations, promote_option1_candidate_to_ga_seed,
    seed_option2_candidates,
};
pub use application::search::pool::EvaluatedPool;
pub use constraints::{BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA, P_ATM_PA};
pub use delivery::{
    load_top5_report_json, save_blueprint_schematic_svg, save_json_pretty, save_top5_report_json,
};
pub use design::build_milestone12_blueprint_candidate_space;
pub use domain::{BlueprintCandidate, OperatingPoint, OptimizationGoal, PatientContext};
pub use error::OptimError;
pub use metrics::{
    compute_blueprint_safety_metrics, compute_blueprint_separation_metrics,
    compute_blueprint_venturi_metrics, compute_residence_metrics, evaluate_blueprint_candidate,
    giersiepen_hi, solve_blueprint_candidate, BlueprintEvaluation, BlueprintSafetyMetrics,
    BlueprintSeparationMetrics, BlueprintSolveSample, BlueprintSolveSummary,
    BlueprintVenturiMetrics, ChannelHemolysis, ResidenceMetrics, SdtMetrics,
    StageBlueprintSeparationSummary, VenturiPlacementMetrics,
};
pub use reporting::{
    audit_goal_candidates, compute_blueprint_report_metrics, is_milestone12_lineage_topology,
    milestone12_lineage_key, pct_diff,
    run_milestone12_validation as run_milestone12_validation_solver,
    shortlist_report_designs, sort_report_designs, validate_milestone12_candidate,
    write_goal_audit_report, write_milestone12_narrative_report, write_milestone12_results,
    GoalAuditArtifacts, GoalAuditEntry, GoalAuditStatus, Milestone12LineageKey,
    Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage, ParetoPoint, ParetoTag,
    ValidationRow,
};
pub use scoring::{score_candidate, score_description, OptimMode, SdtWeights};
