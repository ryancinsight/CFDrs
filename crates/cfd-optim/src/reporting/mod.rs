//! Report generation: ranking helpers and canonical markdown writer.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`ranking`]  | Deterministic tie-break sort, oncology priority score, shortlisting |
//! [`markdown`] | `write_milestone12_results`, `ValidationRow` |
//! [`narrative`] | Milestone 12 narrative templating + figure manifest generation |
//! [`svg_primitives`] | Low-level SVG rendering primitives and generic bar-chart builder |

mod design_record;
mod figures;
mod guardrails;
mod markdown;
mod milestone12_audit;
mod narrative;
mod ranking;
mod report_math;
mod report_metrics;
mod validation_runner;

pub use design_record::{
    compute_blueprint_report_metrics, Milestone12ReportDesign, ParetoPoint, ParetoTag,
};
pub use guardrails::{
    is_milestone12_lineage_topology, milestone12_lineage_key, validate_milestone12_candidate,
    Milestone12LineageKey, Milestone12Stage,
};
pub use markdown::{write_milestone12_results, ValidationRow};
pub use milestone12_audit::{
    audit_goal_candidates, write_goal_audit_report, GoalAuditArtifacts, GoalAuditEntry,
    GoalAuditStatus,
};
pub use narrative::{write_milestone12_narrative_report, Milestone12NarrativeInput};
pub use ranking::{pct_diff, shortlist_report_designs, sort_report_designs};
pub use validation_runner::run_milestone12_validation;
