//! Audit, lineage validation, and guardrail checks for Milestone 12 candidates.

mod goal_audit;
mod guardrails;

pub use goal_audit::{
    audit_goal_candidates, write_goal_audit_report, GoalAuditArtifacts, GoalAuditEntry,
    GoalAuditStatus,
};
pub use guardrails::{
    is_milestone12_lineage_topology, milestone12_lineage_key, validate_milestone12_candidate,
    Milestone12LineageKey, Milestone12Stage,
};
