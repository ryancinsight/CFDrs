//! Report generation: ranking helpers and canonical markdown writer.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`ranking`]  | Deterministic tie-break sort, oncology priority score, shortlisting |
//! [`markdown`] | `write_milestone12_results`, `ValidationRow` |
//! [`narrative`] | Milestone 12 narrative templating + figure manifest generation |

mod contract_trace;
mod design_record;
mod figures;
mod figures_svg;
mod guardrails;
mod markdown;
mod narrative;
mod narrative_sections;
mod ranking;
mod report_metrics;
mod template;

pub use design_record::{compute_blueprint_report_metrics, Milestone12ReportDesign};
pub use guardrails::{
    is_milestone12_lineage_topology, milestone12_lineage_key, validate_milestone12_candidate,
    Milestone12LineageKey, Milestone12Stage,
};
pub use markdown::{write_milestone12_results, ValidationRow};
pub use narrative::{write_milestone12_narrative_report, Milestone12NarrativeInput};
pub use ranking::{pct_diff, shortlist_report_designs, sort_report_designs};
