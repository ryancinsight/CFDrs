//! Report generation: ranking helpers and canonical markdown writer.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`ranking`]  | Deterministic tie-break sort, oncology priority score, shortlisting |
//! [`markdown`] | `write_milestone12_results`, `ValidationRow` |
//! [`narrative`] | Milestone 12 narrative templating + figure manifest generation |

pub mod contract_trace;
pub mod figures;
mod figures_svg;
pub mod markdown;
pub mod narrative;
mod narrative_sections;
pub mod ranking;
pub mod template;

pub use markdown::{write_milestone12_results, ValidationRow};
pub use narrative::{
    write_milestone12_narrative_report, M12Metadata, Milestone12NarrativeArtifacts,
    Milestone12NarrativeInput,
};
pub use ranking::{oncology_priority_score, pct_diff, shortlist_report, sort_by_report_priority};
