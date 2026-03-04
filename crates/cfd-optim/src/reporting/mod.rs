//! Report generation: ranking helpers and canonical markdown writer.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`ranking`]  | Deterministic tie-break sort, oncology priority score, shortlisting |
//! [`markdown`] | `write_milestone12_results`, `ValidationRow` |

pub mod markdown;
pub mod ranking;

pub use markdown::{write_milestone12_results, ValidationRow};
pub use ranking::{oncology_priority_score, pct_diff, shortlist_report, sort_by_report_priority};
