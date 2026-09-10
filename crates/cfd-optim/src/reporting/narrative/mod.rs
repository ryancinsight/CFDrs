//! Milestone 12 narrative report: templating, section builders, addenda, and contract tracing.

mod addenda;
mod contract;
mod sections;
pub(crate) mod template;
mod writer;

pub use writer::{
    Milestone12GaRankingAuditEntry, Milestone12NarrativeInput, write_milestone12_narrative_report,
};
