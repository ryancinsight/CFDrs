//! Shared helper functions for Milestone 12 example pipelines.

use std::path::Path;

use crate::domain::BlueprintCandidate;
use crate::metrics::SdtMetrics;
use crate::reporting::{is_milestone12_lineage_topology, milestone12_lineage_key};

pub use crate::reporting::Milestone12LineageKey;

/// Check whether a candidate belongs to the Milestone 12 selective-routing
/// topology family.
pub fn is_selective_report_topology(candidate: &BlueprintCandidate) -> bool {
    is_milestone12_lineage_topology(candidate)
}

/// Extract the lineage key for a candidate, if it belongs to the Milestone 12
/// selective-routing family.
pub fn blueprint_lineage_key(candidate: &BlueprintCandidate) -> Option<Milestone12LineageKey> {
    milestone12_lineage_key(candidate)
}

/// Relaxed lineage match for GA results: same `PrimitiveSelectiveTree` sequence
/// as the selected lineage, but geometry parameters may differ (the GA mutates
/// channel width, serpentine counts, etc.).
pub fn ga_matches_lineage_sequence(
    candidate: &BlueprintCandidate,
    selected_key: &Milestone12LineageKey,
) -> bool {
    candidate
        .topology_spec()
        .is_ok_and(|spec| spec.stage_sequence_label() == selected_key.stage_sequence_label())
}

/// Check whether a venturi candidate's metrics satisfy the minimum eligibility
/// criteria for Option 2 oncology scoring.
pub fn report_eligible_venturi_oncology(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.therapy_channel_fraction > 0.0
        && metrics.cavitation_number < 1.0
}

/// Save a blueprint schematic SVG with logging.
pub fn save_figure(
    blueprint: &cfd_schematics::NetworkBlueprint,
    path: &Path,
    label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    crate::delivery::save_blueprint_schematic_svg(blueprint, path)?;
    tracing::info!("      \u{2713} {label}  \u{2192}  {}", path.display());
    Ok(())
}
