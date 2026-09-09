//! Canonical context and selected-schematic figure writers for Milestone 12.
//!
//! # Theorem
//! Figures 2–6 are reproducible from workspace SSOT: the selected-design
//! schematics are exported directly from their ranked geometry-authored
//! blueprints, and the concept figures are exported from canonical
//! `cfd-schematics` Milestone 12 topology requests.
//!
//! **Proof sketch**
//! Each figure is produced from one authoritative blueprint source and written
//! afresh to the report figure directory. The selected-design figures consume
//! the blueprints already stored in the ranked report designs, while the concept
//! figures are rebuilt from the canonical Milestone 12 topology catalog using a
//! root split-kind filter. Because no cached loose SVG is accepted as input,
//! the generated figure must reflect the current blueprint SSOT.

use std::path::Path;

use cfd_schematics::{
    build_milestone12_blueprint, enumerate_milestone12_topologies, SplitKind,
    TreatmentActuationMode,
};

use crate::delivery::save_blueprint_schematic_svg;
use crate::reporting::Milestone12ReportDesign;

pub(super) fn write_context_concept_figures(
    figures_dir: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    write_concept_schematic_figure(&figures_dir.join("treatment_zone_plate.svg"), 2)?;
    write_concept_schematic_figure(
        &figures_dir.join("treatment_zone_plate_trifurcation.svg"),
        3,
    )?;
    Ok(())
}

pub(super) fn write_selected_schematic_figure(
    path: &Path,
    design: &Milestone12ReportDesign,
) -> Result<(), Box<dyn std::error::Error>> {
    save_blueprint_schematic_svg(design.candidate.blueprint(), path)
}

fn write_concept_schematic_figure(
    path: &Path,
    root_arity: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let base_request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| {
            matches!(
                request.split_kinds.first(),
                Some(SplitKind::NFurcation(arity)) if *arity == root_arity
            )
        })
        .ok_or_else(|| format!("no Milestone 12 topology found for root arity {root_arity}"))?;

    let request = cfd_schematics::Milestone12TopologyRequest {
        treatment_mode: TreatmentActuationMode::UltrasoundOnly,
        venturi_throat_count: 0,
        venturi_target_channel_ids: Vec::new(),
        ..base_request
    };
    let blueprint = build_milestone12_blueprint(&request)?;
    save_blueprint_schematic_svg(&blueprint, path)
}
