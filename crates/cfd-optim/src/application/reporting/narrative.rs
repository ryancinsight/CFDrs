use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ComponentAuditEntry {
    pub component_id: String,
    pub component_type: String,
    pub governing_physics: String,
    pub safety_relevance: String,
}

pub fn render_goal_narrative(evaluation: &BlueprintObjectiveEvaluation) -> String {
    match evaluation.goal {
        OptimizationGoal::AsymmetricSplitResidenceSeparation => format!(
            "Option 1 evaluates mirrored asymmetric split families as a residence-time and selective-routing study. The selected geometry sustains {:.4} s of treatment-zone dwell while maintaining an RBC peripheral fraction of {:.4}, indicating that the imposed daughter-width asymmetry continues to bias cell-rich core flow toward the treatment path without collapsing the healthy-cell bypass.",
            evaluation.residence.treatment_residence_time_s,
            evaluation.separation.rbc_peripheral_fraction,
        ),
        OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity => format!(
            "Option 2 derives from the same mirrored split lineage and adds venturi throats on treatment-zone channels to test how split-partitioned flow alters cavitation feasibility. The strongest screened cavitation state reached sigma = {:.4} while RBC and WBC exposure fractions were {:.4} and {:.4}, respectively.",
            evaluation
                .venturi
                .placements
                .iter()
                .map(|placement| placement.cavitation_number)
                .fold(f64::INFINITY, f64::min),
            evaluation.venturi.rbc_exposure_fraction,
            evaluation.venturi.wbc_exposure_fraction,
        ),
        OptimizationGoal::InPlaceDeanSerpentineRefinement => format!(
            "The GA stage preserves Option 2 lineage and applies in-place treatment-zone mutations, including serpentine insertion, venturi retargeting, and mirrored split-merge refinement. The resulting refinement score was {:.4}, with a peak Dean number of {:.4}.",
            evaluation.score_or_zero(),
            evaluation
                .venturi
                .placements
                .iter()
                .map(|placement| placement.dean_number)
                .fold(0.0_f64, f64::max),
        ),
    }
}

pub fn build_component_audit(candidate: &BlueprintCandidate) -> Vec<ComponentAuditEntry> {
    let mut audit = Vec::new();
    if let Some(render_hints) = candidate.blueprint.render_hints() {
        audit.push(ComponentAuditEntry {
            component_id: "milestone12_catalog_entry".to_string(),
            component_type: "CanonicalTopologyCatalog".to_string(),
            governing_physics: format!(
                "Stage sequence '{}' with mirror state (mirror_x={}, mirror_y={}) generated exclusively through create_geometry-authored schematics.",
                render_hints.stage_sequence, render_hints.mirror_x, render_hints.mirror_y
            ),
            safety_relevance: "Fixes the geometric context in which split partitioning, venturi cavitation, and remerge transport are interpreted across the report lineage.".to_string(),
        });
    }
    if let Some(topology) = candidate.blueprint.topology_spec() {
        for stage in &topology.split_stages {
            audit.push(ComponentAuditEntry {
                component_id: stage.stage_id.clone(),
                component_type: format!("{:?}", stage.split_kind),
                governing_physics: "Hydraulic resistance partitioning, Zweifach-Fung skimming, inertial focusing, and steric exclusion through unequal daughter widths.".to_string(),
                safety_relevance: "Controls which blood sub-populations are routed toward the treatment lane versus healthy-cell bypass lanes.".to_string(),
            });
            for branch in &stage.branches {
                audit.push(ComponentAuditEntry {
                    component_id: format!("{}:{}", stage.stage_id, branch.label),
                    component_type: format!("{:?}", branch.role),
                    governing_physics: format!(
                        "Rectangular microchannel flow with width {:.6e} m, height {:.6e} m, and length {:.6e} m.",
                        branch.route.width_m, branch.route.height_m, branch.route.length_m
                    ),
                    safety_relevance: if branch.treatment_path {
                        "Defines the treatment-lane residence time, cavitation exposure, and CTC enrichment path.".to_string()
                    } else {
                        "Provides an RBC/WBC-protective bypass path that limits cavitating exposure.".to_string()
                    },
                });
            }
        }
        for placement in &topology.venturi_placements {
            audit.push(ComponentAuditEntry {
                component_id: placement.placement_id.clone(),
                component_type: "VenturiPlacement".to_string(),
                governing_physics: format!(
                    "Vena-contracta screening with throat width {:.6e} m, length {:.6e} m, and placement mode {:?}.",
                    placement.throat_geometry.throat_width_m,
                    placement.throat_geometry.throat_length_m,
                    placement.placement_mode,
                ),
                safety_relevance: "Sets cavitation inception, throat shear, and the selective exposure burden experienced by RBC and WBC populations.".to_string(),
            });
        }
    }
    audit
}
