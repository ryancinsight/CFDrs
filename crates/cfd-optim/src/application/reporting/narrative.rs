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
        OptimizationGoal::SelectiveAcousticResidenceSeparation => format!(
            "Option 1 uses asymmetric bifurcation and trifurcation widths to bias hydraulic resistance toward the treatment lane, yielding a treatment-zone residence time of {:.4} s while maintaining an RBC peripheral fraction of {:.4}.",
            evaluation.residence.treatment_residence_time_s,
            evaluation.separation.rbc_peripheral_fraction,
        ),
        OptimizationGoal::SelectiveVenturiCavitation => format!(
            "Option 2 augments the asymmetric cascade with venturi throats placed according to vena-contracta and Dean-number screening. The strongest screened cavitation state reached sigma = {:.4} while RBC and WBC exposure fractions were {:.4} and {:.4}, respectively.",
            evaluation
                .venturi
                .placements
                .iter()
                .map(|placement| placement.cavitation_number)
                .fold(f64::INFINITY, f64::min),
            evaluation.venturi.rbc_exposure_fraction,
            evaluation.venturi.wbc_exposure_fraction,
        ),
        OptimizationGoal::BlueprintGeneticRefinement => format!(
            "The GA stage preserves Option 2 lineage and refines local curvature, serpentine placement, and throat geometry in place. The resulting refinement score was {:.4}, with a peak Dean number of {:.4}.",
            evaluation.score,
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
    if let Some(topology) = candidate.blueprint.topology_spec() {
        for stage in &topology.split_stages {
            audit.push(ComponentAuditEntry {
                component_id: stage.stage_id.clone(),
                component_type: format!("{:?}", stage.split_kind),
                governing_physics: "Hydraulic resistance partitioning, Zweifach-Fung skimming, and inertial focusing through unequal daughter widths.".to_string(),
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
