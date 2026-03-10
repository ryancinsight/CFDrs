use cfd_1d::{mixed_cascade_separation_kappa_aware, CascadeStage};
use serde::{Deserialize, Serialize};

use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageBlueprintSeparationSummary {
    pub stage_id: String,
    pub n_arms: u8,
    pub treatment_flow_fraction: f64,
    pub wbc_collection_flow_fraction: f64,
    pub rbc_bypass_flow_fraction: f64,
    pub treatment_width_m: f64,
    pub total_branch_width_m: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintSeparationMetrics {
    pub stage_summaries: Vec<StageBlueprintSeparationSummary>,
    pub cancer_center_fraction: f64,
    pub wbc_center_fraction: f64,
    pub rbc_peripheral_fraction: f64,
    pub separation_efficiency: f64,
    pub center_hematocrit_ratio: f64,
}

pub fn compute_blueprint_separation_metrics(
    candidate: &BlueprintCandidate,
) -> Result<BlueprintSeparationMetrics, OptimError> {
    let topology = candidate.topology_spec()?;
    let mut cascade_stages = Vec::with_capacity(topology.split_stages.len());
    let mut stage_summaries = Vec::with_capacity(topology.split_stages.len());

    for stage in &topology.split_stages {
        let branch_widths = stage
            .branches
            .iter()
            .map(|branch| branch.route.width_m)
            .collect::<Vec<_>>();
        let q_fractions =
            conductance_flow_fractions(&branch_widths, stage.branches[0].route.height_m);

        let mut treatment_flow_fraction = 0.0_f64;
        let mut wbc_collection_flow_fraction = 0.0_f64;
        let mut rbc_bypass_flow_fraction = 0.0_f64;
        let mut treatment_width_m = 0.0_f64;
        let total_branch_width_m = branch_widths.iter().sum::<f64>();
        let mut arm_q_fracs = [0.0_f64; 3];

        for (branch, q_fraction) in stage.branches.iter().zip(q_fractions.iter().copied()) {
            match branch.role {
                cfd_schematics::BranchRole::Treatment => {
                    treatment_flow_fraction += q_fraction;
                    treatment_width_m += branch.route.width_m;
                    arm_q_fracs[0] += q_fraction;
                }
                cfd_schematics::BranchRole::WbcCollection => {
                    wbc_collection_flow_fraction += q_fraction;
                    arm_q_fracs[1] += q_fraction;
                }
                cfd_schematics::BranchRole::RbcBypass => {
                    rbc_bypass_flow_fraction += q_fraction;
                    arm_q_fracs[2] += q_fraction;
                }
                cfd_schematics::BranchRole::Neutral => {
                    if arm_q_fracs[1] <= arm_q_fracs[2] {
                        arm_q_fracs[1] += q_fraction;
                        wbc_collection_flow_fraction += q_fraction;
                    } else {
                        arm_q_fracs[2] += q_fraction;
                        rbc_bypass_flow_fraction += q_fraction;
                    }
                }
            }
        }

        let treatment_dh_m = 2.0 * treatment_width_m * stage.branches[0].route.height_m
            / (treatment_width_m + stage.branches[0].route.height_m).max(1.0e-18);
        cascade_stages.push(CascadeStage {
            arm_q_fracs,
            n_arms: stage.branches.len() as u8,
            treatment_dh_m,
        });
        stage_summaries.push(StageBlueprintSeparationSummary {
            stage_id: stage.stage_id.clone(),
            n_arms: stage.branches.len() as u8,
            treatment_flow_fraction,
            wbc_collection_flow_fraction,
            rbc_bypass_flow_fraction,
            treatment_width_m,
            total_branch_width_m,
        });
    }

    let cascade = mixed_cascade_separation_kappa_aware(&cascade_stages);
    Ok(BlueprintSeparationMetrics {
        stage_summaries,
        cancer_center_fraction: cascade.cancer_center_fraction,
        wbc_center_fraction: cascade.wbc_center_fraction,
        rbc_peripheral_fraction: cascade.rbc_peripheral_fraction,
        separation_efficiency: cascade.separation_efficiency,
        center_hematocrit_ratio: cascade.center_hematocrit_ratio,
    })
}

fn conductance_flow_fractions(branch_widths: &[f64], height_m: f64) -> Vec<f64> {
    let h = height_m.max(1.0e-9);
    let conductances = branch_widths
        .iter()
        .map(|width| {
            let width = width.max(1.0e-9);
            let resistance = 1.0 / (width.powi(3) * h);
            let minor_loss_factor = 1.3 * (width / h).sqrt();
            1.0 / (resistance * (1.0 + minor_loss_factor * width / h))
        })
        .collect::<Vec<_>>();
    let total = conductances.iter().sum::<f64>().max(1.0e-18);
    conductances
        .into_iter()
        .map(|conductance| conductance / total)
        .collect()
}
