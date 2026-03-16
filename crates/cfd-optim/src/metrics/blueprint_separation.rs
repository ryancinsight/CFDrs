use cfd_1d::{mixed_cascade_separation_kappa_aware, CascadeStage, PeripheralRecovery};
use serde::{Deserialize, Serialize};

use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageBlueprintSeparationSummary {
    pub stage_id: String,
    pub n_arms: u8,
    pub treatment_flow_fraction: f64,
    pub wbc_exclusion_flow_fraction: f64,
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
        let mut wbc_exclusion_flow_fraction = 0.0_f64;
        let mut rbc_bypass_flow_fraction = 0.0_f64;
        let mut treatment_width_m = 0.0_f64;
        let total_branch_width_m = branch_widths.iter().sum::<f64>();
        let mut arm_q_fracs = [0.0_f64; 5];

        let mut arm_idx = 0;
        for (branch, q_fraction) in stage.branches.iter().zip(q_fractions.iter().copied()) {
            match branch.role {
                cfd_schematics::BranchRole::Treatment => {
                    treatment_flow_fraction += q_fraction;
                    treatment_width_m += branch.route.width_m;
                    arm_q_fracs[0] += q_fraction;
                }
                cfd_schematics::BranchRole::WbcCollection => {
                    wbc_exclusion_flow_fraction += q_fraction;
                    arm_idx += 1;
                    if arm_idx < 5 {
                        arm_q_fracs[arm_idx] += q_fraction;
                    }
                }
                cfd_schematics::BranchRole::RbcBypass => {
                    rbc_bypass_flow_fraction += q_fraction;
                    arm_idx += 1;
                    if arm_idx < 5 {
                        arm_q_fracs[arm_idx] += q_fraction;
                    }
                }
                cfd_schematics::BranchRole::Neutral => {
                    wbc_exclusion_flow_fraction += q_fraction;
                    arm_idx += 1;
                    if arm_idx < 5 {
                        arm_q_fracs[arm_idx] += q_fraction;
                    }
                }
            }
        }

        let treatment_dh_m = 2.0 * treatment_width_m * stage.branches[0].route.height_m
            / (treatment_width_m + stage.branches[0].route.height_m).max(1.0e-18);

        // Build peripheral recovery entries from recovery_sub_split specs.
        let mut peripheral_recoveries: [Option<PeripheralRecovery>; 4] = [None; 4];
        let mut n_recoveries = 0_u8;
        let mut recovery_arm_idx = 0_usize;
        for (branch_idx, branch) in stage.branches.iter().enumerate() {
            if branch.treatment_path {
                continue;
            }
            recovery_arm_idx += 1;
            if let Some(ref sub_split) = branch.recovery_sub_split {
                if n_recoveries < 4 && !sub_split.sub_branches.is_empty() {
                    let sub_widths: Vec<f64> =
                        sub_split.sub_branches.iter().map(|sb| sb.width_m).collect();
                    let sub_height = sub_split
                        .sub_branches
                        .first()
                        .map_or(branch.route.height_m, |sb| sb.height_m);
                    let sub_q = conductance_flow_fractions(&sub_widths, sub_height);
                    let mut sub_arm_q_fracs = [0.0_f64; 5];
                    for (i, &q) in sub_q.iter().enumerate().take(5) {
                        sub_arm_q_fracs[i] = q;
                    }
                    let recovery_w = sub_split
                        .sub_branches
                        .get(sub_split.recovery_arm_index)
                        .map_or(sub_widths[0], |sb| sb.width_m);
                    let recovery_dh_m =
                        2.0 * recovery_w * sub_height / (recovery_w + sub_height).max(1e-18);
                    // source_arm_idx uses the same index as in arm_q_fracs
                    let source_idx = if branch_idx == 0 { 0 } else { recovery_arm_idx.min(4) };
                    peripheral_recoveries[n_recoveries as usize] = Some(PeripheralRecovery {
                        source_arm_idx: source_idx,
                        sub_arm_q_fracs,
                        n_sub_arms: sub_q.len().clamp(2, 5) as u8,
                        recovery_arm_idx: sub_split.recovery_arm_index,
                        recovery_dh_m,
                    });
                    n_recoveries += 1;
                }
            }
        }

        cascade_stages.push(CascadeStage {
            arm_q_fracs,
            n_arms: stage.branches.len() as u8,
            treatment_dh_m,
            peripheral_recoveries,
            n_recoveries,
        });
        stage_summaries.push(StageBlueprintSeparationSummary {
            stage_id: stage.stage_id.clone(),
            n_arms: stage.branches.len() as u8,
            treatment_flow_fraction,
            wbc_exclusion_flow_fraction,
            rbc_bypass_flow_fraction,
            treatment_width_m,
            total_branch_width_m,
        });
    }

    let cascade = mixed_cascade_separation_kappa_aware(&cascade_stages);

    // Amini (2014) confinement-dependent lift correction.
    // For narrow treatment channels (D_h < 200 µm), cancer cells experience
    // enhanced focusing due to wall confinement (κ > κ_ref = 0.1).
    // This multiplier increases the effective cancer_center_fraction.
    let cancer_center_corrected = if let Some(last_stage) = cascade_stages.last() {
        let treatment_dh_m = last_stage.treatment_dh_m;
        if treatment_dh_m > 0.0 {
            let ctc_diameter = 17.5e-6; // MCF7 breast cancer CTC diameter [m]
            // Clamp kappa to physical range [0, 0.5]. Amini (2014) validated
            // for κ ∈ [0.07, 0.3]; beyond 0.5 the particle blocks the channel.
            let kappa = (ctc_diameter / (2.0 * treatment_dh_m).max(1e-9)).clamp(0.0, 0.5);
            let amini_factor = cfd_1d::amini_confinement_correction(kappa);
            (cascade.cancer_center_fraction * amini_factor).clamp(0.0, 1.0)
        } else {
            cascade.cancer_center_fraction
        }
    } else {
        cascade.cancer_center_fraction
    };

    Ok(BlueprintSeparationMetrics {
        stage_summaries,
        cancer_center_fraction: cancer_center_corrected,
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
