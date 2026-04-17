use cfd_1d::{
    acoustic_contrast_factor, mixed_cascade_separation_kappa_aware,
    parallel_channel_flow_fractions, CascadeStage, PeripheralRecovery, KAPPA_CTC, KAPPA_PLASMA,
    RHO_CTC, RHO_PLASMA,
};
use cfd_schematics::topology::TreatmentActuationMode;
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
        let branch_dimensions = stage
            .branches
            .iter()
            .map(|branch| (branch.route.width_m, branch.route.height_m))
            .collect::<Vec<_>>();
        let branch_widths = branch_dimensions
            .iter()
            .map(|(width_m, _)| *width_m)
            .collect::<Vec<_>>();
        let q_fractions = conductance_flow_fractions(&branch_dimensions);

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
                    let sub_dimensions: Vec<(f64, f64)> = sub_split
                        .sub_branches
                        .iter()
                        .map(|sub_branch| (sub_branch.width_m, sub_branch.height_m))
                        .collect();
                    let sub_height = sub_dimensions
                        .first()
                        .map_or(branch.route.height_m, |(_, height_m)| *height_m);
                    let sub_q = conductance_flow_fractions(&sub_dimensions);
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
                    let source_idx = if branch_idx == 0 {
                        0
                    } else {
                        recovery_arm_idx.min(4)
                    };
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

        // Compute inflow velocity. Total flow for the chip = 500 mL/min = 8.333e-6 m3/s
        let chip_flow_rate_m3_s = 500.0 / 60.0 / 1e6;
        let parent_area_m2 = total_branch_width_m * stage.branches[0].route.height_m;
        // In a real network this would be scaled by upstream splits, but for the linear blueprint proxy
        // we assume the full chip flow (per parallel sequence) passes through this stage's parent.
        let parent_v_in_m_s = if parent_area_m2 > 1e-12 {
            chip_flow_rate_m3_s / parent_area_m2
        } else {
            0.0
        };

        cascade_stages.push(CascadeStage {
            arm_q_fracs,
            n_arms: stage.branches.len() as u8,
            treatment_dh_m,
            parent_v_in_m_s,
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

    let is_acoustic = topology.treatment_mode == TreatmentActuationMode::UltrasoundOnly
        || topology.treatment_mode == TreatmentActuationMode::VenturiCavitation;

    let cancer_center_corrected = if is_acoustic {
        let mut total_acoustic_yield = 0.0_f64;
        let ctc_radius = 17.5e-6_f64 / 2.0;
        let pressure_amp_pa = 500_000.0_f64; // proxy for 30Vpp equivalent
        let e_ac = cfd_1d::acoustic_energy_density(
            pressure_amp_pa,
            RHO_PLASMA,
            cfd_1d::physics::hemolysis::acoustic_radiation::SPEED_OF_SOUND_PLASMA,
        );
        let phi_ctc = acoustic_contrast_factor(RHO_CTC, RHO_PLASMA, KAPPA_CTC, KAPPA_PLASMA);
        let mu = 0.0035_f64;

        for stage in &cascade_stages {
            if stage.treatment_dh_m > 0.0 && stage.parent_v_in_m_s > 0.0 {
                // Approximate treatment length per stage
                let stage_length_m = 0.04;
                let t_transit = stage_length_m / stage.parent_v_in_m_s;
                let w = stage.treatment_dh_m;
                let k = std::f64::consts::PI / w;

                // Theorem: Gor'kov Acoustic Radiation Force Drift.
                // Transverse migration velocity v_drift = F_RAD / (6 * pi * mu * R)
                // where F_RAD_max = 4 * pi / 3 * R^3 * E_ac * k * Phi
                let f_rad_max =
                    (4.0 * std::f64::consts::PI / 3.0) * ctc_radius.powi(3) * k * e_ac * phi_ctc;
                let v_drift = f_rad_max / (6.0 * std::f64::consts::PI * mu * ctc_radius);

                // Migration distance to centerline is w/4 on average
                let t_mig = (w / 4.0) / v_drift.max(1.0e-15);
                // Exponential yield model for drift in uniform flow
                let yield_frac = 1.0 - (-t_transit / t_mig).exp();
                total_acoustic_yield += yield_frac * (1.0 / cascade_stages.len() as f64);
            }
        }
        (cascade.cancer_center_fraction
            + (1.0 - cascade.cancer_center_fraction) * total_acoustic_yield)
            .clamp(0.0, 1.0)
    } else {
        // Amini (2014) confinement-dependent lift correction.
        // For narrow treatment channels (D_h < 200 µm), cancer cells experience
        // enhanced focusing due to wall confinement (κ > κ_ref = 0.1).
        if let Some(last_stage) = cascade_stages.last() {
            let treatment_dh_m = last_stage.treatment_dh_m;
            if treatment_dh_m > 0.0 {
                let ctc_diameter = 17.5e-6; // MCF7 breast cancer CTC diameter [m]
                let kappa = (ctc_diameter / (2.0 * treatment_dh_m).max(1e-9)).clamp(0.0, 0.5);
                let amini_factor = cfd_1d::amini_confinement_correction(kappa);
                (cascade.cancer_center_fraction * amini_factor).clamp(0.0, 1.0)
            } else {
                cascade.cancer_center_fraction
            }
        } else {
            cascade.cancer_center_fraction
        }
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

fn conductance_flow_fractions(branch_dimensions: &[(f64, f64)]) -> Vec<f64> {
    parallel_channel_flow_fractions(branch_dimensions)
}

#[cfg(test)]
mod tests {
    use super::conductance_flow_fractions;

    #[test]
    fn conductance_flow_fractions_match_shared_rectangular_helper_behavior() {
        let fractions = conductance_flow_fractions(&[(2.0e-3, 1.0e-3), (1.0e-3, 1.0e-3)]);
        assert_eq!(fractions.len(), 2);
        assert!(fractions[0] > fractions[1]);
        assert!((fractions.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn height_sensitive_conductance_flow_fractions_track_rectangular_conductance() {
        let equal_height = conductance_flow_fractions(&[(1.0e-3, 1.0e-3), (1.0e-3, 1.0e-3)]);
        let unequal_height = conductance_flow_fractions(&[(1.0e-3, 2.0e-3), (1.0e-3, 1.0e-3)]);

        assert!((equal_height[0] - 0.5).abs() < 1.0e-12);
        assert!(
            unequal_height[0] > unequal_height[1],
            "taller equal-width branch should draw more flow under rectangular conductance"
        );
    }
}
