//! Outlet collection tree accumulation.
//!
//! Mirror of the inlet distribution tree.  Placed **after** the venturi
//! section in `compute_metrics` so that the cavitation-number calculation
//! only sees inlet-side pressure drop.

use cfd_core::physics::fluid::blood::CassonBlood;

use super::flow_path::FlowAccumulator;
use crate::constraints::TREATMENT_HEIGHT_MM;
use crate::design::{DesignCandidate, DesignTopology};

/// Accumulate pressure drop, shear, haemolysis, path length, and residence
/// time for the outlet-side collection tree (mirror of the inlet distribution
/// tree).
///
/// Placed **after** the venturi section so that the σ calculation (which uses
/// `total_dp` at the venturi) only sees inlet-side pressure — adding the
/// outlet tree beforehand would incorrectly reduce the apparent available
/// pressure and understate cavitation potential.
pub(super) fn accumulate_outlet_tree(
    acc: &mut FlowAccumulator,
    blood: &CassonBlood<f64>,
    candidate: &DesignCandidate,
    q: f64,
    w: f64,
    h: f64,
) {
    match candidate.topology {
        // No symmetric outlet merge tree for these topologies:
        DesignTopology::SingleVenturi
        | DesignTopology::SerialDoubleVenturi
        | DesignTopology::VenturiSerpentine
        | DesignTopology::SerpentineGrid
        | DesignTopology::CellSeparationVenturi
        | DesignTopology::WbcCancerSeparationVenturi => {}

        DesignTopology::BifurcationVenturi => {
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q / 2.0, w, h, branch_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TrifurcationVenturi => {
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.50e-3;
            acc.add_rect(blood, q / 3.0, w, h, branch_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::DoubleBifurcationVenturi => {
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            acc.add_rect(blood, q / 4.0, w, h, branch2_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TripleBifurcationVenturi => {
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            acc.add_rect(blood, q / 8.0, w, h, branch3_len);
            acc.add_rect(blood, q / 4.0, w, h, branch2_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::DoubleTrifurcationVenturi => {
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            acc.add_rect(blood, q / 9.0, w, h, branch2_len);
            acc.add_rect(blood, q / 3.0, w, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::BifurcationTrifurcationVenturi => {
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q / 6.0, w, h, branch2_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TrifurcationBifurcationVenturi => {
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q / 6.0, w, h, branch2_len);
            acc.add_rect(blood, q / 3.0, w, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TripleTrifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let w1 = w * frac;
            let w2 = w1 * frac;
            let w3 = w2 * frac;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.07e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.09e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            acc.add_rect(blood, q * q_frac.powi(3), w3, h, branch3_len);
            acc.add_rect(blood, q * q_frac.powi(2), w2, h, branch2_len);
            acc.add_rect(blood, q * q_frac, w1, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TrifurcationBifurcationBifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let w1 = w * frac;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.11e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.18e-3;
            acc.add_rect(blood, q * q_frac / 4.0, w1 * 0.25, h, branch3_len);
            acc.add_rect(blood, q * q_frac / 2.0, w1 * 0.5, h, branch2_len);
            acc.add_rect(blood, q * q_frac, w1, h, branch1_len);
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::QuadTrifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let level_lens = [
                TREATMENT_HEIGHT_MM * 0.05e-3,
                TREATMENT_HEIGHT_MM * 0.06e-3,
                TREATMENT_HEIGHT_MM * 0.08e-3,
                TREATMENT_HEIGHT_MM * 0.10e-3,
                TREATMENT_HEIGHT_MM * 0.12e-3,
            ];
            let qn4 = q * q_frac.powi(4);
            let qn3 = q * q_frac.powi(3);
            let qn2 = q * q_frac.powi(2);
            let qn1 = q * q_frac;
            let wn4 = w * frac.powi(4);
            let wn3 = w * frac.powi(3);
            let wn2 = w * frac.powi(2);
            let wn1 = w * frac;
            acc.add_rect(blood, qn4, wn4, h, level_lens[0]);
            acc.add_rect(blood, qn3, wn3, h, level_lens[1]);
            acc.add_rect(blood, qn2, wn2, h, level_lens[2]);
            acc.add_rect(blood, qn1, wn1, h, level_lens[3]);
            acc.add_rect(blood, q, w, h, level_lens[4]);
        }

        DesignTopology::AsymmetricTrifurcationVenturi => {
            let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let q_frac = cfd_1d::tri_center_q_frac(candidate.trifurcation_center_frac);
            acc.add_rect(
                blood,
                q * q_frac,
                w * candidate.trifurcation_center_frac,
                h,
                branch_len,
            );
            acc.add_rect(blood, q, w, h, trunk_len);
        }
        DesignTopology::PrimitiveSelectiveTree { sequence } => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let stage_len = TREATMENT_HEIGHT_MM * 0.14e-3;
            let mut stages: Vec<(f64, f64)> = Vec::with_capacity(sequence.levels() as usize);
            let mut qn = q;
            let mut wn = w;
            for stage_idx in 0..sequence.levels() {
                let is_tri = ((sequence.split_types() >> stage_idx) & 1) == 1;
                if is_tri {
                    let center_frac = primitive_tri_stage_fraction(candidate, sequence, stage_idx);
                    qn *= cfd_1d::tri_center_q_frac(center_frac);
                    wn *= center_frac;
                } else {
                    qn *= candidate.cif_terminal_bi_treat_frac();
                    wn *= candidate.cif_terminal_bi_treat_frac();
                }
                stages.push((qn, wn));
            }
            for (stage_q, stage_w) in stages.into_iter().rev() {
                acc.add_rect(blood, stage_q, stage_w, h, stage_len);
            }
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::BifurcationSerpentine => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::TrifurcationSerpentine => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.33e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        DesignTopology::AsymmetricBifurcationSerpentine => {
            let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 / 6.0;
            acc.add_rect(blood, q, w, h, trunk_len);
        }

        // Leukapheresis topologies: straight single-path — no symmetric outlet merge tree.
        DesignTopology::ConstrictionExpansionArray { .. }
        | DesignTopology::SpiralSerpentine { .. }
        | DesignTopology::ParallelMicrochannelArray { .. } => {}

        DesignTopology::AdaptiveTree {
            levels,
            split_types,
        } => {
            if levels > 0 {
                let level_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (f64::from(levels) + 1.0);
                let mut total_fan = 1usize;
                for i in 0..levels as usize {
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    total_fan *= fan;
                }
                let mut divisor = total_fan;
                for i in (0..levels as usize).rev() {
                    acc.add_rect(blood, q / divisor as f64, w, h, level_len);
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    divisor /= fan;
                }
                acc.add_rect(blood, q, w, h, level_len);
            }
        }
    }
}

fn primitive_tri_stage_fraction(
    candidate: &DesignCandidate,
    sequence: crate::design::PrimitiveSplitSequence,
    stage_idx: u8,
) -> f64 {
    let tri_indices: Vec<u8> = (0..sequence.levels())
        .filter(|idx| ((sequence.split_types() >> idx) & 1) == 1)
        .collect();
    if tri_indices.last().copied() == Some(stage_idx) {
        candidate.cif_terminal_tri_center_frac()
    } else {
        candidate.cif_pretri_center_frac()
    }
}

