//! Flow-path accumulation helpers for `compute_metrics`.
//!
//! The inlet distribution tree and outlet collection tree each contain a
//! large topology match.  Rather than keep these inline in one 850-line
//! function, they live here as free functions operating on a shared
//! [`FlowAccumulator`].

use cfd_core::physics::fluid::blood::CassonBlood;

use super::compute::giersiepen_hi;
use crate::constraints::TREATMENT_HEIGHT_MM;
use crate::design::{DesignCandidate, DesignTopology};

// ── Rectangular-duct helper (Shah-London Poiseuille) ─────────────────────────

/// Returns `(ΔP [Pa], wall_shear_stress [Pa])` for a straight rectangular
/// channel carrying flow `q_ch` with dimensions `w_ch × h_ch` and length `len`.
///
/// Wall shear rate at narrow wall: `γ ≈ 6 Q / (w h²)`.
/// ΔP = `12 μ L Q / (w h³) × f_SL`  where `f_SL = 1 − 0.63/(w/h)`.
pub(super) fn rect_metrics(
    blood: &CassonBlood<f64>,
    q_ch: f64,
    w_ch: f64,
    h_ch: f64,
    len: f64,
) -> (f64, f64) {
    let v = q_ch / (w_ch * h_ch);
    let gamma = 6.0 * v / h_ch;
    let mu = blood.apparent_viscosity(gamma);
    let aspect = (w_ch / h_ch).max(1.0);
    let f_sl = 1.0 - 0.63 / aspect;
    let dp = (12.0 * mu * len * q_ch / (w_ch * h_ch.powi(3)) * f_sl).max(0.0);
    let shear = (mu * gamma).max(0.0);
    (dp, shear)
}

// ── Flow accumulator ─────────────────────────────────────────────────────────

/// Mutable accumulation state threaded through the inlet and outlet tree
/// match blocks.
pub(super) struct FlowAccumulator {
    pub total_dp: f64,
    pub max_main_shear_pa: f64,
    pub total_hi: f64,
    pub total_path_len_m: f64,
    pub residence_time_s: f64,
}

impl FlowAccumulator {
    pub fn new() -> Self {
        Self {
            total_dp: 0.0,
            max_main_shear_pa: 0.0,
            total_hi: 0.0,
            total_path_len_m: 0.0,
            residence_time_s: 0.0,
        }
    }

    /// Accumulate one straight rectangular channel segment.
    pub fn add_rect(
        &mut self,
        blood: &CassonBlood<f64>,
        q_ch: f64,
        w_ch: f64,
        h_ch: f64,
        len: f64,
    ) {
        let (dp, shear) = rect_metrics(blood, q_ch, w_ch, h_ch, len);
        self.total_dp += dp;
        self.max_main_shear_pa = self.max_main_shear_pa.max(shear);
        let t = len * w_ch * h_ch / q_ch.max(1e-12);
        self.total_hi += giersiepen_hi(shear, t);
        self.total_path_len_m += len;
        self.residence_time_s += t;
    }
}

// ── Inlet distribution tree ──────────────────────────────────────────────────

/// Accumulate pressure drop, shear, haemolysis, path length, and residence
/// time for the inlet-side distribution tree.  Each topology arm mirrors
/// the channel layout from inlet trunk down to the venturi.
pub(super) fn accumulate_inlet_tree(
    acc: &mut FlowAccumulator,
    blood: &CassonBlood<f64>,
    candidate: &DesignCandidate,
    q: f64,
    w: f64,
    h: f64,
) {
    match candidate.topology {
        DesignTopology::SingleVenturi | DesignTopology::SerialDoubleVenturi => {
            // Inlet → venturi(s) → outlet.  No distribution tree.
        }
        DesignTopology::BifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 2.0, w, h, branch_len);
        }
        DesignTopology::TrifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.5e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 3.0, w, h, branch_len);
        }
        DesignTopology::VenturiSerpentine | DesignTopology::SerpentineGrid => {
            // Venturi and serpentine handled separately; no extra trunk.
        }
        DesignTopology::CellSeparationVenturi | DesignTopology::WbcCancerSeparationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.5e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
        }
        DesignTopology::DoubleBifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q / 4.0, w, h, branch2_len);
        }
        DesignTopology::TripleBifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q / 4.0, w, h, branch2_len);
            acc.add_rect(blood, q / 8.0, w, h, branch3_len);
        }
        DesignTopology::DoubleTrifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 3.0, w, h, branch1_len);
            acc.add_rect(blood, q / 9.0, w, h, branch2_len);
        }
        DesignTopology::BifurcationTrifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 2.0, w, h, branch1_len);
            acc.add_rect(blood, q / 6.0, w, h, branch2_len);
        }
        DesignTopology::TrifurcationBifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q / 3.0, w, h, branch1_len);
            acc.add_rect(blood, q / 6.0, w, h, branch2_len);
        }
        DesignTopology::TripleTrifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let w1 = w * frac;
            let w2 = w1 * frac;
            let w3 = w2 * frac;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.09e-3;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.07e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q * q_frac, w1, h, branch1_len);
            acc.add_rect(blood, q * q_frac.powi(2), w2, h, branch2_len);
            acc.add_rect(blood, q * q_frac.powi(3), w3, h, branch3_len);
        }
        DesignTopology::TrifurcationBifurcationBifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let w1 = w * frac;
            let trunk_len = TREATMENT_HEIGHT_MM * 0.18e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.11e-3;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(blood, q * q_frac, w1, h, branch1_len);
            acc.add_rect(blood, q * q_frac / 2.0, w1 * 0.5, h, branch2_len);
            acc.add_rect(blood, q * q_frac / 4.0, w1 * 0.25, h, branch3_len);
        }
        DesignTopology::QuadTrifurcationVenturi => {
            let frac = candidate.trifurcation_center_frac;
            let q_frac = cfd_1d::tri_center_q_frac(frac);
            let level_lens = [
                TREATMENT_HEIGHT_MM * 0.12e-3,
                TREATMENT_HEIGHT_MM * 0.10e-3,
                TREATMENT_HEIGHT_MM * 0.08e-3,
                TREATMENT_HEIGHT_MM * 0.06e-3,
                TREATMENT_HEIGHT_MM * 0.05e-3,
            ];
            let mut qn = q;
            let mut wn = w;
            acc.add_rect(blood, qn, wn, h, level_lens[0]);
            for i in 1..=4 {
                qn *= q_frac;
                wn *= frac;
                acc.add_rect(blood, qn, wn, h, level_lens[i]);
            }
        }
        DesignTopology::AsymmetricTrifurcationVenturi => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let q_frac = cfd_1d::tri_center_q_frac(candidate.trifurcation_center_frac);
            acc.add_rect(blood, q, w, h, trunk_len);
            acc.add_rect(
                blood,
                q * q_frac,
                w * candidate.trifurcation_center_frac,
                h,
                branch_len,
            );
        }
        DesignTopology::PrimitiveSelectiveTree { sequence } => {
            let trunk_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let stage_len = TREATMENT_HEIGHT_MM * 0.14e-3;
            acc.add_rect(blood, q, w, h, trunk_len);

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
                acc.add_rect(blood, qn, wn, h, stage_len);
            }
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
        DesignTopology::ConstrictionExpansionArray { n_cycles } => {
            let wide_w = w;
            let narrow_w = w * 0.40;
            let seg_len = TREATMENT_HEIGHT_MM * 1e-3 / (n_cycles as f64 * 2.0);
            let narrow_len = seg_len * 0.5;
            for _ in 0..n_cycles {
                acc.add_rect(blood, q, wide_w, h, seg_len);
                acc.add_rect(blood, q, narrow_w, h, narrow_len);
            }
        }
        DesignTopology::SpiralSerpentine { n_turns } => {
            let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
            for _ in 0..n_turns {
                acc.add_rect(blood, q, w, h, turn_len);
            }
        }
        DesignTopology::ParallelMicrochannelArray { n_channels } => {
            let q_per = q / n_channels as f64;
            let ch_len = TREATMENT_HEIGHT_MM * 1e-3;
            acc.add_rect(blood, q_per, w, h, ch_len);
        }
        DesignTopology::AdaptiveTree {
            levels,
            split_types,
        } => {
            if levels > 0 {
                let level_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (f64::from(levels) + 1.0);
                let mut divisor = 1usize;
                acc.add_rect(blood, q, w, h, level_len);
                for i in 0..levels as usize {
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    divisor *= fan;
                    acc.add_rect(blood, q / divisor as f64, w, h, level_len);
                }
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

