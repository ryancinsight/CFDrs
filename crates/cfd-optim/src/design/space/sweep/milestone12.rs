//! Milestone-12 restricted candidate space (tri-first PST only).

use super::super::builder::primitive_selective_candidate;
use super::super::dimensions::{pst_frac_slices, Milestone12Dimensions};
use crate::constraints::{CHANNEL_HEIGHT_M, THROAT_LENGTH_FACTORS};
use crate::design::space::params::TreatmentZoneMode;
use crate::design::{primitive_sequence_metadata, TRI_FIRST_PRIMITIVE_SELECTIVE_SEQUENCES};
use crate::domain::BlueprintCandidate;

#[cfg(any(test, debug_assertions))]
const MILESTONE12_CAPACITY: usize = 8_192;
#[cfg(not(any(test, debug_assertions)))]
const MILESTONE12_CAPACITY: usize = 80_000;

/// Build only the tri-first primitive selective candidates needed by the
/// Milestone 12 report pipeline.
///
/// This avoids allocating and scanning unrelated topology families when the
/// report only consumes the selective-routing lineage:
/// `Option 1 base -> Option 2 venturi -> GA refinement`.
#[must_use]
pub fn build_milestone12_candidate_space() -> Vec<BlueprintCandidate> {
    let dims = Milestone12Dimensions::active();
    let cav = &dims.cavitation;
    let mut candidates = Vec::with_capacity(MILESTONE12_CAPACITY);
    let mut idx: u32 = 0;

    for &sequence in &TRI_FIRST_PRIMITIVE_SELECTIVE_SEQUENCES {
        let seq_tag = sequence.label().replace('→', "");
        let (has_intermediate_tri, has_any_tri, has_any_bi) = primitive_sequence_metadata(sequence);
        let (pretri_fracs, tri_center_fracs, bi_treat_fracs) =
            pst_frac_slices(has_intermediate_tri, has_any_tri, has_any_bi);

        for &pretri_center_frac in pretri_fracs {
            for &terminal_tri_center_frac in tri_center_fracs {
                for &bi_treat_frac in bi_treat_fracs {
                    // Acoustic tier (ultrasound-only)
                    for &q in dims.flows {
                        for &gauge in dims.gauges {
                            for &w_ch in dims.widths {
                                for &n_segs in &[1_usize, 5, 7] {
                                    idx += 1;
                                    let acoustic_id = format!(
                                        "{:04}-PST-{}-pcf{}-tcf{}-btf{}-uo-q{:.0}ml-g{:.0}kPa-w{:.0}um-h{}-n{}",
                                        idx,
                                        seq_tag,
                                        (pretri_center_frac * 1000.0).round() as u32,
                                        (terminal_tri_center_frac * 1000.0).round() as u32,
                                        (bi_treat_frac * 1000.0).round() as u32,
                                        q * 6e7,
                                        gauge * 1e-3,
                                        w_ch * 1e6,
                                        (CHANNEL_HEIGHT_M * 1e6) as u32,
                                        n_segs,
                                    );
                                    candidates.push(primitive_selective_candidate(
                                        acoustic_id,
                                        sequence,
                                        q,
                                        gauge,
                                        0.0,
                                        0.0,
                                        w_ch,
                                        n_segs,
                                        pretri_center_frac,
                                        terminal_tri_center_frac,
                                        bi_treat_frac,
                                        TreatmentZoneMode::UltrasoundOnly,
                                        0,
                                    ));
                                }
                            }
                        }
                    }

                    // Venturi + cavitation tiers (with throat count sweep)
                    for &vt_count in dims.venturi_counts {
                        // Standard venturi tier
                        for &q in dims.flows {
                            for &gauge in dims.gauges {
                                for &d_throat in dims.throats {
                                    for &tl_factor in &THROAT_LENGTH_FACTORS {
                                        for &w_ch in dims.widths {
                                            for &n_segs in &[1_usize, 5, 7] {
                                                idx += 1;
                                                let throat_len = d_throat * tl_factor;
                                                let venturi_id = format!(
                                                    "{:04}-PST-{}-pcf{}-tcf{}-btf{}-vt{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}-n{}",
                                                    idx,
                                                    seq_tag,
                                                    (pretri_center_frac * 1000.0).round() as u32,
                                                    (terminal_tri_center_frac * 1000.0).round() as u32,
                                                    (bi_treat_frac * 1000.0).round() as u32,
                                                    vt_count,
                                                    q * 6e7,
                                                    gauge * 1e-3,
                                                    d_throat * 1e6,
                                                    tl_factor as u32,
                                                    w_ch * 1e6,
                                                    (CHANNEL_HEIGHT_M * 1e6) as u32,
                                                    n_segs,
                                                );
                                                candidates.push(primitive_selective_candidate(
                                                    venturi_id,
                                                    sequence,
                                                    q,
                                                    gauge,
                                                    d_throat,
                                                    throat_len,
                                                    w_ch,
                                                    n_segs,
                                                    pretri_center_frac,
                                                    terminal_tri_center_frac,
                                                    bi_treat_frac,
                                                    TreatmentZoneMode::VenturiThroats,
                                                    vt_count,
                                                ));
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        // Cavitation tier
                        for &q in cav.flows {
                            for &gauge in cav.gauges {
                                for &d_throat in cav.throats {
                                    for &tl_factor in &THROAT_LENGTH_FACTORS {
                                        for &w_ch in dims.widths {
                                            for &n_segs in &[1_usize, 5] {
                                                idx += 1;
                                                let throat_len = d_throat * tl_factor;
                                                let cavitation_id = format!(
                                                    "{:04}-PST-{}-pcf{}-tcf{}-btf{}-vt{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}-n{}",
                                                    idx,
                                                    seq_tag,
                                                    (pretri_center_frac * 1000.0).round() as u32,
                                                    (terminal_tri_center_frac * 1000.0).round() as u32,
                                                    (bi_treat_frac * 1000.0).round() as u32,
                                                    vt_count,
                                                    q * 6e7,
                                                    gauge * 1e-3,
                                                    d_throat * 1e6,
                                                    tl_factor as u32,
                                                    w_ch * 1e6,
                                                    (CHANNEL_HEIGHT_M * 1e6) as u32,
                                                    n_segs,
                                                );
                                                candidates.push(primitive_selective_candidate(
                                                    cavitation_id,
                                                    sequence,
                                                    q,
                                                    gauge,
                                                    d_throat,
                                                    throat_len,
                                                    w_ch,
                                                    n_segs,
                                                    pretri_center_frac,
                                                    terminal_tri_center_frac,
                                                    bi_treat_frac,
                                                    TreatmentZoneMode::VenturiThroats,
                                                    vt_count,
                                                ));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    candidates
}
