use std::sync::Arc;

use cfd_schematics::topology::presets::enumerate_milestone12_topologies;
use cfd_schematics::TreatmentActuationMode;

use super::parameters::CandidateParams;
use crate::constraints::THROAT_LENGTH_FACTORS;
use crate::design::space::dimensions::{pst_frac_slices, Milestone12Dimensions};

#[cfg(any(test, debug_assertions))]
const MILESTONE12_CAPACITY: usize = 8_192;
#[cfg(not(any(test, debug_assertions)))]
const MILESTONE12_CAPACITY: usize = 500_000;

#[must_use]
pub fn generate_milestone12_candidate_params() -> Vec<CandidateParams> {
    let dims = Milestone12Dimensions::active();
    let cav = &dims.cavitation;
    let mut candidates = Vec::with_capacity(MILESTONE12_CAPACITY);
    let mut idx: u32 = 0;

    let topologies: Vec<Arc<_>> = enumerate_milestone12_topologies()
        .into_iter()
        .map(Arc::new)
        .collect();

    for request in &topologies {
        let split_kinds = &request.split_kinds;
        let has_any_tri = split_kinds
            .iter()
            .any(|split_kind| matches!(split_kind, cfd_schematics::SplitKind::NFurcation(3)));
        let has_intermediate_tri = split_kinds
            .iter()
            .skip(1)
            .any(|split_kind| matches!(split_kind, cfd_schematics::SplitKind::NFurcation(3)));
        let has_any_bi = split_kinds
            .iter()
            .any(|split_kind| matches!(split_kind, cfd_schematics::SplitKind::NFurcation(2)));
        let (pretri_fracs, tri_center_fracs, bi_treat_fracs) =
            pst_frac_slices(has_intermediate_tri, has_any_tri, has_any_bi);
        let supports_venturi_materialization = split_kinds
            .iter()
            .all(|split_kind| matches!(split_kind, cfd_schematics::SplitKind::NFurcation(2..=5)));

        for &pretri_center_frac in pretri_fracs {
            for &terminal_tri_center_frac in tri_center_fracs {
                for &bi_treat_frac in bi_treat_fracs {
                    // Acoustic tier (ultrasound-only)
                    for &q in dims.flows {
                        for &gauge in dims.gauges {
                            for &w_ch in dims.widths {
                                for &n_segs in &[1_usize, 5, 7] {
                                    idx += 1;
                                    candidates.push(CandidateParams {
                                        idx,
                                        request: Arc::clone(request),
                                        q,
                                        gauge,
                                        d_throat: 0.0,
                                        throat_len: 0.0,
                                        w_ch,
                                        n_segs,
                                        pretri_center_frac,
                                        terminal_tri_center_frac,
                                        bi_treat_frac,
                                        treatment_actuation_mode:
                                            TreatmentActuationMode::UltrasoundOnly,
                                        vt_count: 0,
                                    });
                                }
                            }
                        }
                    }

                    if supports_venturi_materialization {
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
                                                    candidates.push(CandidateParams {
                                                        idx,
                                                        request: Arc::clone(request),
                                                        q,
                                                        gauge,
                                                        d_throat,
                                                        throat_len: d_throat * tl_factor,
                                                        w_ch,
                                                        n_segs,
                                                        pretri_center_frac,
                                                        terminal_tri_center_frac,
                                                        bi_treat_frac,
                                                        treatment_actuation_mode: TreatmentActuationMode::VenturiCavitation,
                                                        vt_count,
                                                    });
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
                                                    candidates.push(CandidateParams {
                                                        idx,
                                                        request: Arc::clone(request),
                                                        q,
                                                        gauge,
                                                        d_throat,
                                                        throat_len: d_throat * tl_factor,
                                                        w_ch,
                                                        n_segs,
                                                        pretri_center_frac,
                                                        terminal_tri_center_frac,
                                                        bi_treat_frac,
                                                        treatment_actuation_mode: TreatmentActuationMode::VenturiCavitation,
                                                        vt_count,
                                                    });
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
    }

    candidates
}
