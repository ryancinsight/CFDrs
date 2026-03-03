//! Parametric sweep and random sampling of `DesignCandidate` objects.

use super::candidate::{CrossSectionShape, DesignCandidate};
use super::topology::DesignTopology;
use crate::constraints::{
    BIFURCATION_ARM_RATIOS, CHANNEL_HEIGHTS_M, CHANNEL_HEIGHT_M, CHANNEL_WIDTHS_M, FLOW_RATES_M3_S,
    INLET_GAUGES_PA, LEUKA_CHANNEL_HEIGHT_M, LEUKA_CHANNEL_WIDTHS_M, LEUKA_FLOW_RATES_M3_S,
    SERPENTINE_BEND_RADIUS_M, SERPENTINE_SEGMENT_COUNTS, THROAT_DIAMETERS_M, THROAT_LENGTH_FACTORS,
    TREATMENT_WIDTH_MM, TRIFURCATION_CENTER_FRACS, TRIFURCATION_LEFT_FRACS, VENTURI_INLET_DIAM_M,
};
use rand::Rng;

/// Build the complete parametric sweep over all topology families.
///
/// Approximately **300,000+** candidates are generated before physics-based
/// filtering.  The exact count depends on which dimensions each topology uses:
/// topologies without a venturi skip the throat-diameter sweep; topologies
/// without a serpentine skip the segment-count sweep.  New deep-trifurcation
/// topologies (T1–T5) add an extra loop over `TRIFURCATION_CENTER_FRACS`.
#[must_use]
pub fn build_candidate_space() -> Vec<DesignCandidate> {
    // AsymmetricBifurcationSerpentine is handled separately below (arm ratio sweep).
    let topologies = [
        DesignTopology::SingleVenturi,
        DesignTopology::BifurcationVenturi,
        DesignTopology::TrifurcationVenturi,
        DesignTopology::VenturiSerpentine,
        DesignTopology::SerpentineGrid,
        DesignTopology::CellSeparationVenturi,
        DesignTopology::WbcCancerSeparationVenturi,
        DesignTopology::DoubleBifurcationVenturi,
        DesignTopology::TripleBifurcationVenturi,
        DesignTopology::DoubleTrifurcationVenturi,
        DesignTopology::BifurcationTrifurcationVenturi,
        DesignTopology::SerialDoubleVenturi,
        DesignTopology::BifurcationSerpentine,
        DesignTopology::TrifurcationSerpentine,
    ];

    // throat diameter = 0 sentinel for topologies without venturi
    let throat_none: [f64; 1] = [0.0];
    // throat_length_factor = 1 sentinel for topologies without venturi
    let tl_factor_none: [f64; 1] = [2.0];

    let mut candidates = Vec::with_capacity(2048);
    let mut idx: u32 = 0;

    for &topology in &topologies {
        let throat_iter: &[f64] = if topology.has_venturi() {
            &THROAT_DIAMETERS_M
        } else {
            &throat_none
        };
        let tl_factors: &[f64] = if topology.has_venturi() {
            &THROAT_LENGTH_FACTORS
        } else {
            &tl_factor_none
        };

        let seg_counts: &[usize] = if topology.has_serpentine() {
            &SERPENTINE_SEGMENT_COUNTS
        } else {
            &[1_usize] // unused but keeps the loop structure uniform
        };

        for &q in &FLOW_RATES_M3_S {
            for &gauge in &INLET_GAUGES_PA {
                for &d_throat in throat_iter {
                    for &tl_factor in tl_factors {
                        for &w_ch in &CHANNEL_WIDTHS_M {
                            for &h_ch in &CHANNEL_HEIGHTS_M {
                                for &n_segs in seg_counts {
                                    idx += 1;
                                    let throat_len = if d_throat > 0.0 {
                                        d_throat * tl_factor
                                    } else {
                                        0.0
                                    };
                                    let tl_tag = tl_factor as u32;
                                    let h_tag = (h_ch * 1e6) as u32; // µm
                                                                     // seg_length always spans the full treatment width
                                    let seg_len = TREATMENT_WIDTH_MM * 1e-3;

                                    let id = format!(
                                        "{:04}-{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}-n{}",
                                        idx,
                                        topology.short(),
                                        q * 6e7,
                                        gauge * 1e-3,
                                        d_throat * 1e6,
                                        tl_tag,
                                        w_ch * 1e6,
                                        h_tag,
                                        n_segs,
                                    );

                                    candidates.push(DesignCandidate {
                                        id,
                                        topology,
                                        flow_rate_m3_s: q,
                                        inlet_gauge_pa: gauge,
                                        throat_diameter_m: d_throat,
                                        inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                        throat_length_m: throat_len,
                                        channel_width_m: w_ch,
                                        channel_height_m: h_ch,
                                        serpentine_segments: n_segs,
                                        segment_length_m: seg_len,
                                        bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                        feed_hematocrit: 0.45,
                                        trifurcation_center_frac: 1.0 / 3.0,
                                        cif_pretri_center_frac: 1.0 / 3.0,
                                        cif_terminal_tri_center_frac: 1.0 / 3.0,
                                        cif_terminal_bi_treat_frac: 0.68,
                                        asymmetric_narrow_frac: 0.5,
                                        trifurcation_left_frac: 1.0 / 3.0,
                                        cross_section_shape: CrossSectionShape::Rectangular,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ── AsymmetricBifurcationSerpentine with arm ratio sweep (E1) ────────────
    for &narrow_frac in &BIFURCATION_ARM_RATIOS {
        let nar_tag = (narrow_frac * 1000.0).round() as u32;
        for &q in &FLOW_RATES_M3_S {
            for &gauge in &INLET_GAUGES_PA {
                for &w_ch in &CHANNEL_WIDTHS_M {
                    for &h_ch in &CHANNEL_HEIGHTS_M {
                        for &n_segs in &SERPENTINE_SEGMENT_COUNTS {
                            idx += 1;
                            let seg_len = TREATMENT_WIDTH_MM * 1e-3;
                            let id = format!(
                                "{:04}-AB-nar{}-q{:.0}ml-g{:.0}kPa-w{:.0}um-h{}-n{}",
                                idx,
                                nar_tag,
                                q * 6e7,
                                gauge * 1e-3,
                                w_ch * 1e6,
                                (h_ch * 1e6) as u32,
                                n_segs,
                            );
                            candidates.push(DesignCandidate {
                                id,
                                topology: DesignTopology::AsymmetricBifurcationSerpentine,
                                flow_rate_m3_s: q,
                                inlet_gauge_pa: gauge,
                                throat_diameter_m: 0.0,
                                inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                throat_length_m: 0.0,
                                channel_width_m: w_ch,
                                channel_height_m: h_ch,
                                serpentine_segments: n_segs,
                                segment_length_m: seg_len,
                                bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                feed_hematocrit: 0.45,
                                trifurcation_center_frac: 1.0 / 3.0,
                                cif_pretri_center_frac: 1.0 / 3.0,
                                cif_terminal_tri_center_frac: 1.0 / 3.0,
                                cif_terminal_bi_treat_frac: 0.68,
                                asymmetric_narrow_frac: narrow_frac,
                                trifurcation_left_frac: 1.0 / 3.0,
                                cross_section_shape: CrossSectionShape::Rectangular,
                            });
                        }
                    }
                }
            }
        }
    }

    // ── Leukapheresis topologies (micro-scale) ──────────────────────────────
    //
    // These three topologies require micro-scale dimensions (D_h < 170 µm) for
    // inertial focusing.  At millifluidic scale (CHANNEL_WIDTHS_M: 2–6 mm),
    // κ_WBC = a/D_h ≈ 0.003–0.009 — far below the 0.07 threshold — so WBC
    // recovery is zero for all millifluidic candidates.
    //
    // Use:
    //   LEUKA_CHANNEL_WIDTHS_M  (100, 200, 400 µm) → D_h = 75–96 µm
    //   LEUKA_CHANNEL_HEIGHT_M  (60 µm)
    //   LEUKA_FLOW_RATES_M3_S   (2.5, 5.0, 15 µL/s per chip = 150–900 µL/min)
    //   feed_hematocrit = 0.04  (4% diluted blood — leukapheresis pre-dilution)
    //
    // NOTE: SpiralSerpentine uses Dean flow, but the 1D Dean drag model predicts
    // outer-wall focusing only.  Real spiral devices (Nivedita 2017) achieve
    // inner-wall WBC focusing (x̃ ≈ 0.1–0.2).  This 1D limitation means
    // SpiralSerpentine will score near-zero WBC recovery here; the GA will
    // prefer ConstrictionExpansionArray and ParallelMicrochannelArray (straight
    // channels where the inertial lift model is accurate).
    // Expanded leukapheresis sweep:
    // - CE: low / mid / high cycle count
    // - SP: low / mid / high turn count
    // - PM: low / mid / high channel count
    // This materially widens scenario-space exploration for milestone analysis.
    let mut push_leuka = |topology: DesignTopology, variant_tag: String| {
        for &q in &LEUKA_FLOW_RATES_M3_S {
            for &gauge in &INLET_GAUGES_PA {
                for &w_ch in &LEUKA_CHANNEL_WIDTHS_M {
                    idx += 1;
                    let id = format!(
                        "{:04}-{}-LK-{}-q{:.0}ulm-g{:.0}kPa-w{:.0}um",
                        idx,
                        topology.short(),
                        variant_tag,
                        q * 6e10,     // m³/s → µL/min (2.5e-9 → 150, 1.5e-8 → 900)
                        gauge * 1e-3, // Pa → kPa
                        w_ch * 1e6,   // m → µm
                    );

                    candidates.push(DesignCandidate {
                        id,
                        topology,
                        flow_rate_m3_s: q,
                        inlet_gauge_pa: gauge,
                        throat_diameter_m: 0.0, // no venturi throat
                        inlet_diameter_m: VENTURI_INLET_DIAM_M,
                        throat_length_m: 0.0,
                        channel_width_m: w_ch,
                        channel_height_m: LEUKA_CHANNEL_HEIGHT_M, // 60 µm fixed
                        serpentine_segments: 1, // topology carries n_turns/n_cycles
                        segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                        bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                        feed_hematocrit: 0.04, // 4% diluted for leukapheresis
                        trifurcation_center_frac: 1.0 / 3.0,
                        cif_pretri_center_frac: 1.0 / 3.0,
                        cif_terminal_tri_center_frac: 1.0 / 3.0,
                        cif_terminal_bi_treat_frac: 0.68,
                        asymmetric_narrow_frac: 0.5,
                        trifurcation_left_frac: 1.0 / 3.0,
                        cross_section_shape: CrossSectionShape::Rectangular,
                    });
                }
            }
        }
    };

    for &n_cycles in &[6_usize, 10, 14] {
        push_leuka(
            DesignTopology::ConstrictionExpansionArray { n_cycles },
            format!("nc{n_cycles}"),
        );
    }
    for &n_turns in &[4_usize, 8, 12] {
        push_leuka(
            DesignTopology::SpiralSerpentine { n_turns },
            format!("nt{n_turns}"),
        );
    }
    for &n_channels in &[64_usize, 100, 160] {
        push_leuka(
            DesignTopology::ParallelMicrochannelArray { n_channels },
            format!("nch{n_channels}"),
        );
    }

    // ── New deep-trifurcation topologies (millifluidic scale, width-scaled) ──
    //
    // T1: TrifurcationBifurcationVenturi — fix orphaned schematics preset; symmetric
    // T2: TripleTrifurcationVenturi       — 27 outlets; swept over center_frac
    // T3: TrifurcationBifurcationBifurcation — 12 outlets; swept over center_frac
    // T4: QuadTrifurcationVenturi         — 81 outlets; swept over center_frac
    // T5: CascadeCenterTrifurcationSeparator — 1/2/3-level CCT; swept over center_frac
    //
    // All use millifluidic channel dimensions (CHANNEL_WIDTHS_M, CHANNEL_HEIGHT_M)
    // and whole-blood feed hematocrit (0.45).  T2/T3/T4/T5 additionally loop over
    // TRIFURCATION_CENTER_FRACS to vary the width-fraction sweep.

    // T1: symmetric (no center_frac variation)
    for &q in &FLOW_RATES_M3_S {
        for &gauge in &INLET_GAUGES_PA {
            for &d_throat in &THROAT_DIAMETERS_M {
                for &tl_factor in &THROAT_LENGTH_FACTORS {
                    for &w_ch in &CHANNEL_WIDTHS_M {
                        for &h_ch in &CHANNEL_HEIGHTS_M {
                            idx += 1;
                            let throat_len = d_throat * tl_factor;
                            let id = format!(
                                "{:04}-TB-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}",
                                idx,
                                q * 6e7,
                                gauge * 1e-3,
                                d_throat * 1e6,
                                tl_factor as u32,
                                w_ch * 1e6,
                                (h_ch * 1e6) as u32,
                            );
                            candidates.push(DesignCandidate {
                                id,
                                topology: DesignTopology::TrifurcationBifurcationVenturi,
                                flow_rate_m3_s: q,
                                inlet_gauge_pa: gauge,
                                throat_diameter_m: d_throat,
                                inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                throat_length_m: throat_len,
                                channel_width_m: w_ch,
                                channel_height_m: h_ch,
                                serpentine_segments: 1,
                                segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                feed_hematocrit: 0.45,
                                trifurcation_center_frac: 1.0 / 3.0,
                                cif_pretri_center_frac: 1.0 / 3.0,
                                cif_terminal_tri_center_frac: 1.0 / 3.0,
                                cif_terminal_bi_treat_frac: 0.68,
                                asymmetric_narrow_frac: 0.5,
                                trifurcation_left_frac: 1.0 / 3.0,
                                cross_section_shape: CrossSectionShape::Rectangular,
                            });
                        }
                    }
                }
            }
        }
    }

    // T2/T3/T4: swept over center_frac
    let scaled_topologies: &[DesignTopology] = &[
        DesignTopology::TripleTrifurcationVenturi,
        DesignTopology::TrifurcationBifurcationBifurcationVenturi,
        DesignTopology::QuadTrifurcationVenturi,
    ];

    for &topology in scaled_topologies {
        for &center_frac in &TRIFURCATION_CENTER_FRACS {
            let frac_tag = (center_frac * 1000.0).round() as u32;
            for &q in &FLOW_RATES_M3_S {
                for &gauge in &INLET_GAUGES_PA {
                    for &d_throat in &THROAT_DIAMETERS_M {
                        for &tl_factor in &THROAT_LENGTH_FACTORS {
                            for &w_ch in &CHANNEL_WIDTHS_M {
                                for &h_ch in &CHANNEL_HEIGHTS_M {
                                    idx += 1;
                                    let throat_len = d_throat * tl_factor;
                                    let id = format!(
                                        "{:04}-{}-cf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}",
                                        idx,
                                        topology.short(),
                                        frac_tag,
                                        q * 6e7,
                                        gauge * 1e-3,
                                        d_throat * 1e6,
                                        tl_factor as u32,
                                        w_ch * 1e6,
                                        (h_ch * 1e6) as u32,
                                    );
                                    candidates.push(DesignCandidate {
                                        id,
                                        topology,
                                        flow_rate_m3_s: q,
                                        inlet_gauge_pa: gauge,
                                        throat_diameter_m: d_throat,
                                        inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                        throat_length_m: throat_len,
                                        channel_width_m: w_ch,
                                        channel_height_m: h_ch,
                                        serpentine_segments: 1,
                                        segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                        bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                        feed_hematocrit: 0.45,
                                        trifurcation_center_frac: center_frac,
                                        cif_pretri_center_frac: 1.0 / 3.0,
                                        cif_terminal_tri_center_frac: 1.0 / 3.0,
                                        cif_terminal_bi_treat_frac: 0.68,
                                        asymmetric_narrow_frac: 0.5,
                                        trifurcation_left_frac: 1.0 / 3.0,
                                        cross_section_shape: CrossSectionShape::Rectangular,
                                    });
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Dedicated CCT/CIF deep sweep (bounded runtime grid).
    const CCT_CIF_FLOWS: [f64; 3] = [1.667e-6, 2.500e-6, 3.333e-6];
    const CCT_CIF_GAUGES: [f64; 3] = [100_000.0, 200_000.0, 300_000.0];
    const CCT_CIF_THROATS: [f64; 3] = [50e-6, 75e-6, 100e-6];
    const CCT_CIF_WIDTHS: [f64; 3] = [2.0e-3, 4.0e-3, 6.0e-3];
    const CIF_CENTER_FRACS: [f64; 3] = [0.33, 0.45, 0.55];
    const CIF_BI_TREAT_FRACS: [f64; 3] = [0.60, 0.68, 0.76];

    // CCT grid (n_levels 1–3).
    for &n_levels in &[1u8, 2, 3] {
        for &center_frac in &CIF_CENTER_FRACS {
            let cf_tag = (center_frac * 1000.0).round() as u32;
            for &q in &CCT_CIF_FLOWS {
                for &gauge in &CCT_CIF_GAUGES {
                    for &d_throat in &CCT_CIF_THROATS {
                        for &tl_factor in &THROAT_LENGTH_FACTORS {
                            for &w_ch in &CCT_CIF_WIDTHS {
                                idx += 1;
                                let throat_len = d_throat * tl_factor;
                                let id = format!(
                                    "{:04}-CCT-lv{}-cf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um",
                                    idx,
                                    n_levels,
                                    cf_tag,
                                    q * 6e7,
                                    gauge * 1e-3,
                                    d_throat * 1e6,
                                    tl_factor as u32,
                                    w_ch * 1e6,
                                );
                                candidates.push(DesignCandidate {
                                    id,
                                    topology: DesignTopology::CascadeCenterTrifurcationSeparator {
                                        n_levels,
                                    },
                                    flow_rate_m3_s: q,
                                    inlet_gauge_pa: gauge,
                                    throat_diameter_m: d_throat,
                                    inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                    throat_length_m: throat_len,
                                    channel_width_m: w_ch,
                                    channel_height_m: CHANNEL_HEIGHT_M,
                                    serpentine_segments: 1,
                                    segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                    bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                    feed_hematocrit: 0.45,
                                    trifurcation_center_frac: center_frac,
                                    cif_pretri_center_frac: 1.0 / 3.0,
                                    cif_terminal_tri_center_frac: 1.0 / 3.0,
                                    cif_terminal_bi_treat_frac: 0.68,
                                    asymmetric_narrow_frac: 0.5,
                                    trifurcation_left_frac: 1.0 / 3.0,
                                    cross_section_shape: CrossSectionShape::Rectangular,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    // CIF grid.
    for &n_pretri in &[1u8, 2, 3] {
        for &pretri_center_frac in &CIF_CENTER_FRACS {
            for &terminal_tri_center_frac in &CIF_CENTER_FRACS {
                for &bi_treat_frac in &CIF_BI_TREAT_FRACS {
                    let pcf_tag = (pretri_center_frac * 1000.0).round() as u32;
                    let tcf_tag = (terminal_tri_center_frac * 1000.0).round() as u32;
                    let btf_tag = (bi_treat_frac * 1000.0).round() as u32;
                    for &q in &CCT_CIF_FLOWS {
                        for &gauge in &CCT_CIF_GAUGES {
                            for &d_throat in &CCT_CIF_THROATS {
                                for &tl_factor in &THROAT_LENGTH_FACTORS {
                                    for &w_ch in &CCT_CIF_WIDTHS {
                                        idx += 1;
                                        let throat_len = d_throat * tl_factor;
                                        let id = format!(
                                            "{:04}-CIF-pt{}-pcf{}-tcf{}-btf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um",
                                            idx,
                                            n_pretri,
                                            pcf_tag,
                                            tcf_tag,
                                            btf_tag,
                                            q * 6e7,
                                            gauge * 1e-3,
                                            d_throat * 1e6,
                                            tl_factor as u32,
                                            w_ch * 1e6,
                                        );
                                        candidates.push(DesignCandidate {
                                            id,
                                            topology:
                                                DesignTopology::IncrementalFiltrationTriBiSeparator {
                                                    n_pretri,
                                                },
                                            flow_rate_m3_s: q,
                                            inlet_gauge_pa: gauge,
                                            throat_diameter_m: d_throat,
                                            inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                            throat_length_m: throat_len,
                                            channel_width_m: w_ch,
                                            channel_height_m: CHANNEL_HEIGHT_M,
                                            serpentine_segments: 1,
                                            segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                            bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                            feed_hematocrit: 0.45,
                                            trifurcation_center_frac: pretri_center_frac,
                                            cif_pretri_center_frac: pretri_center_frac,
                                            cif_terminal_tri_center_frac: terminal_tri_center_frac,
                                            cif_terminal_bi_treat_frac: bi_treat_frac,
                                            asymmetric_narrow_frac: 0.5,
                                            trifurcation_left_frac: 1.0 / 3.0,
                                            cross_section_shape: CrossSectionShape::Rectangular,
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

    // Oncology-focused CIF refinement sweep.
    //
    // Purpose: add a focused high-performing band around SDT-selective operating
    // points identified by prior runs (high cancer-targeted cavitation with low
    // RBC venturi exposure) while still covering lower-flow hemolysis-safe space.
    const ONCO_FLOWS: [f64; 6] = [1.667e-6, 2.500e-6, 3.333e-6, 4.167e-6, 5.000e-6, 5.833e-6]; // 100..350 mL/min
    const ONCO_GAUGES: [f64; 4] = [200_000.0, 300_000.0, 400_000.0, 500_000.0];
    const ONCO_THROATS: [f64; 4] = [40e-6, 45e-6, 50e-6, 55e-6];
    const ONCO_WIDTHS: [f64; 2] = [5.0e-3, 6.0e-3];
    const ONCO_PRETRI_FRACS: [f64; 3] = [0.50, 0.54, 0.58];
    const ONCO_TERM_TRI_FRACS: [f64; 4] = [0.50, 0.53, 0.56, 0.60];
    const ONCO_BI_TREAT_FRACS: [f64; 4] = [0.76, 0.80, 0.84, 0.85];
    for &n_pretri in &[2u8, 3] {
        for &pretri_center_frac in &ONCO_PRETRI_FRACS {
            for &terminal_tri_center_frac in &ONCO_TERM_TRI_FRACS {
                for &bi_treat_frac in &ONCO_BI_TREAT_FRACS {
                    let pcf_tag = (pretri_center_frac * 1000.0).round() as u32;
                    let tcf_tag = (terminal_tri_center_frac * 1000.0).round() as u32;
                    let btf_tag = (bi_treat_frac * 1000.0).round() as u32;
                    for &q in &ONCO_FLOWS {
                        for &gauge in &ONCO_GAUGES {
                            for &d_throat in &ONCO_THROATS {
                                for &tl_factor in &THROAT_LENGTH_FACTORS {
                                    for &w_ch in &ONCO_WIDTHS {
                                        idx += 1;
                                        let throat_len = d_throat * tl_factor;
                                        let id = format!(
                                            "{:04}-CIFX-pt{}-pcf{}-tcf{}-btf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um",
                                            idx,
                                            n_pretri,
                                            pcf_tag,
                                            tcf_tag,
                                            btf_tag,
                                            q * 6e7,
                                            gauge * 1e-3,
                                            d_throat * 1e6,
                                            tl_factor as u32,
                                            w_ch * 1e6,
                                        );
                                        candidates.push(DesignCandidate {
                                            id,
                                            topology:
                                                DesignTopology::IncrementalFiltrationTriBiSeparator {
                                                    n_pretri,
                                                },
                                            flow_rate_m3_s: q,
                                            inlet_gauge_pa: gauge,
                                            throat_diameter_m: d_throat,
                                            inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                            throat_length_m: throat_len,
                                            channel_width_m: w_ch,
                                            channel_height_m: CHANNEL_HEIGHT_M,
                                            serpentine_segments: 1,
                                            segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                            bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                            feed_hematocrit: 0.45,
                                            trifurcation_center_frac: pretri_center_frac,
                                            cif_pretri_center_frac: pretri_center_frac,
                                            cif_terminal_tri_center_frac: terminal_tri_center_frac,
                                            cif_terminal_bi_treat_frac: bi_treat_frac,
                                            asymmetric_narrow_frac: 0.5,
                                            trifurcation_left_frac: 1.0 / 3.0,
                                            cross_section_shape: CrossSectionShape::Rectangular,
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

    // ── AsymmetricTrifurcationVenturi grid (F1) ───────────────────────────────
    // Sweep center_frac × left_frac (valid pairs only: center + left ≤ 0.85).
    for &center_frac in &TRIFURCATION_CENTER_FRACS {
        for &left_frac in &TRIFURCATION_LEFT_FRACS {
            if center_frac + left_frac > 0.85 {
                continue; // right arm would be < 0.0
            }
            let cf_tag = (center_frac * 1000.0).round() as u32;
            let lf_tag = (left_frac * 1000.0).round() as u32;
            for &q in &CCT_CIF_FLOWS {
                for &gauge in &CCT_CIF_GAUGES {
                    for &d_throat in &CCT_CIF_THROATS {
                        for &tl_factor in &THROAT_LENGTH_FACTORS {
                            for &w_ch in &CCT_CIF_WIDTHS {
                                idx += 1;
                                let throat_len = d_throat * tl_factor;
                                let id = format!(
                                    "{:04}-ATV-cf{}-lf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um",
                                    idx,
                                    cf_tag,
                                    lf_tag,
                                    q * 6e7,
                                    gauge * 1e-3,
                                    d_throat * 1e6,
                                    tl_factor as u32,
                                    w_ch * 1e6,
                                );
                                candidates.push(DesignCandidate {
                                    id,
                                    topology: DesignTopology::AsymmetricTrifurcationVenturi,
                                    flow_rate_m3_s: q,
                                    inlet_gauge_pa: gauge,
                                    throat_diameter_m: d_throat,
                                    inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                    throat_length_m: throat_len,
                                    channel_width_m: w_ch,
                                    channel_height_m: CHANNEL_HEIGHT_M,
                                    serpentine_segments: 1,
                                    segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                    bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                    feed_hematocrit: 0.45,
                                    trifurcation_center_frac: center_frac,
                                    cif_pretri_center_frac: 1.0 / 3.0,
                                    cif_terminal_tri_center_frac: 1.0 / 3.0,
                                    cif_terminal_bi_treat_frac: 0.68,
                                    asymmetric_narrow_frac: 0.5,
                                    trifurcation_left_frac: left_frac,
                                    cross_section_shape: CrossSectionShape::Rectangular,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    // ── TriBiTriSelectiveVenturi grid (F2c) ───────────────────────────────────
    // Sweep trifurcation_center_frac × cif_terminal_bi_treat_frac.
    for &center_frac in &CIF_CENTER_FRACS {
        for &bi_treat_frac in &CIF_BI_TREAT_FRACS {
            let cf_tag = (center_frac * 1000.0).round() as u32;
            let btf_tag = (bi_treat_frac * 1000.0).round() as u32;
            for &q in &CCT_CIF_FLOWS {
                for &gauge in &CCT_CIF_GAUGES {
                    for &d_throat in &CCT_CIF_THROATS {
                        for &tl_factor in &THROAT_LENGTH_FACTORS {
                            for &w_ch in &CCT_CIF_WIDTHS {
                                idx += 1;
                                let throat_len = d_throat * tl_factor;
                                let id = format!(
                                    "{:04}-TBT-cf{}-btf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um",
                                    idx,
                                    cf_tag,
                                    btf_tag,
                                    q * 6e7,
                                    gauge * 1e-3,
                                    d_throat * 1e6,
                                    tl_factor as u32,
                                    w_ch * 1e6,
                                );
                                candidates.push(DesignCandidate {
                                    id,
                                    topology: DesignTopology::TriBiTriSelectiveVenturi,
                                    flow_rate_m3_s: q,
                                    inlet_gauge_pa: gauge,
                                    throat_diameter_m: d_throat,
                                    inlet_diameter_m: VENTURI_INLET_DIAM_M,
                                    throat_length_m: throat_len,
                                    channel_width_m: w_ch,
                                    channel_height_m: CHANNEL_HEIGHT_M,
                                    serpentine_segments: 1,
                                    segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                    bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                    feed_hematocrit: 0.45,
                                    trifurcation_center_frac: center_frac,
                                    cif_pretri_center_frac: 1.0 / 3.0,
                                    cif_terminal_tri_center_frac: 1.0 / 3.0,
                                    cif_terminal_bi_treat_frac: bi_treat_frac,
                                    asymmetric_narrow_frac: 0.5,
                                    trifurcation_left_frac: 1.0 / 3.0,
                                    cross_section_shape: CrossSectionShape::Rectangular,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    candidates
}

/// Generate `n` random [`DesignCandidate`] objects by sampling uniformly from
/// the continuous design space.
///
/// Unlike [`build_candidate_space`], this does **not** use a fixed parameter
/// grid — each candidate is independently sampled, giving broad, non-repeating
/// coverage of the topology parameter space.  Useful for:
///
/// - Pre-populating a GA with a diverse initial population.
/// - Running a random-search baseline against the genetic algorithm.
/// - Generating large datasets for analysis or machine-learning.
///
/// # Parameter sampling ranges
///
/// | Parameter | Range | Units |
/// |-----------|-------|-------|
/// | topology  | uniform 0–23 | — |
/// | flow_rate | \[`FLOW_RATES_M3_S` min, max\] | m³/s |
/// | inlet_gauge_pa | \[50 000, 300 000\] | Pa |
/// | throat_diameter | \[30 µm, 500 µm\] (venturi only) | m |
/// | channel_width | \[500 µm, 6 mm\] millifluidic / \[100 µm, 400 µm\] leuka | m |
/// | trifurcation_center_frac | \[0.25, 0.65\] | — |
pub fn sample_random_candidates(n: usize, rng: &mut impl Rng) -> Vec<DesignCandidate> {
    // All topology templates (same order as ALL_EVO_TOPOLOGIES in evo.rs)
    const N_TOPOS: usize = 27;
    let topo_templates: [DesignTopology; N_TOPOS] = [
        DesignTopology::SingleVenturi,
        DesignTopology::BifurcationVenturi,
        DesignTopology::TrifurcationVenturi,
        DesignTopology::VenturiSerpentine,
        DesignTopology::SerpentineGrid,
        DesignTopology::CellSeparationVenturi,
        DesignTopology::WbcCancerSeparationVenturi,
        DesignTopology::DoubleBifurcationVenturi,
        DesignTopology::TripleBifurcationVenturi,
        DesignTopology::DoubleTrifurcationVenturi,
        DesignTopology::BifurcationTrifurcationVenturi,
        DesignTopology::SerialDoubleVenturi,
        DesignTopology::BifurcationSerpentine,
        DesignTopology::TrifurcationSerpentine,
        DesignTopology::AsymmetricBifurcationSerpentine,
        DesignTopology::ConstrictionExpansionArray { n_cycles: 10 },
        DesignTopology::SpiralSerpentine { n_turns: 8 },
        DesignTopology::ParallelMicrochannelArray { n_channels: 100 },
        DesignTopology::TrifurcationBifurcationVenturi,
        DesignTopology::TripleTrifurcationVenturi,
        DesignTopology::TrifurcationBifurcationBifurcationVenturi,
        DesignTopology::QuadTrifurcationVenturi,
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels: 2 },
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 2 },
        DesignTopology::AsymmetricTrifurcationVenturi,
        DesignTopology::TriBiTriSelectiveVenturi,
        DesignTopology::AdaptiveTree {
            levels: 2,
            split_types: 0b0101,
        },
    ];
    (0..n)
        .map(|i| {
            // Pick topology uniformly
            let topo_base = topo_templates[rng.gen_range(0..N_TOPOS)];
            // Instantiate variant-specific discrete parameters
            let topology = match topo_base {
                DesignTopology::ConstrictionExpansionArray { .. } => {
                    DesignTopology::ConstrictionExpansionArray {
                        n_cycles: rng.gen_range(2..=20),
                    }
                }
                DesignTopology::SpiralSerpentine { .. } => DesignTopology::SpiralSerpentine {
                    n_turns: rng.gen_range(2..=20),
                },
                DesignTopology::ParallelMicrochannelArray { .. } => {
                    DesignTopology::ParallelMicrochannelArray {
                        n_channels: rng.gen_range(10..=200),
                    }
                }
                DesignTopology::CascadeCenterTrifurcationSeparator { .. } => {
                    DesignTopology::CascadeCenterTrifurcationSeparator {
                        n_levels: rng.gen_range(1u8..=3u8),
                    }
                }
                DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => {
                    DesignTopology::IncrementalFiltrationTriBiSeparator {
                        n_pretri: rng.gen_range(1u8..=3u8),
                    }
                }
                DesignTopology::AdaptiveTree { .. } => DesignTopology::AdaptiveTree {
                    levels: rng.gen_range(1u8..=3u8),
                    split_types: rng.gen::<u8>() & 0x0F,
                },
                other => other,
            };
            let is_leuka = matches!(
                topology,
                DesignTopology::ConstrictionExpansionArray { .. }
                    | DesignTopology::SpiralSerpentine { .. }
                    | DesignTopology::ParallelMicrochannelArray { .. }
            );
            // Flow rate — log-linear between array min and max
            let q = if is_leuka {
                let lo = LEUKA_FLOW_RATES_M3_S[0];
                let hi = LEUKA_FLOW_RATES_M3_S[LEUKA_FLOW_RATES_M3_S.len() - 1];
                lo * (hi / lo).powf(rng.gen::<f64>())
            } else {
                let lo = FLOW_RATES_M3_S[0];
                let hi = FLOW_RATES_M3_S[FLOW_RATES_M3_S.len() - 1];
                lo * (hi / lo).powf(rng.gen::<f64>())
            };
            // Inlet gauge pressure: 50–300 kPa
            let gauge = 50_000.0 + rng.gen::<f64>() * 250_000.0;
            // Throat diameter: 30–500 µm (venturi only)
            let d_throat = if topology.has_venturi() {
                30e-6 + rng.gen::<f64>() * 470e-6
            } else {
                0.0
            };
            // Throat length: random factor [1.5, 15.0] × diameter
            let tl_factor = if topology.has_venturi() {
                1.5 + rng.gen::<f64>() * 13.5
            } else {
                2.0
            };
            let throat_len = d_throat * tl_factor;
            // Channel dimensions
            let (w_ch, h_ch) = if is_leuka {
                let w = LEUKA_CHANNEL_WIDTHS_M[0]
                    + rng.gen::<f64>()
                        * (LEUKA_CHANNEL_WIDTHS_M[LEUKA_CHANNEL_WIDTHS_M.len() - 1]
                            - LEUKA_CHANNEL_WIDTHS_M[0]);
                (w, LEUKA_CHANNEL_HEIGHT_M)
            } else {
                // 500 µm – 6 mm wide; height log-linear [0.3 mm, 3.0 mm]
                let h = 0.3e-3 * (10.0_f64.powf(rng.gen::<f64>()));
                (500e-6 + rng.gen::<f64>() * 5.5e-3, h)
            };
            // Serpentine segments: 2–12
            let n_segs = if topology.has_serpentine() {
                rng.gen_range(2..=12)
            } else {
                1
            };
            let seg_len = (0.5 + rng.gen::<f64>()) * TREATMENT_WIDTH_MM * 1e-3;
            let bend_r = (0.05 + rng.gen::<f64>() * 0.20) * seg_len;
            // Trifurcation center fraction: 0.25–0.65
            let trifurcation_center_frac = if matches!(
                topology,
                DesignTopology::TrifurcationVenturi
                    | DesignTopology::DoubleTrifurcationVenturi
                    | DesignTopology::BifurcationTrifurcationVenturi
                    | DesignTopology::TrifurcationSerpentine
                    | DesignTopology::TrifurcationBifurcationVenturi
                    | DesignTopology::TripleTrifurcationVenturi
                    | DesignTopology::TrifurcationBifurcationBifurcationVenturi
                    | DesignTopology::QuadTrifurcationVenturi
                    | DesignTopology::CascadeCenterTrifurcationSeparator { .. }
                    | DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
                    | DesignTopology::AsymmetricTrifurcationVenturi
                    | DesignTopology::TriBiTriSelectiveVenturi
            ) {
                0.25 + rng.gen::<f64>() * 0.40
            } else {
                1.0 / 3.0
            };
            let (cif_pretri_center_frac, cif_terminal_tri_center_frac, cif_terminal_bi_treat_frac) =
                if matches!(
                    topology,
                    DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
                        | DesignTopology::TriBiTriSelectiveVenturi
                ) {
                    (
                        trifurcation_center_frac,
                        0.25 + rng.gen::<f64>() * 0.40,
                        0.50 + rng.gen::<f64>() * 0.35,
                    )
                } else {
                    (1.0 / 3.0, 1.0 / 3.0, 0.68)
                };
            let asymmetric_narrow_frac =
                if matches!(topology, DesignTopology::AsymmetricBifurcationSerpentine) {
                    0.20 + rng.gen::<f64>() * 0.50 // [0.20, 0.70]
                } else {
                    0.5
                };
            let trifurcation_left_frac =
                if matches!(topology, DesignTopology::AsymmetricTrifurcationVenturi) {
                    // ensure center + left ≤ 0.85
                    let max_left = (0.85 - trifurcation_center_frac).max(0.08);
                    0.08 + rng.gen::<f64>() * (max_left - 0.08)
                } else {
                    1.0 / 3.0
                };
            let feed_hematocrit = if is_leuka { 0.04 } else { 0.45 };
            let id = if let DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } =
                topology
            {
                format!(
                    "RND{:05}-CIF-pt{}-pcf{}-tcf{}-btf{}-q{:.0}-g{:.0}-d{:.0}-w{:.0}",
                    i,
                    n_pretri,
                    (cif_pretri_center_frac * 1000.0).round() as u32,
                    (cif_terminal_tri_center_frac * 1000.0).round() as u32,
                    (cif_terminal_bi_treat_frac * 1000.0).round() as u32,
                    q * 6e7,
                    gauge * 1e-3,
                    d_throat * 1e6,
                    w_ch * 1e6,
                )
            } else {
                format!(
                    "RND{:05}-{}-q{:.0}-g{:.0}-d{:.0}-w{:.0}",
                    i,
                    topology.short(),
                    q * 6e7,
                    gauge * 1e-3,
                    d_throat * 1e6,
                    w_ch * 1e6,
                )
            };
            DesignCandidate {
                id,
                topology,
                flow_rate_m3_s: q,
                inlet_gauge_pa: gauge,
                throat_diameter_m: d_throat,
                inlet_diameter_m: VENTURI_INLET_DIAM_M,
                throat_length_m: throat_len,
                channel_width_m: w_ch,
                channel_height_m: h_ch,
                serpentine_segments: n_segs,
                segment_length_m: seg_len,
                bend_radius_m: bend_r,
                feed_hematocrit,
                trifurcation_center_frac,
                cif_pretri_center_frac,
                cif_terminal_tri_center_frac,
                cif_terminal_bi_treat_frac,
                asymmetric_narrow_frac,
                trifurcation_left_frac,
                cross_section_shape: CrossSectionShape::Rectangular,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oncology_focused_cif_refinement_band_present() {
        let candidates = build_candidate_space();
        let has_focused = candidates.iter().any(|c| {
            c.id.contains("-CIFX-")
                && matches!(
                    c.topology,
                    DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
                )
                && c.inlet_gauge_pa >= 400_000.0
                && c.flow_rate_m3_s >= 4.167e-6
                && c.throat_diameter_m <= 55e-6
                && c.channel_width_m >= 5.0e-3
        });
        assert!(
            has_focused,
            "candidate space should include oncology-focused CIF refinement points"
        );
    }
}
