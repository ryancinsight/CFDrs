//! Genome → [`DesignCandidate`] decoding.

use crate::{
    constraints::{
        FLOW_RATES_M3_S, INLET_GAUGES_PA, LEUKA_CHANNEL_HEIGHT_M, LEUKA_CHANNEL_WIDTHS_M,
        LEUKA_FLOW_RATES_M3_S, PLATE_HEIGHT_MM, THROAT_DIAMETERS_M, TREATMENT_WIDTH_MM,
        VENTURI_INLET_DIAM_M,
    },
    design::{CrossSectionShape, DesignCandidate, DesignTopology, TreatmentZoneMode},
};

use super::{MillifluidicGenome, ALL_EVO_TOPOLOGIES, MIN_ADAPTIVE_CH_M};

/// Decode a normalised genome into a [`DesignCandidate`].
///
/// Continuous gene 0 is rounded to select the topology.
/// Gene 5 is rounded to select the segment count.
/// Gene 8 encodes topology-specific discrete parameters for leukapheresis and
/// staged-separation variants.
#[must_use]
pub fn decode_genome(g: &MillifluidicGenome, id_prefix: &str) -> DesignCandidate {
    let genes = &g.genes;

    // Gene 0: topology index
    let topo_idx = (genes[0] * (ALL_EVO_TOPOLOGIES.len() as f64 - 1e-9))
        .floor()
        .clamp(0.0, (ALL_EVO_TOPOLOGIES.len() - 1) as f64) as usize;

    let topo_template = ALL_EVO_TOPOLOGIES[topo_idx];

    let is_leuka = matches!(
        topo_template,
        DesignTopology::ConstrictionExpansionArray { .. }
            | DesignTopology::SpiralSerpentine { .. }
            | DesignTopology::ParallelMicrochannelArray { .. }
    );

    // Gene 4: channel width — micro-scale for leuka, millifluidic for others.
    let w_ch = if is_leuka {
        LEUKA_CHANNEL_WIDTHS_M[0]
            + genes[4]
                * (LEUKA_CHANNEL_WIDTHS_M[LEUKA_CHANNEL_WIDTHS_M.len() - 1]
                    - LEUKA_CHANNEL_WIDTHS_M[0])
    } else {
        2.0e-3 + genes[4] * 4.0e-3
    };

    // Decode topology variant from gene 8 + genes 9–12 when needed.
    let topology = match topo_template {
        DesignTopology::AdaptiveTree { .. } => {
            let max_depth = {
                let mut d = 0u8;
                let mut w = w_ch;
                while d < 4 {
                    w /= 2.0;
                    if w < MIN_ADAPTIVE_CH_M {
                        break;
                    }
                    d += 1;
                }
                d
            };
            let depth = if max_depth == 0 {
                0u8
            } else {
                ((genes[8] * (f64::from(max_depth) + 1.0)).floor() as u8).min(max_depth)
            };
            let mut split_types = 0u8;
            for i in 0..4usize {
                if genes[9 + i] > 0.5 {
                    split_types |= 1u8 << i;
                }
            }
            let mask = if depth == 0 {
                0u8
            } else {
                (1u8 << depth).wrapping_sub(1)
            };
            split_types &= mask;
            DesignTopology::AdaptiveTree {
                levels: depth,
                split_types,
            }
        }
        DesignTopology::ConstrictionExpansionArray { .. } => {
            let n_cycles = (2.0 + genes[8] * 18.0).round() as usize;
            DesignTopology::ConstrictionExpansionArray { n_cycles }
        }
        DesignTopology::SpiralSerpentine { .. } => {
            let n_turns = (2.0 + genes[8] * 18.0).round() as usize;
            DesignTopology::SpiralSerpentine { n_turns }
        }
        DesignTopology::ParallelMicrochannelArray { .. } => {
            let raw_n = (10.0 + genes[8] * 490.0).round() as usize;
            let pitch_m = w_ch * 2.0;
            let max_n = ((PLATE_HEIGHT_MM - 10.0) * 1e-3 / pitch_m).floor() as usize;
            let n_channels = raw_n.min(max_n).max(1);
            DesignTopology::ParallelMicrochannelArray { n_channels }
        }
        DesignTopology::CascadeCenterTrifurcationSeparator { .. } => {
            let n_levels = (1.0 + genes[8] * 2.0).round().clamp(1.0, 3.0) as u8;
            DesignTopology::CascadeCenterTrifurcationSeparator { n_levels }
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => {
            let n_pretri = (1.0 + genes[8] * 2.0).round().clamp(1.0, 3.0) as u8;
            DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri }
        }
        DesignTopology::DoubleTrifurcationCIFVenturi { .. } => {
            let center_throat_count = (1.0 + (genes[8] * 3.999).floor()) as u8;
            DesignTopology::DoubleTrifurcationCIFVenturi {
                center_throat_count,
            }
        }
        other => other,
    };

    // Feed hematocrit: leukapheresis topologies use diluted blood.
    let feed_hematocrit = match topology {
        DesignTopology::ConstrictionExpansionArray { .. }
        | DesignTopology::SpiralSerpentine { .. }
        | DesignTopology::ParallelMicrochannelArray { .. } => 0.04,
        _ => 0.45,
    };

    // Gene 1: flow rate (log-linear).
    let q = if is_leuka {
        let q_min = LEUKA_FLOW_RATES_M3_S[0];
        let q_max = LEUKA_FLOW_RATES_M3_S[LEUKA_FLOW_RATES_M3_S.len() - 1];
        q_min * (q_max / q_min).powf(genes[1])
    } else {
        let q_min = FLOW_RATES_M3_S[0];
        let q_max = FLOW_RATES_M3_S[FLOW_RATES_M3_S.len() - 1];
        q_min * (q_max / q_min).powf(genes[1])
    };

    // Gene 2: inlet gauge pressure (linear)
    let p_min = INLET_GAUGES_PA[0];
    let p_max = INLET_GAUGES_PA[INLET_GAUGES_PA.len() - 1];
    let gauge = p_min + genes[2] * (p_max - p_min);

    // Gene 3: throat diameter
    let d_min = THROAT_DIAMETERS_M[0];
    let d_max = THROAT_DIAMETERS_M[THROAT_DIAMETERS_M.len() - 1];
    let d_throat = if topology.has_venturi() {
        d_min + genes[3] * (d_max - d_min)
    } else {
        0.0
    };

    // Gene 5: serpentine segment count (2 – 12)
    let n_segs = if topology.has_serpentine() {
        2 + (genes[5] * 10.0).round() as usize
    } else {
        1
    };

    // Gene 6: segment length
    let seg_len = (0.5 + genes[6]) * TREATMENT_WIDTH_MM * 1e-3;

    // Gene 7: bend radius
    let bend_r = (0.05 + genes[7] * 0.20) * seg_len;

    // Gene 17: throat length factor [1.5, 15.0] × diameter.
    let tl_factor = if topology.has_venturi() {
        1.5 + genes[17] * 13.5
    } else {
        2.0
    };
    let throat_len = if d_throat > 0.0 {
        d_throat * tl_factor
    } else {
        0.0
    };

    // Gene 13: trifurcation center-arm width fraction.
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
            | DesignTopology::DoubleTrifurcationCIFVenturi { .. }
            | DesignTopology::AsymmetricTrifurcationVenturi
            | DesignTopology::TriBiTriSelectiveVenturi
    ) {
        0.25 + genes[13] * 0.40
    } else {
        1.0 / 3.0
    };
    // Genes 14-15: staged CIF / TBT controls (terminal trifurcation + terminal bifurcation).
    let (cif_pretri_center_frac, cif_terminal_tri_center_frac, cif_terminal_bi_treat_frac) = if matches!(
        topology,
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
            | DesignTopology::DoubleTrifurcationCIFVenturi { .. }
            | DesignTopology::TriBiTriSelectiveVenturi
    ) {
        (
            trifurcation_center_frac,
            0.25 + genes[14] * 0.40,
            0.50 + genes[15] * 0.35,
        )
    } else {
        (1.0 / 3.0, 1.0 / 3.0, 0.68)
    };

    // Gene 16: asymmetric bifurcation narrow-arm fraction [0.20, 0.70].
    let asymmetric_narrow_frac =
        if matches!(topology, DesignTopology::AsymmetricBifurcationSerpentine) {
            0.20 + genes[16] * 0.50
        } else {
            0.5
        };

    // Gene 18: channel height — millifluidic log-linear [0.3 mm, 3.0 mm].
    let h_ch = if is_leuka {
        LEUKA_CHANNEL_HEIGHT_M
    } else {
        // log-linear: 0.3e-3 × 10^gene[18]  → [0.3 mm, 3.0 mm]
        0.3e-3 * (10.0_f64.powf(genes[18]))
    };

    // Gene 19: asymmetric trifurcation left-arm fraction [0.08, clamped by center_frac].
    let trifurcation_left_frac =
        if matches!(topology, DesignTopology::AsymmetricTrifurcationVenturi) {
            let max_left = (0.85 - trifurcation_center_frac).max(0.08);
            (0.08 + genes[19] * (max_left - 0.08)).clamp(0.08, max_left)
        } else {
            1.0 / 3.0
        };

    // Gene 12: centerline venturi throat count for selective centerline topologies.
    // Reuses one AdaptiveTree split-type gene for non-Adaptive topologies.
    let centerline_venturi_throat_count = match topology {
        DesignTopology::DoubleTrifurcationCIFVenturi {
            center_throat_count,
        } => center_throat_count,
        DesignTopology::CascadeCenterTrifurcationSeparator { .. }
        | DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
        | DesignTopology::AsymmetricTrifurcationVenturi
        | DesignTopology::TriBiTriSelectiveVenturi
        | DesignTopology::CellSeparationVenturi
        | DesignTopology::WbcCancerSeparationVenturi => {
            (1.0 + (genes[12] * 3.999).floor()) as u8
        }
        _ => 1,
    };

    let treatment_zone_mode = if topology.has_venturi() {
        TreatmentZoneMode::VenturiThroats
    } else {
        TreatmentZoneMode::UltrasoundOnly
    };

    // Encode topology identity tag
    let topo_tag = match topology {
        DesignTopology::AdaptiveTree {
            levels,
            split_types,
        } => {
            format!("AT-d{}-s{:04b}", levels, split_types & 0x0F)
        }
        DesignTopology::ConstrictionExpansionArray { n_cycles } => {
            format!("CE-c{n_cycles}")
        }
        DesignTopology::SpiralSerpentine { n_turns } => {
            format!("SP-t{n_turns}")
        }
        DesignTopology::ParallelMicrochannelArray { n_channels } => {
            format!("PM-n{n_channels}")
        }
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
            format!("CCT-lv{n_levels}-vt{centerline_venturi_throat_count}")
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            format!(
                "CIF-pt{}-pcf{}-tcf{}-btf{}-vt{}",
                n_pretri,
                (cif_pretri_center_frac * 1000.0).round() as u32,
                (cif_terminal_tri_center_frac * 1000.0).round() as u32,
                (cif_terminal_bi_treat_frac * 1000.0).round() as u32,
                centerline_venturi_throat_count,
            )
        }
        DesignTopology::DoubleTrifurcationCIFVenturi {
            center_throat_count,
        } => {
            format!(
                "DTCV-s1cf{}-s2cf{}-ct{}",
                (cif_pretri_center_frac * 1000.0).round() as u32,
                (cif_terminal_tri_center_frac * 1000.0).round() as u32,
                center_throat_count,
            )
        }
        DesignTopology::AsymmetricTrifurcationVenturi => {
            format!(
                "ATV-cf{}-lf{}-vt{}",
                (trifurcation_center_frac * 1000.0).round() as u32,
                (trifurcation_left_frac * 1000.0).round() as u32,
                centerline_venturi_throat_count,
            )
        }
        DesignTopology::TriBiTriSelectiveVenturi => {
            format!(
                "TBT-cf{}-btf{}-vt{}",
                (trifurcation_center_frac * 1000.0).round() as u32,
                (cif_terminal_bi_treat_frac * 1000.0).round() as u32,
                centerline_venturi_throat_count,
            )
        }
        _ => format!("t{topo_idx}"),
    };

    let id = format!(
        "{}-EVO-{}-q{:.0}-g{:.0}-d{:.0}-w{:.0}-h{:.0}-n{}",
        id_prefix,
        topo_tag,
        q * 6e7,
        gauge * 1e-3,
        d_throat * 1e6,
        w_ch * 1e6,
        h_ch * 1e6,
        n_segs,
    );

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
        treatment_zone_mode,
        centerline_venturi_throat_count,
    }
}

/// Encode a [`DesignCandidate`] into a [`MillifluidicGenome`] suitable for
/// warm-starting the genetic algorithm from known-good parametric candidates.
///
/// This is the approximate inverse of [`decode_genome`].  The encoding is
/// intentionally lossy for discrete parameters (topology variant, segment count)
/// since the genome is continuous; when decoded again it will recover the same
/// discrete values via rounding / floor operations.
#[must_use]
pub fn candidate_to_genome(c: &DesignCandidate) -> MillifluidicGenome {
    let n_topos = ALL_EVO_TOPOLOGIES.len();

    // Gene 0: topology index (find index in ALL_EVO_TOPOLOGIES; use midpoint of slot)
    let topo_idx = ALL_EVO_TOPOLOGIES
        .iter()
        .position(|t| topology_matches(t, &c.topology))
        .unwrap_or(0);
    let gene0 = (topo_idx as f64 + 0.5) / n_topos as f64;

    let is_leuka = matches!(
        c.topology,
        DesignTopology::ConstrictionExpansionArray { .. }
            | DesignTopology::SpiralSerpentine { .. }
            | DesignTopology::ParallelMicrochannelArray { .. }
    );

    // Gene 1: flow rate (log-linear inverse)
    let (q_min, q_max) = if is_leuka {
        (
            LEUKA_FLOW_RATES_M3_S[0],
            LEUKA_FLOW_RATES_M3_S[LEUKA_FLOW_RATES_M3_S.len() - 1],
        )
    } else {
        (
            FLOW_RATES_M3_S[0],
            FLOW_RATES_M3_S[FLOW_RATES_M3_S.len() - 1],
        )
    };
    let gene1 = if q_max > q_min && c.flow_rate_m3_s > 0.0 {
        ((c.flow_rate_m3_s / q_min).ln() / (q_max / q_min).ln()).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Gene 2: inlet gauge pressure (linear inverse)
    let p_min = INLET_GAUGES_PA[0];
    let p_max = INLET_GAUGES_PA[INLET_GAUGES_PA.len() - 1];
    let gene2 = ((c.inlet_gauge_pa - p_min) / (p_max - p_min)).clamp(0.0, 1.0);

    // Gene 3: throat diameter (linear inverse over [d_min, d_max])
    let d_min = THROAT_DIAMETERS_M[0];
    let d_max = THROAT_DIAMETERS_M[THROAT_DIAMETERS_M.len() - 1];
    let gene3 = if c.throat_diameter_m > 0.0 {
        ((c.throat_diameter_m - d_min) / (d_max - d_min)).clamp(0.0, 1.0)
    } else {
        0.5
    };

    // Gene 4: channel width (linear inverse for millifluidic [2mm, 6mm])
    let gene4 = if is_leuka {
        let w_min = LEUKA_CHANNEL_WIDTHS_M[0];
        let w_max = LEUKA_CHANNEL_WIDTHS_M[LEUKA_CHANNEL_WIDTHS_M.len() - 1];
        ((c.channel_width_m - w_min) / (w_max - w_min)).clamp(0.0, 1.0)
    } else {
        ((c.channel_width_m - 2.0e-3) / 4.0e-3).clamp(0.0, 1.0)
    };

    // Gene 5: serpentine segment count (linear inverse [2, 12])
    let gene5 = if c.serpentine_segments >= 2 {
        ((c.serpentine_segments - 2) as f64 / 10.0).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Gene 6: segment length fraction
    let gene6 = (c.segment_length_m / (TREATMENT_WIDTH_MM * 1e-3) - 0.5).clamp(0.0, 1.0);

    // Gene 7: bend radius fraction
    let gene7 = if c.segment_length_m > 0.0 {
        ((c.bend_radius_m / c.segment_length_m) - 0.05).clamp(0.0, 0.20) / 0.20
    } else {
        0.5
    };

    // Gene 8: discrete topology parameter (n_levels / n_pretri / n_cycles / etc.)
    let gene8 = match c.topology {
        DesignTopology::ConstrictionExpansionArray { n_cycles } => {
            ((n_cycles as f64 - 2.0) / 18.0).clamp(0.0, 1.0)
        }
        DesignTopology::SpiralSerpentine { n_turns } => {
            ((n_turns as f64 - 2.0) / 18.0).clamp(0.0, 1.0)
        }
        DesignTopology::ParallelMicrochannelArray { n_channels } => {
            ((n_channels as f64 - 10.0) / 490.0).clamp(0.0, 1.0)
        }
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
            ((n_levels as f64 - 1.0) / 2.0).clamp(0.0, 1.0)
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            ((n_pretri as f64 - 1.0) / 2.0).clamp(0.0, 1.0)
        }
        DesignTopology::DoubleTrifurcationCIFVenturi {
            center_throat_count,
        } => ((center_throat_count as f64 - 1.0) / 3.0).clamp(0.0, 1.0),
        _ => 0.5,
    };

    // Genes 9-12: AdaptiveTree per-level split types.
    // For non-Adaptive topologies, gene12 also stores centerline venturi count.
    let mut genes_9_12 = [0.0f64; 4];
    if let DesignTopology::AdaptiveTree {
        levels,
        split_types,
    } = c.topology
    {
        for i in 0..usize::from(levels.min(4)) {
            genes_9_12[i] = if ((split_types >> i) & 1) == 1 {
                1.0
            } else {
                0.0
            };
        }
    } else if matches!(
        c.topology,
        DesignTopology::CascadeCenterTrifurcationSeparator { .. }
            | DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
            | DesignTopology::AsymmetricTrifurcationVenturi
            | DesignTopology::TriBiTriSelectiveVenturi
            | DesignTopology::CellSeparationVenturi
            | DesignTopology::WbcCancerSeparationVenturi
            | DesignTopology::DoubleTrifurcationCIFVenturi { .. }
    ) {
        let stages = c.centerline_venturi_throat_count.clamp(1, 4);
        genes_9_12[3] = (f64::from(stages) - 1.0) / 3.0;
    }

    // Gene 13: trifurcation center frac [0.25, 0.65]
    let gene13 = ((c.trifurcation_center_frac - 0.25) / 0.40).clamp(0.0, 1.0);

    // Gene 14: CIF terminal tri center frac [0.25, 0.65]
    let gene14 = ((c.cif_terminal_tri_center_frac - 0.25) / 0.40).clamp(0.0, 1.0);

    // Gene 15: CIF terminal bi treat frac [0.50, 0.85]
    let gene15 = ((c.cif_terminal_bi_treat_frac - 0.50) / 0.35).clamp(0.0, 1.0);

    // Gene 16: asymmetric narrow frac [0.20, 0.70]
    let gene16 = ((c.asymmetric_narrow_frac - 0.20) / 0.50).clamp(0.0, 1.0);

    // Gene 17: throat length factor [1.5, 15.0]
    let tl_factor = if c.throat_diameter_m > 0.0 {
        c.throat_length_m / c.throat_diameter_m
    } else {
        2.0
    };
    let gene17 = ((tl_factor - 1.5) / 13.5).clamp(0.0, 1.0);

    // Gene 18: channel height log-linear [0.3mm, 3.0mm] — log10 normalised
    let gene18 = (c.channel_height_m / 0.3e-3).log10().clamp(0.0, 1.0);

    // Gene 19: trifurcation left frac [0.08, max_left]
    let gene19 = {
        let max_left = (0.85 - c.trifurcation_center_frac).max(0.08);
        ((c.trifurcation_left_frac - 0.08) / (max_left - 0.08)).clamp(0.0, 1.0)
    };

    let genes = [
        gene0,
        gene1,
        gene2,
        gene3,
        gene4,
        gene5,
        gene6,
        gene7,
        gene8,
        genes_9_12[0],
        genes_9_12[1],
        genes_9_12[2],
        genes_9_12[3],
        gene13,
        gene14,
        gene15,
        gene16,
        gene17,
        gene18,
        gene19,
    ];

    MillifluidicGenome { genes }
}

/// Check if a static topology slot matches a candidate's topology (ignoring variant params).
fn topology_matches(slot: &DesignTopology, candidate: &DesignTopology) -> bool {
    match (slot, candidate) {
        (
            DesignTopology::CascadeCenterTrifurcationSeparator { .. },
            DesignTopology::CascadeCenterTrifurcationSeparator { .. },
        ) => true,
        (
            DesignTopology::IncrementalFiltrationTriBiSeparator { .. },
            DesignTopology::IncrementalFiltrationTriBiSeparator { .. },
        ) => true,
        (
            DesignTopology::DoubleTrifurcationCIFVenturi { .. },
            DesignTopology::DoubleTrifurcationCIFVenturi { .. },
        ) => true,
        (
            DesignTopology::ConstrictionExpansionArray { .. },
            DesignTopology::ConstrictionExpansionArray { .. },
        ) => true,
        (DesignTopology::SpiralSerpentine { .. }, DesignTopology::SpiralSerpentine { .. }) => true,
        (
            DesignTopology::ParallelMicrochannelArray { .. },
            DesignTopology::ParallelMicrochannelArray { .. },
        ) => true,
        (DesignTopology::AdaptiveTree { .. }, DesignTopology::AdaptiveTree { .. }) => true,
        _ => slot == candidate,
    }
}
