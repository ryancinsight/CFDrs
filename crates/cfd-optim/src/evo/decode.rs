//! Genome → [`DesignCandidate`] decoding.

use crate::{
    constraints::{
        FLOW_RATES_M3_S, INLET_GAUGES_PA, LEUKA_CHANNEL_HEIGHT_M, LEUKA_CHANNEL_WIDTHS_M,
        LEUKA_FLOW_RATES_M3_S, PLATE_HEIGHT_MM, THROAT_DIAMETERS_M, TREATMENT_WIDTH_MM,
        VENTURI_INLET_DIAM_M,
    },
    design::{CrossSectionShape, DesignCandidate, DesignTopology},
};

use super::{MillifluidicGenome, ALL_EVO_TOPOLOGIES, MIN_ADAPTIVE_CH_M};

/// Decode a normalised genome into a [`DesignCandidate`].
///
/// Continuous gene 0 is rounded to select the topology.
/// Gene 5 is rounded to select the segment count.
/// Gene 8 encodes topology-specific discrete parameters for leukapheresis and
/// staged-separation variants (indices 15–17, 22–24).
#[must_use]
pub fn decode_genome(g: &MillifluidicGenome, id_prefix: &str) -> DesignCandidate {
    let genes = &g.genes;

    // Gene 0: topology index
    let topo_idx = (genes[0] * (ALL_EVO_TOPOLOGIES.len() as f64 - 1e-9))
        .floor()
        .clamp(0.0, (ALL_EVO_TOPOLOGIES.len() - 1) as f64) as usize;

    let n_topos = ALL_EVO_TOPOLOGIES.len();

    // Leukapheresis topologies (indices 15–17) require micro-scale channels.
    let is_leuka = topo_idx == 15 || topo_idx == 16 || topo_idx == 17;

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
    let topology = if topo_idx == n_topos - 1 {
        // ── AdaptiveTree: gene 8 → depth, genes 9–12 → split_types ───────
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
    } else if topo_idx == 15 {
        let n_cycles = (2.0 + genes[8] * 18.0).round() as usize;
        DesignTopology::ConstrictionExpansionArray { n_cycles }
    } else if topo_idx == 16 {
        let n_turns = (2.0 + genes[8] * 18.0).round() as usize;
        DesignTopology::SpiralSerpentine { n_turns }
    } else if topo_idx == 17 {
        let raw_n = (10.0 + genes[8] * 490.0).round() as usize;
        let pitch_m = w_ch * 2.0;
        let max_n = ((PLATE_HEIGHT_MM - 10.0) * 1e-3 / pitch_m).floor() as usize;
        let n_channels = raw_n.min(max_n).max(1);
        DesignTopology::ParallelMicrochannelArray { n_channels }
    } else if topo_idx == 22 {
        let n_levels = (1.0 + genes[8] * 2.0).round().clamp(1.0, 3.0) as u8;
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels }
    } else if topo_idx == 23 {
        let n_pretri = (1.0 + genes[8] * 2.0).round().clamp(1.0, 3.0) as u8;
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri }
    } else {
        ALL_EVO_TOPOLOGIES[topo_idx]
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
            format!("CCT-lv{n_levels}")
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            format!(
                "CIF-pt{}-pcf{}-tcf{}-btf{}",
                n_pretri,
                (cif_pretri_center_frac * 1000.0).round() as u32,
                (cif_terminal_tri_center_frac * 1000.0).round() as u32,
                (cif_terminal_bi_treat_frac * 1000.0).round() as u32,
            )
        }
        DesignTopology::AsymmetricTrifurcationVenturi => {
            format!(
                "ATV-cf{}-lf{}",
                (trifurcation_center_frac * 1000.0).round() as u32,
                (trifurcation_left_frac * 1000.0).round() as u32,
            )
        }
        DesignTopology::TriBiTriSelectiveVenturi => {
            format!(
                "TBT-cf{}-btf{}",
                (trifurcation_center_frac * 1000.0).round() as u32,
                (cif_terminal_bi_treat_frac * 1000.0).round() as u32,
            )
        }
        _ => format!("t{topo_idx}"),
    };

    let id = format!(
        "{}-EVO-{}-q{:.0}-g{:.0}-d{:.0}-w{:.0}-n{}",
        id_prefix,
        topo_tag,
        q * 6e7,
        gauge * 1e-3,
        d_throat * 1e6,
        w_ch * 1e6,
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
    }
}
