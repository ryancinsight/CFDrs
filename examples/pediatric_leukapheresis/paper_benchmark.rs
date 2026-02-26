//! Reference design candidates from leukapheresis literature.
//!
//! Two devices from the PMC10645346 review are replicated as [`DesignCandidate`]
//! instances so that Phase 1 can pass them through `compute_metrics()` and
//! compare predictions with the authors' reported values.
//!
//! ## Designs
//!
//! | Reference | Topology | D_h | Q | HCT |
//! |---|---|---|---|---|
//! | Nivedita 2017 | Inertial spiral | 133 µm | 1.8 mL/min | 2 % |
//! | Wu Z. 2019 | Constriction-expansion | 75 µm | 150 µL/min | 0.1 % |

use cfd_optim::{DesignCandidate, DesignTopology};

// ── Nivedita et al. (2017) — inertial spiral ─────────────────────────────────

/// Channel width of the Nivedita spiral. Units: m.
pub const NIVEDITA_WIDTH_M: f64 = 400e-6;
/// Channel height of the Nivedita spiral. Units: m.
pub const NIVEDITA_HEIGHT_M: f64 = 80e-6;
/// Number of turns in the Nivedita spiral.
pub const NIVEDITA_N_TURNS: usize = 14;
/// Mean bend radius of the Nivedita spiral. Units: m.
pub const NIVEDITA_BEND_M: f64 = 2.5e-3;
/// Volumetric flow rate per chip (1.8 mL/min). Units: m³/s.
pub const NIVEDITA_FLOW_M3S: f64 = 1.8e-6 / 60.0;
/// Feed HCT for the Nivedita experiment (~1:9 dilution of whole blood).
pub const NIVEDITA_HCT: f64 = 0.02;

/// Build a [`DesignCandidate`] that replicates the Nivedita et al. (2017) spiral chip.
///
/// Hydraulic diameter D_h = 2 × 400 × 80 / (400 + 80) ≈ 133 µm.
/// Confinement ratio κ_WBC ≈ 12 µm / 133 µm ≈ 0.090 (marginal inertial focusing).
pub fn nivedita_spiral() -> DesignCandidate {
    let w = NIVEDITA_WIDTH_M;
    let h = NIVEDITA_HEIGHT_M;
    DesignCandidate {
        id: "nivedita_spiral_2017".to_owned(),
        topology: DesignTopology::SpiralSerpentine { n_turns: NIVEDITA_N_TURNS },
        inlet_gauge_pa:         50_000.0,   // ~0.5 bar, typical syringe-pump outlet
        flow_rate_m3_s:         NIVEDITA_FLOW_M3S,
        throat_diameter_m:      w,          // not a venturi; set to channel width
        inlet_diameter_m:       w,
        throat_length_m:        1e-6,       // unused; set to negligible positive value
        channel_width_m:        w,
        channel_height_m:       h,
        serpentine_segments:    NIVEDITA_N_TURNS,
        // Arc length of one 360° turn at the mean bend radius.
        segment_length_m:       2.0 * std::f64::consts::PI * NIVEDITA_BEND_M,
        bend_radius_m:          NIVEDITA_BEND_M,
        feed_hematocrit:        NIVEDITA_HCT,
    }
}

// ── Wu Z. et al. (2019) — constriction-expansion ─────────────────────────────

/// Wide-section channel width for the Wu design. Units: m.
pub const WU_WIDE_WIDTH_M: f64 = 300e-6;
/// Channel height for the Wu design. Units: m.
pub const WU_HEIGHT_M: f64 = 60e-6;
/// Number of constriction-expansion cycles in the Wu design.
pub const WU_N_CYCLES: usize = 20;
/// Volumetric flow rate per chip (150 µL/min). Units: m³/s.
pub const WU_FLOW_M3S: f64 = 150e-9 / 60.0;
/// Feed HCT for the Wu experiment (heavily diluted, ~0.1 %). Dimensionless.
pub const WU_HCT: f64 = 0.001;

/// Build a [`DesignCandidate`] that replicates the Wu Z. et al. (2019) chip.
///
/// Narrow constriction width = 150 µm (50 % of 300 µm wide section).
/// D_h_narrow = 2 × 150 × 60 / (150 + 60) ≈ 86 µm.
pub fn wu_constriction() -> DesignCandidate {
    let w = WU_WIDE_WIDTH_M;
    let h = WU_HEIGHT_M;
    DesignCandidate {
        id: "wu_constriction_2019".to_owned(),
        topology: DesignTopology::ConstrictionExpansionArray { n_cycles: WU_N_CYCLES },
        inlet_gauge_pa:         30_000.0,   // ~0.3 bar, syringe pump
        flow_rate_m3_s:         WU_FLOW_M3S,
        throat_diameter_m:      w * 0.50,   // narrow constriction width
        inlet_diameter_m:       w,
        throat_length_m:        250e-6,     // approximate narrow-section length
        channel_width_m:        w,
        channel_height_m:       h,
        serpentine_segments:    WU_N_CYCLES,
        // Approximate straight-segment length: total device length / (2 × n_cycles).
        segment_length_m:       5e-3,       // ~5 mm per half-cycle
        bend_radius_m:          0.0,        // straight channel; no bends
        feed_hematocrit:        WU_HCT,
    }
}
