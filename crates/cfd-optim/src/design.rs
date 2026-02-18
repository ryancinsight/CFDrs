//! Design topology families and candidate parameter space for SDT optimization.
//!
//! Five topology families are evaluated:
//!
//! | ID | Topology | Cavitation sites | 6×6 coverage |
//! |----|----------|-----------------|--------------|
//! | 1  | `SingleVenturi`          | 1  | point (~1 well) |
//! | 2  | `BifurcationVenturi`     | 4  | 4-quadrant grid |
//! | 3  | `TrifurcationVenturi`    | 3  | 3-column zones  |
//! | 4  | `VenturiSerpentine`      | 1  | full 6×6 grid   |
//! | 5  | `SerpentineGrid`         | 0  | full 6×6 grid   |

use crate::constraints::*;
use serde::{Deserialize, Serialize};
use std::fmt;

// ── Topology families ────────────────────────────────────────────────────────

/// The five millifluidic topology families for SDT / exposure optimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DesignTopology {
    /// Inlet → single venturi throat → outlet.
    ///
    /// Maximum cavitation intensity at one point; minimal spatial coverage.
    /// Best when a single treatment well is targeted.
    SingleVenturi,

    /// Inlet → 2-stage symmetric binary bifurcation (1 → 4 branches) → venturi
    /// throat in each branch → outlets.
    ///
    /// Four cavitation sites arranged over four quadrants of the 6 × 6 grid.
    /// Balanced cavitation + spatial coverage.
    BifurcationVenturi,

    /// Inlet → symmetric trifurcation (1 → 3 branches) → venturi throat in
    /// each branch → outlets.
    ///
    /// Three cavitation sites spanning the treatment zone columns.
    TrifurcationVenturi,

    /// Inlet → venturi throat (cavitation generation) → serpentine channels
    /// sweeping all 36 wells of the 6 × 6 treatment zone.
    ///
    /// Combines a single cavitation stage with full-grid exposure and
    /// extended residence time.
    VenturiSerpentine,

    /// Inlet → pure serpentine grid covering all 36 wells of the 6 × 6 zone.
    ///
    /// No dedicated cavitation stage; maximises uniform exposure, residence
    /// time, and light / ultrasound path length.  FDA compliance is easier to
    /// achieve because the flow is distributed over a large channel area.
    SerpentineGrid,
}

impl DesignTopology {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            Self::SingleVenturi => "Single Venturi",
            Self::BifurcationVenturi => "Bifurcation + Venturi ×4",
            Self::TrifurcationVenturi => "Trifurcation + Venturi ×3",
            Self::VenturiSerpentine => "Venturi → Serpentine Grid",
            Self::SerpentineGrid => "Serpentine Grid (full 6×6)",
        }
    }

    /// Short 3-character code used in candidate IDs.
    pub fn short(self) -> &'static str {
        match self {
            Self::SingleVenturi => "SV",
            Self::BifurcationVenturi => "BV",
            Self::TrifurcationVenturi => "TV",
            Self::VenturiSerpentine => "VS",
            Self::SerpentineGrid => "SG",
        }
    }

    /// Returns `true` if this topology includes at least one venturi throat.
    pub fn has_venturi(self) -> bool {
        matches!(
            self,
            Self::SingleVenturi
                | Self::BifurcationVenturi
                | Self::TrifurcationVenturi
                | Self::VenturiSerpentine
        )
    }

    /// Returns `true` if this topology includes a downstream serpentine or
    /// bifurcation distribution network covering multiple wells.
    pub fn has_distribution(self) -> bool {
        matches!(
            self,
            Self::BifurcationVenturi
                | Self::TrifurcationVenturi
                | Self::VenturiSerpentine
                | Self::SerpentineGrid
        )
    }

    /// Returns `true` if the distribution stage is serpentine-based.
    pub fn has_serpentine(self) -> bool {
        matches!(self, Self::VenturiSerpentine | Self::SerpentineGrid)
    }

    /// Number of venturi stages operating in parallel.
    pub fn venturi_count(self) -> usize {
        match self {
            Self::SingleVenturi => 1,
            Self::BifurcationVenturi => 4,
            Self::TrifurcationVenturi => 3,
            Self::VenturiSerpentine => 1,
            Self::SerpentineGrid => 0,
        }
    }

    /// Number of parallel outlet branches.
    pub fn outlet_count(self) -> usize {
        match self {
            Self::SingleVenturi => 1,
            Self::BifurcationVenturi => 4,
            Self::TrifurcationVenturi => 3,
            Self::VenturiSerpentine => 1,
            Self::SerpentineGrid => 1,
        }
    }

    /// Fraction of the 36 treatment wells covered by this topology's channels.
    ///
    /// Serpentine designs pass over every well row; venturi-branch designs
    /// cover quadrants / columns; single venturi covers roughly one well.
    pub fn nominal_well_coverage(self) -> f64 {
        match self {
            // Single venturi outlet covers roughly the volume around 1 well
            Self::SingleVenturi => 1.0 / TREATMENT_WELL_COUNT as f64,
            // 4 venturis arranged in 2×2 pattern → each quadrant (9 wells)
            Self::BifurcationVenturi => 1.0,
            // 3 venturis each covering ≈12 wells across 3 column bands
            Self::TrifurcationVenturi => 1.0,
            // Serpentine sweeps all 6 rows → all 36 wells
            Self::VenturiSerpentine | Self::SerpentineGrid => 1.0,
        }
    }
}

impl fmt::Display for DesignTopology {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}

// ── Design candidate ────────────────────────────────────────────────────────

/// A single fully-specified design candidate in the parameter sweep.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesignCandidate {
    /// Unique sequential identifier (e.g. `"0001-SV"`).
    pub id: String,
    /// Topology family.
    pub topology: DesignTopology,
    /// Total volumetric inlet flow rate [m³/s].
    pub flow_rate_m3_s: f64,
    /// Gauge inlet pressure above atmospheric [Pa].
    /// Absolute inlet pressure = `P_ATM_PA + inlet_gauge_pa`.
    pub inlet_gauge_pa: f64,
    /// Venturi throat diameter [m].  Zero / ignored if `!topology.has_venturi()`.
    pub throat_diameter_m: f64,
    /// Venturi inlet (upstream) channel diameter [m].
    pub inlet_diameter_m: f64,
    /// Venturi throat length [m]  (≈ 2× throat diameter).
    pub throat_length_m: f64,
    /// Width of rectangular main / serpentine channels [m].
    pub channel_width_m: f64,
    /// Height of rectangular main / serpentine channels [m].
    pub channel_height_m: f64,
    /// Number of straight segments in the serpentine section.
    pub serpentine_segments: usize,
    /// Length of each serpentine straight segment [m].
    pub segment_length_m: f64,
    /// Radius of serpentine 180° bends [m].
    pub bend_radius_m: f64,
}

impl DesignCandidate {
    /// Absolute inlet pressure [Pa].
    #[inline]
    pub fn inlet_pressure_pa(&self) -> f64 {
        P_ATM_PA + self.inlet_gauge_pa
    }

    /// Cross-sectional area of the circular venturi inlet channel [m²].
    #[inline]
    pub fn inlet_area_m2(&self) -> f64 {
        std::f64::consts::FRAC_PI_4 * self.inlet_diameter_m * self.inlet_diameter_m
    }

    /// Mean velocity at the venturi inlet [m/s] for the *per-venturi* share of flow.
    #[inline]
    pub fn inlet_velocity_m_s(&self) -> f64 {
        let q_per = self.per_venturi_flow();
        q_per / self.inlet_area_m2()
    }

    /// Flow rate allocated to **one** venturi stage [m³/s].
    #[inline]
    pub fn per_venturi_flow(&self) -> f64 {
        let n = self.topology.venturi_count().max(1);
        self.flow_rate_m3_s / n as f64
    }

    /// Mean velocity in the rectangular main channel [m/s].
    #[inline]
    pub fn channel_velocity_m_s(&self) -> f64 {
        let q_per = self.flow_rate_m3_s / self.topology.outlet_count().max(1) as f64;
        q_per / (self.channel_width_m * self.channel_height_m)
    }

    /// Total approximate channel path length [mm].
    pub fn total_path_length_mm(&self) -> f64 {
        let mut len = 0.0_f64;
        // Venturi section (inlet cone + throat + diffuser, approximated as 10× d_inlet)
        if self.topology.has_venturi() {
            len += self.inlet_diameter_m * 10.0 * self.topology.venturi_count() as f64 * 1000.0;
        }
        // Serpentine section
        if self.topology.has_serpentine() {
            let n = self.serpentine_segments as f64;
            let bends = (self.serpentine_segments.saturating_sub(1)) as f64;
            len += (n * self.segment_length_m + bends * std::f64::consts::PI * self.bend_radius_m)
                * 1000.0;
        }
        len
    }
}

// ── Parameter-sweep candidate space ─────────────────────────────────────────

/// Build the complete parametric sweep over all five topology families.
///
/// Approximately **5 × 3 × 3 × 3 × 3 × 2 ≈ 810** candidates are generated
/// before any physics-based filtering.  Topologies without a venturi stage
/// skip the throat-diameter dimension; topologies without a distribution
/// stage skip the segment-count dimension.
pub fn build_candidate_space() -> Vec<DesignCandidate> {
    let topologies = [
        DesignTopology::SingleVenturi,
        DesignTopology::BifurcationVenturi,
        DesignTopology::TrifurcationVenturi,
        DesignTopology::VenturiSerpentine,
        DesignTopology::SerpentineGrid,
    ];

    // throat diameter = 0 placeholder for topologies without venturi
    let throat_none: [f64; 1] = [0.0];

    let mut candidates = Vec::with_capacity(512);
    let mut idx: u32 = 0;

    for &topology in &topologies {
        let throat_iter: &[f64] = if topology.has_venturi() {
            &THROAT_DIAMETERS_M
        } else {
            &throat_none
        };

        let seg_counts: &[usize] = if topology.has_serpentine() {
            &SERPENTINE_SEGMENT_COUNTS
        } else {
            &[1_usize] // unused but keeps the loop structure uniform
        };

        for &q in &FLOW_RATES_M3_S {
            for &gauge in &INLET_GAUGES_PA {
                for &d_throat in throat_iter {
                    for &w_ch in &CHANNEL_WIDTHS_M {
                        for &n_segs in seg_counts {
                            idx += 1;
                            let throat_len = if d_throat > 0.0 {
                                d_throat * 2.0 // 2 diameters long
                            } else {
                                0.0
                            };
                            // seg_length always spans the full treatment width
                            let seg_len = TREATMENT_WIDTH_MM * 1e-3;

                            let id = format!(
                                "{:04}-{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-w{:.0}um-n{}",
                                idx,
                                topology.short(),
                                q * 6e7,        // → mL/min × 10 (avoids decimals)
                                gauge * 1e-3,   // → kPa
                                d_throat * 1e6, // → μm
                                w_ch * 1e6,     // → μm
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
                                channel_height_m: CHANNEL_HEIGHT_M,
                                serpentine_segments: n_segs,
                                segment_length_m: seg_len,
                                bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                            });
                        }
                    }
                }
            }
        }
    }

    candidates
}
