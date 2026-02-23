//! Design topology families and candidate parameter space for SDT optimization.
//!
//! Fifteen topology families are evaluated (14 fixed + 1 GA-only adaptive):
//!
//! | ID | Topology | Cavitation sites | 6×6 coverage |
//! |----|----------|-----------------|--------------|
//! | 1  | `SingleVenturi`                    |  1 | point (~1 well)   |
//! | 2  | `BifurcationVenturi`               |  2 | 2-zone split      |
//! | 3  | `TrifurcationVenturi`              |  3 | 3-column zones    |
//! | 4  | `VenturiSerpentine`                |  1 | full 6×6 grid     |
//! | 5  | `SerpentineGrid`                   |  0 | full 6×6 grid     |
//! | 6  | `CellSeparationVenturi`            |  1 | 50%               |
//! | 7  | `WbcCancerSeparationVenturi`       |  1 | 50%               |
//! | 8  | `DoubleBifurcationVenturi`         |  4 | 2×2 quad grid     |
//! | 9  | `TripleBifurcationVenturi`         |  8 | 2×2×2 octant grid |
//! | 10 | `DoubleTrifurcationVenturi`        |  9 | 3×3 grid          |
//! | 11 | `BifurcationTrifurcationVenturi`   |  6 | 2×3 grid          |
//! | 12 | `SerialDoubleVenturi`              |  2 | point (serial)    |
//! | 13 | `BifurcationSerpentine`            |  0 | full 6×6 (2 arms) |
//! | 14 | `TrifurcationSerpentine`           |  0 | full 6×6 (3 arms) |
//! | GA | `AdaptiveTree`                     | N | full 6×6 grid     |
//!
//! `AdaptiveTree` is GA-only: depth (0–4) and per-level split type (Bi / Tri)
//! are encoded in genome genes 8–12 and constrained so the leaf channel width
//! stays ≥ 150 µm.

use cfd_schematics::{
    interface::presets::{bifurcation_rect, serpentine_rect, trifurcation_rect, venturi_rect},
    NetworkBlueprint,
};

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

    /// Inlet → curved section (Dean flow margination) → venturi throat
    /// (inertial focusing of cancer cells to center) → split outlet:
    /// center channel (cancer cells + cavitation) + peripheral channels
    /// (healthy cells / RBCs pushed to walls by Dean flow).
    ///
    /// Designed for **cancer cell isolation + sonodynamic therapy**:
    /// - Large, stiff cancer cells (MCF-7, κ ≈ 0.35) focus to channel center.
    /// - Small, deformable RBCs (κ ≈ 0.14) focus to channel walls.
    /// - Venturi throat in the center channel generates cavitation for SDT.
    /// - Peripheral channels return healthy cells without cavitation exposure.
    ///
    /// Physics: inertial lift (Di Carlo 2009) + Dean drag (Gossett 2009).
    CellSeparationVenturi,

    /// Inlet → curved serpentine section (Dean flow margination) → venturi
    /// throat → three-way split outlet.
    ///
    /// Designed for **three-population cell separation + SDT**:
    /// - **Cancer cells (MCF-7, 17.5 µm, DI=0.15)** → channel center (x̃ < 0.3)
    /// - **White blood cells (WBC, 10 µm, DI=0.45)** → channel center (x̃ < 0.3)
    /// - **Red blood cells (RBC, 7 µm, DI=0.85)** → channel periphery (x̃ > 0.5)
    ///
    /// The goal is to co-focus WBCs and cancer cells at the channel center for
    /// joint sonodynamic cavitation therapy while routing RBCs to peripheral
    /// channels to minimise haemolysis exposure.
    ///
    /// Physics: inertial lift (Di Carlo 2009) + Dean drag (Gossett 2009).
    /// Smaller channel dimensions than `CellSeparationVenturi` ensure
    /// κ_WBC > 0.07 (WBC inertial focusing threshold).
    WbcCancerSeparationVenturi,

    // ── Multi-level bifurcation / trifurcation trees ─────────────────────────

    /// 2-level symmetric bifurcation → 4 parallel venturi throats (2×2 grid).
    ///
    /// Provides 4 cavitation sites distributed over the four quadrants of the
    /// 6 × 6 treatment zone.  Balanced pressure and flow at all four outlets
    /// by construction.
    DoubleBifurcationVenturi,

    /// 3-level symmetric bifurcation → 8 parallel venturi throats (2×2×2 grid).
    ///
    /// Dense cavitation grid covering all 8 octant zones.  High cavitation
    /// site density; moderate per-site intensity.
    TripleBifurcationVenturi,

    /// 2-level symmetric trifurcation → 9 parallel venturi throats (3×3 grid).
    ///
    /// Provides the most uniform coverage for the full 6 × 6 treatment zone.
    /// Nine cavitation sites arranged in a 3×3 matrix.
    DoubleTrifurcationVenturi,

    /// 1 bifurcation + 1 trifurcation → 6 parallel venturi throats (2×3 grid).
    ///
    /// Combines the flexibility of an asymmetric tree with 6-column-zone
    /// coverage.  Useful when plate geometry aligns with 6 column groups.
    BifurcationTrifurcationVenturi,

    // ── Serial venturis ──────────────────────────────────────────────────────

    /// 2 venturi throats in series on the same flow path.
    ///
    /// First throat nucleates cavitation; second collapses and re-nucleates.
    /// Lower cavitation onset threshold than a single venturi at the same
    /// inlet pressure, because the first venturi pre-conditions the fluid.
    /// Single outlet; maximised cavitation intensity at the cost of coverage.
    SerialDoubleVenturi,

    // ── Bifurcation / trifurcation → serpentine arms ─────────────────────────

    /// Symmetric bifurcation → 2 parallel serpentine arms.
    ///
    /// Each arm covers one half of the 6 × 6 treatment zone.
    /// Maximises residence time while halving the shear rate per arm
    /// (Q/2 per branch).  Ideal for photopheresis-like uniform exposure
    /// with conservative FDA margins.
    BifurcationSerpentine,

    /// Symmetric trifurcation → 3 parallel serpentine arms.
    ///
    /// Each arm covers one third (2-column band) of the treatment zone.
    /// Q/3 per branch gives the lowest shear of any topology; best for
    /// gentle exposure with minimal haemolysis.
    TrifurcationSerpentine,

    // ── GA-only: variable-depth adaptive split tree ───────────────────────────

    /// Variable-depth split-merge tree with venturi throats at every leaf node.
    ///
    /// **Used by the GA only** — not included in the parametric sweep.
    ///
    /// Each level can independently be a binary (Bi) or ternary (Tri) split,
    /// encoded compactly in `split_types`:
    ///
    /// - Bit `i` of `split_types` = `0` → Bifurcation at level `i`
    /// - Bit `i` of `split_types` = `1` → Trifurcation at level `i`
    ///
    /// The maximum depth is constrained at decode time so that the leaf
    /// channel width `w_ch / ∏ fan_i ≥ 150 µm` (minimum fabricatable width).
    ///
    /// # Examples
    ///
    /// | `levels` | `split_types` | Leaf count | Structure          |
    /// |----------|--------------|------------|--------------------|
    /// | 0        | 0            | 1          | Single venturi     |
    /// | 1        | 0b_0         | 2          | [Bi] → 2 venturis  |
    /// | 1        | 0b_1         | 3          | [Tri] → 3 venturis |
    /// | 2        | 0b_00        | 4          | [Bi,Bi] → 4        |
    /// | 2        | 0b_01        | 6          | [Tri,Bi] → 6       |
    /// | 2        | 0b_11        | 9          | [Tri,Tri] → 9      |
    AdaptiveTree {
        /// Number of split levels (0 = single channel, max 4).
        levels: u8,
        /// Packed split-type bits: bit `i` = 0 → Bi, 1 → Tri for level `i`.
        split_types: u8,
    },
}

impl DesignTopology {
    /// Human-readable name.
    pub fn name(self) -> &'static str {
        match self {
            Self::SingleVenturi                  => "Single Venturi",
            Self::BifurcationVenturi             => "Bifurcation + Venturi ×2",
            Self::TrifurcationVenturi            => "Trifurcation + Venturi ×3",
            Self::VenturiSerpentine              => "Venturi → Serpentine Grid",
            Self::SerpentineGrid                 => "Serpentine Grid (full 6×6)",
            Self::CellSeparationVenturi          => "Cell Separation + Venturi (cancer focus)",
            Self::WbcCancerSeparationVenturi     => "WBC+Cancer→Center, RBC→Periphery + Venturi",
            Self::DoubleBifurcationVenturi       => "[Bi,Bi] → 4 Parallel Venturis (2×2 grid)",
            Self::TripleBifurcationVenturi       => "[Bi,Bi,Bi] → 8 Parallel Venturis (octants)",
            Self::DoubleTrifurcationVenturi      => "[Tri,Tri] → 9 Parallel Venturis (3×3 grid)",
            Self::BifurcationTrifurcationVenturi => "[Bi,Tri] → 6 Parallel Venturis (2×3 grid)",
            Self::SerialDoubleVenturi            => "2× Venturi in Series (double cavitation)",
            Self::BifurcationSerpentine          => "Bifurcation → 2 Serpentine Arms",
            Self::TrifurcationSerpentine         => "Trifurcation → 3 Serpentine Arms",
            Self::AdaptiveTree { .. }            => "Adaptive Split Tree",
        }
    }

    /// Short code used in candidate IDs and schematic file names.
    pub fn short(self) -> &'static str {
        match self {
            Self::SingleVenturi                  => "SV",
            Self::BifurcationVenturi             => "BV",
            Self::TrifurcationVenturi            => "TV",
            Self::VenturiSerpentine              => "VS",
            Self::SerpentineGrid                 => "SG",
            Self::CellSeparationVenturi          => "CS",
            Self::WbcCancerSeparationVenturi     => "WC",
            Self::DoubleBifurcationVenturi       => "D2",
            Self::TripleBifurcationVenturi       => "D3",
            Self::DoubleTrifurcationVenturi      => "T9",
            Self::BifurcationTrifurcationVenturi => "B6",
            Self::SerialDoubleVenturi            => "S2",
            Self::BifurcationSerpentine          => "BS",
            Self::TrifurcationSerpentine         => "TS",
            Self::AdaptiveTree { .. }            => "AT",
        }
    }

    /// Returns `true` if this topology includes at least one venturi throat.
    ///
    /// `AdaptiveTree` always includes venturi throats at the leaf nodes.
    pub fn has_venturi(self) -> bool {
        matches!(
            self,
            Self::SingleVenturi
                | Self::BifurcationVenturi
                | Self::TrifurcationVenturi
                | Self::VenturiSerpentine
                | Self::CellSeparationVenturi
                | Self::WbcCancerSeparationVenturi
                | Self::DoubleBifurcationVenturi
                | Self::TripleBifurcationVenturi
                | Self::DoubleTrifurcationVenturi
                | Self::BifurcationTrifurcationVenturi
                | Self::SerialDoubleVenturi
                | Self::AdaptiveTree { .. }
        )
    }

    /// Returns `true` if this topology includes a distribution network
    /// (bifurcation tree, trifurcation tree, or serpentine).
    ///
    /// For `AdaptiveTree`, returns `true` only when `levels > 0` (depth-0 is
    /// equivalent to a single straight channel with a venturi).
    pub fn has_distribution(self) -> bool {
        match self {
            Self::BifurcationVenturi
            | Self::TrifurcationVenturi
            | Self::VenturiSerpentine
            | Self::SerpentineGrid
            | Self::CellSeparationVenturi
            | Self::WbcCancerSeparationVenturi
            | Self::DoubleBifurcationVenturi
            | Self::TripleBifurcationVenturi
            | Self::DoubleTrifurcationVenturi
            | Self::BifurcationTrifurcationVenturi
            | Self::BifurcationSerpentine
            | Self::TrifurcationSerpentine => true,
            Self::AdaptiveTree { levels, .. } => levels > 0,
            _ => false,
        }
    }

    /// Returns `true` if the distribution stage is serpentine-based.
    pub fn has_serpentine(self) -> bool {
        matches!(
            self,
            Self::VenturiSerpentine
                | Self::SerpentineGrid
                | Self::BifurcationSerpentine
                | Self::TrifurcationSerpentine
        )
    }

    /// Total number of venturi stages (parallel + serial).
    pub fn venturi_count(self) -> usize {
        match self {
            Self::SingleVenturi                  => 1,
            Self::BifurcationVenturi             => 2,
            Self::TrifurcationVenturi            => 3,
            Self::VenturiSerpentine              => 1,
            Self::SerpentineGrid                 => 0,
            Self::CellSeparationVenturi          => 1,
            Self::WbcCancerSeparationVenturi     => 1,
            Self::DoubleBifurcationVenturi       => 4,
            Self::TripleBifurcationVenturi       => 8,
            Self::DoubleTrifurcationVenturi      => 9,
            Self::BifurcationTrifurcationVenturi => 6,
            // Serial: 2 in series on one path
            Self::SerialDoubleVenturi            => 2,
            Self::BifurcationSerpentine          => 0,
            Self::TrifurcationSerpentine         => 0,
            // Adaptive: multiply fan factors across all levels
            Self::AdaptiveTree { levels, split_types } => {
                let mut count = 1usize;
                for i in 0..levels as usize {
                    count *= if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                }
                count
            }
        }
    }

    /// Number of venturi stages running in PARALLEL (for flow-per-venturi calc).
    ///
    /// For series topologies (`SerialDoubleVenturi`) both venturis see the full
    /// flow Q, so this returns 1.  For tree topologies the flow splits, so this
    /// equals [`venturi_count`].
    pub fn parallel_venturi_count(self) -> usize {
        match self {
            Self::SerialDoubleVenturi => 1, // full Q through each stage in series
            _ => self.venturi_count(),
        }
    }

    /// Number of serial venturi stages on one flow path.
    ///
    /// `1` for all topologies except `SerialDoubleVenturi` which has `2`.
    pub fn serial_venturi_stages(self) -> usize {
        match self {
            Self::SerialDoubleVenturi => 2,
            _ if self.has_venturi() => 1,
            _ => 0,
        }
    }

    /// Number of parallel outlet branches.
    pub fn outlet_count(self) -> usize {
        match self {
            Self::SingleVenturi                  => 1,
            Self::BifurcationVenturi             => 2,   // 1-level [Bi] → 2 outlets
            Self::TrifurcationVenturi            => 3,
            Self::VenturiSerpentine              => 1,
            Self::SerpentineGrid                 => 1,
            // 2 outlets: center (cancer cells) + peripheral (healthy cells)
            Self::CellSeparationVenturi          => 2,
            // 2 outlets: center (WBC+cancer) + peripheral (RBC)
            Self::WbcCancerSeparationVenturi     => 2,
            Self::DoubleBifurcationVenturi       => 4,   // [Bi,Bi]     → 4 outlets
            Self::TripleBifurcationVenturi       => 8,   // [Bi,Bi,Bi]  → 8 outlets
            Self::DoubleTrifurcationVenturi      => 9,   // [Tri,Tri]   → 9 outlets
            Self::BifurcationTrifurcationVenturi => 6,   // [Bi,Tri]    → 6 outlets
            Self::SerialDoubleVenturi            => 1,   // serial path → single outlet
            Self::BifurcationSerpentine          => 2,   // [Bi] → 2 serpentine arms
            Self::TrifurcationSerpentine         => 3,   // [Tri] → 3 serpentine arms
            // Adaptive tree: same as venturi_count (one outlet per leaf)
            Self::AdaptiveTree { levels, split_types } => {
                let mut count = 1usize;
                for i in 0..levels as usize {
                    count *= if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                }
                count
            }
        }
    }

    /// Fraction of the 36 treatment wells covered by this topology's channels.
    ///
    /// Serpentine designs pass over every well row; venturi-branch designs
    /// cover quadrants / columns; single venturi covers roughly one well.
    pub fn nominal_well_coverage(self) -> f64 {
        match self {
            // Single-point cavitation: covers ~1 well
            Self::SingleVenturi | Self::SerialDoubleVenturi => {
                1.0 / TREATMENT_WELL_COUNT as f64
            }
            // 2 venturis spanning 2 half-plate zones → full 6×6
            Self::BifurcationVenturi => 1.0,
            // 3 venturis each covering ~12 wells → full 6×6
            Self::TrifurcationVenturi => 1.0,
            // Serpentine sweeps all 6 rows → all 36 wells
            Self::VenturiSerpentine
            | Self::SerpentineGrid
            | Self::BifurcationSerpentine
            | Self::TrifurcationSerpentine => 1.0,
            // Center channel covers cancer cell treatment wells (≈ half the plate)
            Self::CellSeparationVenturi => 0.5,
            // Center channel covers WBC+cancer treatment wells (≈ half the plate)
            Self::WbcCancerSeparationVenturi => 0.5,
            // Multi-level tree topologies distribute over the full 6×6 zone
            Self::DoubleBifurcationVenturi
            | Self::TripleBifurcationVenturi
            | Self::DoubleTrifurcationVenturi
            | Self::BifurcationTrifurcationVenturi => 1.0,
            // Adaptive tree: distributes venturis across the full zone
            Self::AdaptiveTree { .. } => 1.0,
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
    ///
    /// For series topologies the full inlet flow passes through each stage;
    /// for tree (parallel) topologies the flow splits across branches.
    #[inline]
    pub fn per_venturi_flow(&self) -> f64 {
        let n = self.topology.parallel_venturi_count().max(1);
        self.flow_rate_m3_s / n as f64
    }

    /// Mean velocity in the rectangular main channel [m/s].
    #[inline]
    pub fn channel_velocity_m_s(&self) -> f64 {
        let q_per = self.flow_rate_m3_s / self.topology.outlet_count().max(1) as f64;
        q_per / (self.channel_width_m * self.channel_height_m)
    }

    /// Build a [`NetworkBlueprint`] for this candidate using `cfd-schematics` presets.
    ///
    /// The blueprint encodes the channel cross-sections in `ChannelSpec.cross_section`
    /// so that downstream consumers (e.g. `compute_metrics`) can read geometry via
    /// `CrossSectionSpec::area()` / `hydraulic_diameter()` rather than recomputing it
    /// from raw parameter fields.
    ///
    /// **Naming convention** (relied on by `compute_metrics`):
    /// - `"throat_section"` — venturi throat (`throat_diameter_m × channel_height_m`)
    /// - `"inlet_section"` — venturi inlet / main channel approach
    /// - `"segment_1"`, `"segment_2"`, … — serpentine straight segments
    /// - `"parent"` — bifurcation / trifurcation trunk
    #[must_use]
    pub fn to_blueprint(&self) -> NetworkBlueprint {
        let w  = self.channel_width_m;
        let h  = self.channel_height_m;
        let dt = self.throat_diameter_m;
        let tl = self.throat_length_m.max(dt * 2.0);
        let n  = self.serpentine_segments;
        let sl = self.segment_length_m;

        match self.topology {
            // ── Single venturi (and cell-separation topologies using a venturi) ──
            DesignTopology::SingleVenturi
            | DesignTopology::CellSeparationVenturi
            | DesignTopology::WbcCancerSeparationVenturi => {
                venturi_rect(&self.id, w, dt, h, tl)
            }

            // ── Pure serpentine grid ──
            DesignTopology::SerpentineGrid => serpentine_rect(&self.id, n, sl, w, h),

            // ── Venturi + downstream serpentine ──
            DesignTopology::VenturiSerpentine => {
                // Start with the venturi (gives inlet_section, throat_section, diffuser_section)
                let mut bp = venturi_rect(&self.id, w, dt, h, tl);
                // Append serpentine channels (segment_1, segment_2, …)
                let ser = serpentine_rect(&self.id, n, sl, w, h);
                for ch in ser.channels {
                    bp.channels.push(ch);
                }
                bp
            }

            // ── Bifurcation tree → venturi in each leaf ──
            DesignTopology::BifurcationVenturi => {
                // Distribution tree: trunk + two first-level branches
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                let mut bp = bifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                // Append representative venturi channels so the throat section is discoverable
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels {
                    bp.channels.push(ch);
                }
                bp
            }

            // ── Trifurcation tree → venturi in each leaf ──
            DesignTopology::TrifurcationVenturi => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.50e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                let mut bp = trifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels {
                    bp.channels.push(ch);
                }
                bp
            }

            // ── DoubleBifurcation [Bi,Bi] → 4 parallel venturis ──
            DesignTopology::DoubleBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let mut bp = bifurcation_rect(&self.id, trunk_len, branch1_len, w, w, h);
                let bp2 = bifurcation_rect(
                    "_l2",
                    TREATMENT_HEIGHT_MM * 0.12e-3,
                    TREATMENT_HEIGHT_MM * 0.10e-3,
                    w, w, h,
                );
                for ch in bp2.channels { bp.channels.push(ch); }
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── TripleBifurcation [Bi,Bi,Bi] → 8 parallel venturis ──
            DesignTopology::TripleBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch_len  = TREATMENT_HEIGHT_MM * 0.12e-3;
                let mut bp = bifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── DoubleTrifurcation [Tri,Tri] → 9 parallel venturis ──
            DesignTopology::DoubleTrifurcationVenturi => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let mut bp = trifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── BifurcationTrifurcation [Bi,Tri] → 6 parallel venturis ──
            DesignTopology::BifurcationTrifurcationVenturi => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let mut bp = bifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let vbp = venturi_rect("_v", w, dt, h, tl);
                for ch in vbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── SerialDoubleVenturi (2 venturis in series on one path) ──
            DesignTopology::SerialDoubleVenturi => {
                // Blueprint shows representative single venturi; the serial
                // second stage is reflected in metrics but not in this schematic.
                venturi_rect(&self.id, w, dt, h, tl)
            }

            // ── BifurcationSerpentine [Bi] → 2 serpentine arms ──
            DesignTopology::BifurcationSerpentine => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                let mut bp = bifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let sbp = serpentine_rect("_s", n, sl, w, h);
                for ch in sbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── TrifurcationSerpentine [Tri] → 3 serpentine arms ──
            DesignTopology::TrifurcationSerpentine => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.33e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                let mut bp = trifurcation_rect(&self.id, trunk_len, branch_len, w, w, h);
                let sbp = serpentine_rect("_s", n, sl, w, h);
                for ch in sbp.channels { bp.channels.push(ch); }
                bp
            }

            // ── AdaptiveTree: variable-depth split tree → venturi at every leaf ──
            DesignTopology::AdaptiveTree { levels, split_types } => {
                if levels == 0 {
                    // Depth-0 acts like SingleVenturi
                    venturi_rect(&self.id, w, dt, h, tl)
                } else {
                    // Represent first split level + venturi at leaves.
                    // Full multi-level tree is rendered via to_channel_system().
                    let level_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (levels as f64 + 1.0);
                    let fan0 = if split_types & 1 == 0 { 2usize } else { 3 };
                    let mut bp = if fan0 == 2 {
                        bifurcation_rect(&self.id, level_len, level_len, w, w, h)
                    } else {
                        trifurcation_rect(&self.id, level_len, level_len, w, w, h)
                    };
                    // Append representative venturi channels
                    let vbp = venturi_rect("_v", w, dt, h, tl);
                    for ch in vbp.channels { bp.channels.push(ch); }
                    bp
                }
            }
        }
    }

    /// Build a [`cfd_schematics::geometry::ChannelSystem`] for 2D schematic
    /// visualisation using `cfd-schematics`.
    ///
    /// All physical units are converted from metres (candidate fields) to
    /// millimetres (the schematics API).  The generated system mirrors the
    /// x/y symmetry shown in existing `cfd-schematics` examples.
    ///
    /// Returns a fully-populated `ChannelSystem` that can be passed directly
    /// to [`cfd_schematics::plot_geometry`].
    #[must_use]
    pub fn to_channel_system(&self) -> cfd_schematics::geometry::ChannelSystem {
        use cfd_schematics::{
            config::{
                ArcConfig, ChannelTypeConfig, FrustumConfig, GeometryConfig, SerpentineConfig,
                TaperProfile,
            },
            geometry::{create_geometry, SplitType},
        };

        let w_mm = self.channel_width_m * 1e3;  // e.g. 2.0 mm
        let h_mm = self.channel_height_m * 1e3; // e.g. 0.5 mm
        // Clamp throat to valid FrustumConfig range: strictly inside (0, w_mm)
        let dt_mm = (self.throat_diameter_m * 1e3)
            .max(0.06)
            .min(w_mm * 0.9);

        let box_dims = (TREATMENT_WIDTH_MM, TREATMENT_HEIGHT_MM);

        let gc = GeometryConfig {
            wall_clearance: 2.0,
            channel_width:  w_mm,
            channel_height: h_mm,
            ..Default::default()
        };

        let frustum = || FrustumConfig {
            inlet_width:     w_mm,
            throat_width:    dt_mm,
            outlet_width:    w_mm,
            taper_profile:   TaperProfile::Smooth,
            smoothness:      50,
            throat_position: 0.5,
        };

        match self.topology {
            DesignTopology::SingleVenturi
            | DesignTopology::CellSeparationVenturi
            | DesignTopology::WbcCancerSeparationVenturi => {
                create_geometry(box_dims, &[], &gc, &ChannelTypeConfig::AllFrustum(frustum()))
            }
            DesignTopology::SerpentineGrid => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
            ),
            DesignTopology::VenturiSerpentine => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::Adaptive {
                    serpentine_config: SerpentineConfig::default(),
                    arc_config: ArcConfig::default(),
                    frustum_config: frustum(),
                },
            ),
            DesignTopology::BifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::TrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::DoubleBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::TripleBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Bifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::DoubleTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::BifurcationTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            // Serial double venturi: single straight path — same schematic as SingleVenturi
            DesignTopology::SerialDoubleVenturi => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllFrustum(frustum()),
            ),
            DesignTopology::BifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
            ),
            DesignTopology::TrifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(SerpentineConfig::default()),
            ),
            // Adaptive tree: decode split sequence from packed bits, use frustum at every leaf
            DesignTopology::AdaptiveTree { levels, split_types } => {
                let splits: Vec<SplitType> = (0..levels as usize)
                    .map(|i| {
                        if (split_types >> i) & 1 == 0 {
                            SplitType::Bifurcation
                        } else {
                            SplitType::Trifurcation
                        }
                    })
                    .collect();
                create_geometry(box_dims, &splits, &gc, &ChannelTypeConfig::AllFrustum(frustum()))
            }
        }
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

/// Build the complete parametric sweep over all 14 topology families.
///
/// Approximately **1134** candidates are generated before physics-based
/// filtering (exact count depends on which dimensions each topology uses:
/// topologies without a venturi skip the throat-diameter sweep; topologies
/// without a serpentine skip the segment-count sweep).
pub fn build_candidate_space() -> Vec<DesignCandidate> {
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
