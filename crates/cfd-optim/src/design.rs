//! Design topology families and candidate parameter space for SDT optimization.
//!
//! Eighteen topology families are evaluated (15 fixed + 3 leukapheresis + 1 GA-only adaptive):
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
//! | 15 | `AsymmetricBifurcationSerpentine`  |  0 | full 6×6 (2 arms) |
//! | 16 | `ConstrictionExpansionArray`       |  0 | full (N cycles)   |
//! | 17 | `SpiralSerpentine`                 |  0 | full (N turns)    |
//! | 18 | `ParallelMicrochannelArray`        |  0 | full (N channels) |
//! | GA | `AdaptiveTree`                     | N | full 6×6 grid     |
//!
//! `AdaptiveTree` is GA-only: depth (0–4) and per-level split type (Bi / Tri)
//! are encoded in genome genes 8–12 and constrained so the leaf channel width
//! stays ≥ 150 µm.
//!
//! Topologies 16–18 (`ConstrictionExpansionArray`, `SpiralSerpentine`,
//! `ParallelMicrochannelArray`) are designed for leukapheresis / cell separation
//! and are included in the GA evolutionary search (`ALL_EVO_TOPOLOGIES`).

use cfd_schematics::{
    interface::presets::{
        asymmetric_bifurcation_serpentine_rect, bifurcation_serpentine_rect,
        bifurcation_trifurcation_venturi_rect, bifurcation_venturi_rect,
        cascade_center_trifurcation_rect, cell_separation_rect,
        constriction_expansion_array_rect, double_bifurcation_venturi_rect,
        double_trifurcation_venturi_rect, parallel_microchannel_array_rect,
        quad_trifurcation_venturi_rect, serial_double_venturi_rect, serpentine_rect,
        spiral_channel_rect, trifurcation_bifurcation_bifurcation_venturi_rect,
        trifurcation_bifurcation_venturi_rect, trifurcation_serpentine_rect,
        trifurcation_venturi_rect, triple_bifurcation_venturi_rect,
        triple_trifurcation_venturi_rect, venturi_rect, venturi_serpentine_rect,
    },
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

    /// 1 trifurcation + 1 bifurcation → 6 parallel venturi throats (3×2 grid).
    ///
    /// Complement to `BifurcationTrifurcationVenturi`: trifurcation at the first
    /// level distributes across three column-bands; each band then bifurcates.
    TrifurcationBifurcationVenturi,

    /// 3-level trifurcation tree → 27 parallel venturi throats (3×3×3 grid).
    ///
    /// Channel widths **scale** at every level: center arm = `trifurcation_center_frac`
    /// × parent width; each peripheral arm = `(1 − frac) / 2` × parent width.
    /// At millifluidic inlet widths with `frac ≈ 1/3`, the level-3 terminal channels
    /// reach D_h < 200 µm — inertial-focusing range for WBCs and cancer cells.
    TripleTrifurcationVenturi,

    /// Trifurcation → Bifurcation → Bifurcation → 12 parallel venturi throats (3×2×2 grid).
    ///
    /// Mixed 3-level tree; channel widths scaled at each split using
    /// `trifurcation_center_frac` at the trifurcation stage and equal halving
    /// at each bifurcation stage.
    TrifurcationBifurcationBifurcationVenturi,

    /// 4-level trifurcation tree → 81 parallel venturi throats (3^4 grid).
    ///
    /// Maximum treatment-zone coverage of any fixed topology.  Width scaling at
    /// all four levels produces extreme terminal-channel narrowing at high `frac`.
    QuadTrifurcationVenturi,

    /// CIF-inspired cascade: only the center arm is re-trifurcated at each stage.
    ///
    /// Level 1: inlet → Tri(asym) → center_arm + L1_bypass + R1_bypass
    /// Level 2: center_arm → Tri(asym) → center_center + L2_bypass + R2_bypass
    /// Level N (optional): repeat on center_center
    ///
    /// After the cascade the deepest center arm passes through a **venturi throat**
    /// (cancer-enriched, cavitation treatment), then all arms — venturi stream and
    /// all bypass channels — **converge at a single outlet**.
    ///
    /// Zweifach-Fung junction routing concentrates large/stiff cells (cancer, WBC)
    /// into the center arm with each cascade level; peripheral bypass channels carry
    /// the RBC-enriched fraction through wide, low-shear paths, minimising haemolysis.
    ///
    /// Hemolysis and PAI are computed using the **local hematocrit** at the venturi
    /// (lower than feed hematocrit because RBCs are preferentially diverted to bypass),
    /// rewarding designs that maximally protect RBCs while targeting cancer cells.
    CascadeCenterTrifurcationSeparator {
        /// Number of cascade trifurcation levels (1–3).
        n_levels: u8,
    },

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

    /// Asymmetric bifurcation → wide serpentine arm (cancer/WBC) + narrow serpentine bypass (RBC).
    ///
    /// Exploits the **Zweifach-Fung bifurcation law**: at a T-junction, large stiff
    /// cells (cancer cells, WBCs) preferentially enter the arm carrying the greater
    /// volumetric flow (the wider, lower-resistance arm).  Small deformable RBCs
    /// distribute more uniformly and partly enter the narrow bypass arm.
    ///
    /// Wide arm: `channel_width_m`; narrow arm: `channel_width_m × 0.5`.
    /// Both arms are serpentine wave channels (Dean-flow focusing within each arm).
    AsymmetricBifurcationSerpentine,

    // ── Leukapheresis / inertial cell separation ─────────────────────────────

    /// N alternating wide→narrow constriction-expansion cycles.
    ///
    /// Each constriction forces cells to re-equilibrate laterally.
    /// Cumulative WBC margination enhancement: `E_N = 1 − (1−E_1)^n_cycles`.
    ///
    /// Physics: repeated inertial re-equilibration (Wu Z. et al., 2019, *Lab Chip*).
    /// Effective for WBC separation at moderate HCT; each cycle adds margination gain.
    ConstrictionExpansionArray {
        /// Number of wide→narrow cycles (2–20 typical range).
        n_cycles: usize,
    },

    /// Tight spiral channel (N full 360° turns, high Dean number).
    ///
    /// Continuously curved path creates stable secondary Dean-flow circulation
    /// that separates cells by size and deformability.  High Dean number
    /// (De = Re·√(D_h / 2R)) separates large stiff cells from small deformable RBCs.
    ///
    /// Physics: Nivedita et al. (2017, *Sci. Rep.*) inertial spiral design.
    SpiralSerpentine {
        /// Number of full spiral turns (2–20 typical range).
        n_turns: usize,
    },

    /// N identical parallel microchannels for clinical throughput scaling.
    ///
    /// Achieves micro-scale inertial focusing (`D_h < 150 µm`, `κ_WBC > 0.15`)
    /// at clinical flow rates (≥ 10 mL/min) by replicating N identical units in
    /// the millifluidic chip shell.  Each channel carries `Q_total / N`.
    ///
    /// Manufacturing: N rectangular channels etched side-by-side.
    ParallelMicrochannelArray {
        /// Number of parallel microchannels (10–500 typical range).
        n_channels: usize,
    },

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
            Self::BifurcationTrifurcationVenturi            => "[Bi,Tri] → 6 Parallel Venturis (2×3 grid)",
            Self::TrifurcationBifurcationVenturi            => "[Tri,Bi] → 6 Parallel Venturis (3×2 grid)",
            Self::TripleTrifurcationVenturi                 => "[Tri,Tri,Tri] → 27 Scaled-Width Venturis",
            Self::TrifurcationBifurcationBifurcationVenturi => "[Tri,Bi,Bi] → 12 Scaled-Width Venturis",
            Self::QuadTrifurcationVenturi                   => "[Tri,Tri,Tri,Tri] → 81 Scaled-Width Venturis",
            Self::CascadeCenterTrifurcationSeparator { .. } => "Cascade Center Trifurcation (CIF-inspired, venturi centre)",
            Self::SerialDoubleVenturi            => "2× Venturi in Series (double cavitation)",
            Self::BifurcationSerpentine              => "Bifurcation → 2 Serpentine Arms",
            Self::TrifurcationSerpentine             => "Trifurcation → 3 Serpentine Arms",
            Self::AsymmetricBifurcationSerpentine    => "Asymmetric Bifurcation → Wide+Narrow Serpentine (Zweifach-Fung)",
            Self::ConstrictionExpansionArray { .. }  => "Constriction-Expansion Array (WBC Margination)",
            Self::SpiralSerpentine { .. }            => "Spiral Serpentine Channel (Dean-Flow Separation)",
            Self::ParallelMicrochannelArray { .. }   => "Parallel Microchannel Array (Clinical Throughput)",
            Self::AdaptiveTree { .. }                => "Adaptive Split Tree",
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
            Self::BifurcationTrifurcationVenturi            => "B6",
            Self::TrifurcationBifurcationVenturi            => "TB",
            Self::TripleTrifurcationVenturi                 => "T27",
            Self::TrifurcationBifurcationBifurcationVenturi => "TBB",
            Self::QuadTrifurcationVenturi                   => "T81",
            Self::CascadeCenterTrifurcationSeparator { .. } => "CCT",
            Self::SerialDoubleVenturi            => "S2",
            Self::BifurcationSerpentine              => "BS",
            Self::TrifurcationSerpentine             => "TS",
            Self::AsymmetricBifurcationSerpentine    => "AB",
            Self::ConstrictionExpansionArray { .. }  => "CE",
            Self::SpiralSerpentine { .. }            => "SP",
            Self::ParallelMicrochannelArray { .. }   => "PM",
            Self::AdaptiveTree { .. }                => "AT",
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
                | Self::TrifurcationBifurcationVenturi
                | Self::TripleTrifurcationVenturi
                | Self::TrifurcationBifurcationBifurcationVenturi
                | Self::QuadTrifurcationVenturi
                | Self::CascadeCenterTrifurcationSeparator { .. }
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
            | Self::TrifurcationBifurcationVenturi
            | Self::TripleTrifurcationVenturi
            | Self::TrifurcationBifurcationBifurcationVenturi
            | Self::QuadTrifurcationVenturi
            | Self::CascadeCenterTrifurcationSeparator { .. }
            | Self::BifurcationSerpentine
            | Self::TrifurcationSerpentine
            | Self::ConstrictionExpansionArray { .. }
            | Self::SpiralSerpentine { .. }
            | Self::ParallelMicrochannelArray { .. } => true,
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
                | Self::AsymmetricBifurcationSerpentine
                | Self::ConstrictionExpansionArray { .. }
                | Self::SpiralSerpentine { .. }
                | Self::ParallelMicrochannelArray { .. }
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
            Self::BifurcationTrifurcationVenturi            => 6,
            Self::TrifurcationBifurcationVenturi            => 6,
            Self::TripleTrifurcationVenturi                 => 27,
            Self::TrifurcationBifurcationBifurcationVenturi => 12,
            Self::QuadTrifurcationVenturi                   => 81,
            Self::CascadeCenterTrifurcationSeparator { .. } => 1, // single centre-arm venturi
            // Serial: 2 in series on one path
            Self::SerialDoubleVenturi            => 2,
            Self::BifurcationSerpentine              => 0,
            Self::TrifurcationSerpentine             => 0,
            Self::AsymmetricBifurcationSerpentine    => 0,
            Self::ConstrictionExpansionArray { .. }  => 0,
            Self::SpiralSerpentine { .. }            => 0,
            Self::ParallelMicrochannelArray { .. }   => 0,
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

    /// Number of physical outlet ports on the device.
    ///
    /// All topologies produced by the composite closed-loop presets converge to
    /// **exactly one outlet**.  Splits at bifurcation / trifurcation junctions
    /// are balanced by corresponding merge junctions before the single outlet.
    pub fn outlet_count(self) -> usize {
        1
    }

    /// Number of parallel serpentine arms (for per-arm flow-rate calculations).
    ///
    /// `1` for all topologies except `BifurcationSerpentine` (2) and
    /// `TrifurcationSerpentine` (3), where the flow splits into equal parallel
    /// serpentine arms.
    pub fn serpentine_arm_count(self) -> usize {
        match self {
            Self::BifurcationSerpentine              => 2,
            Self::TrifurcationSerpentine             => 3,
            Self::AsymmetricBifurcationSerpentine    => 2, // wide + narrow arms
            Self::ParallelMicrochannelArray { n_channels } => n_channels,
            _                                        => 1,
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
            | Self::TrifurcationSerpentine
            | Self::AsymmetricBifurcationSerpentine => 1.0,
            // Center channel covers cancer cell treatment wells (≈ half the plate)
            Self::CellSeparationVenturi => 0.5,
            // Center channel covers WBC+cancer treatment wells (≈ half the plate)
            Self::WbcCancerSeparationVenturi => 0.5,
            // Multi-level tree topologies distribute over the full 6×6 zone
            Self::DoubleBifurcationVenturi
            | Self::TripleBifurcationVenturi
            | Self::DoubleTrifurcationVenturi
            | Self::BifurcationTrifurcationVenturi
            | Self::TrifurcationBifurcationVenturi
            | Self::TripleTrifurcationVenturi
            | Self::TrifurcationBifurcationBifurcationVenturi
            | Self::QuadTrifurcationVenturi => 1.0,
            // CCT: center arm covers ~half the treatment zone; bypass arms cover the other half
            Self::CascadeCenterTrifurcationSeparator { .. } => 1.0,
            // Leukapheresis topologies: cover full channel length
            Self::ConstrictionExpansionArray { .. }
            | Self::SpiralSerpentine { .. }
            | Self::ParallelMicrochannelArray { .. } => 1.0,
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
    /// Feed hematocrit (RBC volume fraction) at device inlet [0.0–0.45].
    ///
    /// Used by the cell-cell interaction model for leukapheresis topologies
    /// (`ConstrictionExpansionArray`, `SpiralSerpentine`, `ParallelMicrochannelArray`).
    /// Default for SDT / cancer separation: `0.45` (whole blood).
    /// Post-dilution leukapheresis: typically `0.01–0.10`.
    pub feed_hematocrit: f64,

    /// Fraction of the parent-arm width allocated to the **center arm** at each
    /// trifurcation split.  Each peripheral arm receives `(1 − frac) / 2`.
    ///
    /// - `1.0 / 3.0` (default) → symmetric split, all arms equal width.
    /// - `> 1/3` → center-biased: larger center flow fraction for stronger
    ///   Zweifach-Fung routing of large/stiff cells (cancer, WBC) toward center.
    ///
    /// Applies to `TripleTrifurcationVenturi`, `TrifurcationBifurcationBifurcationVenturi`,
    /// `QuadTrifurcationVenturi`, `CascadeCenterTrifurcationSeparator`, and all existing
    /// trifurcation-bearing topologies.  Ignored by bifurcation-only and serpentine topologies.
    #[serde(default = "default_tri_frac")]
    pub trifurcation_center_frac: f64,
}

fn default_tri_frac() -> f64 {
    1.0 / 3.0
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
    ///
    /// For serpentine topologies this is the per-arm velocity (each arm carries
    /// `Q / serpentine_arm_count()`).  For all other topologies it is the inlet
    /// trunk velocity.
    #[inline]
    pub fn channel_velocity_m_s(&self) -> f64 {
        let n_arms = self.topology.serpentine_arm_count().max(1);
        let q_per  = self.flow_rate_m3_s / n_arms as f64;
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
            // ── Single venturi ──
            DesignTopology::SingleVenturi => venturi_rect(&self.id, w, dt, h, tl),

            // ── Pure serpentine grid ──
            DesignTopology::SerpentineGrid => serpentine_rect(&self.id, n, sl, w, h),

            // ── Venturi → serpentine (series, closed-loop) ──
            DesignTopology::VenturiSerpentine => {
                venturi_serpentine_rect(&self.id, w, dt, h, tl, n, sl)
            }

            // ── Bifurcation + venturi in each branch (closed-loop) ──
            DesignTopology::BifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                bifurcation_venturi_rect(&self.id, trunk_len, w, dt, h, tl)
            }

            // ── Trifurcation + venturi in each branch (closed-loop) ──
            DesignTopology::TrifurcationVenturi => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.50e-3;
                trifurcation_venturi_rect(&self.id, trunk_len, w, dt, h, tl)
            }

            // ── Cell separation: center venturi + peripheral bypass (closed-loop) ──
            DesignTopology::CellSeparationVenturi
            | DesignTopology::WbcCancerSeparationVenturi => {
                let approach_len = TREATMENT_HEIGHT_MM * 0.5e-3;
                cell_separation_rect(&self.id, approach_len, w, dt, h, tl, approach_len)
            }

            // ── DoubleBifurcation [Bi,Bi] → 4 parallel venturis (closed-loop) ──
            DesignTopology::DoubleBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                double_bifurcation_venturi_rect(&self.id, trunk_len, branch1_len, w, dt, h, tl)
            }

            // ── TripleBifurcation [Bi,Bi,Bi] → 8 parallel venturis (closed-loop) ──
            DesignTopology::TripleBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
                triple_bifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, branch2_len, w, dt, h, tl,
                )
            }

            // ── DoubleTrifurcation [Tri,Tri] → 9 parallel venturis (closed-loop) ──
            DesignTopology::DoubleTrifurcationVenturi => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                double_trifurcation_venturi_rect(&self.id, trunk_len, branch_len, w, dt, h, tl)
            }

            // ── BifurcationTrifurcation [Bi,Tri] → 6 parallel venturis (closed-loop) ──
            DesignTopology::BifurcationTrifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                bifurcation_trifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, w, dt, h, tl,
                )
            }

            // ── TrifurcationBifurcation [Tri,Bi] → 6 parallel venturis (closed-loop) ──
            DesignTopology::TrifurcationBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.25e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
                trifurcation_bifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, w, dt, h, tl,
                )
            }

            // ── TripleTrifurcation [Tri,Tri,Tri] → 27 scaled-width venturis ──
            DesignTopology::TripleTrifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.09e-3;
                let frac = self.trifurcation_center_frac;
                triple_trifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, branch2_len, w, frac, dt, h, tl,
                )
            }

            // ── TrifurcationBifurcationBifurcation [Tri,Bi,Bi] → 12 scaled-width venturis ──
            DesignTopology::TrifurcationBifurcationBifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.18e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.11e-3;
                let frac = self.trifurcation_center_frac;
                trifurcation_bifurcation_bifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, branch2_len, w, frac, dt, h, tl,
                )
            }

            // ── QuadTrifurcation [Tri^4] → 81 scaled-width venturis ──
            DesignTopology::QuadTrifurcationVenturi => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.12e-3;
                let branch1_len = TREATMENT_HEIGHT_MM * 0.10e-3;
                let branch2_len = TREATMENT_HEIGHT_MM * 0.08e-3;
                let branch3_len = TREATMENT_HEIGHT_MM * 0.06e-3;
                let frac = self.trifurcation_center_frac;
                quad_trifurcation_venturi_rect(
                    &self.id, trunk_len, branch1_len, branch2_len, branch3_len, w, frac, dt, h, tl,
                )
            }

            // ── CascadeCenterTrifurcationSeparator → single venturi on center arm ──
            DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
                let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
                let branch_len  = TREATMENT_HEIGHT_MM * 0.15e-3;
                let frac = self.trifurcation_center_frac;
                cascade_center_trifurcation_rect(
                    &self.id, trunk_len, branch_len, n_levels, w, frac, dt, tl, h,
                )
            }

            // ── Two venturis in series on one path (closed-loop) ──
            DesignTopology::SerialDoubleVenturi => {
                let inter_len = TREATMENT_HEIGHT_MM * 0.40e-3; // 18 mm inter-venturi
                serial_double_venturi_rect(&self.id, w, dt, h, tl, inter_len)
            }

            // ── Bifurcation + serpentine in each arm (closed-loop) ──
            DesignTopology::BifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
                bifurcation_serpentine_rect(&self.id, trunk_len, n, sl, w, h)
            }

            // ── Trifurcation + serpentine in each arm (closed-loop) ──
            DesignTopology::TrifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 0.33e-3;
                trifurcation_serpentine_rect(&self.id, trunk_len, n, sl, w, h)
            }

            // ── Asymmetric bifurcation → wide serpentine (cancer/WBC) + narrow bypass (RBC) ──
            DesignTopology::AsymmetricBifurcationSerpentine => {
                let trunk_len = TREATMENT_HEIGHT_MM * 1e-3 / 6.0;
                asymmetric_bifurcation_serpentine_rect(&self.id, trunk_len, n, sl, w, h)
            }

            // ── Constriction-expansion array (wide→narrow cycles, WBC margination) ──
            DesignTopology::ConstrictionExpansionArray { n_cycles } => {
                let wide_w   = w;
                let narrow_w = w * 0.40;
                // Each cycle has a wide + narrow segment; split treatment length evenly
                let seg_len  = TREATMENT_HEIGHT_MM * 1e-3 / (n_cycles as f64 * 2.0);
                constriction_expansion_array_rect(
                    &self.id, n_cycles, seg_len, seg_len * 0.5, wide_w, narrow_w, h,
                )
            }

            // ── Spiral serpentine (N turns, high Dean number, WBC/RBC separation) ──
            DesignTopology::SpiralSerpentine { n_turns } => {
                let turn_len = TREATMENT_HEIGHT_MM * 1e-3 / n_turns as f64;
                spiral_channel_rect(&self.id, n_turns, turn_len, w, h)
            }

            // ── Parallel microchannel array (N identical channels, clinical throughput) ──
            DesignTopology::ParallelMicrochannelArray { n_channels } => {
                let ch_len = TREATMENT_HEIGHT_MM * 1e-3;
                parallel_microchannel_array_rect(&self.id, n_channels, ch_len, w, h)
            }

            // ── AdaptiveTree: variable-depth split tree → venturi at every leaf ──
            DesignTopology::AdaptiveTree { levels, split_types } => {
                let trunk_len  = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (levels as f64 + 1.0);
                let branch_len = trunk_len * 0.8;
                let branch2_len = branch_len * 0.8; // for 3-level trees
                match levels {
                    0 => venturi_rect(&self.id, w, dt, h, tl),
                    1 => {
                        let level_len = TREATMENT_HEIGHT_MM * 1e-3 / 4.0;
                        if split_types & 1 == 0 {
                            bifurcation_venturi_rect(&self.id, level_len, w, dt, h, tl)
                        } else {
                            trifurcation_venturi_rect(&self.id, level_len, w, dt, h, tl)
                        }
                    }
                    2 => match split_types & 0b11 {
                        // Bi→Bi = 4 leaves
                        0b00 => double_bifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Bi→Tri = 6 leaves
                        0b10 => bifurcation_trifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Tri→Bi = 6 leaves (new preset)
                        0b01 => trifurcation_bifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                        // Tri→Tri = 9 leaves
                        _ => double_trifurcation_venturi_rect(
                            &self.id, trunk_len, branch_len, w, dt, h, tl,
                        ),
                    },
                    // levels 3–4: approximate with best-fit 3-level preset
                    _ => {
                        if split_types & 1 == 0 {
                            triple_bifurcation_venturi_rect(
                                &self.id, trunk_len, branch_len, branch2_len, w, dt, h, tl,
                            )
                        } else {
                            double_trifurcation_venturi_rect(
                                &self.id, trunk_len, branch_len, w, dt, h, tl,
                            )
                        }
                    }
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

        // ── Wave-shape configs ───────────────────────────────────────────
        // Sine wave — smooth Dean-flow visualisation for tree / venturi topologies.
        let sine_wave = SerpentineConfig {
            fill_factor: 0.75,
            wave_density_factor: 2.5,
            gaussian_width_factor: 4.0,
            ..SerpentineConfig::default()
        };

        // Square wave — high-density sharp wave for cell-size separation designs.
        let square_wave = SerpentineConfig {
            fill_factor: 0.80,
            wave_density_factor: 4.0,
            gaussian_width_factor: 4.0,
            ..SerpentineConfig::default()
        }
        .with_square_wave();

        match self.topology {
            // ── Single venturi: sine wave, no splits ──
            DesignTopology::SingleVenturi | DesignTopology::SerialDoubleVenturi => {
                create_geometry(
                    box_dims,
                    &[],
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(sine_wave),
                )
            }

            // ── Cell separation topologies: square wave + asymmetric bifurcation ──
            DesignTopology::CellSeparationVenturi | DesignTopology::WbcCancerSeparationVenturi => {
                create_geometry(
                    box_dims,
                    &[SplitType::AsymmetricBifurcation { ratio: 0.67 }],
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(square_wave),
                )
            }

            // ── Pure serpentine grid ──
            DesignTopology::SerpentineGrid => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Venturi + serpentine (series): adaptive (sine + frustum) ──
            DesignTopology::VenturiSerpentine => {
                let frustum_cfg = FrustumConfig {
                    inlet_width: w_mm,
                    throat_width: dt_mm,
                    outlet_width: w_mm,
                    taper_profile: TaperProfile::Smooth,
                    smoothness: 50,
                    throat_position: 0.5,
                };
                create_geometry(
                    box_dims,
                    &[],
                    &gc,
                    &ChannelTypeConfig::Adaptive {
                        serpentine_config: sine_wave,
                        arc_config: ArcConfig::default(),
                        frustum_config: frustum_cfg,
                    },
                )
            }

            // ── Binary venturi trees → sine wave ──
            DesignTopology::BifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::DoubleBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TripleBifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Bifurcation, SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::DoubleTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Trifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::BifurcationTrifurcationVenturi => create_geometry(
                box_dims,
                &[SplitType::Bifurcation, SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Symmetric serpentine arms → sine wave ──
            DesignTopology::BifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Bifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),
            DesignTopology::TrifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::Trifurcation],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Asymmetric bifurcation → square wave (cell-size separation) ──
            DesignTopology::AsymmetricBifurcationSerpentine => create_geometry(
                box_dims,
                &[SplitType::AsymmetricBifurcation { ratio: 0.67 }],
                &gc,
                &ChannelTypeConfig::AllSerpentine(square_wave),
            ),

            // ── Constriction-expansion array: square wave (WBC margination) ──
            DesignTopology::ConstrictionExpansionArray { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(square_wave),
            ),

            // ── Spiral serpentine: sine wave (Dean-flow separation) ──
            DesignTopology::SpiralSerpentine { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Parallel microchannel array: sine wave ──
            DesignTopology::ParallelMicrochannelArray { .. } => create_geometry(
                box_dims,
                &[],
                &gc,
                &ChannelTypeConfig::AllSerpentine(sine_wave),
            ),

            // ── Adaptive tree: decode split sequence from packed bits → sine wave ──
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
                create_geometry(
                    box_dims,
                    &splits,
                    &gc,
                    &ChannelTypeConfig::AllSerpentine(sine_wave),
                )
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

/// Build the complete parametric sweep over all topology families.
///
/// Approximately **1500+** candidates are generated before physics-based
/// filtering.  The exact count depends on which dimensions each topology uses:
/// topologies without a venturi skip the throat-diameter sweep; topologies
/// without a serpentine skip the segment-count sweep.  New deep-trifurcation
/// topologies (T1–T5) add an extra loop over `TRIFURCATION_CENTER_FRACS`.
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
        DesignTopology::AsymmetricBifurcationSerpentine,
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
                                feed_hematocrit: 0.45, // whole blood default
                                trifurcation_center_frac: 1.0 / 3.0, // symmetric default
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
    let leuka_topologies = [
        DesignTopology::ConstrictionExpansionArray { n_cycles: 10 },
        DesignTopology::SpiralSerpentine { n_turns: 8 },
        DesignTopology::ParallelMicrochannelArray { n_channels: 100 },
    ];

    for &topology in &leuka_topologies {
        for &q in &LEUKA_FLOW_RATES_M3_S {
            for &gauge in &INLET_GAUGES_PA {
                for &w_ch in &LEUKA_CHANNEL_WIDTHS_M {
                    idx += 1;
                    let id = format!(
                        "{:04}-{}-LK-q{:.0}ulm-g{:.0}kPa-w{:.0}um",
                        idx,
                        topology.short(),
                        q * 6e10,       // m³/s → µL/min (2.5e-9 → 150, 1.5e-8 → 900)
                        gauge * 1e-3,   // Pa → kPa
                        w_ch * 1e6,     // m → µm
                    );

                    candidates.push(DesignCandidate {
                        id,
                        topology,
                        flow_rate_m3_s: q,
                        inlet_gauge_pa: gauge,
                        throat_diameter_m: 0.0,         // no venturi throat
                        inlet_diameter_m: VENTURI_INLET_DIAM_M,
                        throat_length_m: 0.0,
                        channel_width_m: w_ch,
                        channel_height_m: LEUKA_CHANNEL_HEIGHT_M, // 60 µm fixed
                        serpentine_segments: 1,  // topology carries n_turns/n_cycles
                        segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                        bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                        feed_hematocrit: 0.04,   // 4% diluted for leukapheresis
                        trifurcation_center_frac: 1.0 / 3.0, // symmetric default
                    });
                }
            }
        }
    }

    // ── New deep-trifurcation topologies (millifluidic scale, width-scaled) ──
    //
    // T1: TrifurcationBifurcationVenturi — fix orphaned schematics preset; symmetric
    // T2: TripleTrifurcationVenturi       — 27 outlets; swept over center_frac
    // T3: TrifurcationBifurcationBifurcation — 12 outlets; swept over center_frac
    // T4: QuadTrifurcationVenturi         — 81 outlets; swept over center_frac
    // T5: CascadeCenterTrifurcationSeparator — 2-level CCT; swept over center_frac
    //
    // All use millifluidic channel dimensions (CHANNEL_WIDTHS_M, CHANNEL_HEIGHT_M)
    // and whole-blood feed hematocrit (0.45).  T2/T3/T4/T5 additionally loop over
    // TRIFURCATION_CENTER_FRACS to vary the width-fraction sweep.

    // T1: symmetric (no center_frac variation)
    for &q in &FLOW_RATES_M3_S {
        for &gauge in &INLET_GAUGES_PA {
            for &d_throat in &THROAT_DIAMETERS_M {
                for &w_ch in &CHANNEL_WIDTHS_M {
                    idx += 1;
                    let throat_len = d_throat * 2.0;
                    let id = format!(
                        "{:04}-TB-q{:.0}ml-g{:.0}kPa-dt{:.0}um-w{:.0}um",
                        idx, q * 6e7, gauge * 1e-3, d_throat * 1e6, w_ch * 1e6,
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
                        channel_height_m: CHANNEL_HEIGHT_M,
                        serpentine_segments: 1,
                        segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                        bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                        feed_hematocrit: 0.45,
                        trifurcation_center_frac: 1.0 / 3.0,
                    });
                }
            }
        }
    }

    // T2/T3/T4/T5: swept over center_frac
    let scaled_topologies: &[(DesignTopology, bool)] = &[
        (DesignTopology::TripleTrifurcationVenturi,                 true),  // has venturi
        (DesignTopology::TrifurcationBifurcationBifurcationVenturi, true),
        (DesignTopology::QuadTrifurcationVenturi,                   true),
        (DesignTopology::CascadeCenterTrifurcationSeparator { n_levels: 2 }, true),
    ];

    for &(topology, _) in scaled_topologies {
        for &center_frac in &TRIFURCATION_CENTER_FRACS {
            let frac_tag = (center_frac * 1000.0).round() as u32; // e.g. 333, 450, 550
            for &q in &FLOW_RATES_M3_S {
                for &gauge in &INLET_GAUGES_PA {
                    for &d_throat in &THROAT_DIAMETERS_M {
                        for &w_ch in &CHANNEL_WIDTHS_M {
                            idx += 1;
                            let throat_len = d_throat * 2.0;
                            let id = format!(
                                "{:04}-{}-cf{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-w{:.0}um",
                                idx,
                                topology.short(),
                                frac_tag,
                                q * 6e7,
                                gauge * 1e-3,
                                d_throat * 1e6,
                                w_ch * 1e6,
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
                                serpentine_segments: 1,
                                segment_length_m: TREATMENT_WIDTH_MM * 1e-3,
                                bend_radius_m: SERPENTINE_BEND_RADIUS_M,
                                feed_hematocrit: 0.45,
                                trifurcation_center_frac: center_frac,
                            });
                        }
                    }
                }
            }
        }
    }

    candidates
}
