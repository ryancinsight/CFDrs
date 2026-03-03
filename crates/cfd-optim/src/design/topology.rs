//! `DesignTopology` enum — millifluidic topology families.

use crate::constraints::TREATMENT_WELL_COUNT;
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
    /// `κ_WBC` > 0.07 (WBC inertial focusing threshold).
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
    /// reach `D_h` < 200 µm — inertial-focusing range for WBCs and cancer cells.
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
    /// Level 1: inlet → Tri(asym) → `center_arm` + `L1_bypass` + `R1_bypass`
    /// Level 2: `center_arm` → Tri(asym) → `center_center` + `L2_bypass` + `R2_bypass`
    /// Level N (optional): repeat on `center_center`
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

    /// Controlled incremental filtration (CIF): pre-trifurcation cascade,
    /// then terminal trifurcation + bifurcation skimming before one venturi.
    ///
    /// Sequence:
    /// 1. `n_pretri` center-only trifurcations for progressive RBC skimming.
    /// 2. Terminal trifurcation to intensify center-path enrichment.
    /// 3. Terminal bifurcation to route the high-flow treatment arm to venturi.
    ///
    /// All bypass arms merge back internally and the device still exposes a
    /// single external outlet.
    IncrementalFiltrationTriBiSeparator {
        /// Number of initial center-cascade trifurcation levels (1–3).
        n_pretri: u8,
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
    /// Wide arm: `channel_width_m`; narrow arm: `channel_width_m × asymmetric_narrow_frac`.
    /// Both arms are serpentine wave channels (Dean-flow focusing within each arm).
    AsymmetricBifurcationSerpentine,

    /// Asymmetric 3-stream trifurcation with center-only venturi for selective SDT.
    ///
    /// Three independent arm widths:
    /// - **Center arm** (cancer-enriched, `trifurcation_center_frac` × parent width)
    ///   → venturi cavitation treatment
    /// - **Left arm** (WBC collection, `trifurcation_left_frac` × parent width)
    ///   → HealthyBypass, collects WBCs
    /// - **Right arm** (`1 - center_frac - left_frac` × parent width)
    ///   → HealthyBypass, RBC waste stream (narrow → CFL margination)
    ///
    /// Independent left/right fractions enable 3-stream sorting: cancer treatment,
    /// WBC recovery, and RBC waste in a single-pass millifluidic chip.
    AsymmetricTrifurcationVenturi,

    /// Tri → Bi → Tri cascade with center-only venturi for ultra-selective SDT.
    ///
    /// Three progressive focusing stages ensure only the most cancer-enriched
    /// center stream reaches the venturi cavitation zone:
    /// - Stage 1: Trifurcation (L/R → HealthyBypass, center → stage 2)
    /// - Stage 2: Bifurcation (bypass → HealthyBypass, treatment → stage 3)
    /// - Stage 3: Trifurcation (L/R → HealthyBypass, center → venturi)
    ///
    /// `trifurcation_center_frac` controls stage-1 and stage-3 splits;
    /// `cif_terminal_bi_treat_frac` controls the stage-2 bifurcation.
    TriBiTriSelectiveVenturi,

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
    /// (De = `Re·√(D_h` / 2R)) separates large stiff cells from small deformable RBCs.
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
    /// Human-readable topology name.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::SingleVenturi => "Single Venturi",
            Self::BifurcationVenturi => "Bifurcation + Venturi",
            Self::TrifurcationVenturi => "Trifurcation + Venturi",
            Self::VenturiSerpentine => "Venturi → Serpentine",
            Self::SerpentineGrid => "Serpentine Grid",
            Self::CellSeparationVenturi => "Cell Separation (Cancer-targeted SDT)",
            Self::WbcCancerSeparationVenturi => "WBC+Cancer Separation + SDT",
            Self::DoubleBifurcationVenturi => "Double Bifurcation [Bi,Bi] → 4× Venturi",
            Self::TripleBifurcationVenturi => "Triple Bifurcation [Bi,Bi,Bi] → 8× Venturi",
            Self::DoubleTrifurcationVenturi => "Double Trifurcation [Tri,Tri] → 9× Venturi",
            Self::BifurcationTrifurcationVenturi => "Bifurcation → Trifurcation → 6× Venturi",
            Self::TrifurcationBifurcationVenturi => "Trifurcation → Bifurcation → 6× Venturi",
            Self::TripleTrifurcationVenturi => "Triple Trifurcation [Tri,Tri,Tri] → 27× Venturi",
            Self::TrifurcationBifurcationBifurcationVenturi => {
                "Trifurcation → Bifurcation → Bifurcation → 12× Venturi"
            }
            Self::QuadTrifurcationVenturi => "Quad Trifurcation [Tri^4] → 81× Venturi",
            Self::CascadeCenterTrifurcationSeparator { .. } => {
                "Cascade Center Trifurcation Separator"
            }
            Self::IncrementalFiltrationTriBiSeparator { .. } => {
                "Incremental Filtration Tri→Tri/Bi Separator"
            }
            Self::SerialDoubleVenturi => "Serial Double Venturi",
            Self::BifurcationSerpentine => "Bifurcation → 2× Serpentine",
            Self::TrifurcationSerpentine => "Trifurcation → 3× Serpentine",
            Self::AsymmetricBifurcationSerpentine => "Asymmetric Bifurcation → Serpentine",
            Self::AsymmetricTrifurcationVenturi => {
                "Asymmetric 3-Stream Trifurcation + Center Venturi"
            }
            Self::TriBiTriSelectiveVenturi => "Tri→Bi→Tri Selective Center Venturi",
            Self::ConstrictionExpansionArray { .. } => "Constriction-Expansion Array",
            Self::SpiralSerpentine { .. } => "Spiral Serpentine",
            Self::ParallelMicrochannelArray { .. } => "Parallel Microchannel Array",
            Self::AdaptiveTree { .. } => "Adaptive Tree",
        }
    }

    /// Short (≤3 char) topology code for candidate IDs.
    #[must_use]
    pub fn short(self) -> &'static str {
        match self {
            Self::SingleVenturi => "SV",
            Self::BifurcationVenturi => "BV",
            Self::TrifurcationVenturi => "TV",
            Self::VenturiSerpentine => "VS",
            Self::SerpentineGrid => "SG",
            Self::CellSeparationVenturi => "CS",
            Self::WbcCancerSeparationVenturi => "WC",
            Self::DoubleBifurcationVenturi => "D2",
            Self::TripleBifurcationVenturi => "D3",
            Self::DoubleTrifurcationVenturi => "T9",
            Self::BifurcationTrifurcationVenturi => "B6",
            Self::TrifurcationBifurcationVenturi => "TB",
            Self::TripleTrifurcationVenturi => "T27",
            Self::TrifurcationBifurcationBifurcationVenturi => "TBB",
            Self::QuadTrifurcationVenturi => "T81",
            Self::CascadeCenterTrifurcationSeparator { .. } => "CCT",
            Self::IncrementalFiltrationTriBiSeparator { .. } => "CIF",
            Self::SerialDoubleVenturi => "S2",
            Self::BifurcationSerpentine => "BS",
            Self::TrifurcationSerpentine => "TS",
            Self::AsymmetricBifurcationSerpentine => "AB",
            Self::AsymmetricTrifurcationVenturi => "ATV",
            Self::TriBiTriSelectiveVenturi => "TBT",
            Self::ConstrictionExpansionArray { .. } => "CE",
            Self::SpiralSerpentine { .. } => "SP",
            Self::ParallelMicrochannelArray { .. } => "PM",
            Self::AdaptiveTree { .. } => "AT",
        }
    }

    /// Returns `true` if this topology includes at least one venturi throat.
    ///
    /// `AdaptiveTree` always includes venturi throats at the leaf nodes.
    #[must_use]
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
                | Self::IncrementalFiltrationTriBiSeparator { .. }
                | Self::AsymmetricTrifurcationVenturi
                | Self::TriBiTriSelectiveVenturi
                | Self::SerialDoubleVenturi
                | Self::AdaptiveTree { .. }
        )
    }

    /// Returns `true` if this topology includes a distribution network
    /// (bifurcation tree, trifurcation tree, or serpentine).
    ///
    /// For `AdaptiveTree`, returns `true` only when `levels > 0` (depth-0 is
    /// equivalent to a single straight channel with a venturi).
    #[must_use]
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
            | Self::IncrementalFiltrationTriBiSeparator { .. }
            | Self::AsymmetricTrifurcationVenturi
            | Self::TriBiTriSelectiveVenturi
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
    #[must_use]
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
    #[must_use]
    pub fn venturi_count(self) -> usize {
        match self {
            Self::SingleVenturi => 1,
            Self::BifurcationVenturi => 2,
            Self::TrifurcationVenturi => 3,
            Self::VenturiSerpentine => 1,
            Self::SerpentineGrid => 0,
            Self::CellSeparationVenturi => 1,
            Self::WbcCancerSeparationVenturi => 1,
            Self::DoubleBifurcationVenturi => 4,
            Self::TripleBifurcationVenturi => 8,
            Self::DoubleTrifurcationVenturi => 9,
            Self::BifurcationTrifurcationVenturi => 6,
            Self::TrifurcationBifurcationVenturi => 6,
            Self::TripleTrifurcationVenturi => 27,
            Self::TrifurcationBifurcationBifurcationVenturi => 12,
            Self::QuadTrifurcationVenturi => 81,
            Self::CascadeCenterTrifurcationSeparator { .. } => 1, // single centre-arm venturi
            Self::IncrementalFiltrationTriBiSeparator { .. } => 1, // single treatment-arm venturi
            Self::AsymmetricTrifurcationVenturi => 1,             // single center-arm venturi
            Self::TriBiTriSelectiveVenturi => 1,                  // single center-path venturi
            // Serial: 2 in series on one path
            Self::SerialDoubleVenturi => 2,
            Self::BifurcationSerpentine => 0,
            Self::TrifurcationSerpentine => 0,
            Self::AsymmetricBifurcationSerpentine => 0,
            Self::ConstrictionExpansionArray { .. } => 0,
            Self::SpiralSerpentine { .. } => 0,
            Self::ParallelMicrochannelArray { .. } => 0,
            // Adaptive: multiply fan factors across all levels
            Self::AdaptiveTree {
                levels,
                split_types,
            } => {
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
    #[must_use]
    pub fn parallel_venturi_count(self) -> usize {
        match self {
            Self::SerialDoubleVenturi => 1, // full Q through each stage in series
            _ => self.venturi_count(),
        }
    }

    /// Number of serial venturi stages on one flow path.
    ///
    /// `1` for all topologies except `SerialDoubleVenturi` which has `2`.
    #[must_use]
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
    #[must_use]
    pub fn outlet_count(self) -> usize {
        1
    }

    /// Number of parallel serpentine arms (for per-arm flow-rate calculations).
    ///
    /// `1` for all topologies except `BifurcationSerpentine` (2) and
    /// `TrifurcationSerpentine` (3), where the flow splits into equal parallel
    /// serpentine arms.
    #[must_use]
    pub fn serpentine_arm_count(self) -> usize {
        match self {
            Self::BifurcationSerpentine => 2,
            Self::TrifurcationSerpentine => 3,
            Self::AsymmetricBifurcationSerpentine => 2, // wide + narrow arms
            Self::ParallelMicrochannelArray { n_channels } => n_channels,
            _ => 1,
        }
    }

    /// Number of terminal (leaf) branches in one half of the split tree.
    ///
    /// This is the product of all split factors.  The geometry generator
    /// mirrors the tree so the total parallel channels in the device is
    /// `2 × terminal_branch_count()`.  Used by [`adaptive_box_dims`] to
    /// compute足 the vertical footprint.
    #[must_use]
    pub fn terminal_branch_count(self) -> usize {
        match self {
            Self::SingleVenturi
            | Self::SerialDoubleVenturi
            | Self::VenturiSerpentine
            | Self::SerpentineGrid
            | Self::ConstrictionExpansionArray { .. }
            | Self::SpiralSerpentine { .. } => 1,

            Self::BifurcationVenturi
            | Self::BifurcationSerpentine
            | Self::AsymmetricBifurcationSerpentine
            | Self::CellSeparationVenturi
            | Self::WbcCancerSeparationVenturi => 2,

            Self::TrifurcationVenturi
            | Self::TrifurcationSerpentine
            | Self::AsymmetricTrifurcationVenturi => 3,

            Self::DoubleBifurcationVenturi => 4,
            Self::BifurcationTrifurcationVenturi
            | Self::TrifurcationBifurcationVenturi => 6,
            Self::TripleBifurcationVenturi => 8,
            Self::DoubleTrifurcationVenturi => 9,
            Self::TrifurcationBifurcationBifurcationVenturi => 12,
            Self::TriBiTriSelectiveVenturi => 18,  // 3 × 2 × 3
            Self::TripleTrifurcationVenturi => 27,
            Self::QuadTrifurcationVenturi => 81,

            Self::ParallelMicrochannelArray { n_channels } => n_channels,

            Self::CascadeCenterTrifurcationSeparator { n_levels } => {
                3usize.pow(n_levels as u32)
            }
            Self::IncrementalFiltrationTriBiSeparator { n_pretri } => {
                3usize.pow(n_pretri as u32 + 1) * 2
            }
            Self::AdaptiveTree {
                levels,
                split_types,
            } => {
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
    #[must_use]
    pub fn nominal_well_coverage(self) -> f64 {
        match self {
            // Single-point cavitation: covers ~1 well
            Self::SingleVenturi | Self::SerialDoubleVenturi => 1.0 / TREATMENT_WELL_COUNT as f64,
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
            // CIF tri→tri/bi: treatment + bypass paths span the full channelized zone
            Self::IncrementalFiltrationTriBiSeparator { .. } => 1.0,
            // Asymmetric trifurcation: 3-stream split across the full plate
            Self::AsymmetricTrifurcationVenturi => 1.0,
            // Tri→Bi→Tri: three-stage cascade spans the full channelized zone
            Self::TriBiTriSelectiveVenturi => 1.0,
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
