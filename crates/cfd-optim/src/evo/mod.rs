//! Evolutionary / genetic algorithm optimizer for millifluidic SDT devices.
//!
//! [`GeneticOptimizer`] searches the design space using a real-coded genetic
//! algorithm with tournament selection, Simulated Binary Crossover (SBX),
//! and polynomial mutation.
//!
//! Each individual is a [`MillifluidicGenome`] — a 20-gene real-valued vector
//! encoding one fully-specified [`DesignCandidate`].  The optimizer covers
//! **all 26 fixed topology families** plus the GA-only `AdaptiveTree`
//! (variable-depth split tree), giving **27 searchable families**.
//!
//! # Submodules
//!
//! | Module       | Purpose                                           |
//! |--------------|---------------------------------------------------|
//! | `decode`     | Genome → `DesignCandidate` mapping                |
//! | `operators`  | Tournament selection, SBX crossover, mutation      |
//! | `optimizer`  | `GeneticOptimizer` struct and `run()` entry point  |
//!
//! # References
//! - Deb, K. & Agrawal, R. B. (1995). Simulated binary crossover for continuous
//!   search space. *Complex Systems*, 9, 115–148.
//! - Deb, K. & Deb, D. (2012). Analysing mutation schemes for real-parameter
//!   genetic algorithms. *Int. J. Artif. Intell. Soft Comput.*, 4, 1–28.

mod decode;
mod operators;
mod optimizer;

pub use decode::decode_genome;
pub use optimizer::{EvolutionResult, GeneticOptimizer};

use crate::design::DesignTopology;

// ── All topology families available to the genetic search ───────────────────

pub(super) const ALL_EVO_TOPOLOGIES: [DesignTopology; 27] = [
    DesignTopology::SingleVenturi,                                      // 0
    DesignTopology::BifurcationVenturi,                                 // 1
    DesignTopology::TrifurcationVenturi,                                // 2
    DesignTopology::VenturiSerpentine,                                  // 3
    DesignTopology::SerpentineGrid,                                     // 4
    DesignTopology::CellSeparationVenturi,                              // 5
    DesignTopology::WbcCancerSeparationVenturi,                         // 6
    DesignTopology::DoubleBifurcationVenturi,                           // 7
    DesignTopology::TripleBifurcationVenturi,                           // 8
    DesignTopology::DoubleTrifurcationVenturi,                          // 9
    DesignTopology::BifurcationTrifurcationVenturi,                     // 10
    DesignTopology::SerialDoubleVenturi,                                // 11
    DesignTopology::BifurcationSerpentine,                              // 12
    DesignTopology::TrifurcationSerpentine,                             // 13
    DesignTopology::AsymmetricBifurcationSerpentine,                    // 14
    DesignTopology::ConstrictionExpansionArray { n_cycles: 10 },        // 15
    DesignTopology::SpiralSerpentine { n_turns: 8 },                    // 16
    DesignTopology::ParallelMicrochannelArray { n_channels: 100 },      // 17
    DesignTopology::TrifurcationBifurcationVenturi,                     // 18
    DesignTopology::TripleTrifurcationVenturi,                          // 19
    DesignTopology::TrifurcationBifurcationBifurcationVenturi,          // 20
    DesignTopology::QuadTrifurcationVenturi,                            // 21
    DesignTopology::CascadeCenterTrifurcationSeparator { n_levels: 2 }, // 22
    DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 2 }, // 23
    DesignTopology::AsymmetricTrifurcationVenturi,                      // 24
    DesignTopology::TriBiTriSelectiveVenturi,                           // 25
    DesignTopology::AdaptiveTree {
        levels: 0,
        split_types: 0,
    }, // 26
];

// ── Genome definition ─────────────────────────────────────────────────────────

/// Real-coded genome for a millifluidic design.
///
/// Gene layout (all normalised to `[0, 1]`):
///
/// | Index | Meaning |
/// |-------|---------|
/// | 0  | Topology index (continuous, rounded to int) |
/// | 1  | Flow rate |
/// | 2  | Inlet gauge pressure |
/// | 3  | Throat diameter (venturi only) |
/// | 4  | Channel width |
/// | 5  | Serpentine segment count |
/// | 6  | Segment length fraction |
/// | 7  | Bend radius fraction |
/// | 8  | Discrete topology parameter (indices 15–17, 22–24) |
/// | 9–12 | AdaptiveTree per-level split type (0 = Bi, 1 = Tri) |
/// | 13 | `trifurcation_center_frac` / CIF pretri center frac ∈ [0.25, 0.65] |
/// | 14 | CIF terminal tri center frac ∈ [0.25, 0.65] |
/// | 15 | CIF terminal bifurcation treat frac / TBT bi frac ∈ [0.50, 0.85] |
/// | 16 | `asymmetric_narrow_frac` ∈ [0.20, 0.70] (AsymmetricBifurcationSerpentine only) |
/// | 17 | Throat length factor ∈ [1.5, 15.0] × throat diameter (venturi only) |
/// | 18 | Channel height log-linear ∈ [0.3 mm, 3.0 mm] (millifluidic only) |
/// | 19 | `trifurcation_left_frac` ∈ [0.08, 0.50] (AsymmetricTrifurcationVenturi only) |
#[derive(Debug, Clone)]
pub struct MillifluidicGenome {
    /// Normalised gene values ∈ [0, 1].
    pub genes: Vec<f64>,
}

/// Number of genes per individual.
pub(super) const N_GENES: usize = 20;

/// Minimum leaf channel width for `AdaptiveTree` depth constraint [m].
pub(super) const MIN_ADAPTIVE_CH_M: f64 = 150e-6;
