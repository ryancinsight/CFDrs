//! Physics-based metric computation for SDT millifluidic design candidates.
//!
//! Each [`DesignCandidate`] is evaluated against the 1-D physics models from
//! `cfd-1d` (venturi + serpentine) and the Casson blood rheology from
//! `cfd-core`.  The resulting [`SdtMetrics`] struct gathers all values needed
//! for multi-objective scoring.
//!
//! ## Sub-modules
//!
//! | Module | Responsibility |
//! |--------|----------------|
//! | [`sdt_metrics`] | `SdtMetrics` struct definition |
//! | [`network_solve`] | Full 1D solved-network extraction helpers |
//! | [`separation`] | 3-population + leukapheresis separation models |
//! | [`compute`] | `compute_metrics()` entry point + `giersiepen_hi()` |

mod compute;
mod network_solve;
mod sdt_metrics;
mod separation;

pub use compute::{compute_metrics, giersiepen_hi};
pub use sdt_metrics::SdtMetrics;
