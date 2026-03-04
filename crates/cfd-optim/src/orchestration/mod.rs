//! High-level orchestration: parametric sweep + ranking engine.
//!
//! [`optimizer`] owns `SdtOptimizer`, the entry point for running a full
//! parametric sweep, ranking candidates, and optionally performing robustness
//! sweeps.  It is separate from `evo/optimizer.rs` (`GeneticOptimizer`) to
//! enforce SRP: parametric ranking ≠ evolutionary search.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`optimizer`] | `SdtOptimizer`, `RankedDesign`, `RobustSweepConfig`, `OptimStats` |

pub mod optimizer;

pub use optimizer::{OptimStats, RankedDesign, RobustScoreStats, RobustSweepConfig, SdtOptimizer};
