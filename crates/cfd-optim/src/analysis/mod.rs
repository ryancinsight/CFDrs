//! Post-processing analysis for SDT design optimization results.
//!
//! Provides multi-objective Pareto front computation ([`pareto`]) and
//! parametric robustness sweeps ([`robustness`]).
//!
//! ## Quick start
//!
//! ```ignore
//! use cfd_optim::analysis::{compute_sdt_pareto_front, robustness_sweep, STANDARD_PERTURBATIONS};
//! use cfd_optim::scoring::{OptimMode, SdtWeights};
//!
//! // Build Pareto front from all feasible RankedDesigns
//! let front = compute_sdt_pareto_front(&feasible_designs);
//! println!("Pareto front size: {}", front.len());
//!
//! // Evaluate robustness of top candidate
//! let report = robustness_sweep(&top.candidate, OptimMode::RbcProtectedSdt,
//!                               &SdtWeights::default(), &STANDARD_PERTURBATIONS);
//! println!("Is robust: {}", report.is_robust);
//! ```

pub mod pareto;
pub mod robustness;

pub use pareto::{compute_sdt_pareto_front, SdtParetoFront};
pub use robustness::{robustness_sweep, RobustnessReport, STANDARD_PERTURBATIONS};
