//! Convenience helpers for chip-level clinical scaling and parallel assembly.
//!
//! [`scale_to_clinical`] computes how many parallel chip replicas are needed to
//! reach a target flow rate and verifies the combined ECV against the paediatric
//! safety budget.
//!
//! [`make_parallel_candidate`] constructs a logically-scaled [`DesignCandidate`]
//! that represents N identical chips operating in parallel — useful for passing
//! to `compute_metrics()` when the `ParallelMicrochannelArray` topology is used.

use cfd_optim::design::DesignCandidate;

// ── Clinical scaling result ───────────────────────────────────────────────────

/// Summary of a parallel-chip clinical feasibility assessment.
#[derive(Debug, Clone)]
pub struct ClinicalScaleResult {
    /// Number of parallel chip instances required to reach the target flow.
    pub n_parallel: usize,
    /// Achieved total flow rate across all parallel chips. Units: mL / min.
    pub total_flow_ml_min: f64,
    /// Estimated combined ECV of the parallel chip assembly. Units: mL.
    pub total_ecv_ml: f64,
    /// Paediatric ECV budget (10 % of patient TBV). Units: mL.
    pub ecv_budget_ml: f64,
    /// Whether the combined ECV satisfies the budget constraint.
    pub ecv_ok: bool,
}

impl ClinicalScaleResult {
    /// Print a human-readable feasibility summary.
    pub fn print(&self, target_ml_min: f64) {
        println!("  Parallel chips : {}", self.n_parallel);
        println!(
            "  Total flow     : {:.2} mL/min  (target ≥ {:.0} mL/min)",
            self.total_flow_ml_min, target_ml_min,
        );
        println!(
            "  ECV            : {:.2} mL  (budget {:.1} mL  →  {})",
            self.total_ecv_ml,
            self.ecv_budget_ml,
            if self.ecv_ok { "PASS" } else { "FAIL" },
        );
    }
}

// ── Scaling function ──────────────────────────────────────────────────────────

/// Compute the number of parallel chip replicas needed to reach `target_ml_min`
/// and assess ECV feasibility for the given patient.
///
/// # Arguments
/// * `single_chip_flow_m3s` — volumetric flow rate of one chip (m³/s)
/// * `single_chip_ecv_ml`   — channel-volume ECV of one chip (mL)
/// * `ecv_budget_ml`        — maximum ECV for the patient (mL)
/// * `target_ml_min`        — required total throughput (mL/min)
pub fn scale_to_clinical(
    single_chip_flow_m3s: f64,
    single_chip_ecv_ml: f64,
    ecv_budget_ml: f64,
    target_ml_min: f64,
) -> ClinicalScaleResult {
    let single_flow_ml_min = single_chip_flow_m3s * 6.0e7; // m³/s → mL/min
    let n_parallel = if single_flow_ml_min > 0.0 {
        ((target_ml_min / single_flow_ml_min).ceil() as usize).max(1)
    } else {
        1
    };
    let total_flow = single_flow_ml_min * n_parallel as f64;
    let total_ecv = single_chip_ecv_ml * n_parallel as f64;
    ClinicalScaleResult {
        n_parallel,
        total_flow_ml_min: total_flow,
        total_ecv_ml: total_ecv,
        ecv_budget_ml,
        ecv_ok: total_ecv <= ecv_budget_ml,
    }
}

// ── Parallel candidate builder ────────────────────────────────────────────────

/// Build a logically-scaled version of `base` that represents `n_parallel`
/// identical chips operating in parallel.
///
/// The flow rate is multiplied by `n_parallel` so that `compute_metrics()`
/// sees the aggregate throughput while the channel geometry remains unchanged.
/// This is meaningful when the topology is `ParallelMicrochannelArray`; for
/// other topologies it serves as a throughput-annotation only.
pub fn make_parallel_candidate(base: &DesignCandidate, n_parallel: usize) -> DesignCandidate {
    let mut scaled = base.clone();
    scaled.id = format!("{}_x{n_parallel}", base.id);
    scaled.flow_rate_m3_s = base.flow_rate_m3_s * n_parallel as f64;
    scaled
}
