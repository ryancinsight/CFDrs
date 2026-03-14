mod design_record;
mod validation;

pub use design_record::{
    compute_blueprint_report_metrics, Milestone12ReportDesign, ParetoPoint, ParetoTag,
};
pub use validation::ValidationRow;

/// Format a boolean gate as a PASS/FAIL string for markdown tables.
pub(in crate::reporting) fn pass_fail(value: bool) -> &'static str {
    if value { "PASS" } else { "FAIL" }
}

/// Pediatric 3 kg neonatal ECV as a percentage of the 25.5 mL limit.
pub(in crate::reporting) fn pediatric_limit_pct(ecv_ml: f64) -> f64 {
    100.0 * ecv_ml / 25.5_f64
}
