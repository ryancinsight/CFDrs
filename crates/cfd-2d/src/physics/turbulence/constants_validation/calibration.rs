//! Turbulence-constant calibration against DNS via Coeus least-squares.
//!
//! The sensitivity analysis in this module answers "how wrong is each
//! constant?" by perturbing it ±10%. This submodule answers the complementary
//! question — "which constants best reproduce the DNS reference?" — by
//! nonlinear least-squares. The fit is solved by Coeus's
//! Levenberg-Marquardt solver ([`coeus_optim::least_squares`]), keeping Coeus
//! the single source of truth for first-order model calibration in CFDrs.
//!
//! # Identifiability of the k-ε constants
//!
//! In the equilibrium log-law region the standard k-ε model reduces to
//!
//! ```text
//! κ² = √C_μ · (C_ε2 − C_ε1) · σ_ε
//! ```
//!
//! so the mean-velocity profile `u⁺ = (1/κ) ln y⁺ + B` identifies only the
//! *combination* κ — not the four constants separately (Pope 2000, §10.4).
//! This module therefore calibrates the identifiable quantity (κ, together
//! with the log-law intercept `B`) and reports how the fitted κ compares to
//! the value the standard constants imply, rather than over-claiming a
//! per-constant fit that the equilibrium data cannot support.

use coeus_optim::{
    levenberg_marquardt, LeastSquaresProblem, LeastSquaresReport, LevenbergMarquardtConfig,
    ProblemError, Termination,
};

use super::dns_database::DnsChannelFlowDatabase;
use crate::physics::turbulence::constants::{C1_EPSILON, C2_EPSILON, C_MU, SIGMA_EPSILON};

/// Lower bound of the equilibrium log-law region in wall units.
const LOG_LAW_Y_MIN: f64 = 30.0;

/// Fraction of `Re_τ` that bounds the log-law region from above.
const LOG_LAW_Y_MAX_FRACTION: f64 = 0.3;

/// The equilibrium log-law problem: fit `u⁺ = (1/κ) ln y⁺ + B` to DNS.
///
/// Free parameters are `[κ, B]` (the von Karman constant and the log-law
/// intercept). The DNS samples are restricted to `30 ≤ y⁺ ≤ 0.3·Re_τ`, where
/// the equilibrium law-of-the-wall holds.
pub struct LogLawProblem {
    /// DNS mean-velocity samples in the log-law region, as `(y⁺, u⁺)` pairs.
    points: Vec<(f64, f64)>,
}

impl LogLawProblem {
    /// Build the problem from the DNS database, selecting the log-law region.
    #[must_use]
    pub fn from_database(db: &DnsChannelFlowDatabase) -> Self {
        let y_max = LOG_LAW_Y_MAX_FRACTION * db.re_tau;
        let points = db
            .mean_velocity_profile
            .iter()
            .copied()
            .filter(|(y, _)| *y >= LOG_LAW_Y_MIN && *y <= y_max)
            .collect();
        Self { points }
    }

    /// The number of free parameters (κ and B).
    const PARAMETER_COUNT: usize = 2;

    /// The log-law-region samples this problem fits against.
    #[must_use]
    pub fn points(&self) -> &[(f64, f64)] {
        &self.points
    }
}

impl LeastSquaresProblem<f64> for LogLawProblem {
    fn residual_count(&self) -> usize {
        self.points.len()
    }

    fn parameter_count(&self) -> usize {
        Self::PARAMETER_COUNT
    }

    fn residuals(&self, parameters: &[f64], residuals: &mut [f64]) -> Result<(), ProblemError> {
        let kappa = parameters[0];
        let intercept = parameters[1];
        if kappa <= 0.0 {
            return Err(ProblemError::Domain {
                reason: format!("non-positive von Karman constant κ = {kappa}"),
            });
        }
        for (slot, (y, u_dns)) in residuals.iter_mut().zip(&self.points) {
            *slot = y.ln() / kappa + intercept - u_dns;
        }
        Ok(())
    }

    fn jacobian(&self, parameters: &[f64], jacobian: &mut [f64]) -> Result<(), ProblemError> {
        let kappa = parameters[0];
        if kappa <= 0.0 {
            return Err(ProblemError::Domain {
                reason: format!("non-positive von Karman constant κ = {kappa}"),
            });
        }
        let kappa_sq = kappa * kappa;
        for (row, (y, _)) in self.points.iter().enumerate() {
            // ∂u/∂κ = −ln y / κ², ∂u/∂B = 1.
            jacobian[row * Self::PARAMETER_COUNT] = -y.ln() / kappa_sq;
            jacobian[row * Self::PARAMETER_COUNT + 1] = 1.0;
        }
        Ok(())
    }
}

/// Outcome of a log-law calibration.
#[derive(Debug, Clone)]
pub struct LogLawCalibration {
    /// Calibrated von Karman constant κ.
    pub kappa: f64,
    /// Calibrated log-law intercept B.
    pub intercept: f64,
    /// `0.5‖r‖²` at the calibrated point.
    pub cost: f64,
    /// Iterations executed by the solver.
    pub iterations: usize,
    /// Why the solver stopped.
    pub termination: Termination,
}

impl LogLawCalibration {
    /// Root-mean-square residual of the calibrated fit against the DNS.
    #[must_use]
    pub fn rms_residual(&self, residual_count: usize) -> f64 {
        if residual_count == 0 {
            return 0.0;
        }
        (2.0 * self.cost / residual_count as f64).sqrt()
    }
}

/// Calibrate the log-law `(κ, B)` against the DNS mean velocity.
///
/// The solver is started from the textbook values `κ = 0.41`, `B = 5.0` and
/// refines them to the DNS samples in the log-law region.
///
/// # Panics
///
/// Panics if the solver returns an error. The problem is well-posed
/// (two parameters, at least six residuals), so the only error paths are
/// programmer error in this module rather than data-dependent failure.
#[must_use]
pub fn calibrate_log_law(db: &DnsChannelFlowDatabase) -> LogLawCalibration {
    let problem = LogLawProblem::from_database(db);
    let initial = [0.41_f64, 5.0_f64];
    let config = LevenbergMarquardtConfig::default();
    let report: LeastSquaresReport<f64> =
        levenberg_marquardt(&problem, &initial, &config).expect("log-law fit is well-posed");
    LogLawCalibration {
        kappa: report.parameters[0],
        intercept: report.parameters[1],
        cost: report.cost,
        iterations: report.iterations,
        termination: report.termination,
    }
}

/// The von Karman constant implied by the standard k-ε constants.
///
/// This is the equilibrium relation `κ² = √C_μ · (C_ε2 − C_ε1) · σ_ε`.
#[must_use]
pub fn kappa_from_k_epsilon_constants(
    c_mu: f64,
    c1_epsilon: f64,
    c2_epsilon: f64,
    sigma_epsilon: f64,
) -> f64 {
    (c_mu.sqrt() * (c2_epsilon - c1_epsilon) * sigma_epsilon).sqrt()
}

/// The κ implied by CFDrs's standard k-ε constants.
#[must_use]
pub fn standard_kappa() -> f64 {
    kappa_from_k_epsilon_constants(C_MU, C1_EPSILON, C2_EPSILON, SIGMA_EPSILON)
}

/// Calibrate the k-ε log-law against DNS and compare to the standard constants.
#[derive(Debug, Clone)]
pub struct KepsilonCalibrationReport {
    /// The raw log-law fit.
    pub log_law: LogLawCalibration,
    /// κ implied by the standard k-ε constants.
    pub standard_kappa: f64,
    /// Whether the calibrated κ stays within 10% of the standard value.
    pub within_standard_band: bool,
}

impl KepsilonCalibrationReport {
    /// Run the calibration and produce the comparison report.
    #[must_use]
    pub fn run(db: &DnsChannelFlowDatabase) -> Self {
        let log_law = calibrate_log_law(db);
        let standard_kappa = standard_kappa();
        let within_standard_band = (log_law.kappa - standard_kappa).abs() <= 0.1 * standard_kappa;
        Self {
            log_law,
            standard_kappa,
            within_standard_band,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn database() -> DnsChannelFlowDatabase {
        DnsChannelFlowDatabase::moser_1999_re590()
    }

    #[test]
    fn standard_constants_imply_kappa_near_0_41() {
        let kappa = standard_kappa();
        assert!((kappa - 0.41).abs() < 0.05, "κ = {kappa}");
    }

    #[test]
    fn log_law_region_contains_only_wall_units_between_30_and_177() {
        let db = database();
        let problem = LogLawProblem::from_database(&db);
        assert!(problem.residual_count() >= problem.parameter_count());
        for (y, _) in problem.points() {
            assert!(*y >= LOG_LAW_Y_MIN && *y <= LOG_LAW_Y_MAX_FRACTION * db.re_tau);
        }
    }

    #[test]
    fn jacobian_matches_finite_differences() {
        let db = database();
        let problem = LogLawProblem::from_database(&db);
        let parameters = [0.41_f64, 5.0_f64];
        let n = problem.residual_count();
        let p = problem.parameter_count();

        let mut analytic = vec![0.0; n * p];
        problem
            .jacobian(&parameters, &mut analytic)
            .expect("jacobian evaluation");

        let h = 1e-6;
        for column in 0..p {
            let mut plus = parameters;
            let mut minus = parameters;
            plus[column] += h;
            minus[column] -= h;
            let mut r_plus = vec![0.0; n];
            let mut r_minus = vec![0.0; n];
            problem
                .residuals(&plus, &mut r_plus)
                .expect("residual evaluation");
            problem
                .residuals(&minus, &mut r_minus)
                .expect("residual evaluation");
            for row in 0..n {
                let finite = (r_plus[row] - r_minus[row]) / (2.0 * h);
                let analytic_value = analytic[row * p + column];
                let tol = 1e-5 * (1.0 + analytic_value.abs());
                assert!(
                    (finite - analytic_value).abs() < tol,
                    "row {row} col {column}: finite {finite} vs analytic {analytic_value}"
                );
            }
        }
    }

    #[test]
    fn calibration_converges_and_reduces_error() {
        let db = database();
        let report = KepsilonCalibrationReport::run(&db);
        assert!(report.log_law.termination.is_converged());
        assert!(report.log_law.kappa > 0.2 && report.log_law.kappa < 0.6);
        assert!(report.log_law.cost.is_finite());
        // The fit must do better than the standard-constants profile at the
        // same intercept, otherwise calibration bought nothing.
        let problem = LogLawProblem::from_database(&db);
        let standard = [standard_kappa(), 5.0_f64];
        let mut residuals = vec![0.0; problem.residual_count()];
        problem
            .residuals(&standard, &mut residuals)
            .expect("residual evaluation");
        let standard_cost = 0.5 * residuals.iter().map(|r| r * r).sum::<f64>();
        assert!(
            report.log_law.cost < standard_cost,
            "calibrated {:.6} should beat standard {:.6}",
            report.log_law.cost,
            standard_cost
        );
    }
}
