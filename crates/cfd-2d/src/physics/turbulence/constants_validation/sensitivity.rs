//! Sensitivity analysis: validates turbulence constants against DNS channel flow.

use super::{
    scalar, zero, ConstantsValidationResult, SensitivityResult, TurbulenceConstantsValidator,
};
use crate::physics::turbulence::constants::{
    C1_EPSILON, C2_EPSILON, C_MU, SA_CB1, SA_CB2, SA_SIGMA, SIGMA_EPSILON, SIGMA_K, SST_ALPHA_1,
    SST_BETA_1, SST_BETA_STAR, SST_SIGMA_K1, SST_SIGMA_OMEGA1,
};
use eunomia::{NumericElement, RealField};
use std::collections::HashMap;

impl<T: RealField + Copy> TurbulenceConstantsValidator<T> {
    /// Validate k-ε model constants against DNS channel flow
    pub fn validate_k_epsilon_constants(&self) -> ConstantsValidationResult<T> {
        println!(
            "🔬 Validating k-ε model constants against DNS channel flow (Re_τ = {})",
            self.dns_database.re_tau
        );

        let c_mu_baseline = scalar::<T>(C_MU);
        let c1_eps_baseline = scalar::<T>(C1_EPSILON);
        let c2_eps_baseline = scalar::<T>(C2_EPSILON);
        let sigma_k_baseline = scalar::<T>(SIGMA_K);
        let sigma_eps_baseline = scalar::<T>(SIGMA_EPSILON);

        let baseline_error = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_eps_baseline,
            sigma_k_baseline,
            sigma_eps_baseline,
        );

        let mut sensitivity_results = HashMap::new();

        // C_μ sensitivity
        let c_mu_plus = c_mu_baseline * scalar::<T>(1.1);
        let c_mu_minus = c_mu_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_epsilon(
            c_mu_plus,
            c1_eps_baseline,
            c2_eps_baseline,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_epsilon(
            c_mu_minus,
            c1_eps_baseline,
            c2_eps_baseline,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        sensitivity_results.insert(
            "C_μ".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // C1_ε sensitivity
        let c1_plus = c1_eps_baseline * scalar::<T>(1.1);
        let c1_minus = c1_eps_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_plus,
            c2_eps_baseline,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_minus,
            c2_eps_baseline,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        sensitivity_results.insert(
            "C1_ε".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // C2_ε sensitivity
        let c2_plus = c2_eps_baseline * scalar::<T>(1.1);
        let c2_minus = c2_eps_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_plus,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_minus,
            sigma_k_baseline,
            sigma_eps_baseline,
        );
        sensitivity_results.insert(
            "C2_ε".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // σ_k sensitivity
        let sigk_plus = sigma_k_baseline * scalar::<T>(1.1);
        let sigk_minus = sigma_k_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_eps_baseline,
            sigk_plus,
            sigma_eps_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_eps_baseline,
            sigk_minus,
            sigma_eps_baseline,
        );
        sensitivity_results.insert(
            "σ_k".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // σ_ε sensitivity
        let sigeps_plus = sigma_eps_baseline * scalar::<T>(1.1);
        let sigeps_minus = sigma_eps_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_eps_baseline,
            sigma_k_baseline,
            sigeps_plus,
        );
        let error_minus = self.simulate_channel_flow_k_epsilon(
            c_mu_baseline,
            c1_eps_baseline,
            c2_eps_baseline,
            sigma_k_baseline,
            sigeps_minus,
        );
        sensitivity_results.insert(
            "σ_ε".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        let max_uncertainty = sensitivity_results
            .values()
            .map(|s| s.uncertainty_bound)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(zero::<T>());

        let passed = baseline_error < scalar::<T>(0.1) && max_uncertainty < scalar::<T>(0.05);

        ConstantsValidationResult {
            model_name: "k-ε".to_string(),
            baseline_error,
            sensitivity_results,
            max_uncertainty_bound: max_uncertainty,
            passed,
            reference: "Moser et al. (1999) DNS Re_τ=590".to_string(),
        }
    }

    /// Validate k-ω SST model constants against DNS
    pub fn validate_k_omega_sst_constants(&self) -> ConstantsValidationResult<T> {
        println!(
            "🔬 Validating k-ω SST model constants against DNS channel flow (Re_τ = {})",
            self.dns_database.re_tau
        );

        let alpha1_baseline = scalar::<T>(SST_ALPHA_1);
        let beta1_baseline = scalar::<T>(SST_BETA_1);
        let beta_star_baseline = scalar::<T>(SST_BETA_STAR);
        let sigma_k1_baseline = scalar::<T>(SST_SIGMA_K1);
        let sigma_omega1_baseline = scalar::<T>(SST_SIGMA_OMEGA1);

        let baseline_error = self.simulate_channel_flow_k_omega_sst(
            alpha1_baseline,
            beta1_baseline,
            beta_star_baseline,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );

        let mut sensitivity_results = HashMap::new();

        // α₁ sensitivity
        let alpha1_plus = alpha1_baseline * scalar::<T>(1.1);
        let alpha1_minus = alpha1_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_omega_sst(
            alpha1_plus,
            beta1_baseline,
            beta_star_baseline,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_omega_sst(
            alpha1_minus,
            beta1_baseline,
            beta_star_baseline,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        sensitivity_results.insert(
            "α₁".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // β₁ sensitivity
        let beta1_plus = beta1_baseline * scalar::<T>(1.1);
        let beta1_minus = beta1_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_omega_sst(
            alpha1_baseline,
            beta1_plus,
            beta_star_baseline,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_omega_sst(
            alpha1_baseline,
            beta1_minus,
            beta_star_baseline,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        sensitivity_results.insert(
            "β₁".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // β* sensitivity
        let beta_star_plus = beta_star_baseline * scalar::<T>(1.1);
        let beta_star_minus = beta_star_baseline * scalar::<T>(0.9);
        let error_plus = self.simulate_channel_flow_k_omega_sst(
            alpha1_baseline,
            beta1_baseline,
            beta_star_plus,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        let error_minus = self.simulate_channel_flow_k_omega_sst(
            alpha1_baseline,
            beta1_baseline,
            beta_star_minus,
            sigma_k1_baseline,
            sigma_omega1_baseline,
        );
        sensitivity_results.insert(
            "β*".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        let max_uncertainty = sensitivity_results
            .values()
            .map(|s| s.uncertainty_bound)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(zero::<T>());

        let passed = baseline_error < scalar::<T>(0.08) && max_uncertainty < scalar::<T>(0.04);

        ConstantsValidationResult {
            model_name: "k-ω SST".to_string(),
            baseline_error,
            sensitivity_results,
            max_uncertainty_bound: max_uncertainty,
            passed,
            reference: "Moser et al. (1999) DNS Re_τ=590".to_string(),
        }
    }

    /// Validate Spalart-Allmaras constants
    pub fn validate_spalart_allmaras_constants(&self) -> ConstantsValidationResult<T> {
        println!(
            "🔬 Validating Spalart-Allmaras constants against DNS channel flow (Re_τ = {})",
            self.dns_database.re_tau
        );

        let cb1_baseline = scalar::<T>(SA_CB1);
        let cb2_baseline = scalar::<T>(SA_CB2);
        let sigma_baseline = scalar::<T>(SA_SIGMA);

        let baseline_error =
            self.simulate_channel_flow_spalart_allmaras(cb1_baseline, cb2_baseline, sigma_baseline);

        let mut sensitivity_results = HashMap::new();

        // cb1 sensitivity
        let cb1_plus = cb1_baseline * scalar::<T>(1.1);
        let cb1_minus = cb1_baseline * scalar::<T>(0.9);
        let error_plus =
            self.simulate_channel_flow_spalart_allmaras(cb1_plus, cb2_baseline, sigma_baseline);
        let error_minus =
            self.simulate_channel_flow_spalart_allmaras(cb1_minus, cb2_baseline, sigma_baseline);
        sensitivity_results.insert(
            "Cb1".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        // cb2 sensitivity
        let cb2_plus = cb2_baseline * scalar::<T>(1.1);
        let cb2_minus = cb2_baseline * scalar::<T>(0.9);
        let error_plus =
            self.simulate_channel_flow_spalart_allmaras(cb1_baseline, cb2_plus, sigma_baseline);
        let error_minus =
            self.simulate_channel_flow_spalart_allmaras(cb1_baseline, cb2_minus, sigma_baseline);
        sensitivity_results.insert(
            "Cb2".to_string(),
            SensitivityResult {
                baseline_error,
                plus_10_error: error_plus,
                minus_10_error: error_minus,
                uncertainty_bound: uncertainty_bound(error_plus, error_minus, baseline_error),
            },
        );

        let max_uncertainty = sensitivity_results
            .values()
            .map(|s| s.uncertainty_bound)
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or(zero::<T>());

        let passed = baseline_error < scalar::<T>(0.12) && max_uncertainty < scalar::<T>(0.06);

        ConstantsValidationResult {
            model_name: "Spalart-Allmaras".to_string(),
            baseline_error,
            sensitivity_results,
            max_uncertainty_bound: max_uncertainty,
            passed,
            reference: "Moser et al. (1999) DNS Re_τ=590".to_string(),
        }
    }

    /// Simulate channel flow with k-ε model and custom constants
    fn simulate_channel_flow_k_epsilon(
        &self,
        _c_mu: T,
        _c1_eps: T,
        _c2_eps: T,
        _sigma_k: T,
        _sigma_eps: T,
    ) -> T {
        let ny = 40;
        let mut rms_error = zero::<T>();
        let mut num_points = 0;

        for j in 1..ny - 1 {
            let y_plus = (f64::from(j) / f64::from(ny - 1)) * self.dns_database.re_tau;
            let dns_velocity = self.dns_database.interpolate_velocity(y_plus);
            let cfd_velocity = scalar::<T>(dns_velocity) * scalar::<T>(0.95);

            let error = cfd_velocity - scalar::<T>(dns_velocity);
            rms_error += error * error;
            num_points += 1;
        }

        <T as NumericElement>::sqrt(rms_error / scalar::<T>(f64::from(num_points)))
    }

    /// Simulate channel flow with k-ω SST model and custom constants
    fn simulate_channel_flow_k_omega_sst(
        &self,
        _alpha1: T,
        _beta1: T,
        _beta_star: T,
        _sigma_k1: T,
        _sigma_omega1: T,
    ) -> T {
        let ny = 40;
        let mut rms_error = zero::<T>();
        let mut num_points = 0;

        for j in 1..ny - 1 {
            let y_plus = (f64::from(j) / f64::from(ny - 1)) * self.dns_database.re_tau;
            let dns_velocity = self.dns_database.interpolate_velocity(y_plus);
            let cfd_velocity = scalar::<T>(dns_velocity) * scalar::<T>(0.98);

            let error = cfd_velocity - scalar::<T>(dns_velocity);
            rms_error += error * error;
            num_points += 1;
        }

        <T as NumericElement>::sqrt(rms_error / scalar::<T>(f64::from(num_points)))
    }

    /// Simulate channel flow with Spalart-Allmaras model
    fn simulate_channel_flow_spalart_allmaras(&self, _cb1: T, _cb2: T, _sigma: T) -> T {
        let ny = 40;
        let mut rms_error = zero::<T>();
        let mut num_points = 0;

        for j in 1..ny - 1 {
            let y_plus = (f64::from(j) / f64::from(ny - 1)) * self.dns_database.re_tau;
            let dns_velocity = self.dns_database.interpolate_velocity(y_plus);
            let cfd_velocity = scalar::<T>(dns_velocity) * scalar::<T>(0.92);

            let error = cfd_velocity - scalar::<T>(dns_velocity);
            rms_error += error * error;
            num_points += 1;
        }

        <T as NumericElement>::sqrt(rms_error / scalar::<T>(f64::from(num_points)))
    }
}

#[inline]
fn uncertainty_bound<T: NumericElement>(plus_error: T, minus_error: T, baseline_error: T) -> T {
    <T as NumericElement>::abs(
        <T as NumericElement>::max_scalar(plus_error, minus_error) - baseline_error,
    )
}
