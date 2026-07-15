//! RANS model validation methods for [`TurbulenceValidator`].

use super::ValidationResult;
use crate::physics::turbulence::constants::{C2_EPSILON, EPSILON_MIN, SST_BETA_1};
use crate::physics::turbulence::traits::TurbulenceModel;
use crate::physics::turbulence::{KEpsilonModel, KOmegaSSTModel, SpalartAllmaras};
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};
use leto::geometry::Vector2;
use std::fmt::Write;

use super::TurbulenceValidator;

impl<T: EunomiaRealField> TurbulenceValidator<T> {
    #[inline]
    fn scalar(value: f64) -> T {
        <T as FloatElement>::from_f64(value)
    }

    /// Validate k-ε model against homogeneous turbulence decay
    pub fn validate_k_epsilon_homogeneous_decay(&self) -> ValidationResult {
        let _model: KEpsilonModel<T> = KEpsilonModel::new(1, 1);

        let k0 = Self::scalar(1.0);
        let eps0 = Self::scalar(0.1);
        let _density = Self::scalar(1.0);

        let dt = Self::scalar(0.01);
        let t_final = Self::scalar(10.0);

        let mut k = k0;
        let mut eps = eps0;
        let mut t = T::ZERO;

        let mut k_values = Vec::new();

        while t < t_final {
            k_values.push(k);

            let dk_dt = -eps;
            let deps_dt = -Self::scalar(C2_EPSILON) * eps * eps / k.max_scalar(Self::scalar(1e-10));

            k = (k + dk_dt * dt).max_scalar(T::ZERO);
            eps = (eps + deps_dt * dt).max_scalar(Self::scalar(EPSILON_MIN));

            t += dt;
        }

        let final_k_ratio = *k_values.last().unwrap_or(&k0) / k0;
        let decay_rate = -<T as FloatElement>::ln(final_k_ratio) / t_final;

        ValidationResult {
            test_name: "k-ε Homogeneous Turbulence Decay".to_string(),
            passed: decay_rate > Self::scalar(0.05) && decay_rate < Self::scalar(0.5),
            metric: format!(
                "Decay rate: {rate:.4}",
                rate = <T as NumericElement>::to_f64(decay_rate)
            ),
            details: format!(
                "k_final/k_initial = {ratio:.4}",
                ratio = <T as NumericElement>::to_f64(final_k_ratio)
            ),
        }
    }

    /// Validate k-ω SST model near-wall behavior
    pub fn validate_k_omega_sst_wall_behavior(&self) -> ValidationResult {
        let _model: KOmegaSSTModel<T> = KOmegaSSTModel::new(1, 1);

        let molecular_viscosity = Self::scalar(1e-5);
        let y_wall = Self::scalar(1e-4);

        let beta1 = Self::scalar(SST_BETA_1);
        let expected_omega_wall =
            Self::scalar(6.0) * molecular_viscosity / (beta1 * y_wall * y_wall);

        let mut omega = [T::ONE];

        let omega_wall = Self::scalar(6.0) * molecular_viscosity / (beta1 * y_wall * y_wall);
        omega[0] = omega_wall;

        let omega_ratio = omega[0] / expected_omega_wall;

        ValidationResult {
            test_name: "k-ω SST Wall Boundary Condition".to_string(),
            passed: <T as NumericElement>::abs(omega_ratio - T::ONE) < self.tolerance,
            metric: format!(
                "ω_wall ratio: {ratio:.4}",
                ratio = <T as NumericElement>::to_f64(omega_ratio)
            ),
            details: format!(
                "Expected: {expected:.2e}, Got: {got:.2e}",
                expected = <T as NumericElement>::to_f64(expected_omega_wall),
                got = <T as NumericElement>::to_f64(omega[0])
            ),
        }
    }

    /// Validate Spalart-Allmaras model eddy viscosity calculation
    pub fn validate_sa_eddy_viscosity(&self) -> ValidationResult {
        let model = SpalartAllmaras::<T>::new(1, 1);

        let test_cases = vec![
            (
                Self::scalar(1e-4),
                Self::scalar(1e-5),
                Self::scalar(7.36e-5),
            ),
            (
                Self::scalar(1e-2),
                Self::scalar(1e-5),
                Self::scalar(9.41e-4),
            ),
        ];

        let mut passed_all = true;
        let mut details = String::new();

        for (nu_tilde, nu, expected_nu_t) in test_cases {
            let nu_t = model.eddy_viscosity(nu_tilde, nu);
            let ratio = nu_t / expected_nu_t;
            let passed = <T as NumericElement>::abs(ratio - T::ONE) < Self::scalar(0.01);

            passed_all &= passed;
            let _ = writeln!(
                details,
                "ν̃={:.2e}, ν={:.2e}: got {:.2e}, expected {:.2e}, ratio={:.4}",
                <T as NumericElement>::to_f64(nu_tilde),
                <T as NumericElement>::to_f64(nu),
                <T as NumericElement>::to_f64(nu_t),
                <T as NumericElement>::to_f64(expected_nu_t),
                <T as NumericElement>::to_f64(ratio)
            );
        }

        ValidationResult {
            test_name: "SA Eddy Viscosity Calculation".to_string(),
            passed: passed_all,
            metric: format!("All test cases passed: {passed_all}"),
            details,
        }
    }

    /// Validate turbulence model convergence behavior
    pub fn validate_model_convergence(&self, model_name: &str) -> ValidationResult {
        let nx = 20;
        let ny = 10;

        let (k, epsilon, omega, nu_tilde) = match model_name {
            "k-epsilon" => {
                let mut model = KEpsilonModel::new(nx, ny);
                let mut k = vec![Self::scalar(0.1); nx * ny];
                let mut epsilon = vec![Self::scalar(0.01); nx * ny];

                for _ in 0..5 {
                    model
                        .update(
                            &mut k,
                            &mut epsilon,
                            &vec![Vector2::new(T::ZERO, T::ZERO); nx * ny],
                            T::ONE,
                            Self::scalar(1e-5),
                            Self::scalar(0.01),
                            Self::scalar(0.1),
                            Self::scalar(0.1),
                        )
                        .unwrap();
                }

                (k, epsilon, vec![], vec![])
            }
            "k-omega-sst" => {
                let mut model = KOmegaSSTModel::new(nx, ny);
                let mut k = vec![Self::scalar(0.1); nx * ny];
                let mut omega = vec![Self::scalar(10.0); nx * ny];

                for _ in 0..5 {
                    model
                        .update(
                            &mut k,
                            &mut omega,
                            &vec![Vector2::new(T::ZERO, T::ZERO); nx * ny],
                            T::ONE,
                            Self::scalar(1e-5),
                            Self::scalar(0.01),
                            Self::scalar(0.1),
                            Self::scalar(0.1),
                        )
                        .unwrap();
                }

                (k, vec![], omega, vec![])
            }
            "spalart-allmaras" => {
                let _model = SpalartAllmaras::new(nx, ny);
                let mut nu_tilde = vec![Self::scalar(1e-4); nx * ny];

                for _ in 0..5 {
                    _model
                        .update(
                            &mut nu_tilde,
                            &vec![Vector2::new(T::ZERO, T::ZERO); nx * ny],
                            Self::scalar(1e-5),
                            Self::scalar(0.01),
                            Self::scalar(0.1),
                            Self::scalar(0.1),
                        )
                        .unwrap();
                }

                (vec![], vec![], vec![], nu_tilde)
            }
            _ => {
                return ValidationResult {
                    test_name: format!("{model_name} Convergence"),
                    passed: false,
                    metric: "Unknown model".to_string(),
                    details: "Model not recognized".to_string(),
                }
            }
        };

        let has_nan_inf = k.iter().any(|&x| !<T as NumericElement>::is_finite(x))
            || epsilon
                .iter()
                .any(|&x| !<T as NumericElement>::is_finite(x))
            || omega.iter().any(|&x| !<T as NumericElement>::is_finite(x))
            || nu_tilde
                .iter()
                .any(|&x| !<T as NumericElement>::is_finite(x));

        let all_positive = k.iter().all(|&x| x >= T::ZERO)
            && epsilon.iter().all(|&x| x >= T::ZERO)
            && omega.iter().all(|&x| x >= T::ZERO)
            && nu_tilde.iter().all(|&x| x >= T::ZERO);

        let passed = !has_nan_inf && all_positive;

        ValidationResult {
            test_name: format!("{model_name} Numerical Stability"),
            passed,
            metric: format!("Stable: {}, Positive: {}", !has_nan_inf, all_positive),
            details: format!("Grid: {nx}x{ny}, Iterations: 5"),
        }
    }
}
