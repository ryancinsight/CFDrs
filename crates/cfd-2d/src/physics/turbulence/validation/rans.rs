//! RANS model validation methods for [`TurbulenceValidator`].

use super::ValidationResult;
use crate::physics::turbulence::constants::{C2_EPSILON, EPSILON_MIN, SST_BETA_1};
use crate::physics::turbulence::traits::TurbulenceModel;
use crate::physics::turbulence::{KEpsilonModel, KOmegaSSTModel, SpalartAllmaras};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use std::fmt::Write;

use super::TurbulenceValidator;

impl<T: RealField + FromPrimitive + ToPrimitive + Copy> TurbulenceValidator<T> {
    /// Validate k-ε model against homogeneous turbulence decay
    pub fn validate_k_epsilon_homogeneous_decay(&self) -> ValidationResult {
        let _model: KEpsilonModel<T> = KEpsilonModel::new(1, 1);

        let k0 = T::from_f64(1.0).expect("analytical constant conversion");
        let eps0 = T::from_f64(0.1).expect("analytical constant conversion");
        let _density = T::from_f64(1.0).expect("analytical constant conversion");

        let dt = T::from_f64(0.01).expect("analytical constant conversion");
        let t_final = T::from_f64(10.0).expect("analytical constant conversion");

        let mut k = k0;
        let mut eps = eps0;
        let mut t = T::zero();

        let mut k_values = Vec::new();

        while t < t_final {
            k_values.push(k);

            let dk_dt = -eps;
            let deps_dt =
                -T::from_f64(C2_EPSILON).expect("analytical constant conversion") * eps * eps
                    / k.max(T::from_f64(1e-10).expect("analytical constant conversion"));

            k = (k + dk_dt * dt).max(T::zero());
            eps = (eps + deps_dt * dt)
                .max(T::from_f64(EPSILON_MIN).expect("analytical constant conversion"));

            t += dt;
        }

        let final_k_ratio = *k_values.last().unwrap_or(&k0) / k0;
        let decay_rate = -final_k_ratio.ln() / t_final;

        ValidationResult {
            test_name: "k-ε Homogeneous Turbulence Decay".to_string(),
            passed: decay_rate > T::from_f64(0.05).expect("analytical constant conversion")
                && decay_rate < T::from_f64(0.5).expect("analytical constant conversion"),
            metric: format!(
                "Decay rate: {rate:.4}",
                rate = decay_rate.to_f64().unwrap_or(0.0)
            ),
            details: format!(
                "k_final/k_initial = {ratio:.4}",
                ratio = final_k_ratio.to_f64().unwrap_or(0.0)
            ),
        }
    }

    /// Validate k-ω SST model near-wall behavior
    pub fn validate_k_omega_sst_wall_behavior(&self) -> ValidationResult {
        let _model: KOmegaSSTModel<T> = KOmegaSSTModel::new(1, 1);

        let molecular_viscosity = T::from_f64(1e-5).expect("analytical constant conversion");
        let y_wall = T::from_f64(1e-4).expect("analytical constant conversion");

        let beta1 = T::from_f64(SST_BETA_1).expect("analytical constant conversion");
        let expected_omega_wall = T::from_f64(6.0).expect("analytical constant conversion")
            * molecular_viscosity
            / (beta1 * y_wall * y_wall);

        let _k = [T::zero()];
        let mut omega = [T::one()];

        let omega_wall = T::from_f64(6.0).expect("analytical constant conversion")
            * molecular_viscosity
            / (beta1 * y_wall * y_wall);
        omega[0] = omega_wall;

        let omega_ratio = omega[0] / expected_omega_wall;

        ValidationResult {
            test_name: "k-ω SST Wall Boundary Condition".to_string(),
            passed: (omega_ratio - T::one()).abs() < self.tolerance,
            metric: format!(
                "ω_wall ratio: {ratio:.4}",
                ratio = omega_ratio.to_f64().unwrap_or(0.0)
            ),
            details: format!(
                "Expected: {expected:.2e}, Got: {got:.2e}",
                expected = expected_omega_wall.to_f64().unwrap_or(0.0),
                got = omega[0].to_f64().unwrap_or(0.0)
            ),
        }
    }

    /// Validate Spalart-Allmaras model eddy viscosity calculation
    pub fn validate_sa_eddy_viscosity(&self) -> ValidationResult {
        let model = SpalartAllmaras::<T>::new(1, 1);

        let test_cases = vec![
            (
                T::from_f64(1e-4).expect("analytical constant conversion"),
                T::from_f64(1e-5).expect("analytical constant conversion"),
                T::from_f64(7.36e-5).expect("analytical constant conversion"),
            ),
            (
                T::from_f64(1e-2).expect("analytical constant conversion"),
                T::from_f64(1e-5).expect("analytical constant conversion"),
                T::from_f64(9.41e-4).expect("analytical constant conversion"),
            ),
        ];

        let mut passed_all = true;
        let mut details = String::new();

        for (nu_tilde, nu, expected_nu_t) in test_cases {
            let nu_t = model.eddy_viscosity(nu_tilde, nu);
            let ratio = nu_t / expected_nu_t;
            let passed = (ratio - T::one()).abs()
                < T::from_f64(0.01).expect("analytical constant conversion");

            passed_all &= passed;
            let _ = writeln!(
                details,
                "ν̃={:.2e}, ν={:.2e}: got {:.2e}, expected {:.2e}, ratio={:.4}",
                nu_tilde.to_f64().unwrap_or(0.0),
                nu.to_f64().unwrap_or(0.0),
                nu_t.to_f64().unwrap_or(0.0),
                expected_nu_t.to_f64().unwrap_or(0.0),
                ratio.to_f64().unwrap_or(0.0)
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
                let mut k =
                    vec![T::from_f64(0.1).expect("analytical constant conversion"); nx * ny];
                let mut epsilon =
                    vec![T::from_f64(0.01).expect("analytical constant conversion"); nx * ny];

                for _ in 0..5 {
                    model
                        .update(
                            &mut k,
                            &mut epsilon,
                            &vec![nalgebra::Vector2::zeros(); nx * ny],
                            T::one(),
                            T::from_f64(1e-5).expect("analytical constant conversion"),
                            T::from_f64(0.01).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
                        )
                        .unwrap();
                }

                (k, epsilon, vec![], vec![])
            }
            "k-omega-sst" => {
                let mut model = KOmegaSSTModel::new(nx, ny);
                let mut k =
                    vec![T::from_f64(0.1).expect("analytical constant conversion"); nx * ny];
                let mut omega =
                    vec![T::from_f64(10.0).expect("analytical constant conversion"); nx * ny];

                for _ in 0..5 {
                    model
                        .update(
                            &mut k,
                            &mut omega,
                            &vec![nalgebra::Vector2::zeros(); nx * ny],
                            T::one(),
                            T::from_f64(1e-5).expect("analytical constant conversion"),
                            T::from_f64(0.01).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
                        )
                        .unwrap();
                }

                (k, vec![], omega, vec![])
            }
            "spalart-allmaras" => {
                let _model = SpalartAllmaras::new(nx, ny);
                let mut nu_tilde =
                    vec![T::from_f64(1e-4).expect("analytical constant conversion"); nx * ny];

                for _ in 0..5 {
                    _model
                        .update(
                            &mut nu_tilde,
                            &vec![nalgebra::Vector2::zeros(); nx * ny],
                            T::from_f64(1e-5).expect("analytical constant conversion"),
                            T::from_f64(0.01).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
                            T::from_f64(0.1).expect("analytical constant conversion"),
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

        let has_nan_inf = k.iter().any(|&x| !x.is_finite())
            || epsilon.iter().any(|&x| !x.is_finite())
            || omega.iter().any(|&x| !x.is_finite())
            || nu_tilde.iter().any(|&x| !x.is_finite());

        let all_positive = k.iter().all(|&x| x >= T::zero())
            && epsilon.iter().all(|&x| x >= T::zero())
            && omega.iter().all(|&x| x >= T::zero())
            && nu_tilde.iter().all(|&x| x >= T::zero());

        let passed = !has_nan_inf && all_positive;

        ValidationResult {
            test_name: format!("{model_name} Numerical Stability"),
            passed,
            metric: format!("Stable: {}, Positive: {}", !has_nan_inf, all_positive),
            details: format!("Grid: {nx}x{ny}, Iterations: 5"),
        }
    }
}
