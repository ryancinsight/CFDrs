use super::*;
use super::spalding::spalding_u_plus;
use crate::physics::turbulence::constants::{KAPPA, E_WALL_FUNCTION, C_MU, Y_PLUS_VISCOUS_SUBLAYER};
use approx::assert_relative_eq;

#[test]
fn test_new_standard_wall_treatment() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    assert_relative_eq!(treatment.kappa, KAPPA, epsilon = 1e-10);
    assert_relative_eq!(treatment.e_wall, E_WALL_FUNCTION, epsilon = 1e-10);
}

#[test]
fn test_new_blended_wall_treatment() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
    assert_relative_eq!(treatment.kappa, KAPPA, epsilon = 1e-10);
}

#[test]
fn test_low_reynolds_u_plus() {
    let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
    assert_relative_eq!(treatment.u_plus(1.0), 1.0, epsilon = 1e-10);
    assert_relative_eq!(treatment.u_plus(5.0), 5.0, epsilon = 1e-10);
    assert_relative_eq!(treatment.u_plus(10.0), 10.0, epsilon = 1e-10);
}

#[test]
fn test_standard_wall_function_viscous_sublayer() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus = 3.0;
    let u_plus = treatment.u_plus(y_plus);
    assert_relative_eq!(u_plus, y_plus, epsilon = 1e-10);
}

#[test]
fn test_standard_wall_function_log_law() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus = 100.0;
    let u_plus = treatment.u_plus(y_plus);
    let expected = (y_plus.ln() / KAPPA) + 5.5;
    assert_relative_eq!(u_plus, expected, epsilon = 1e-6);
}

#[test]
fn test_standard_wall_function_buffer_layer() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus = 15.0;
    let u_plus = treatment.u_plus(y_plus);
    assert!(u_plus > 5.0);
    assert!(u_plus < (30.0_f64.ln() / KAPPA) + 5.5);
}

#[test]
fn test_blended_wall_function_smooth() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
    let y_plus_values = vec![1.0, 5.0, 10.0, 30.0, 100.0, 1000.0];

    for &y_plus in &y_plus_values {
        let u_plus = treatment.u_plus(y_plus);
        assert!(u_plus > 0.0);
        assert!(u_plus < y_plus * 2.0);
    }
}

#[test]
fn test_blended_wall_function_asymptotic_behavior() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Blended);

    let y_plus_low = 1.0;
    let u_plus_low = treatment.u_plus(y_plus_low);
    assert_relative_eq!(u_plus_low, y_plus_low, epsilon = 0.1);

    let y_plus_high = 1000.0;
    let u_plus_high = treatment.u_plus(y_plus_high);
    let expected_high = (y_plus_high.ln() / KAPPA) + 5.5;
    assert_relative_eq!(u_plus_high, expected_high, epsilon = 1.0);
}

#[test]
fn test_calculate_y_plus_positive() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus = treatment.calculate_y_plus(0.001, 10.0, 1.0, 1e-5);
    assert!(y_plus > 0.0);
    assert!(y_plus.is_finite());
}

#[test]
fn test_calculate_y_plus_scales_with_velocity() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus_1 = treatment.calculate_y_plus(0.001, 10.0, 1.0, 1e-5);
    let y_plus_2 = treatment.calculate_y_plus(0.001, 20.0, 1.0, 1e-5);
    assert!(y_plus_2 > y_plus_1);
}

#[test]
fn test_wall_shear_stress_positive() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let tau_w = treatment.wall_shear_stress(10.0, 0.001, 1.0, 1e-5);
    assert!(tau_w > 0.0);
    assert!(tau_w.is_finite());
}

#[test]
fn test_wall_shear_stress_scales_with_velocity() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let tau_1 = treatment.wall_shear_stress(10.0, 0.001, 1.0, 1e-5);
    let tau_2 = treatment.wall_shear_stress(20.0, 0.001, 1.0, 1e-5);
    assert!(tau_2 > tau_1);
}

#[test]
fn test_wall_k_positive() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let k = treatment.wall_k(0.5, C_MU);
    assert!(k > 0.0);
    assert!(k.is_finite());
}

#[test]
fn test_wall_k_scales_with_u_tau_squared() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let k_1 = treatment.wall_k(1.0, C_MU);
    let k_2 = treatment.wall_k(2.0, C_MU);
    assert_relative_eq!(k_2 / k_1, 4.0, epsilon = 1e-10);
}

#[test]
fn test_wall_epsilon_positive() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let epsilon = treatment.wall_epsilon(0.5, 0.001);
    assert!(epsilon > 0.0);
    assert!(epsilon.is_finite());
}

#[test]
fn test_wall_epsilon_scales_with_u_tau_cubed() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let eps_1 = treatment.wall_epsilon(1.0, 0.001);
    let eps_2 = treatment.wall_epsilon(2.0, 0.001);
    assert_relative_eq!(eps_2 / eps_1, 8.0, epsilon = 1e-10);
}

#[test]
fn test_wall_omega_low_reynolds() {
    let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
    let omega = treatment.wall_omega(0.001, 1e-5, 1.0);
    assert!(omega > 0.0);
    assert!(omega.is_finite());
}

#[test]
fn test_wall_omega_standard() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let omega = treatment.wall_omega(0.001, 1e-5, 1.0);
    assert!(omega > 0.0);
    assert!(omega.is_finite());
}

#[test]
fn test_wall_function_type_standard() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    matches!(treatment.wall_function_type(), WallFunction::Standard);
}

#[test]
fn test_wall_function_type_blended() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Blended);
    matches!(treatment.wall_function_type(), WallFunction::Blended);
}

#[test]
fn test_wall_function_type_low_reynolds() {
    let treatment = WallTreatment::<f64>::new(WallFunction::LowReynolds);
    matches!(treatment.wall_function_type(), WallFunction::LowReynolds);
}

#[test]
fn test_standard_wall_function_continuity() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_trans_1 = Y_PLUS_VISCOUS_SUBLAYER - 0.1;
    let y_trans_2 = Y_PLUS_VISCOUS_SUBLAYER + 0.1;
    let u_1 = treatment.u_plus(y_trans_1);
    let u_2 = treatment.u_plus(y_trans_2);
    assert!((u_2 - u_1).abs() < 1.0);
}

#[test]
fn test_u_plus_monotonically_increasing() {
    let treatment = WallTreatment::<f64>::new(WallFunction::Standard);
    let y_plus_values = vec![1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0];
    let mut prev_u_plus = 0.0;

    for &y_plus in &y_plus_values {
        let u_plus = treatment.u_plus(y_plus);
        assert!(u_plus > prev_u_plus);
        prev_u_plus = u_plus;
    }
}

// --- Spalding's universal law of the wall tests ---

#[test]
fn test_spalding_viscous_sublayer() {
    let u_plus = spalding_u_plus(1.0);
    assert_relative_eq!(u_plus, 1.0, epsilon = 0.01);
}

#[test]
fn test_spalding_buffer_layer() {
    let u_plus = spalding_u_plus(10.0);
    let log_law_value = 10.0_f64.ln() / 0.41 + 5.0;
    assert!(u_plus > 0.0, "u+ must be positive");
    assert!(u_plus < 10.0, "u+ must be less than y+ in buffer layer");
    assert!(
        u_plus < log_law_value + 1.0,
        "u+ should be near log-law value"
    );
}

#[test]
fn test_spalding_log_law() {
    let u_plus = spalding_u_plus(100.0);
    let expected = 100.0_f64.ln() / 0.41 + 5.0;
    let rel_error = (u_plus - expected).abs() / expected;
    assert!(
        rel_error < 0.02,
        "Spalding at y+=100: u+={u_plus}, expected~{expected}, rel_error={rel_error}"
    );
}

#[test]
fn test_spalding_c1_continuity() {
    let eps = 0.01;

    for &y_plus in &[5.0, 30.0] {
        let u_minus = spalding_u_plus(y_plus - eps);
        let u_center = spalding_u_plus(y_plus);
        let u_plus_val = spalding_u_plus(y_plus + eps);

        let deriv_left = (u_center - u_minus) / eps;
        let deriv_right = (u_plus_val - u_center) / eps;

        let deriv_diff = (deriv_right - deriv_left).abs();
        assert!(
            deriv_diff < 0.01,
            "Derivative discontinuity at y+={y_plus}: left={deriv_left}, right={deriv_right}, diff={deriv_diff}"
        );
    }
}
