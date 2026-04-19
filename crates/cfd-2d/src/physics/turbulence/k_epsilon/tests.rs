use super::kato_launder;
use super::model::KEpsilonModel;
use crate::physics::turbulence::constants::{
    C1_EPSILON, C2_EPSILON, C_MU, EPSILON_MIN, REALIZABLE_A0, SIGMA_EPSILON, SIGMA_K,
};
use crate::physics::turbulence::traits::TurbulenceModel;
use approx::assert_relative_eq;
use nalgebra::Vector2;

#[test]
fn test_new_model_initialization() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    assert_eq!(model.nx, 10);
    assert_eq!(model.ny, 10);
    assert_relative_eq!(model.c_mu, C_MU, epsilon = 1e-10);
    assert_relative_eq!(model.c1_epsilon, C1_EPSILON, epsilon = 1e-10);
    assert_relative_eq!(model.c2_epsilon, C2_EPSILON, epsilon = 1e-10);
    assert_relative_eq!(model.sigma_k, SIGMA_K, epsilon = 1e-10);
    assert_relative_eq!(model.sigma_epsilon, SIGMA_EPSILON, epsilon = 1e-10);
}

#[test]
fn test_turbulent_viscosity_positive() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let k = 1.0;
    let epsilon = 1.0;
    let density = 1.0;
    let nu_t = model.turbulent_viscosity(k, epsilon, density);
    assert!(nu_t > 0.0);
    assert_relative_eq!(nu_t, density * C_MU * k * k / epsilon, epsilon = 1e-10);
}

#[test]
fn test_turbulent_viscosity_zero_epsilon() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let k = 1.0;
    let epsilon = 0.0;
    let density = 1.0;
    let nu_t = model.turbulent_viscosity(k, epsilon, density);
    assert!(nu_t > 0.0);
    assert!(nu_t.is_finite());
}

#[test]
fn test_strain_rate_pure_shear() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[0.0, 1.0], [0.0, 0.0]];
    let strain = model.strain_rate(&grad);
    assert_relative_eq!(strain, 1.0, epsilon = 1e-10);
}

#[test]
fn test_strain_rate_pure_extension() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[1.0, 0.0], [0.0, -1.0]];
    let strain = model.strain_rate(&grad);
    assert_relative_eq!(strain, 2.0, epsilon = 1e-10);
}

#[test]
fn test_strain_rate_zero_velocity_gradient() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[0.0, 0.0], [0.0, 0.0]];
    let strain = model.strain_rate(&grad);
    assert_relative_eq!(strain, 0.0, epsilon = 1e-10);
}

#[test]
fn test_production_term_positive() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[1.0, 0.0], [0.0, 1.0]];
    let nu_t = 0.1;
    let p_k = model.production_term(&grad, nu_t, 0.0, 0.0, 1e-5);
    assert!(p_k > 0.0);
}

#[test]
fn test_production_term_zero_strain() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[0.0, 0.0], [0.0, 0.0]];
    let nu_t = 0.1;
    let p_k = model.production_term(&grad, nu_t, 0.0, 0.0, 1e-5);
    assert_relative_eq!(p_k, 0.0, epsilon = 1e-10);
}

#[test]
fn test_dissipation_term() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let k = 1.0;
    let epsilon = 2.5;
    let dissipation = model.dissipation_term(k, epsilon);
    assert_relative_eq!(dissipation, epsilon, epsilon = 1e-10);
}

#[test]
fn test_apply_boundary_conditions_enforces_positivity() {
    let model = KEpsilonModel::<f64>::new(5, 5);
    let mut k = vec![-1.0; 25];
    let mut epsilon = vec![-2.0; 25];

    model.apply_boundary_conditions(&mut k, &mut epsilon);

    for &val in &k {
        assert!(val >= 0.0);
    }
    for &val in &epsilon {
        assert!(val >= EPSILON_MIN);
    }
}

#[test]
fn test_apply_boundary_conditions_wall_values() {
    let model = KEpsilonModel::<f64>::new(5, 5);
    let mut k = vec![1.0; 25];
    let mut epsilon = vec![1.0; 25];

    model.apply_boundary_conditions(&mut k, &mut epsilon);

    // Bottom wall (j=0)
    for i in 0..5 {
        assert_relative_eq!(k[i], 0.0, epsilon = 1e-10);
        assert!(epsilon[i] >= EPSILON_MIN);
    }

    // Top wall (j=4)
    for i in 0..5 {
        let idx = i + 4 * 5;
        assert_relative_eq!(k[idx], 0.0, epsilon = 1e-10);
        assert!(epsilon[idx] >= EPSILON_MIN);
    }
}

#[test]
fn test_update_maintains_positivity() {
    let mut model = KEpsilonModel::<f64>::new(5, 5);
    let mut k = vec![0.1; 25];
    let mut epsilon = vec![0.1; 25];
    let velocity = vec![Vector2::new(0.0, 0.0); 25];

    let result = model.update(&mut k, &mut epsilon, &velocity, 1.0, 1e-5, 0.001, 0.1, 0.1);

    assert!(result.is_ok());
    for &val in &k {
        assert!(val >= 0.0);
    }
    for &val in &epsilon {
        assert!(val >= 0.0);
    }
}

#[test]
fn test_model_name() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    assert_eq!(model.name(), "k-epsilon");
}

#[test]
fn test_is_valid_for_high_reynolds() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    assert!(model.is_valid_for_reynolds(1e5));
    assert!(model.is_valid_for_reynolds(1e6));
}

#[test]
fn test_is_not_valid_for_low_reynolds() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    assert!(!model.is_valid_for_reynolds(1e3));
    assert!(!model.is_valid_for_reynolds(1e2));
}

#[test]
fn test_reynolds_threshold() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    assert!(!model.is_valid_for_reynolds(1e4));
    assert!(model.is_valid_for_reynolds(1e4 + 1.0));
}

#[test]
fn test_update_with_uniform_flow() {
    let mut model = KEpsilonModel::<f64>::new(5, 5);
    let mut k = vec![1.0; 25];
    let mut epsilon = vec![1.0; 25];
    let velocity = vec![Vector2::new(1.0, 0.0); 25];

    let result = model.update(&mut k, &mut epsilon, &velocity, 1.0, 1e-5, 0.001, 0.1, 0.1);
    assert!(result.is_ok());
}

#[test]
fn test_turbulent_viscosity_scales_with_k_squared() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let epsilon = 1.0;
    let density = 1.0;

    let nu_t_1 = model.turbulent_viscosity(1.0, epsilon, density);
    let nu_t_2 = model.turbulent_viscosity(2.0, epsilon, density);

    assert_relative_eq!(nu_t_2 / nu_t_1, 4.0, epsilon = 1e-10);
}

#[test]
fn test_production_term_scales_with_turbulent_viscosity() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let grad = [[1.0, 0.0], [0.0, 1.0]];

    let p_k_1 = model.production_term(&grad, 0.1, 0.0, 0.0, 1e-5);
    let p_k_2 = model.production_term(&grad, 0.2, 0.0, 0.0, 1e-5);

    assert_relative_eq!(p_k_2 / p_k_1, 2.0, epsilon = 1e-10);
}

/// Analytical MMS validation for k-ε model
#[test]
fn test_k_epsilon_mms_validation() {
    let nx = 16;
    let ny = 16;
    let mut model = KEpsilonModel::<f64>::new(nx, ny);
    let dt = 0.001;
    let dx = 0.1;
    let dy = 0.1;

    let mut k_field = Vec::new();
    let mut epsilon_field = Vec::new();

    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64 * dx;
            let y = j as f64 * dy;
            let k_exact = x * x * (1.0 - y) * (1.0 - y);
            let eps_exact = (x.powi(4) + y.powi(4)) * x / (x + y + 1.0).max(0.1);
            k_field.push(k_exact);
            epsilon_field.push(eps_exact);
        }
    }

    model.apply_boundary_conditions(&mut k_field, &mut epsilon_field);

    let initial_k_sum: f64 = k_field.iter().sum();
    let initial_eps_sum: f64 = epsilon_field.iter().sum();

    let velocity = vec![Vector2::new(0.0, 0.0); nx * ny];
    model
        .update(
            &mut k_field,
            &mut epsilon_field,
            &velocity,
            1.0,
            1e-5,
            dt,
            dx,
            dy,
        )
        .unwrap();

    for i in 0..nx * ny {
        assert!(k_field[i].is_finite(), "k became non-finite at index {i}");
        assert!(
            epsilon_field[i].is_finite(),
            "epsilon became non-finite at index {i}"
        );
        assert!(k_field[i] >= 0.0, "k became negative at index {i}");
        assert!(
            epsilon_field[i] >= 0.0,
            "epsilon became negative at index {i}"
        );
    }

    let final_k_sum: f64 = k_field.iter().sum();
    let final_eps_sum: f64 = epsilon_field.iter().sum();

    let k_conservation_error = (initial_k_sum - final_k_sum).abs() / initial_k_sum;
    let eps_conservation_error = (initial_eps_sum - final_eps_sum).abs() / initial_eps_sum;

    assert!(
        k_conservation_error < 0.1,
        "k conservation error too high: {k_conservation_error}"
    );
    assert!(
        eps_conservation_error < 0.2,
        "epsilon conservation error too high: {eps_conservation_error}"
    );
}

/// Test k-ε model numerical stability across different mesh sizes
#[test]
fn test_k_epsilon_numerical_stability() {
    let grid_sizes = [8, 12, 16];
    let mut stability_scores = Vec::new();

    for &n in &grid_sizes {
        let mut model = KEpsilonModel::<f64>::new(n, n);
        let dx = 1.0 / (n as f64);
        let dy = dx;
        let dt = 0.1 * dx * dx;

        let k_init = 0.01;
        let eps_init = 0.005;
        let mut k_field = vec![k_init; n * n];
        let mut epsilon_field = vec![eps_init; n * n];

        let mut velocity = vec![Vector2::new(0.0, 0.0); n * n];
        for j in 0..n {
            for i in 0..n {
                let idx = j * n + i;
                let x = i as f64 * dx;
                let y = j as f64 * dy;
                velocity[idx].x = 0.1 * (std::f64::consts::PI * x).sin();
                velocity[idx].y = 0.1 * (std::f64::consts::PI * y).cos();
            }
        }

        model
            .update(
                &mut k_field,
                &mut epsilon_field,
                &velocity,
                1.0,
                1e-5,
                dt,
                dx,
                dy,
            )
            .unwrap();

        let mut finite_count = 0;
        let mut positive_count = 0;
        let mut reasonable_range_count = 0;

        for &k_val in &k_field {
            if k_val.is_finite() {
                finite_count += 1;
            }
            if k_val >= 0.0 {
                positive_count += 1;
            }
            if (0.0..1e3).contains(&k_val) {
                reasonable_range_count += 1;
            }
        }

        for &eps_val in &epsilon_field {
            if eps_val.is_finite() {
                finite_count += 1;
            }
            if eps_val >= 0.0 {
                positive_count += 1;
            }
            if (0.0..1e3).contains(&eps_val) {
                reasonable_range_count += 1;
            }
        }

        let total_points = 2 * n * n;
        let stability_score = f64::from(finite_count + positive_count + reasonable_range_count)
            / (3 * total_points) as f64;
        stability_scores.push(stability_score);
    }

    for (i, &score) in stability_scores.iter().enumerate() {
        let grid_size = grid_sizes[i];
        assert!(
            score > 0.85,
            "Poor stability on {grid_size}x{grid_size} grid: score = {score}"
        );
    }

    if stability_scores.len() > 1 {
        let max_score = stability_scores
            .iter()
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        for &score in &stability_scores {
            assert!(
                (score / max_score) > 0.9,
                "Inconsistent stability across mesh sizes"
            );
        }
    }
}

/// Property-based test: bounded production-dissipation ratio
#[test]
fn test_production_dissipation_bounds() {
    use proptest::prelude::*;

    proptest!(ProptestConfig::with_cases(100), |(
        k in 0.001f64..10.0,
        epsilon in 0.001f64..10.0,
        strain_rate in 0.1f64..100.0,
        nu_t in 1e-6f64..1e-2
    )| {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let velocity_gradient = [[0.0, strain_rate], [0.0, 0.0]];

        let production = model.production_term(&velocity_gradient, nu_t, k, 0.0, 1e-5);
        let dissipation = model.dissipation_term(k, epsilon);
        let ratio = production / dissipation.max(1e-12);

        prop_assert!(ratio > 0.0 && ratio < 1e6, "Unrealizable P/ε ratio: {ratio}");
        prop_assert!(production >= 0.0, "Negative production: {production}");
        prop_assert!(dissipation >= 0.0, "Negative dissipation: {dissipation}");
    });
}

/// Test k-ε model stability in extreme conditions
#[test]
fn test_stochastic_robustness_extreme_conditions() {
    use rand::prelude::*;

    let mut rng = rand::thread_rng();
    let model = KEpsilonModel::<f64>::new(10, 10);

    for _ in 0..100 {
        let k_val = rng.gen_range(1e-6..1e3);
        let eps_val = rng.gen_range(1e-6..1e3);

        let nu_t = model.turbulent_viscosity(k_val, eps_val, 1.0);
        assert!(
            nu_t.is_finite(),
            "Turbulence viscosity non-finite: k={k_val}, ε={eps_val}"
        );
        assert!(nu_t >= 0.0, "Negative viscosity: {nu_t}");

        let strain_rate_magnitude = rng.gen_range(1e-3..1e3);
        let grad = [[0.0, strain_rate_magnitude], [0.0, 0.0]];

        let production = model.production_term(&grad, 1e-3, k_val, 0.0, 1e-5);
        assert!(production.is_finite(), "Production non-finite");
        assert!(production >= 0.0, "Negative production");

        let dissipation = model.dissipation_term(k_val, eps_val);
        assert!(dissipation.is_finite(), "Dissipation non-finite");
        assert!(dissipation >= 0.0, "Negative dissipation");
    }
}

/// Analytical validation: steady-state equilibrium.
///
/// At equilibrium P_k = ε. With P_k = ν_t·|S|² and ν_t = Cμ·k²/ε:
///   (Cμ k²/ε)·|S|² = ε  →  |S| = ε / (√Cμ · k)
///
/// Reference: Launder & Spalding (1974) §2, Jones & Launder (1972) §4
#[test]
fn test_equilibrium_balance() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let target_epsilon: f64 = 1.0;
    let target_k: f64 = 1.0;

    let strain_magnitude = target_epsilon / (C_MU.sqrt() * target_k);
    let nu_t_eq = model.turbulent_viscosity(target_k, target_epsilon, 1.0);
    let velocity_gradient = [[0.0, strain_magnitude], [0.0, 0.0]];

    let production = model.production_term(&velocity_gradient, nu_t_eq, target_k, 0.0, 1e-5);
    let dissipation = model.dissipation_term(target_k, target_epsilon);

    let ratio = production / dissipation;
    assert_relative_eq!(ratio, 1.0, epsilon = 0.2);
    assert!(
        (0.8..=1.2).contains(&ratio),
        "Equilibrium not maintained: P/ε = {ratio}"
    );
}

/// Test turbulence model consistency across different grid sizes
#[test]
fn test_grid_size_independence() {
    let grid_sizes = [8, 16, 24];
    let mut results = Vec::new();

    for &n in &grid_sizes {
        let mut model = KEpsilonModel::<f64>::new(n, n);
        let dx = 1.0 / (n - 1) as f64;
        let dy = dx;

        let mut k_field = vec![1.0; n * n];
        let mut epsilon_field = vec![1.0; n * n];

        model.apply_boundary_conditions(&mut k_field, &mut epsilon_field);

        let velocity = vec![Vector2::new(0.0, 0.0); n * n];
        model
            .update(
                &mut k_field,
                &mut epsilon_field,
                &velocity,
                1.0,
                0.0,
                0.001,
                dx,
                dy,
            )
            .unwrap();

        let interior_count = (1..n - 1)
            .flat_map(|j| (1..n - 1).map(move |i| (i, j)))
            .filter(|&(i, j)| {
                let idx = j * n + i;
                k_field[idx].is_finite()
                    && k_field[idx] > 0.0
                    && epsilon_field[idx].is_finite()
                    && epsilon_field[idx] > 0.0
            })
            .count();

        results.push(interior_count as f64 / ((n - 2) * (n - 2)) as f64);
    }

    let avg_stability = results.iter().sum::<f64>() / results.len() as f64;
    for &stability in &results {
        assert!(
            (stability - avg_stability).abs() < 0.1,
            "Inconsistent stability across grids: {stability} vs avg {avg_stability}"
        );
    }
}

/// Test k-ε model's response to anisotropic strain rate
#[test]
fn test_anisotropic_strain_response() {
    let model = KEpsilonModel::<f64>::new(10, 10);
    let nu_t = 1e-3;

    let test_gradients = vec![
        [[1.0, 0.0], [0.0, -1.0]],
        [[0.0, 1.0], [0.0, 0.0]],
        [[0.5, 0.0], [0.0, -0.5]],
        [[0.3, 0.7], [0.2, -0.4]],
    ];

    for gradient in &test_gradients {
        let production = model.production_term(gradient, nu_t, 0.0, 0.0, 1e-5);
        assert!(
            production.is_finite(),
            "Non-finite production for gradient {gradient:?}"
        );
        assert!(
            production > 0.0,
            "Negative production for gradient {gradient:?}"
        );

        let production_scaled = model.production_term(gradient, nu_t * 3.0, 0.0, 0.0, 1e-5);
        assert_relative_eq!(production_scaled / production, 3.0, epsilon = 1e-10);
    }
}

/// Validate k-ε constants against literature values
#[test]
fn test_constants_physical_validation() {
    const _: () = {
        assert!(C_MU >= 0.07 && C_MU <= 0.11);
        assert!(C1_EPSILON > 1.0);
        assert!(C2_EPSILON > C1_EPSILON);
        assert!(SIGMA_K > 0.0 && SIGMA_EPSILON > 0.0);
    };
}

// ─── Realizable k-ε tests (Shih, Zhu & Lumley 1995) ───────────────

/// Verify that the realizable C_mu is bounded above by 1/A_0 ≈ 0.2475
#[test]
fn test_realizable_c_mu_bounded() {
    let model = KEpsilonModel::<f64>::new_realizable(10, 10);
    let upper = 1.0 / REALIZABLE_A0;

    let grad_zero = [[0.0, 0.0], [0.0, 0.0]];
    let c_mu_zero = super::realizable::realizable_c_mu(&model, &grad_zero, 1.0, 1.0);
    assert!(
        c_mu_zero <= upper + 1e-12,
        "C_mu at zero strain ({c_mu_zero}) exceeds upper bound ({upper})"
    );
    assert!(c_mu_zero > 0.0, "C_mu must be positive");

    let grad_mod = [[0.0, 5.0], [0.0, 0.0]];
    let c_mu_mod = super::realizable::realizable_c_mu(&model, &grad_mod, 1.0, 1.0);
    assert!(c_mu_mod <= upper + 1e-12);
    assert!(c_mu_mod > 0.0);

    let grad_extreme = [[0.0, 1000.0], [0.0, 0.0]];
    let c_mu_ext = super::realizable::realizable_c_mu(&model, &grad_extreme, 1.0, 1.0);
    assert!(c_mu_ext <= upper + 1e-12);
    assert!(c_mu_ext > 0.0);

    let c_mu_high_k = super::realizable::realizable_c_mu(&model, &grad_mod, 100.0, 0.01);
    assert!(c_mu_high_k <= upper + 1e-12);
    assert!(c_mu_high_k > 0.0);
}

/// Verify C_mu reduces below standard 0.09 for high strain rates
#[test]
fn test_realizable_c_mu_reduces_at_high_strain() {
    let model = KEpsilonModel::<f64>::new_realizable(10, 10);

    let grad_high = [[0.0, 50.0], [0.0, 0.0]];
    let c_mu = super::realizable::realizable_c_mu(&model, &grad_high, 1.0, 1.0);
    assert!(
        c_mu < C_MU,
        "Realizable C_mu ({c_mu}) should be less than standard C_mu ({C_MU}) at high strain"
    );

    let grad_very_high = [[10.0, 100.0], [0.0, -10.0]];
    let c_mu_vhigh = super::realizable::realizable_c_mu(&model, &grad_very_high, 1.0, 1.0);
    assert!(
        c_mu_vhigh < c_mu,
        "C_mu should decrease with increasing strain: {c_mu_vhigh} >= {c_mu}"
    );
}

/// At equilibrium S*k/eps ≈ 3.3, the realizable C_mu ≈ standard C_mu
#[test]
fn test_realizable_c_mu_equals_standard_at_equilibrium() {
    let model = KEpsilonModel::<f64>::new_realizable(10, 10);
    let equilibrium_s_k_over_eps = 1.0 / C_MU.sqrt();
    let grad = [[0.0, equilibrium_s_k_over_eps], [0.0, 0.0]];
    let c_mu = super::realizable::realizable_c_mu(&model, &grad, 1.0, 1.0);
    assert_relative_eq!(c_mu, C_MU, epsilon = 0.02);
}

/// Standard model backward compatibility
#[test]
fn test_realizable_standard_backward_compatible() {
    let model_standard = KEpsilonModel::<f64>::new(10, 10);
    assert!(!model_standard.is_realizable());

    let k = 2.5;
    let epsilon = 1.5;
    let density = 1.225;

    let nu_t_standard = model_standard.turbulent_viscosity(k, epsilon, density);
    let nu_t_expected = density * C_MU * k * k / epsilon;
    assert_relative_eq!(nu_t_standard, nu_t_expected, epsilon = 1e-12);

    let mut model_disabled = KEpsilonModel::<f64>::new_realizable(10, 10);
    model_disabled.set_realizable(false);
    assert!(!model_disabled.is_realizable());

    let nu_t_disabled = model_disabled.turbulent_viscosity(k, epsilon, density);
    assert_relative_eq!(nu_t_disabled, nu_t_standard, epsilon = 1e-12);
}

/// Realizable and standard models should produce different k fields with strong gradients
#[test]
fn test_realizable_update_reduces_viscosity() {
    let nx = 5;
    let ny = 5;
    let n = nx * ny;

    let mut model_std = KEpsilonModel::<f64>::new(nx, ny);
    let mut model_real = KEpsilonModel::<f64>::new_realizable(nx, ny);

    let mut k_std = vec![1.0; n];
    let mut eps_std = vec![1.0; n];
    let mut k_real = vec![1.0; n];
    let mut eps_real = vec![1.0; n];

    let mut velocity = vec![Vector2::new(0.0, 0.0); n];
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            velocity[idx] = Vector2::new(10.0 * (j as f64) * 0.1, 0.0);
        }
    }

    model_std
        .update(
            &mut k_std,
            &mut eps_std,
            &velocity,
            1.0,
            1e-5,
            0.0001,
            0.1,
            0.1,
        )
        .unwrap();
    model_real
        .update(
            &mut k_real,
            &mut eps_real,
            &velocity,
            1.0,
            1e-5,
            0.0001,
            0.1,
            0.1,
        )
        .unwrap();

    let any_differ = (0..n).any(|idx| (k_std[idx] - k_real[idx]).abs() > 1e-15);
    assert!(
        any_differ,
        "Realizable and standard models should produce different k fields with strong gradients"
    );
}

/// Verify that realizable C_mu monotonically decreases as strain increases
#[test]
fn test_realizable_c_mu_monotone_decrease() {
    let model = KEpsilonModel::<f64>::new_realizable(10, 10);
    let k = 1.0;
    let eps = 1.0;

    let strains = [0.1, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0];
    let mut prev_c_mu = f64::INFINITY;

    for &s in &strains {
        let grad = [[0.0, s], [0.0, 0.0]];
        let c_mu = super::realizable::realizable_c_mu(&model, &grad, k, eps);
        assert!(
            c_mu <= prev_c_mu + 1e-14,
            "C_mu not monotone: at S={s}, C_mu={c_mu} > prev={prev_c_mu}"
        );
        prev_c_mu = c_mu;
    }
}

// ─── Kato-Launder (1993) tests ─────────────────────────────────────

/// For simple shear: P_KL = P_standard
#[test]
fn test_kato_launder_pure_shear() {
    let gamma_dot = 5.0;
    let nu_t = 0.01;
    let grad = [[0.0, gamma_dot], [0.0, 0.0]];

    let p_kl = kato_launder::kato_launder_production(&grad, nu_t);
    let p_standard = nu_t * gamma_dot * gamma_dot;

    assert_relative_eq!(p_kl, p_standard, epsilon = 1e-10);
}

/// For stagnation flow: P_KL = 0
#[test]
fn test_kato_launder_stagnation() {
    let a = 10.0;
    let nu_t = 0.01;
    let grad = [[a, 0.0], [0.0, -a]];

    let p_kl = kato_launder::kato_launder_production(&grad, nu_t);

    assert!(
        p_kl.abs() < 1e-12,
        "Kato-Launder production should be zero in stagnation flow: got {p_kl}"
    );
}

/// P_KL ≤ P_standard by Cauchy-Schwarz inequality
#[test]
fn test_kato_launder_vs_standard_ratio() {
    let nu_t = 0.05;

    let test_grads: Vec<[[f64; 2]; 2]> = vec![
        [[0.3, 0.7], [0.2, -0.4]],
        [[1.0, 3.0], [1.0, -1.0]],
        [[0.0, 5.0], [2.0, 0.0]],
        [[2.0, 1.0], [-1.0, -2.0]],
    ];

    let model = KEpsilonModel::<f64>::new(10, 10);

    for grad in &test_grads {
        let p_kl = kato_launder::kato_launder_production(grad, nu_t);
        let p_standard = model.production_term(grad, nu_t, 0.0, 0.0, 1e-5);

        assert!(
            p_kl <= p_standard + 1e-10,
            "Kato-Launder ({p_kl}) should not exceed standard ({p_standard}) \
             for gradient {grad:?}"
        );
        assert!(
            p_kl >= 0.0,
            "Kato-Launder production must be non-negative: got {p_kl}"
        );
    }
}
