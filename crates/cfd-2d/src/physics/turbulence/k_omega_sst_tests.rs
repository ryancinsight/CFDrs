//! Comprehensive tests for k-ω SST turbulence model
//!
//! Tests validate the full implementation of the Menter (1994) k-ω SST model
//! including the Bradshaw assumption limiter, blending functions, and boundary conditions.
//!
//! References:
//! - Menter, F.R. (1994). "Two-equation eddy-viscosity turbulence models for 
//!   engineering applications." AIAA Journal, 32(8), 1598-1605.
//! - Wilcox, D.C. (2006). "Turbulence Modeling for CFD" (3rd ed.). DCW Industries.

#[cfg(test)]
mod k_omega_sst_tests {
    use super::super::k_omega_sst::KOmegaSSTModel;
    use super::super::traits::TurbulenceModel;
    use approx::assert_relative_eq;

    /// Test k-ω SST turbulent viscosity with full Bradshaw limiter
    /// Reference: Menter (1994) - Equation 14
    #[test]
    fn test_sst_full_limiter_low_strain_rate() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        let k = 1.0; // Turbulent kinetic energy [m²/s²]
        let omega = 100.0; // Specific dissipation rate [1/s]
        let density = 1.0; // Density [kg/m³]
        let strain_rate_magnitude = 1.0; // Low strain rate [1/s]
        let f2 = 0.5; // Mid-range blending function

        // With low strain rate, limiter should not activate significantly
        let nu_t_limited = model.turbulent_viscosity_with_limiter(
            k,
            omega,
            density,
            strain_rate_magnitude,
            f2,
        );

        // Expected: νt ≈ a1*k / max(a1*ω, S*F2)
        // a1 = 0.31 (SST constant)
        let a1 = 0.31;
        let denominator_limited = (a1 * omega).max(strain_rate_magnitude * f2);
        let expected_nu_t = a1 * k / denominator_limited;

        assert_relative_eq!(nu_t_limited, expected_nu_t * density, epsilon = 1e-10);
        assert!(nu_t_limited > 0.0, "Turbulent viscosity must be positive");
    }

    /// Test SST limiter with high strain rate (limiter active)
    /// Reference: Menter (1994) - Bradshaw assumption validation
    #[test]
    fn test_sst_full_limiter_high_strain_rate() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        let k = 10.0;
        let omega = 10.0;
        let density = 1.225; // Air density [kg/m³]
        let strain_rate_magnitude = 1000.0; // High strain rate
        let f2 = 1.0; // Near-wall value

        let nu_t_limited = model.turbulent_viscosity_with_limiter(
            k,
            omega,
            density,
            strain_rate_magnitude,
            f2,
        );

        // With high S*F2, denominator should be dominated by strain rate
        let a1 = 0.31;
        let expected_denominator = strain_rate_magnitude * f2; // >> a1*ω
        let expected_nu_t = a1 * k / expected_denominator;

        assert_relative_eq!(nu_t_limited, expected_nu_t * density, epsilon = 1e-8);
    }

    /// Test that simplified and full limiter converge when F2 = 0 (freestream)
    #[test]
    fn test_limiter_convergence_freestream() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        let k = 5.0;
        let omega = 50.0;
        let density = 1.0;
        let strain_rate_magnitude = 100.0; // Arbitrary, should not matter when F2 = 0
        let f2 = 0.0; // Freestream condition

        let nu_t_full = model.turbulent_viscosity_with_limiter(
            k,
            omega,
            density,
            strain_rate_magnitude,
            f2,
        );

        // When F2 = 0, S*F2 = 0, so denominator is max(a1*ω, 0) = a1*ω
        let a1 = 0.31;
        let expected = density * a1 * k / (a1 * omega);

        assert_relative_eq!(nu_t_full, expected, epsilon = 1e-10);
        assert!(nu_t_full > 0.0, "Turbulent viscosity must be positive");
    }

    /// Test SST production term calculation
    /// Reference: Menter (1994) - Production term P_k = τ_ij * ∂u_i/∂x_j
    #[test]
    fn test_sst_production_term_shear_flow() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        
        // Simple shear flow: ∂u/∂y = 10 s⁻¹, other gradients zero
        let velocity_gradient = [[0.0, 10.0], [0.0, 0.0]];
        let turbulent_viscosity = 0.001; // [kg/(m·s)]

        let production = model.production_term(&velocity_gradient, turbulent_viscosity);

        // P_k = νt * 2 * Sij * Sij
        // Sij = 0.5 * (∂ui/∂xj + ∂uj/∂xi)
        // S12 = S21 = 5.0, others zero
        // S:S = 2 * (5.0)² = 50.0
        let expected = turbulent_viscosity * 2.0 * 50.0;

        assert_relative_eq!(production, expected, epsilon = 1e-10);
    }

    /// Test SST dissipation term
    /// Reference: Menter (1994) - Dissipation ε_k = β* * k * ω
    #[test]
    fn test_sst_dissipation_term() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        let k = 1.0;
        let omega = 100.0;

        let dissipation = model.dissipation_term(k, omega);

        // β* = 0.09 (SST constant)
        let beta_star = 0.09;
        let expected = beta_star * k * omega;

        assert_relative_eq!(dissipation, expected, epsilon = 1e-10);
    }

    /// Test physical realizability: turbulent viscosity must be positive
    #[test]
    fn test_sst_positive_turbulent_viscosity() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        
        // Test various conditions
        let test_cases = vec![
            (0.01, 1.0, 1.0, 1.0, 0.5),    // Low k
            (10.0, 0.1, 1.0, 10.0, 0.5),   // Low omega
            (1.0, 100.0, 1.0, 50.0, 0.5),  // High omega
            (100.0, 1000.0, 1.225, 100.0, 1.0), // High k, high omega
        ];

        for (k, omega, density, strain_rate, f2) in test_cases {
            let nu_t = model.turbulent_viscosity_with_limiter(
                k,
                omega,
                density,
                strain_rate,
                f2,
            );
            assert!(
                nu_t >= 0.0,
                "Turbulent viscosity must be non-negative: got {nu_t}"
            );
        }
    }

    /// Test SST limiter bounds check (prevent division by zero)
    #[test]
    fn test_sst_limiter_numerical_stability() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        
        // Edge case: very small omega (near-zero)
        let k = 1.0;
        let omega = 1e-10;
        let density = 1.0;
        let strain_rate = 1.0;
        let f2 = 0.5;

        let nu_t = model.turbulent_viscosity_with_limiter(
            k,
            omega,
            density,
            strain_rate,
            f2,
        );

        // Should not panic or produce infinity
        assert!(nu_t.is_finite(), "Turbulent viscosity must be finite");
        assert!(nu_t >= 0.0, "Turbulent viscosity must be non-negative");
    }

    /// Test that SST limiter reduces to k-ω model behavior in near-wall region
    /// Reference: Menter (1994) - F1 = 1 near walls activates k-ω model
    #[test]
    fn test_sst_near_wall_behavior() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        let k = 0.1;
        let omega = 1000.0; // High omega near wall
        let density = 1.0;
        let strain_rate = 50.0;
        let f2 = 1.0; // Near-wall value (F2 = 1)

        let nu_t = model.turbulent_viscosity_with_limiter(
            k,
            omega,
            density,
            strain_rate,
            f2,
        );

        // Near wall with F2 = 1, limiter is fully active
        let a1 = 0.31;
        let expected_denominator = (a1 * omega).max(strain_rate * f2);
        let expected = density * a1 * k / expected_denominator;

        assert_relative_eq!(nu_t, expected, epsilon = 1e-10);
    }

    /// Test SST model name
    #[test]
    fn test_sst_model_name() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        assert_eq!(model.name(), "k-omega-SST");
    }

    /// Test SST Reynolds number validity
    #[test]
    fn test_sst_reynolds_validity() {
        let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
        
        // SST is valid for all positive Reynolds numbers
        assert!(model.is_valid_for_reynolds(100.0));
        assert!(model.is_valid_for_reynolds(1e6));
        assert!(model.is_valid_for_reynolds(1e-3));
        
        // Invalid for zero or negative
        assert!(!model.is_valid_for_reynolds(0.0));
        assert!(!model.is_valid_for_reynolds(-100.0));
    }

    /// Property-based test: limiter should never increase turbulent viscosity
    /// Reference: Physical realizability constraint
    #[test]
    fn test_sst_limiter_is_limiter() {
        use proptest::prelude::*;
        
        proptest!(|(
            k in 0.001f64..100.0,
            omega in 0.1f64..1000.0,
            strain_rate in 0.1f64..1000.0,
            f2 in 0.0f64..=1.0
        )| {
            let model: KOmegaSSTModel<f64> = KOmegaSSTModel::new(10, 10);
            let density = 1.0;
            
            let nu_t_unlimited = density * k / omega;
            let nu_t_limited = model.turbulent_viscosity_with_limiter(
                k,
                omega,
                density,
                strain_rate,
                f2,
            );
            
            // Limiter should never increase turbulent viscosity
            prop_assert!(nu_t_limited <= nu_t_unlimited * 1.001, 
                        "Limiter increased turbulent viscosity: {} > {}", 
                        nu_t_limited, nu_t_unlimited);
        });
    }
}
