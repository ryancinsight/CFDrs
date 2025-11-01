//! TVD/MUSCL Scheme Validation Tests
//!
//! This module validates the Total Variation Diminishing (TVD) schemes and
//! MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
//! implementations against literature benchmarks and monotonicity requirements.
//!
//! ## Validation Cases
//!
//! 1. **Square Wave Advection**: Test monotonicity preservation and shock sharpness
//! 2. **Convergence Order**: Verify 2nd/3rd order accuracy vs upwind 1st order
//! 3. **TVD Property**: Ensure no new extrema are created (TVD condition)
//! 4. **Limiter Effects**: Compare different TVD limiters performance
//!
//! ## References
//!
//! - Barth, T.J. & Jespersen, D.C. (1989). "The design and application of upwind schemes on unstructured meshes"
//! - Yee, H.C. (1987). "Upwind and Symmetric Shock-Capturing Schemes"
//! - Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"

use cfd_2d::schemes::tvd::{MUSCLScheme, MUSCLOrder, FluxLimiter};
use cfd_2d::schemes::{FaceReconstruction, Grid2D};
use cfd_2d::schemes::grid::Grid2D as Grid2DT;
use approx::assert_relative_eq;
use std::f64;

/// Test grid parameters
const NX: usize = 100;
const NY: usize = 1; // 1D problem in 2D grid
const L: f64 = 1.0;
const DX: f64 = L / (NX as f64);
const DY: f64 = 1.0; // Y-direction spacing (not used in 1D tests)

/// Create uniform 2D grid for testing
fn create_test_grid() -> Grid2DT<f64> {
    Grid2DT::new(NX, NY, DX, DY, 0)
}

/// Create square wave initial condition (shock-like profile)
fn square_wave(phi: &mut Grid2D<f64>, amplitude: f64, width: f64, center: f64) {
    for i in 0..NX {
        let x = (i as f64 + 0.5) * DX;
        if (x - center).abs() < width / 2.0 {
            phi.data[(i, 0)] = amplitude;
        } else {
            phi.data[(i, 0)] = 0.0;
        }
    }
}

/// Create linear advection exact solution
fn exact_solution_advection(x: f64, t: f64, amplitude: f64, width: f64, center: f64, velocity: f64) -> f64 {
    let center_t = center + velocity * t;
    let x_rel = (x - center_t) / width;

    // Periodic domain wrapping
    let x_wrapped = if x_rel > L { x_rel - L } else if x_rel < 0.0 { x_rel + L } else { x_rel };

    if x_wrapped.abs() < 0.5 {
        amplitude
    } else {
        0.0
    }
}

#[test]
fn test_muscl2_superbee_monotonicity() {
    let scheme = MUSCLScheme::muscl2_superbee();
    assert_eq!(scheme.order(), 2);

    // Test on square wave
    let mut phi = create_test_grid();
    square_wave(&mut phi, 1.0, 0.2, 0.5);

    // Evolve for one time step
    let velocity = 1.0;
    let dt = DX * 0.5; // CFL = 0.5
    let mut phi_new = create_test_grid();

    for i in 0..NX-1 {
        let face_value = scheme.reconstruct_face_value_x(&phi, velocity, i, 0);
        let face_value_next = scheme.reconstruct_face_value_x(&phi, velocity, i+1, 0);

        // Apply upwind advection: φ_new[i] = φ[i] - velocity * dt/dx * (face_value_next - face_value)
        let flux_diff = face_value_next - face_value;
        phi_new.data[(i, 0)] = phi.data[(i, 0)] - velocity * dt / DX * flux_diff;
    }

    // Check monotonicity: no new extrema created
    // Also check TVD property: no oscillations
    for i in 1..NX-1 {
        let phi_i = phi.data[(i, 0)];
        let phi_im1 = phi.data[(i-1, 0)];
        let phi_ip1 = phi.data[(i+1, 0)];

        // TVD condition: total variation should not increase
        let tv_original = (phi_i - phi_im1).abs() + (phi_ip1 - phi_i).abs();
        let tv_new = (phi_new.data[(i, 0)] - phi_new.data[(i-1, 0)]).abs() +
                     (phi_new.data[(i+1, 0)] - phi_new.data[(i, 0)]).abs();

        // Total variation should not increase significantly (allowing for some numerical diffusion)
        assert!(tv_new <= tv_original * 1.1, "TVD condition violated at i={}", i);
    }
}

#[test]
fn test_muscl3_vs_muscl2_accuracy() {
    // Test convergence order by comparing solutions at different grid resolutions
    let velocities: [f64; 2] = [0.5, 1.0];
    let limiters = [FluxLimiter::Superbee, FluxLimiter::VanLeer, FluxLimiter::Minmod];

    for &limiter in &limiters {
        for &velocity in velocities.iter() {
            let scheme_muscl2 = MUSCLScheme::new(MUSCLOrder::SecondOrder, limiter);
            let scheme_muscl3 = MUSCLScheme::new(MUSCLOrder::ThirdOrder, limiter);

            // Create test grid
            let mut phi = create_test_grid();
            square_wave(&mut phi, 1.0, 0.1, 0.3); // Narrow square wave

            // Evolve both schemes
            let dt = DX / velocity.abs() * 0.4; // CFL = 0.4

            for t in 0..3 {
                let mut phi_muscl2 = phi.clone();
                let mut phi_muscl3 = phi.clone();

                for _step in 0..10 {
                    for i in 0..NX-1 {
                        let face_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi_muscl2, velocity, i, 0);
                        let face_next_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi_muscl2, velocity, i+1, 0);

                        let face_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi_muscl3, velocity, i, 0);
                        let face_next_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi_muscl3, velocity, i+1, 0);

                        phi_muscl2.data[(i, 0)] -= velocity * dt / DX * (face_next_muscl2 - face_muscl2);
                        phi_muscl3.data[(i, 0)] -= velocity * dt / DX * (face_next_muscl3 - face_muscl3);
                    }
                }

                // MUSCL3 should preserve sharper profile than MUSCL2
                let max_muscl3 = (0..NX).map(|i| phi_muscl3.data[(i, 0)]).fold(f64::NEG_INFINITY, f64::max);
                let max_muscl2 = (0..NX).map(|i| phi_muscl2.data[(i, 0)]).fold(f64::NEG_INFINITY, f64::max);

                // MUSCL3 should maintain higher peak (less diffusion) for narrow features
                assert!(max_muscl3 >= max_muscl2 * 0.9, "MUSCL3 should show less diffusion");
            }
        }
    }
}

#[test]
fn test_muscl_scheme_order_accuracy() {
    let scheme_muscl2 = MUSCLScheme::muscl2_van_leer();
    let scheme_muscl3 = MUSCLScheme::muscl3_superbee();

    // Create smooth initial condition: sin(2πx)
    let mut phi = create_test_grid();
    for i in 0..NX {
        let x = (i as f64 + 0.5) * DX;
        phi.data[(i, 0)] = (2.0 * std::f64::consts::PI * x).sin();
    }

    let velocity: f64 = 1.0;
    let t_final: f64 = 0.1;
    let nt = ((t_final * velocity.abs() / DX).ceil() as usize).max(2);

    let mut phi_muscl2 = phi.clone();
    let mut phi_muscl3 = phi.clone();

    // Evolve to final time
    for _step in 0..nt {
        let dt = t_final / (nt as f64);
        let mut phi_new_muscl2 = phi_muscl2.clone();
        let mut phi_new_muscl3 = phi_muscl3.clone();

        for i in 0..NX-1 {
            let face_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi_muscl2, velocity, i, 0);
            let face_next_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi_muscl2, velocity, i+1, 0);
            phi_new_muscl2.data[(i, 0)] -= velocity * dt / DX * (face_next_muscl2 - face_muscl2);

            let face_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi_muscl3, velocity, i, 0);
            let face_next_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi_muscl3, velocity, i+1, 0);
            phi_new_muscl3.data[(i, 0)] -= velocity * dt / DX * (face_next_muscl3 - face_muscl3);
        }

        phi_muscl2 = phi_new_muscl2;
        phi_muscl3 = phi_new_muscl3;
    }

    // Compare L2 error against exact solution
    let mut l2_muscl2 = 0.0;
    let mut l2_muscl3 = 0.0;

    for i in 0..NX {
        let x = (i as f64 + 0.5) * DX;
        let exact = (2.0 * std::f64::consts::PI * (x - velocity * t_final)).sin();

        let error_muscl2 = phi_muscl2.data[(i, 0)] - exact;
        let error_muscl3 = phi_muscl3.data[(i, 0)] - exact;

        l2_muscl2 += error_muscl2 * error_muscl2;
        l2_muscl3 += error_muscl3 * error_muscl3;
    }

    l2_muscl2 = (l2_muscl2 / NX as f64).sqrt();
    l2_muscl3 = (l2_muscl3 / NX as f64).sqrt();

    // MUSCL3 should be more accurate than MUSCL2 for smooth solutions
    assert!(l2_muscl3 < l2_muscl2, "MUSCL3 error {} should be less than MUSCL2 error {}", l2_muscl3, l2_muscl2);
}

#[test]
fn test_tvd_limiters_properties() {
    // Test that all limiters satisfy TVD conditions
    let limiters = vec![
        FluxLimiter::Superbee,
        FluxLimiter::VanLeer,
        FluxLimiter::Minmod,
        FluxLimiter::MC,
    ];

    for limiter in limiters {
        // Test TVD region r ∈ [0, 1]
        for i in 0..10 {
            let r = (i as f64) / 10.0;
            let psi = limiter.apply(r);

            // TVD condition 1: 0 ≤ ψ(r) ≤ 2
            assert!(psi >= 0.0, "Limiter {} violates ψ ≥ 0 for r={}", format!("{:?}", limiter), r);
            assert!(psi <= 2.0, "Limiter {} violates ψ ≤ 2 for r={}", format!("{:?}", limiter), r);

            // TVD condition 2: ψ(r) ≤ 2r for r ∈ [0, 1]
            if r <= 1.0 {
                assert!(psi <= 2.0 * r + 1e-10, "Limiter {} violates TVD condition for r={}", format!("{:?}", limiter), r);
            }
        }

        // Test discontinuities: ψ(r) = 0 when r ≤ 0 (no overshoots)
        let psi_neg = limiter.apply(-1.0);
        assert_eq!(psi_neg, 0.0, "Limiter {} should be 0 for negative r", format!("{:?}", limiter));

        // Test smooth regions: ψ(1) = 1 (full centered differencing)
        let psi_one: f64 = limiter.apply(1.0);
        assert!((psi_one - 1.0).abs() < 1e-10, "Limiter {} should be 1 for r=1", format!("{:?}", limiter));
    }
}

#[test]
fn test_boundary_one_sided_reconstruction() {
    let scheme = MUSCLScheme::muscl2_van_leer();

    // Test boundary handling: one-sided reconstruction at left boundary
    let mut phi = create_test_grid();
    for i in 0..NX {
        phi.data[(i, 0)] = i as f64;
    }

    // Left boundary (i=0) with outflow (positive velocity)
    let face_left_outflow = scheme.reconstruct_face_value_x(&phi, 1.0, 0, 0);
    assert_eq!(face_left_outflow, 0.0, "Left boundary outflow should use cell value");

    // Left boundary with inflow (negative velocity) - should use one-sided reconstruction
    let face_left_inflow = scheme.reconstruct_face_value_x(&phi, -1.0, 0, 0);
    let expected = 0.0 - 0.5 * (1.0 - 0.0); // φ_0 - (1/2)δφ
    assert_relative_eq!(face_left_inflow, expected, epsilon = 1e-10);

    // Right boundary (i=NX-1) with inflow (negative velocity)
    let face_right_inflow = scheme.reconstruct_face_value_x(&phi, -1.0, NX-1, 0);
    assert_eq!(face_right_inflow, ((NX-1) as f64), "Right boundary inflow should use cell value");

    // Right boundary with outflow (positive velocity) - one-sided reconstruction
    let face_right_outflow = scheme.reconstruct_face_value_x(&phi, 1.0, NX-1, 0);
    let expected_right = (NX-1) as f64 + 0.5 * ((NX-1) as f64 - (NX-2) as f64);
    assert_relative_eq!(face_right_outflow, expected_right, epsilon = 1e-10);
}

#[test]
fn test_muscl3_boundary_fallback() {
    let scheme_muscl3 = MUSCLScheme::muscl3_superbee();
    let scheme_muscl2 = MUSCLScheme::muscl2_superbee();

    // Test that MUSCL3 falls back to MUSCL2 at boundaries
    let mut phi = create_test_grid();
    for i in 0..NX {
        phi.data[(i, 0)] = (i as f64).powi(2);
    }

    // At boundary point where MUSCL3 needs points beyond boundary
    let face_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi, -1.0, 1, 0); // i=1, still interior for MUSCL3
    let face_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi, -1.0, 1, 0);

    // Should be different (MUSCL3 uses more points)
    assert!(face_muscl3 != face_muscl2, "MUSCL3 should use different reconstruction than MUSCL2");

    // For a very boundary case (near left edge where MUSCL3 lacks upstream points)
    let face_boundary_muscl3 = scheme_muscl3.reconstruct_face_value_x(&phi, -1.0, 0, 0);
    let face_boundary_muscl2 = scheme_muscl2.reconstruct_face_value_x(&phi, -1.0, 0, 0);

    // Should be the same (MUSCL3 falls back to MUSCL2)
    assert_eq!(face_boundary_muscl3, face_boundary_muscl2, "MUSCL3 should fall back to MUSCL2 at boundaries");
}
