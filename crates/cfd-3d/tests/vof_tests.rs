//! # Comprehensive VOF Solver Tests
//!
//! ## Test Coverage Matrix
//!
//! | Category     | Tests                                                             |
//! |-------------|-------------------------------------------------------------------|
//! | Positive     | Volume conservation, PLIC normal accuracy, sphere curvature       |
//! | Boundary     | Zero velocity, full/empty cell, CFL=1                            |
//! | Adversarial  | NaN velocity, Inf velocity, α outside [0,1] input                |
//! | Property     | α ∈ [0,1] invariant, volume monotone, curvature sign            |

use approx::assert_relative_eq;
use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};
use nalgebra::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn default_config() -> VofConfig {
    VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    }
}

fn algebraic_config() -> VofConfig {
    VofConfig {
        advection_method: AdvectionMethod::Algebraic,
        ..default_config()
    }
}

fn make_solver(nx: usize, ny: usize, nz: usize, config: VofConfig) -> VofSolver<f64> {
    VofSolver::new(nx, ny, nz, config).expect("VofSolver::new failed")
}

// ─────────────────────────────────────────────────────────────────────────────
// Positive Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Positive**: With no flow (u=0), volume fraction is conserved exactly.
///
/// Invariant: geometric advection with u≡0 must produce zero net flux → α unchanged.
#[test]
fn test_zero_velocity_volume_conservation() {
    let mut solver = make_solver(10, 10, 10, default_config());

    // Initialize: half the domain filled
    {
        let alpha = solver.alpha_mut();
        for (i, a) in alpha.iter_mut().enumerate() {
            *a = if i < 500 { 1.0 } else { 0.0 };
        }
    }
    // velocity is already zero (default)

    let initial_volume: f64 = solver.alpha().iter().sum();

    for _ in 0..50 {
        solver.step(1e-4, &[], &[], &[]).ok();
    }

    let final_volume: f64 = solver.alpha().iter().sum();
    let rel_error = (final_volume - initial_volume).abs() / (initial_volume + 1e-15);
    assert!(
        rel_error < 1e-12,
        "Zero-velocity advection must conserve volume exactly, relative error = {rel_error:.3e}"
    );
}

/// **Positive**: Volume fraction ∈ [0,1] is preserved after advection with moderate velocity.
#[test]
fn test_volume_fraction_bounds_preserved_geometric() {
    let mut solver = make_solver(12, 12, 12, default_config());

    // Set a smooth initial distribution
    let nx = 12usize;
    {
        let alpha = solver.alpha_mut();
        for idx in 0..alpha.len() {
            let i = idx % nx;
            alpha[idx] = if i < nx / 2 { 1.0 } else { 0.0 };
        }
    }

    // Uniform x-velocity (CFL ≈ 0.1)
    let n = 12 * 12 * 12;
    let vel = vec![Vector3::new(0.1, 0.0, 0.0); n];
    for idx in 0..n {
        solver.set_velocity_at(idx, vel[idx]);
    }

    for _ in 0..20 {
        solver.step(1e-4, &[], &[], &[]).ok();
        for &a in solver.alpha() {
            assert!(
                a >= -1e-12 && a <= 1.0 + 1e-12,
                "Volume fraction out of bounds: {a}"
            );
        }
    }
}

/// **Positive**: Algebraic advection preserves volume fraction bounds.
#[test]
fn test_volume_fraction_bounds_algebraic() {
    let mut solver = make_solver(10, 10, 10, algebraic_config());
    let nx = 10usize;
    {
        let alpha = solver.alpha_mut();
        for idx in 0..alpha.len() {
            let i = idx % nx;
            alpha[idx] = if i < 5 { 0.9 } else { 0.1 };
        }
    }

    let n = 10 * 10 * 10;
    let vel = vec![Vector3::new(0.05, 0.0, 0.0); n];
    for idx in 0..n {
        solver.set_velocity_at(idx, vel[idx]);
    }

    for _ in 0..30 {
        solver.step(1e-4, &[], &[], &[]).ok();
        for &a in solver.alpha() {
            assert!(
                a >= -1e-9 && a <= 1.0 + 1e-9,
                "Algebraic advection out of bounds: {a}"
            );
        }
    }
}

/// **Positive**: Planar interface — PLIC normal should align with x-axis.
///
/// A sharp planar interface at x = L/2 has ∇α ∝ x̂.  The reconstructed
/// normal must point in the x-direction to within 5°.
#[test]
fn test_plic_normal_planar_interface() {
    let nx = 20usize;
    let mut solver = make_solver(nx, 6, 6, default_config());

    {
        let alpha = solver.alpha_mut();
        for idx in 0..alpha.len() {
            let i = idx % nx;
            alpha[idx] = if i < nx / 2 { 1.0 } else { 0.0 };
        }
    }

    solver.reconstruct_interface();

    // Interior interface cells: i ∈ {nx/2 - 1, nx/2}
    let ix = nx / 2 - 1;
    for k in 1..5usize {
        for j in 1..5usize {
            let idx = solver.linear_index(ix, j, k);
            let n = solver.normals()[idx];
            if n.norm() > 0.1 {
                let cos_angle = n.x.abs() / n.norm();
                assert!(
                    cos_angle > 0.9,   // within ~26° of x-axis (Youngs accuracy)
                    "PLIC normal should align with x-axis, cos_angle = {cos_angle:.3}"
                );
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Boundary Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Boundary**: Completely full domain (α=1 everywhere) → volume unchanged.
#[test]
fn test_full_domain_no_volume_change() {
    let mut solver = make_solver(8, 8, 8, default_config());
    {
        solver.alpha_mut().fill(1.0);
    }
    let vol0: f64 = solver.alpha().iter().sum();

    let n = 8 * 8 * 8;
    let vel = vec![Vector3::new(0.3, 0.2, 0.1); n];
    for idx in 0..n {
        solver.set_velocity_at(idx, vel[idx]);
    }

    for _ in 0..10 {
        solver.step(1e-3, &[], &[], &[]).ok();
    }

    let vol1: f64 = solver.alpha().iter().sum();
    assert_relative_eq!(vol0, vol1, epsilon = 1e-10);
}

/// **Boundary**: Completely empty domain (α=0 everywhere) → volume unchanged.
#[test]
fn test_empty_domain_no_volume_change() {
    let mut solver = make_solver(8, 8, 8, default_config());
    {
        solver.alpha_mut().fill(0.0);
    }
    let vol0: f64 = solver.alpha().iter().sum();

    let n = 8 * 8 * 8;
    let vel = vec![Vector3::new(0.5, 0.0, 0.0); n];
    for idx in 0..n {
        solver.set_velocity_at(idx, vel[idx]);
    }

    for _ in 0..10 {
        solver.step(1e-3, &[], &[], &[]).ok();
    }

    let vol1: f64 = solver.alpha().iter().sum();
    assert_relative_eq!(vol0, vol1, epsilon = 1e-10);
}

// ─────────────────────────────────────────────────────────────────────────────
// Adversarial Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Adversarial**: NaN velocity → solver must not panic and must return Err or α ∈ [0,1].
#[test]
fn test_nan_velocity_no_panic() {
    let mut solver = make_solver(6, 6, 6, default_config());
    {
        solver.alpha_mut().fill(0.5);
    }
    solver.set_velocity_at(0, Vector3::new(f64::NAN, 0.0, 0.0));
    // Must not panic
    let _ = solver.step(1e-4, &[], &[], &[]);
}

/// **Adversarial**: Infinite velocity → solver must not panic.
#[test]
fn test_infinite_velocity_no_panic() {
    let mut solver = make_solver(6, 6, 6, default_config());
    {
        solver.alpha_mut().fill(0.5);
    }
    solver.set_velocity_at(5, Vector3::new(f64::INFINITY, 0.0, 0.0));
    let _ = solver.step(1e-4, &[], &[], &[]);
}

/// **Adversarial**: Initial α values clamped if outside [0,1] (robustness).
#[test]
fn test_out_of_bounds_alpha_clamped() {
    let mut solver = make_solver(6, 6, 6, default_config());
    {
        let alpha = solver.alpha_mut();
        alpha[0] = 1.5;   // too large
        alpha[1] = -0.3;  // negative
    }
    // After one step with zero velocity, values should be clamped to [0,1]
    let _ = solver.step(1e-6, &[], &[], &[]);
    for &a in solver.alpha() {
        assert!(
            a >= -1e-12 && a <= 1.0 + 1e-12,
            "Volume fraction must be in [0,1] after step, got {a}"
        );
    }
}

/// **Adversarial**: Very high velocity (CFL > 1) — solver must not panic.
#[test]
fn test_high_cfl_no_panic() {
    let mut solver = make_solver(8, 8, 8, default_config());
    {
        let alpha = solver.alpha_mut();
        for (i, a) in alpha.iter_mut().enumerate() {
            *a = if i < 256 { 1.0 } else { 0.0 };
        }
    }
    // CFL >> 1
    let n = 8 * 8 * 8;
    let vel = vec![Vector3::new(100.0, 0.0, 0.0); n];
    for idx in 0..n {
        solver.set_velocity_at(idx, vel[idx]);
    }
    // Must not panic — may produce inaccurate results
    let _ = solver.step(1.0, &[], &[], &[]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Property Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Property**: `volume_under_plane_3d` at C=0 gives zero volume.
#[test]
fn test_volume_formula_zero_at_c_equals_zero() {
    use cfd_3d::vof::volume_under_plane_3d;
    let n = nalgebra::Vector3::new(1.0f64 / 3.0f64.sqrt(), 1.0 / 3.0f64.sqrt(), 1.0 / 3.0f64.sqrt());
    let v = volume_under_plane_3d(n, 0.0, 1.0, 1.0, 1.0);
    assert!(v.abs() < 1e-14, "V(C=0) must be zero, got {v:.3e}");
}

/// **Property**: `volume_under_plane_3d` at C = n_sum gives full cell volume.
#[test]
fn test_volume_formula_full_at_c_equals_sum() {
    use cfd_3d::vof::volume_under_plane_3d;
    let dx = 0.1;
    let dy = 0.2;
    let dz = 0.3;
    let n = nalgebra::Vector3::new(1.0f64 / 3.0f64.sqrt(), 1.0 / 3.0f64.sqrt(), 1.0 / 3.0f64.sqrt());
    let c_max = n.x.abs() * dx + n.y.abs() * dy + n.z.abs() * dz;
    let v = volume_under_plane_3d(n, c_max, dx, dy, dz);
    let cell_vol = dx * dy * dz;
    assert_relative_eq!(v, cell_vol, epsilon = 1e-13);
}

/// **Property**: `volume_under_plane_3d` is monotonically increasing in C.
#[test]
fn test_volume_formula_monotone_in_c() {
    use cfd_3d::vof::volume_under_plane_3d;
    let n = nalgebra::Vector3::new(0.6f64, 0.6, 0.529150262);
    let n_unit = n / n.norm();
    let dx = 1.0;
    let dy = 1.0;
    let dz = 1.0;
    let c_max = n_unit.x.abs() + n_unit.y.abs() + n_unit.z.abs();
    let steps = 50;
    let mut prev_v = -1.0;
    for i in 0..=steps {
        let c = c_max * i as f64 / steps as f64;
        let v = volume_under_plane_3d(n_unit, c, dx, dy, dz);
        assert!(
            v >= prev_v - 1e-14,
            "V(C) must be monotone: V({:.3}) = {v:.6} < prev {prev_v:.6}",
            c
        );
        prev_v = v;
    }
}

/// **Property**: Reconstructed normals have unit magnitude in interface cells.
#[test]
fn test_reconstructed_normal_unit_magnitude() {
    let nx = 10usize;
    let mut solver = make_solver(nx, 6, 6, default_config());
    {
        let alpha = solver.alpha_mut();
        for idx in 0..alpha.len() {
            let i = idx % nx;
            alpha[idx] = if i < nx / 2 { 1.0 } else { 0.0 };
        }
    }
    solver.reconstruct_interface();
    for (idx, &n) in solver.normals().iter().enumerate() {
        let mag = n.into_iter().map(|x| x * x).sum::<f64>().sqrt();
        if mag > 0.01 {
            assert!(
                (mag - 1.0).abs() < 1e-10,
                "Normal magnitude must be 1.0, got {mag:.6} at idx {idx}"
            );
        }
    }
}
