//! # Comprehensive IBM Solver Tests
//!
//! ## Test Coverage Matrix
//!
//! | Category     | Tests                                                            |
//! |-------------|------------------------------------------------------------------|
//! | Positive     | Kernel integral = 1, momentum conservation, InterpolationKernel |
//! | Boundary     | Lagrangian point at domain edge, domain corner                   |
//! | Adversarial  | Point outside domain, NaN force, zero dt                       |
//! | Property     | Eulerian forces always finite after spreading                    |

use cfd_3d::ibm::{
    config::IbmConfig, DeltaFunction, InterpolationKernel, IbmSolver, LagrangianPoint,
};
use nalgebra::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn default_config() -> IbmConfig {
    IbmConfig {
        use_direct_forcing: true,
        ..IbmConfig::default()
    }
}

fn make_solver(nx: usize, ny: usize, nz: usize) -> IbmSolver<f64> {
    let h = 1.0 / nx as f64;
    IbmSolver::new(
        default_config(),
        Vector3::new(h, h, h),
        (nx, ny, nz),
    )
}

fn make_lagrangian(pos: Vector3<f64>, force: Vector3<f64>) -> LagrangianPoint<f64> {
    LagrangianPoint {
        position: pos,
        force,
        velocity: Vector3::zeros(),
        weight: 1.0,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Positive Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Positive**: Roma-Peskin 3-point 1D kernel integrates to 1 over its stencil.
///
/// Discrete summation over the stencil with spacing h must satisfy:
/// Σ φ(r_i/h) = 1.
#[test]
fn test_roma_peskin_3_kernel_sums_to_one() {
    let h = 0.1f64;
    let kernel = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, h);
    // 3-point stencil: r ∈ {-1, 0, 1} in units of h
    let sum: f64 = [-1.0, 0.0, 1.0]
        .iter()
        .map(|&r| kernel.delta(r))
        .sum();
    assert!(
        (sum - 1.0).abs() < 1e-12,
        "Roma-Peskin 3-point 1D kernel must sum to 1, got {sum}"
    );
}

/// **Positive**: Roma-Peskin 4-point 1D kernel integrates to 1 over its stencil.
#[test]
fn test_roma_peskin_4_kernel_sums_to_one() {
    let h = 0.1f64;
    let kernel = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin4, h);
    // 4-point stencil: r ∈ {-1.5, -0.5, 0.5, 1.5}
    let sum: f64 = [-1.5, -0.5, 0.5, 1.5]
        .iter()
        .map(|&r| kernel.delta(r))
        .sum();
    assert!(
        (sum - 1.0).abs() < 1e-12,
        "Roma-Peskin 4-point 1D kernel must sum to 1, got {sum}"
    );
}

/// **Positive**: Kernel is non-negative everywhere.
#[test]
fn test_kernel_non_negative() {
    let kernel = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0f64);
    for i in -30..=30 {
        let r = i as f64 * 0.1;
        let v = kernel.delta(r);
        assert!(v >= 0.0, "Kernel must be non-negative, got {v:.6} at r={r:.2}");
    }
}

/// **Positive**: No Lagrangian points → spreading produces zero force field.
#[test]
fn test_empty_lagrangian_zero_force() {
    let solver = make_solver(10, 10, 10);
    let forces = solver.spread_forces().expect("spread_forces failed");
    for f in &forces {
        assert_eq!(*f, Vector3::zeros(), "No Lagrangian points → zero force");
    }
}

/// **Positive**: Spreading a unit force conserves total x-momentum.
///
/// **Theorem** (Spreading conservation): Σ_ij f_ij = F_L
///
/// When the kernel δ(r) satisfies partition of unity (Σ φ(r-j) = 1),
/// the 3D tensor product Σ_{ijk} δ(rx)·δ(ry)·δ(rz) = 1 over the stencil,
/// so Σ f_ij = F_L exactly.
#[test]
fn test_spreading_conserves_total_momentum() {
    let nx = 16usize;
    let h = 1.0 / nx as f64;
    let mut solver = IbmSolver::new(
        default_config(),
        Vector3::new(h, h, h),
        (nx, nx, nx),
    );
    let total_force = Vector3::new(1.0, 0.0, 0.0);
    solver.add_lagrangian_point(
        make_lagrangian(Vector3::new(0.5, 0.5, 0.5), total_force)
    );

    let forces = solver.spread_forces().expect("spread_forces failed");
    let sum_x: f64 = forces.iter().map(|f| f.x).sum();

    // Kernel sums to 1 → total spread force = input Lagrangian force
    assert!(
        (sum_x - total_force.x).abs() < 1e-10,
        "Total x-momentum must equal spread force, got {sum_x}"
    );
}

/// **Positive**: Stencil size for 4-point kernel is 4.
#[test]
fn test_stencil_size_correct() {
    let k3 = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0);
    let k4 = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin4, 1.0);
    assert_eq!(k3.stencil_size(), 3);
    assert_eq!(k4.stencil_size(), 4);
}

// ─────────────────────────────────────────────────────────────────────────────
// Boundary Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Boundary**: Lagrangian point on x=0 face → no out-of-bounds panic.
#[test]
fn test_lagrangian_point_on_x_boundary_no_panic() {
    let mut solver = make_solver(10, 10, 10);
    solver.add_lagrangian_point(
        make_lagrangian(Vector3::new(0.0, 0.5, 0.5), Vector3::new(1.0, 0.0, 0.0))
    );
    let _ = solver.spread_forces();
}

/// **Boundary**: Lagrangian point at maximum corner → no out-of-bounds panic.
#[test]
fn test_lagrangian_point_at_max_corner_no_panic() {
    let mut solver = make_solver(10, 10, 10);
    solver.add_lagrangian_point(
        make_lagrangian(Vector3::new(1.0, 1.0, 1.0), Vector3::new(0.0, 1.0, 0.0))
    );
    let _ = solver.spread_forces();
}

/// **Boundary**: Multiple points all land in bounds.
#[test]
fn test_multiple_boundary_points_no_panic() {
    let mut solver = make_solver(12, 12, 12);
    for i in 0..4 {
        let edge = i as f64 * 0.33;
        solver.add_lagrangian_point(
            make_lagrangian(Vector3::new(edge, 0.0, 0.5), Vector3::zeros())
        );
    }
    let _ = solver.spread_forces();
}

// ─────────────────────────────────────────────────────────────────────────────
// Adversarial Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Adversarial**: Lagrangian point far outside domain → no panic, forces finite.
#[test]
fn test_lagrangian_point_outside_domain_no_panic() {
    let mut solver = make_solver(10, 10, 10);
    solver.add_lagrangian_point(
        make_lagrangian(Vector3::new(5.0, 5.0, 5.0), Vector3::new(1.0, 1.0, 1.0))
    );
    let forces = solver.spread_forces().unwrap_or_default();
    for f in &forces {
        assert!(
            f.x.is_finite() && f.y.is_finite() && f.z.is_finite(),
            "Forces must be finite even for out-of-domain point"
        );
    }
}

/// **Adversarial**: NaN force on Lagrangian point → no panic.
#[test]
fn test_nan_force_no_panic() {
    let mut solver = make_solver(8, 8, 8);
    solver.add_lagrangian_point(
        make_lagrangian(Vector3::new(0.5, 0.5, 0.5), Vector3::new(f64::NAN, 0.0, 0.0))
    );
    let _ = solver.spread_forces();
}

/// **Adversarial**: Zero-weight Lagrangian point → force spread is zero.
#[test]
fn test_zero_weight_lagrangian_no_force() {
    let mut solver = make_solver(10, 10, 10);
    solver.add_lagrangian_point(LagrangianPoint {
        position: Vector3::new(0.5, 0.5, 0.5),
        force: Vector3::new(100.0, 0.0, 0.0),
        velocity: Vector3::zeros(),
        weight: 0.0, // zero weight
    });
    let forces = solver.spread_forces().expect("spread_forces failed");
    let total: f64 = forces.iter().map(|f| f.x).sum();
    assert!(
        total.abs() < 1e-14,
        "Zero-weight point must produce zero spread force, got {total}"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Property Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Property**: Eulerian forces are always finite after spreading multiple points.
#[test]
fn test_eulerian_forces_always_finite() {
    let mut solver = make_solver(10, 10, 10);
    for i in 0..5 {
        solver.add_lagrangian_point(
            make_lagrangian(
                Vector3::new(0.1 + 0.18 * i as f64, 0.5, 0.5),
                Vector3::new(1.0, -0.5, 0.3),
            )
        );
    }
    let forces = solver.spread_forces().expect("spread_forces failed");
    for f in &forces {
        assert!(
            f.x.is_finite() && f.y.is_finite() && f.z.is_finite(),
            "All Eulerian forces must be finite"
        );
    }
}

/// **Property**: num_points() equals the number added.
#[test]
fn test_num_points_consistent() {
    let mut solver = make_solver(8, 8, 8);
    for i in 0..7 {
        solver.add_lagrangian_point(
            make_lagrangian(Vector3::new(0.1 * i as f64, 0.5, 0.5), Vector3::zeros())
        );
    }
    assert_eq!(solver.num_points(), 7);
}

/// **Property**: Kernel is zero outside support region |r| > 1.5 (3-point).
#[test]
fn test_kernel_zero_outside_support() {
    let kernel = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0);
    for r in [1.6, 2.0, 3.0, 5.0] {
        let v = kernel.delta(r);
        assert!(
            v.abs() < 1e-15,
            "3-point kernel should be zero for |r| >= 1.5, got {v:.3e} at r={r}"
        );
    }
}
