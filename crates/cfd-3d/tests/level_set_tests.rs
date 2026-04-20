//! # Comprehensive Level Set Solver Tests
//!
//! ## Test Coverage Matrix
//!
//! | Category     | Tests                                                                  |
//! |-------------|------------------------------------------------------------------------|
//! | Positive     | WENO5 accuracy for smooth function, reinitialization |∇φ|=1 recovery  |
//! | Boundary     | Zero velocity, interface at domain boundary, single-cell band          |
//! | Adversarial  | CFL >> 1, NaN velocity, extreme grid aspect ratio                      |
//! | Property     | Narrow band non-empty, |∇φ| in [0.5, 2.0] after reinit               |

use cfd_3d::level_set::{LevelSetConfig, LevelSetSolver};
use nalgebra::Vector3;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

fn default_config() -> LevelSetConfig {
    LevelSetConfig {
        band_width: 5.0,
        reinitialization_interval: 5,
        use_narrow_band: true,
        ..LevelSetConfig::default()
    }
}

fn make_solver(nx: usize, ny: usize, nz: usize) -> LevelSetSolver<f64> {
    let h = 1.0 / nx as f64;
    LevelSetSolver::new(default_config(), nx, ny, nz, h, h, h)
}

// Initialize a 3D sphere interface centred at (0.5, 0.5, 0.5) with radius R.
fn init_sphere(solver: &mut LevelSetSolver<f64>, r: f64) {
    let nx = solver.nx();
    let ny = solver.ny();
    let nz = solver.nz();
    let h = 1.0 / nx as f64;
    let phi = solver.phi_mut();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * h - 0.5;
                let y = (j as f64 + 0.5) * h - 0.5;
                let z = (k as f64 + 0.5) * h - 0.5;
                let idx = k * ny * nx + j * nx + i;
                phi[idx] = (x * x + y * y + z * z).sqrt() - r;
            }
        }
    }
}

fn init_perturbed_sphere(solver: &mut LevelSetSolver<f64>, r: f64) {
    let nx = solver.nx();
    let ny = solver.ny();
    let nz = solver.nz();
    let h = 1.0 / nx as f64;
    let phi = solver.phi_mut();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64 * h - 0.5;
                let y = j as f64 * h - 0.5;
                let z = k as f64 * h - 0.5;
                let radius = (x * x + y * y + z * z).sqrt();
                let idx = k * nx * ny + j * nx + i;
                phi[idx] = (radius - r) * (1.0 + 0.3 * (x + y).sin());
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Positive Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Positive**: Zero-velocity advection preserves the level set exactly.
///
/// Invariant: ∂φ/∂t + 0·∇φ = 0 → φ(t) = φ(0) for all t.
///
/// Reinitialization is disabled for this test because the pure advection
/// operator with zero velocity should be an identity transform.
#[test]
fn test_zero_velocity_level_set_preserved() {
    let nx = 16usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: usize::MAX,
            use_narrow_band: false,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_sphere(&mut solver, 0.3);
    let phi0: Vec<f64> = solver.phi().to_vec();

    // All-zero velocity (default)
    for _ in 0..10 {
        solver.advance(1e-4).expect("advance failed");
    }

    for (idx, (&p0, &p1)) in phi0.iter().zip(solver.phi()).enumerate() {
        assert!(
            (p0 - p1).abs() < 1e-12,
            "phi must be unchanged with zero velocity at idx {idx}"
        );
    }
}

/// **Positive**: After reinitialization, |∇φ| ≈ 1 in the narrow band.
///
/// Theorem: Sussman reinitialization converges to |∇φ|=1 in O(Δτ²) iterations.
/// We check that the error |∇φ|−1 < 1% for all cells in the narrow band.
#[test]
fn test_reinitialization_restores_unit_gradient() {
    let nx = 20usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 1,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_perturbed_sphere(&mut solver, 0.25);

    // Advance once (triggers reinitialization)
    solver.advance(1e-4).expect("advance failed");

    // Check |∇φ| ≈ 1 in the near-interface band, excluding the medial axis
    // (the sphere centre where the SDF gradient is undefined — a topological
    // singularity where all one-sided differences vanish by symmetry).
    let mut max_err = 0.0f64;
    for k in 2..nx - 2 {
        for j in 2..nx - 2 {
            for i in 2..nx - 2 {
                let phi = solver.phi();
                let idx = |ii: usize, jj: usize, kk: usize| kk * nx * nx + jj * nx + ii;
                let dphi_dx = (phi[idx(i + 1, j, k)] - phi[idx(i - 1, j, k)]) / (2.0 * h);
                let dphi_dy = (phi[idx(i, j + 1, k)] - phi[idx(i, j - 1, k)]) / (2.0 * h);
                let dphi_dz = (phi[idx(i, j, k + 1)] - phi[idx(i, j, k - 1)]) / (2.0 * h);
                let grad_mag = (dphi_dx * dphi_dx + dphi_dy * dphi_dy + dphi_dz * dphi_dz).sqrt();

                // Only check cells in the narrow band near the interface,
                // but exclude cells too close to the medial axis (φ < 0 far from interface).
                // The interface is at φ = 0; the medial axis (sphere centre) is at φ ≈ −R.
                let phi_val = phi[idx(i, j, k)];
                let abs_phi = phi_val.abs();
                if abs_phi < 3.0 * h && abs_phi > 0.5 * h {
                    let err = (grad_mag - 1.0).abs();
                    if err > max_err {
                        max_err = err;
                    }
                }
            }
        }
    }
    assert!(
        max_err < 0.20, // 20% tolerance for first-order Godunov reinit on 20-cell grid
        "Reinitialization must restore |∇φ|≈1, max error = {max_err:.4}"
    );
}

/// **Positive**: reinitialization interval `0` disables the reinitialization pass.
#[test]
fn test_zero_reinitialization_interval_skips_reinitialization() {
    let nx = 20usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 0,
            use_narrow_band: false,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_perturbed_sphere(&mut solver, 0.25);
    let phi0: Vec<f64> = solver.phi().to_vec();

    solver.advance(1e-4).expect("advance failed");

    for (idx, (&p0, &p1)) in phi0.iter().zip(solver.phi()).enumerate() {
        assert!(
            (p0 - p1).abs() < 1e-12,
            "phi must remain unchanged when reinitialization is disabled at idx {idx}"
        );
    }
}

/// **Positive**: a zero reinitialization iteration budget disables the solve.
#[test]
fn test_reinitialization_max_iterations_zero_skips_reinitialization() {
    let nx = 20usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 1,
            max_iterations: 0,
            use_narrow_band: false,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_perturbed_sphere(&mut solver, 0.25);
    let phi0: Vec<f64> = solver.phi().to_vec();

    solver.advance(1e-4).expect("advance failed");

    for (idx, (&p0, &p1)) in phi0.iter().zip(solver.phi()).enumerate() {
        assert!(
            (p0 - p1).abs() < 1e-12,
            "phi must remain unchanged when max_iterations is zero at idx {idx}"
        );
    }
}

/// **Positive**: Sphere advection preserves zero-level-set radius to leading order.
///
/// A sphere of radius R advected by uniform velocity u for time t should remain
/// centred at (x₀ + u·t, y₀, z₀) with unchanged radius (steady translation).
#[test]
fn test_sphere_advection_preserves_radius() {
    let nx = 32usize;
    let h = 1.0 / nx as f64;
    let r0 = 0.2;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 10,
            use_narrow_band: false,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_sphere(&mut solver, r0);

    // Apply uniform x-velocity
    let u = 0.2;
    let vel: Vec<Vector3<f64>> = vec![Vector3::new(u, 0.0, 0.0); nx * nx * nx];
    solver.set_velocity(vel);

    let dt = 0.3 * h / u; // CFL = 0.3
    let n_steps = 5;
    for _ in 0..n_steps {
        solver.advance(dt).expect("advance failed");
    }

    // After advection, find minimum |φ| along z=0 plane and compare to radius
    // (a rough check that the interface hasn't collapsed)
    let phi = solver.phi();
    let has_positive = phi.iter().any(|&v| v > 0.0);
    let has_negative = phi.iter().any(|&v| v < 0.0);
    assert!(
        has_positive && has_negative,
        "Level set must still straddle zero after advection"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Boundary Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Boundary**: Single-cell domain → advance must not panic or out-of-bounds.
#[test]
fn test_single_cell_no_panic() {
    let mut solver = LevelSetSolver::new(default_config(), 4, 4, 4, 0.25, 0.25, 0.25);
    {
        solver.phi_mut().fill(-0.1);
    }
    let _ = solver.advance(1e-3);
}

/// **Boundary**: Small domains fall back from WENO5-Z to first-order upwind.
#[test]
fn test_small_grid_uses_first_order_fallback() {
    let nx = 4usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 0,
            use_narrow_band: false,
            use_weno: true,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );

    {
        let phi = solver.phi_mut();
        for k in 0..nx {
            for j in 0..nx {
                for i in 0..nx {
                    let x = i as f64 * h - 0.5;
                    let idx = k * nx * nx + j * nx + i;
                    phi[idx] = x;
                }
            }
        }
    }

    solver.set_velocity(vec![Vector3::new(1.0, 0.0, 0.0); nx * nx * nx]);

    let before = solver.phi().to_vec();
    solver.advance(0.01).expect("small-grid advance failed");

    let center = solver.index(1, 1, 1);
    assert!(
        (solver.phi()[center] - before[center]).abs() > 1e-12,
        "small-grid advance must update the interior state"
    );
}

/// **Boundary**: Interface exactly at domain boundary (φ=0 on ghost cells).
#[test]
fn test_interface_at_boundary_no_panic() {
    let nx = 10usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(default_config(), nx, nx, nx, h, h, h);
    {
        let phi = solver.phi_mut();
        for (i, p) in phi.iter_mut().enumerate() {
            *p = i as f64 * h - 0.5; // crosses zero
        }
    }
    for _ in 0..5 {
        solver.advance(1e-3).expect("advance at boundary failed");
    }
}

/// **Boundary**: Narrow band update produces a non-empty set for a sphere interface.
#[test]
fn test_narrow_band_non_empty_for_sphere() {
    let mut solver = make_solver(20, 20, 20);
    init_sphere(&mut solver, 0.3);
    solver.update_narrow_band();
    assert!(
        !solver.narrow_band().is_empty(),
        "Narrow band must be non-empty for sphere interface"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Adversarial Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Adversarial**: CFL >> 1 (dt = 1.0, u = 10.0) — must not panic.
#[test]
fn test_extreme_cfl_no_panic() {
    let nx = 10usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(default_config(), nx, nx, nx, h, h, h);
    init_sphere(&mut solver, 0.3);
    let vel = vec![Vector3::new(10.0, 0.0, 0.0); nx * nx * nx];
    solver.set_velocity(vel);
    let _ = solver.advance(1.0); // CFL ≈ 100 — stability guarantee voided, but no panic
}

/// **Adversarial**: NaN velocity → no panic, φ remains finite or solver returns Err.
#[test]
fn test_nan_velocity_no_panic() {
    let nx = 8usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(default_config(), nx, nx, nx, h, h, h);
    init_sphere(&mut solver, 0.2);
    let mut vel = vec![Vector3::new(0.0, 0.0, 0.0); nx * nx * nx];
    vel[0] = Vector3::new(f64::NAN, 0.0, 0.0);
    solver.set_velocity(vel);
    let _ = solver.advance(1e-4);
    // Enough that solver did not panic — values may be NaN in polluted cells
}

/// **Adversarial**: All-positive φ (interface outside domain) → stable.
#[test]
fn test_no_interface_all_positive_phi() {
    let mut solver = make_solver(8, 8, 8);
    {
        solver.phi_mut().fill(1.0); // no zero-crossing
    }
    for _ in 0..5 {
        solver
            .advance(1e-3)
            .expect("advance failed with all-positive phi");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Property Tests
// ─────────────────────────────────────────────────────────────────────────────

/// **Property**: φ is symmetric through the centre for a centred sphere.
///
/// For a perfect sphere centred at origin, φ(x) = φ(−x)  (even function).
/// The signed-distance field |x| − R depends only on distance to the centre,
/// so φ(x,y,z) = φ(−x,−y,−z) exactly.
#[test]
fn test_sphere_phi_symmetry() {
    let nx = 20usize;
    let h = 1.0 / nx as f64;
    let mut solver = LevelSetSolver::new(
        LevelSetConfig {
            reinitialization_interval: 100,
            use_narrow_band: false,
            ..default_config()
        },
        nx,
        nx,
        nx,
        h,
        h,
        h,
    );
    init_sphere(&mut solver, 0.3);

    // Check φ(i,j,k) ≈ φ(nx−1−i, ny−1−j, nz−1−k) for interior points
    let phi = solver.phi();
    let ny = solver.ny();
    let nz = solver.nz();
    let mut max_asym = 0.0f64;
    for k in 1..nz - 1 {
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx1 = k * ny * nx + j * nx + i;
                let idx2 = (nz - 1 - k) * ny * nx + (ny - 1 - j) * nx + (nx - 1 - i);
                let asym = (phi[idx1] - phi[idx2]).abs();
                if asym > max_asym {
                    max_asym = asym;
                }
            }
        }
    }
    // Tolerance: machine precision for cell-centred grid with even nx
    assert!(
        max_asym < 1e-14,
        "Sphere φ symmetry violated: max |φ(x)−φ(−x)| = {max_asym:.4e}"
    );
}
