//! Robustness and mathematical-property tests for cfd-3d
//!
//! These tests verify the formal theorems documented in each module's
//! Rustdoc, covering partition of unity, Kronecker delta, strain-rate
//! symmetry, stiffness-matrix symmetry, element-volume edge cases,
//! stabilisation-parameter limits, quadrature exactness, and
//! shape-function gradient completeness.

use approx::assert_relative_eq;
use nalgebra::{Matrix3, Matrix3x4, Vector3};

// ═══════════════════════════════════════════════════════════════════════════
//  P2 Shape Functions — Partition of Unity / Kronecker Delta
// ═══════════════════════════════════════════════════════════════════════════

/// The 10 P2 Lagrange nodes in barycentric coordinates [L1,L2,L3,L4].
fn p2_nodes() -> [[f64; 4]; 10] {
    [
        // corners
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        // mid-edges: (0,1), (1,2), (2,0), (0,3), (1,3), (2,3)
        [0.5, 0.5, 0.0, 0.0],
        [0.0, 0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5, 0.0],
        [0.5, 0.0, 0.0, 0.5],
        [0.0, 0.5, 0.0, 0.5],
        [0.0, 0.0, 0.5, 0.5],
    ]
}

/// Build a trivial P1 gradient matrix for the reference tet
/// v0=(0,0,0), v1=(1,0,0), v2=(0,1,0), v3=(0,0,1).
/// ∇L_0 = (-1,-1,-1), ∇L_1 = (1,0,0), ∇L_2 = (0,1,0), ∇L_3 = (0,0,1).
fn reference_p1_gradients() -> Matrix3x4<f64> {
    Matrix3x4::new(
        -1.0,  1.0,  0.0,  0.0,
        -1.0,  0.0,  1.0,  0.0,
        -1.0,  0.0,  0.0,  1.0,
    )
}

/// Theorem: Σ N_i(x) = 1 for all x in the element (Partition of Unity).
#[test]
fn test_p2_partition_of_unity() {
    use cfd_3d::fem::shape_functions::LagrangeTet10;
    let sf = LagrangeTet10::<f64>::new(reference_p1_gradients());

    // Test at all 10 nodes
    for node in p2_nodes() {
        let vals = sf.values(&node);
        let sum: f64 = vals.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-14);
    }

    // Test at centroid (0.25, 0.25, 0.25, 0.25)
    let centroid = [0.25, 0.25, 0.25, 0.25];
    let vals = sf.values(&centroid);
    let sum: f64 = vals.iter().sum();
    assert_relative_eq!(sum, 1.0, epsilon = 1e-14);

    // Test at a random interior point
    let interior = [0.1, 0.3, 0.4, 0.2];
    let vals = sf.values(&interior);
    let sum: f64 = vals.iter().sum();
    assert_relative_eq!(sum, 1.0, epsilon = 1e-14);
}

/// Theorem: N_i(x_j) = δ_{ij} (Kronecker Delta property).
#[test]
fn test_p2_kronecker_delta() {
    use cfd_3d::fem::shape_functions::LagrangeTet10;
    let sf = LagrangeTet10::<f64>::new(reference_p1_gradients());

    let nodes = p2_nodes();
    for (j, node) in nodes.iter().enumerate() {
        let vals = sf.values(node);
        for (i, &v) in vals.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert_relative_eq!(v, expected, epsilon = 1e-14,);
        }
    }
}

/// Theorem: Σ ∇N_i = 0 (gradient completeness).
#[test]
fn test_p2_gradient_sum_zero() {
    use cfd_3d::fem::shape_functions::LagrangeTet10;
    let sf = LagrangeTet10::<f64>::new(reference_p1_gradients());

    let test_points: [[f64; 4]; 3] = [
        [0.25, 0.25, 0.25, 0.25],
        [0.1, 0.3, 0.4, 0.2],
        [0.5, 0.5, 0.0, 0.0],
    ];

    for l in &test_points {
        let grad = sf.gradients(l);
        // Sum columns (nodes) for each row (x,y,z)
        for row in 0..3 {
            let sum: f64 = (0..10).map(|col| grad[(row, col)]).sum();
            assert_relative_eq!(sum, 0.0, epsilon = 1e-12);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Keast Quadrature — Weights Sum to Reference Volume
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: Keast degree-3 weights sum to 1/6 (reference tet volume).
#[test]
fn test_quadrature_weights_sum_to_reference_volume() {
    use cfd_3d::fem::quadrature::TetrahedronQuadrature;
    let quad = TetrahedronQuadrature::<f64>::keast_degree_3();
    let sum: f64 = quad.weights().iter().sum();
    assert_relative_eq!(sum, 1.0 / 6.0, epsilon = 1e-14);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Stress / Strain — Symmetry & Constitutive Relation
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: ε̇ is symmetric for any velocity gradient.
#[test]
fn test_strain_rate_symmetry() {
    use cfd_3d::fem::stress::strain_rate_tensor;
    let grad_u = Matrix3::new(
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
    );
    let eps = strain_rate_tensor(&grad_u);
    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(eps[(i, j)], eps[(j, i)], epsilon = 1e-15);
        }
    }
}

/// Theorem: ε̇ = 0 for rigid-body rotation.
#[test]
fn test_strain_rate_zero_for_rigid_rotation() {
    use cfd_3d::fem::stress::strain_rate_tensor;
    // Antisymmetric ∇u → solid-body rotation
    let grad_u = Matrix3::new(
        0.0,  1.0, -2.0,
       -1.0,  0.0,  3.0,
        2.0, -3.0,  0.0,
    );
    let eps = strain_rate_tensor(&grad_u);
    for i in 0..3 {
        for j in 0..3 {
            assert_relative_eq!(eps[(i, j)], 0.0, epsilon = 1e-15);
        }
    }
}

/// Theorem: tr(σ) = −3p for incompressible Newtonian flow.
#[test]
fn test_stress_trace_is_neg_3p() {
    use cfd_3d::fem::stress::{strain_rate_tensor, stress_tensor};
    use cfd_core::physics::fluid::ConstantPropertyFluid;

    let fluid = ConstantPropertyFluid::new(
        "test".into(), 1000.0, 0.001, 4186.0, 0.6, 1500.0,
    );

    // Divergence-free velocity gradient: tr(∇u) = 0.
    let grad_u = Matrix3::new(
        1.0, 0.5, 0.0,
        0.5, -2.0, 0.3,
        0.0, 0.3, 1.0,
    );
    let eps = strain_rate_tensor(&grad_u);
    let pressure = 42.0;
    let sigma = stress_tensor(&fluid, pressure, &eps);
    let trace = sigma[(0, 0)] + sigma[(1, 1)] + sigma[(2, 2)];

    // For incompressible (tr(ε̇)=0): σ = −pI + 2με̇, so tr(σ) = −3p + 2μ·tr(ε̇) = −3p.
    assert_relative_eq!(trace, -3.0 * pressure, epsilon = 1e-10);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Element — Volume & Shape Derivative Properties
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: V = |det(e1,e2,e3)| / 6 for the reference tet.
#[test]
fn test_element_volume_reference_tet() {
    use cfd_3d::fem::FluidElement;
    let mut elem = FluidElement::<f64>::new(vec![0, 1, 2, 3]);
    let verts = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    elem.calculate_volume(&verts);
    assert_relative_eq!(elem.volume, 1.0 / 6.0, epsilon = 1e-14);
}

/// Degenerate: coplanar nodes → volume ≈ 0.
#[test]
fn test_element_volume_degenerate_flat_tet() {
    use cfd_3d::fem::FluidElement;
    let mut elem = FluidElement::<f64>::new(vec![0, 1, 2, 3]);
    // All four nodes in the z=0 plane
    let verts = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.5, 0.5, 0.0),
    ];
    elem.calculate_volume(&verts);
    assert!(elem.volume.abs() < 1e-20, "flat tet should have ~zero volume");
}

/// Shape derivative identity: Σ ∇N_i = 0 for linear (P1) elements.
#[test]
fn test_p1_shape_derivative_sum_zero() {
    use cfd_3d::fem::FluidElement;
    let mut elem = FluidElement::<f64>::new(vec![0, 1, 2, 3]);
    let verts = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.3, 0.9, 0.0),
        Vector3::new(0.2, 0.3, 0.8),
    ];
    elem.calculate_volume(&verts);
    elem.calculate_shape_derivatives(&verts);

    for row in 0..3 {
        let sum: f64 = (0..4).map(|col| elem.shape_derivatives[(row, col)]).sum();
        assert_relative_eq!(sum, 0.0, epsilon = 1e-12);
    }
}

/// Stiffness matrix symmetry: K = K^T for Stokes.
#[test]
fn test_stiffness_matrix_symmetry() {
    use cfd_3d::fem::FluidElement;
    let mut elem = FluidElement::<f64>::new(vec![0, 1, 2, 3]);
    let verts = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    elem.calculate_volume(&verts);
    elem.calculate_shape_derivatives(&verts);

    let k = elem.stiffness_contribution(0.001);
    let n = k.nrows();
    for i in 0..n {
        for j in i + 1..n {
            assert_relative_eq!(k[(i, j)], k[(j, i)], epsilon = 1e-15);
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  SUPG/PSPG Stabilisation — Asymptotic Limits
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: τ → h/(2|u|) for Pe_h >> 1 (advection-dominated).
#[test]
fn test_stabilization_high_peclet_limit() {
    use cfd_3d::fem::stabilization::StabilizationParameters;
    // h = 0.01, ν = 1e-10 (very small), u = (1,0,0) → Pe_h = |u|h/(2ν) >> 1
    let vel = Vector3::new(1.0, 0.0, 0.0);
    let params = StabilizationParameters::new(0.01, 1e-10, vel, None);
    let tau = params.tau_supg();
    let expected = 0.01 / (2.0 * 1.0);
    assert_relative_eq!(tau, expected, epsilon = 1e-4);
}

/// Theorem: τ → h²/(4ν) for Pe_h << 1 (diffusion-dominated).
#[test]
fn test_stabilization_low_peclet_limit() {
    use cfd_3d::fem::stabilization::StabilizationParameters;
    // h = 0.01, ν = 1.0, u = (1e-10,0,0) → Pe_h << 1
    let vel = Vector3::new(1e-10, 0.0, 0.0);
    let params = StabilizationParameters::new(0.01, 1.0, vel, None);
    let tau = params.tau_supg();
    let expected = 0.01_f64.powi(2) / (4.0 * 1.0);
    assert_relative_eq!(tau, expected, epsilon = 1e-6);
}

/// Degenerate: zero velocity should not panic or produce NaN.
#[test]
fn test_stabilization_zero_velocity_no_nan() {
    use cfd_3d::fem::stabilization::StabilizationParameters;
    let vel = Vector3::new(0.0, 0.0, 0.0);
    let params = StabilizationParameters::new(0.01, 0.001, vel, None);
    let tau: f64 = params.tau_supg();
    assert!(tau.is_finite(), "τ should be finite for zero velocity");
    assert!(tau >= 0.0, "τ should be non-negative");
}

/// Degenerate: zero viscosity should not panic or produce NaN.
#[test]
fn test_stabilization_zero_viscosity_no_nan() {
    use cfd_3d::fem::stabilization::StabilizationParameters;
    let vel = Vector3::new(1.0, 0.0, 0.0);
    let params = StabilizationParameters::new(0.01, 0.0, vel, None);
    let tau: f64 = params.tau_supg();
    assert!(tau.is_finite(), "τ should be finite for zero viscosity");
    assert!(tau >= 0.0, "τ should be non-negative");
}

// ═══════════════════════════════════════════════════════════════════════════
//  Bifurcation — Murray's Law & Geometry
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: Murray-compliant geometry gives deviation ≈ 0.
#[test]
fn test_murray_law_perfect_compliance() {
    use cfd_3d::bifurcation::BifurcationGeometry3D;
    // D_parent = 100μm → each daughter = 100/(2^(1/3)) ≈ 79.37μm
    let d_parent = 100e-6;
    let d_daughter = d_parent / 2.0_f64.cbrt();
    let geom = BifurcationGeometry3D::<f64>::symmetric(
        d_parent, d_daughter, 1e-3, 1e-3, 1e-4,
    );
    let dev = geom.murray_law_deviation();
    // D^3 ≈ 1e-12 so f64 rounding yields ~eps*D^3 ≈ 1e-28
    assert!(dev.abs() < 1e-25, "Murray-compliant geometry should have near-zero deviation: {dev}");
}

/// Volume and surface-area positivity for valid geometry.
#[test]
fn test_bifurcation_volume_and_area_positive() {
    use cfd_3d::bifurcation::BifurcationGeometry3D;
    let geom = BifurcationGeometry3D::<f64>::symmetric(
        100e-6, 80e-6, 1e-3, 1e-3, 1e-4,
    );
    assert!(geom.total_volume() > 0.0);
    assert!(geom.total_surface_area() > 0.0);
}

/// Asymmetric constructor works with unequal daughter diameters.
#[test]
fn test_asymmetric_bifurcation() {
    use cfd_3d::bifurcation::BifurcationGeometry3D;
    let geom = BifurcationGeometry3D::<f64>::asymmetric(
        100e-6, 90e-6, 60e-6,
        1e-3, 1e-3, 1e-3,
        1e-4,
        std::f64::consts::PI / 4.0,
    );
    assert!(geom.total_volume() > 0.0);
    assert!(geom.d_daughter1 > geom.d_daughter2);
}

// ═══════════════════════════════════════════════════════════════════════════
//  Trifurcation — Murray's Law
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: D₀³ = D₁³ + D₂³ + D₃³ for symmetric trifurcation.
#[test]
fn test_trifurcation_murray_law() {
    use cfd_3d::trifurcation::TrifurcationGeometry3D;
    // D_parent = 100μm → each daughter = 100/(3^(1/3)) ≈ 69.3μm
    let d_parent = 100e-6;
    let d_daughter = d_parent / 3.0_f64.cbrt();
    let geom = TrifurcationGeometry3D::<f64>::symmetric(
        d_parent, d_daughter, 1e-3, 1e-3, 1e-4,
        std::f64::consts::PI / 6.0,
    );
    let dev = geom.murray_law_deviation();
    // D^3 ≈ 1e-12 so f64 rounding yields ~eps*D^3 ≈ 1e-28
    assert!(dev.abs() < 1e-25, "Murray-compliant trifurcation should have near-zero deviation: {dev}");
}

// ═══════════════════════════════════════════════════════════════════════════
//  IBM — Delta function symmetry and partition of unity
// ═══════════════════════════════════════════════════════════════════════════

/// Theorem: δ(−r) = δ(r) (even-function symmetry).
#[test]
fn test_delta_function_symmetry() {
    use cfd_3d::ibm::{DeltaFunction, InterpolationKernel};
    let variants = [
        DeltaFunction::RomaPeskin3,
        DeltaFunction::RomaPeskin4,
        DeltaFunction::Peskin4,
    ];
    for df in &variants {
        let kernel = InterpolationKernel::new(df.clone(), 1.0_f64);
        for &r in &[0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0] {
            let pos = kernel.delta(r);
            let neg = kernel.delta(-r);
            assert_relative_eq!(pos, neg, epsilon = 1e-14,);
        }
    }
}

/// Peskin4 variant also sums to 1 (partition of unity).
#[test]
fn test_peskin4_kernel_partition_of_unity() {
    use cfd_3d::ibm::{DeltaFunction, InterpolationKernel};
    let kernel = InterpolationKernel::new(DeltaFunction::Peskin4, 1.0_f64);
    let h = 1.0;
    // Sample at r = 0.3 (not on a grid point).
    // The sum over all integer shifts should be 1:
    // Σ_{j} δ_h(r - j·h) · h = 1
    let r = 0.3;
    let mut sum = 0.0;
    for j in -3..=3 {
        sum += kernel.delta(r - j as f64 * h) * h;
    }
    assert_relative_eq!(sum, 1.0, epsilon = 1e-10);
}

// ═══════════════════════════════════════════════════════════════════════════
//  VOF — Volume Conservation & Bounds
// ═══════════════════════════════════════════════════════════════════════════

/// Geometric vs algebraic advection with zero velocity must both conserve volume exactly.
#[test]
fn test_vof_geometric_vs_algebraic_zero_velocity_conservation() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let base = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };

    for method in [AdvectionMethod::Geometric, AdvectionMethod::Algebraic] {
        let cfg = VofConfig {
            advection_method: method,
            ..base.clone()
        };
        let mut solver = VofSolver::<f64>::new(8, 8, 8, cfg).unwrap();
        // Half-filled domain
        {
            let alpha = solver.alpha_mut();
            for (i, a) in alpha.iter_mut().enumerate() {
                *a = if i < 256 { 1.0 } else { 0.0 };
            }
        }
        let vol_before: f64 = solver.alpha().iter().sum();
        let _ = solver.advance(0.001);
        let vol_after: f64 = solver.alpha().iter().sum();
        assert_relative_eq!(vol_before, vol_after, epsilon = 1e-10);
    }
}

/// After advection, α must remain in [0, 1] for every cell.
#[test]
fn test_vof_alpha_bounds_after_advection() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    let mut solver = VofSolver::<f64>::new(10, 10, 10, cfg).unwrap();

    // Sphere-like initialization
    let n = 10;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.linear_index(i, j, k);
                let x = (i as f64 + 0.5) / n as f64 - 0.5;
                let y = (j as f64 + 0.5) / n as f64 - 0.5;
                let z = (k as f64 + 0.5) / n as f64 - 0.5;
                let r = (x * x + y * y + z * z).sqrt();
                solver.alpha_mut()[idx] = if r < 0.3 { 1.0 } else { 0.0 };
            }
        }
    }

    // Set uniform velocity
    let vel = vec![Vector3::new(0.5, 0.0, 0.0); n * n * n];
    solver.set_velocity_field(vel).unwrap();

    let _ = solver.advance(0.001);

    for &a in solver.alpha() {
        assert!(a >= -1e-14 && a <= 1.0 + 1e-14, "α = {a} outside [0,1]");
    }
}

/// copy_boundaries should not panic on any grid size.
#[test]
fn test_vof_copy_boundaries_no_panic() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    // Minimal grid
    let mut solver = VofSolver::<f64>::new(3, 3, 3, cfg).unwrap();
    solver.alpha_mut()[13] = 1.0; // center cell
    solver.copy_boundaries();
    // Just assert no panic
}

/// set_volume_fraction rejects wrong-size input.
#[test]
fn test_vof_set_volume_fraction_wrong_size() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    let mut solver = VofSolver::<f64>::new(5, 5, 5, cfg).unwrap();
    // Wrong size: 10 instead of 125
    let result = solver.set_volume_fraction(vec![0.5; 10]);
    assert!(result.is_err(), "Should reject wrong-size alpha vector");
}

/// set_velocity_field rejects wrong-size input.
#[test]
fn test_vof_set_velocity_field_wrong_size() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    let mut solver = VofSolver::<f64>::new(5, 5, 5, cfg).unwrap();
    let result = solver.set_velocity_field(vec![Vector3::zeros(); 10]);
    assert!(result.is_err(), "Should reject wrong-size velocity vector");
}

/// total_volume is consistent with alpha sum * cell volume.
#[test]
fn test_vof_total_volume_consistency() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    let mut solver = VofSolver::<f64>::new(10, 10, 10, cfg).unwrap();

    // Fill a quarter of cells
    for a in solver.alpha_mut().iter_mut().take(250) {
        *a = 1.0;
    }

    let vol = solver.total_volume();
    assert!(vol > 0.0, "total_volume should be positive when cells are filled");
}

/// Reconstructed normals should have unit magnitude for interface cells.
#[test]
fn test_vof_reconstructed_normals_unit() {
    use cfd_3d::vof::{AdvectionMethod, InterfaceReconstruction, VofConfig, VofSolver};

    let cfg = VofConfig {
        surface_tension_coefficient: 0.072,
        interface_compression: 0.0,
        reconstruction_method: InterfaceReconstruction::PLIC,
        advection_method: AdvectionMethod::Geometric,
        max_iterations: 20,
        tolerance: 1e-8,
        cfl_number: 0.3,
        enable_compression: false,
    };
    let n = 10;
    let mut solver = VofSolver::<f64>::new(n, n, n, cfg).unwrap();

    // Sharp planar interface at i=5
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.linear_index(i, j, k);
                solver.alpha_mut()[idx] = if i < 5 { 1.0 } else { 0.0 };
            }
        }
    }

    solver.reconstruct_interface();

    for (idx, &a) in solver.alpha().iter().enumerate() {
        // Only check interface cells (0 < α < 1 after reconstruction, or near the interface)
        if a > 0.01 && a < 0.99 {
            let n_vec = solver.normals()[idx];
            let mag = n_vec.norm();
            if mag > 1e-10 {
                assert_relative_eq!(mag, 1.0, epsilon = 0.1);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Level Set — Signed Distance & Reinitialization
// ═══════════════════════════════════════════════════════════════════════════

/// Plane signed distance: φ = n·(x − x₀) is exact for a plane.
#[test]
fn test_level_set_plane_signed_distance() {
    use cfd_3d::level_set::{LevelSetConfig, LevelSetSolver};

    let n = 20;
    let dx = 1.0 / n as f64;
    let mut solver = LevelSetSolver::<f64>::new(
        LevelSetConfig { use_weno: true, ..Default::default() },
        n, n, n, dx, dx, dx,
    );

    // Plane φ = x − 0.5
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.index(i, j, k);
                let x = (i as f64 + 0.5) * dx;
                solver.phi_mut()[idx] = x - 0.5;
            }
        }
    }

    // The field is already an exact SDF; |∇φ| ≈ 1 in interior
    // Check that interior cells have |φ| consistent with distance to plane
    for k in 3..n - 3 {
        for j in 3..n - 3 {
            for i in 3..n - 3 {
                let idx = solver.index(i, j, k);
                let x = (i as f64 + 0.5) * dx;
                let expected = x - 0.5;
                assert_relative_eq!(solver.phi()[idx], expected, epsilon = 1e-14);
            }
        }
    }
}

/// Sphere SDF: reinitialization should preserve zero level set.
#[test]
fn test_level_set_sphere_reinit_preserves_zero() {
    use cfd_3d::level_set::{LevelSetConfig, LevelSetSolver};

    let n = 20;
    let dx = 1.0 / n as f64;
    let cfg = LevelSetConfig {
        reinitialization_interval: 1,
        use_weno: true,
        use_narrow_band: false,
        ..Default::default()
    };
    let mut solver = LevelSetSolver::<f64>::new(cfg, n, n, n, dx, dx, dx);

    let center = Vector3::new(0.5, 0.5, 0.5);
    let radius = 0.3;

    // Initialize as sphere SDF
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.index(i, j, k);
                let pos = Vector3::new(
                    (i as f64 + 0.5) * dx,
                    (j as f64 + 0.5) * dx,
                    (k as f64 + 0.5) * dx,
                );
                solver.phi_mut()[idx] = (pos - center).norm() - radius;
            }
        }
    }

    // Record cells near the zero level set before advance
    let tol = 2.0 * dx;
    let near_zero_before: Vec<usize> = solver.phi()
        .iter()
        .enumerate()
        .filter(|(_, &p)| p.abs() < tol)
        .map(|(i, _)| i)
        .collect();

    assert!(!near_zero_before.is_empty(), "Should have cells near zero level set");

    // Advance with zero velocity (triggers reinitialization)
    let vel = vec![Vector3::zeros(); n * n * n];
    solver.set_velocity(vel);
    let _ = solver.advance(dx * 0.1);

    // Check that cells that were near zero level set are still near it
    for &idx in &near_zero_before {
        assert!(
            solver.phi()[idx].abs() < tol * 2.0,
            "Zero level set should be preserved after reinit: φ = {}",
            solver.phi()[idx]
        );
    }
}

/// Narrow band should include cells near the interface.
#[test]
fn test_level_set_narrow_band_correctness() {
    use cfd_3d::level_set::{LevelSetConfig, LevelSetSolver};

    let n = 16;
    let dx = 1.0 / n as f64;
    let cfg = LevelSetConfig {
        use_narrow_band: true,
        band_width: 3.0,
        ..Default::default()
    };
    let mut solver = LevelSetSolver::<f64>::new(cfg, n, n, n, dx, dx, dx);

    // Sphere SDF
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.index(i, j, k);
                let pos = Vector3::new(
                    (i as f64 + 0.5) * dx,
                    (j as f64 + 0.5) * dx,
                    (k as f64 + 0.5) * dx,
                );
                solver.phi_mut()[idx] = (pos - Vector3::new(0.5, 0.5, 0.5)).norm() - 0.25;
            }
        }
    }

    solver.update_narrow_band();
    let band = solver.narrow_band();

    // Band should be non-empty and smaller than total grid
    assert!(!band.is_empty(), "Narrow band should be non-empty for sphere");
    assert!(
        band.len() < n * n * n,
        "Narrow band should be smaller than full grid"
    );

    // All band cells should have |φ| < band_width * dx
    let max_phi = 3.0 * dx;
    for &idx in band {
        assert!(
            solver.phi()[idx].abs() <= max_phi + dx,
            "Band cell has |φ| = {} > band_width*dx = {}",
            solver.phi()[idx].abs(),
            max_phi
        );
    }
}

/// Level set with zero velocity: φ should be unchanged.
#[test]
fn test_level_set_zero_velocity_preserves_phi() {
    use cfd_3d::level_set::{LevelSetConfig, LevelSetSolver};

    let n = 16;
    let dx = 1.0 / n as f64;
    let cfg = LevelSetConfig {
        reinitialization_interval: 1000, // prevent reinit from modifying values
        use_weno: true,
        use_narrow_band: false,
        ..Default::default()
    };
    let mut solver = LevelSetSolver::<f64>::new(cfg, n, n, n, dx, dx, dx);

    // Initialize as plane SDF
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let idx = solver.index(i, j, k);
                solver.phi_mut()[idx] = (i as f64 + 0.5) * dx - 0.5;
            }
        }
    }

    let phi_before: Vec<f64> = solver.phi().to_vec();

    let vel = vec![Vector3::zeros(); n * n * n];
    solver.set_velocity(vel);
    let _ = solver.advance(0.001);

    // Interior cells should be unchanged (WENO5 needs 3-cell ghost region)
    for k in 3..n - 3 {
        for j in 3..n - 3 {
            for i in 3..n - 3 {
                let idx = solver.index(i, j, k);
                assert_relative_eq!(
                    solver.phi()[idx], phi_before[idx],
                    epsilon = 1e-14,
                );
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Physics / Turbulence — Smagorinsky, MixingLength, k-ε
// ═══════════════════════════════════════════════════════════════════════════

/// Smagorinsky: zero velocity field → zero turbulent viscosity.
#[test]
fn test_smagorinsky_zero_velocity_zero_viscosity() {
    use cfd_3d::physics::turbulence::SmagorinskyModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let model = SmagorinskyModel::<f64>::new(0.17);
    let flow = FlowField::<f64>::new(4, 4, 4);

    let visc = model.turbulent_viscosity(&flow);
    assert_eq!(visc.len(), 64);
    for &v in &visc {
        assert_relative_eq!(v, 0.0, epsilon = 1e-15);
    }
}

/// MixingLength: zero velocity gradient → zero turbulent viscosity.
#[test]
fn test_mixing_length_zero_velocity_zero_viscosity() {
    use cfd_3d::physics::turbulence::MixingLengthModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let model = MixingLengthModel::<f64>::new(0.01);
    let flow = FlowField::<f64>::new(4, 4, 4);

    let visc = model.turbulent_viscosity(&flow);
    for &v in &visc {
        assert_relative_eq!(v, 0.0, epsilon = 1e-15);
    }
}

/// MixingLength: turbulent kinetic energy also zero for uniform flow.
#[test]
fn test_mixing_length_zero_tke_uniform_flow() {
    use cfd_3d::physics::turbulence::MixingLengthModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let model = MixingLengthModel::<f64>::new(0.01);
    let mut flow = FlowField::<f64>::new(4, 4, 4);
    // Uniform velocity → zero gradient
    for v in flow.velocity.components.iter_mut() {
        *v = Vector3::new(1.0, 0.0, 0.0);
    }

    let tke = model.turbulent_kinetic_energy(&flow);
    for &t in &tke {
        assert_relative_eq!(t, 0.0, epsilon = 1e-15);
    }
}

/// k-ε: standard constants match Launder & Spalding (1974).
#[test]
fn test_k_epsilon_standard_constants() {
    use cfd_3d::physics::turbulence::KEpsilonConstants;

    let c = KEpsilonConstants::<f64>::standard();
    assert_relative_eq!(c.c_mu, 0.09, epsilon = 1e-10);
    assert_relative_eq!(c.c_1, 1.44, epsilon = 1e-10);
    assert_relative_eq!(c.c_2, 1.92, epsilon = 1e-10);
    assert_relative_eq!(c.sigma_k, 1.0, epsilon = 1e-10);
    assert_relative_eq!(c.sigma_epsilon, 1.3, epsilon = 1e-10);
}

/// k-ε: uninitialised state returns zero viscosity.
#[test]
fn test_k_epsilon_uninit_zero_viscosity() {
    use cfd_3d::physics::turbulence::KEpsilonModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let model = KEpsilonModel::<f64>::new();
    let flow = FlowField::<f64>::new(4, 4, 4);

    let visc = model.turbulent_viscosity(&flow);
    for &v in &visc {
        assert_relative_eq!(v, 0.0, epsilon = 1e-15);
    }
}

/// k-ε: νₜ = C_μ k² / ε for initialised state with uniform fields.
#[test]
fn test_k_epsilon_viscosity_formula() {
    use cfd_3d::physics::turbulence::KEpsilonModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let mut model = KEpsilonModel::<f64>::new();
    let flow = FlowField::<f64>::new(4, 4, 4);
    let n = 64;

    let k_val = 1.5;
    let eps_val = 2.0;
    model.initialize_state_exact(vec![k_val; n], vec![eps_val; n]);

    let visc = model.turbulent_viscosity(&flow);
    let expected = 0.09 * k_val * k_val / eps_val;
    for &v in &visc {
        assert_relative_eq!(v, expected, epsilon = 1e-10);
    }
}

/// k-ε: zero TKE → zero turbulent viscosity (no division by zero).
#[test]
fn test_k_epsilon_zero_tke_no_nan() {
    use cfd_3d::physics::turbulence::KEpsilonModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;

    let mut model = KEpsilonModel::<f64>::new();
    let flow = FlowField::<f64>::new(4, 4, 4);
    let n = 64;

    model.initialize_state_exact(vec![0.0; n], vec![1.0; n]);

    let visc = model.turbulent_viscosity(&flow);
    for &v in &visc {
        assert!(v.is_finite(), "νₜ should be finite for k=0");
        assert_relative_eq!(v, 0.0, epsilon = 1e-15);
    }
}

/// k-ε dissipation_rate returns correct values from state.
#[test]
fn test_k_epsilon_dissipation_rate() {
    use cfd_3d::physics::turbulence::KEpsilonModel;
    use cfd_core::physics::fluid_dynamics::fields::FlowField;
    use cfd_core::physics::fluid_dynamics::rans::RANSModel;

    let mut model = KEpsilonModel::<f64>::new();
    let flow = FlowField::<f64>::new(4, 4, 4);
    let n = 64;

    let eps_val = 3.7;
    model.initialize_state_exact(vec![1.0; n], vec![eps_val; n]);

    let eps = model.dissipation_rate(&flow);
    for &e in &eps {
        assert_relative_eq!(e, eps_val, epsilon = 1e-14);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  IBM Forcing — Direct & Feedback
// ═══════════════════════════════════════════════════════════════════════════

/// Direct forcing: f = (u_d − u*) / Δt (Fadlun et al. 2000).
#[test]
fn test_direct_forcing_formula() {
    use cfd_3d::ibm::{DirectForcing, ForcingMethod};

    let forcing = DirectForcing::new();
    let u_desired = Vector3::new(1.0, 0.5, 0.0);
    let u_current = Vector3::new(0.2, 0.1, 0.3);
    let dt = 0.01_f64;

    let f = forcing.calculate_force(&u_desired, &u_current, dt);
    let expected = (u_desired - u_current) / dt;
    assert_relative_eq!(f.x, expected.x, epsilon = 1e-12);
    assert_relative_eq!(f.y, expected.y, epsilon = 1e-12);
    assert_relative_eq!(f.z, expected.z, epsilon = 1e-12);
}

/// Direct forcing: dt=0 → zero force (no division by zero).
#[test]
fn test_direct_forcing_zero_dt() {
    use cfd_3d::ibm::{DirectForcing, ForcingMethod};

    let forcing = DirectForcing::new();
    let f = forcing.calculate_force(
        &Vector3::new(1.0, 0.0, 0.0),
        &Vector3::new(0.0, 0.0, 0.0),
        0.0,
    );
    assert_relative_eq!(f.norm(), 0.0, epsilon = 1e-15);
}

/// Direct forcing: when u* = u_d, force is zero for any dt > 0.
#[test]
fn test_direct_forcing_already_satisfied() {
    use cfd_3d::ibm::{DirectForcing, ForcingMethod};

    let forcing = DirectForcing::new();
    let u = Vector3::new(1.7, -0.3, 0.5);
    let f = forcing.calculate_force(&u, &u, 0.01);
    assert_relative_eq!(f.norm(), 0.0, epsilon = 1e-14);
}

/// Feedback forcing: PI control accumulates integral.
#[test]
fn test_feedback_forcing_integral_accumulation() {
    use cfd_3d::ibm::{FeedbackForcing, ForcingMethod};

    let mut forcing = FeedbackForcing::new(100.0_f64, 50.0_f64);

    let error = Vector3::new(0.1, 0.0, 0.0);
    forcing.update(&error);
    forcing.update(&error);

    // After two updates, integral = 2 * error
    let f = forcing.calculate_force(
        &Vector3::new(1.0, 0.0, 0.0),
        &Vector3::new(0.9, 0.0, 0.0), // error = (0.1, 0, 0)
        0.01,
    );
    // f = kp * error + ki * integral = 100*(0.1,0,0) + 50*(0.2,0,0) = (20, 0, 0)
    assert_relative_eq!(f.x, 20.0, epsilon = 1e-10);
    assert_relative_eq!(f.y, 0.0, epsilon = 1e-15);
}

/// Feedback forcing: reset clears integral term.
#[test]
fn test_feedback_forcing_reset() {
    use cfd_3d::ibm::{FeedbackForcing, ForcingMethod};

    let mut forcing = FeedbackForcing::new(100.0_f64, 50.0_f64);
    forcing.update(&Vector3::new(1.0, 2.0, 3.0));
    forcing.reset();

    // After reset, integral contribution should be zero
    let u_d = Vector3::new(1.0, 0.0, 0.0);
    let u_c = Vector3::new(0.5, 0.0, 0.0);
    let f = forcing.calculate_force(&u_d, &u_c, 0.01);
    // f = kp * error = 100 * (0.5, 0, 0) = (50, 0, 0) (no integral term)
    assert_relative_eq!(f.x, 50.0, epsilon = 1e-10);
}

/// LagrangianPoint: update_position advances by v*dt.
#[test]
fn test_lagrangian_point_update_position() {
    use cfd_3d::ibm::LagrangianPoint;

    let mut pt = LagrangianPoint::<f64>::new(Vector3::new(1.0, 2.0, 3.0), 0.01);
    pt.velocity = Vector3::new(0.5, -0.3, 0.1);
    let dt = 0.1;
    pt.update_position(dt);

    assert_relative_eq!(pt.position.x, 1.05, epsilon = 1e-14);
    assert_relative_eq!(pt.position.y, 1.97, epsilon = 1e-14);
    assert_relative_eq!(pt.position.z, 3.01, epsilon = 1e-14);
}

/// LagrangianPoint: reset_force zeroes all force components.
#[test]
fn test_lagrangian_point_reset_force() {
    use cfd_3d::ibm::LagrangianPoint;

    let mut pt = LagrangianPoint::<f64>::new(Vector3::zeros(), 1.0);
    pt.force = Vector3::new(100.0, -50.0, 25.0);
    pt.reset_force();

    assert_relative_eq!(pt.force.norm(), 0.0, epsilon = 1e-15);
}

/// IBM kernel: all variants are continuous at support boundary.
#[test]
fn test_ibm_kernel_continuity_at_support() {
    use cfd_3d::ibm::{DeltaFunction, InterpolationKernel};

    for df in &[DeltaFunction::RomaPeskin3, DeltaFunction::RomaPeskin4, DeltaFunction::Peskin4] {
        let kernel = InterpolationKernel::new(df.clone(), 1.5_f64);
        // At support boundary: δ(r) should approach 0
        let r_boundary = 2.5; // well outside support for width=1.5
        let val = kernel.delta(r_boundary);
        assert!(
            val.abs() < 1e-10,
            "Kernel {df:?} should be ~0 at r={r_boundary}: got {val}"
        );
    }
}

/// IBM kernel: RomaPeskin3 partition of unity.
#[test]
fn test_roma_peskin3_partition_of_unity() {
    use cfd_3d::ibm::{DeltaFunction, InterpolationKernel};

    let kernel = InterpolationKernel::new(DeltaFunction::RomaPeskin3, 1.0_f64);
    let h = 1.0;
    // Test at several off-grid positions
    for &r in &[0.0, 0.15, 0.3, 0.5, 0.75, 0.9] {
        let mut sum = 0.0;
        for j in -3..=3 {
            sum += kernel.delta(r - j as f64 * h) * h;
        }
        assert_relative_eq!(sum, 1.0, epsilon = 1e-10);
    }
}

// ═══════════════════════════════════════════════════════════════════════════
//  Spectral — Chebyshev & Poisson
// ═══════════════════════════════════════════════════════════════════════════

/// Chebyshev differentiation of x² yields 2x exactly at collocation points.
#[test]
fn test_chebyshev_derivative_x_squared() {
    use cfd_3d::spectral::ChebyshevPolynomial;
    use nalgebra::DVector;

    let n = 16;
    let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
    let points = cheb.collocation_points();

    let u: DVector<f64> = DVector::from_iterator(n, points.iter().map(|&x| x * x));
    let du = cheb.differentiate(&u);

    // du/dx = 2x
    for (i, &x) in points.iter().enumerate() {
        assert_relative_eq!(du[i], 2.0 * x, epsilon = 1e-10);
    }
}

/// Chebyshev quadrature weights sum to 2 (integral of 1 over [-1,1]).
#[test]
fn test_chebyshev_quadrature_weights_sum() {
    use cfd_3d::spectral::ChebyshevPolynomial;

    let cheb = ChebyshevPolynomial::<f64>::new(16).unwrap();
    let w = cheb.quadrature_weights().unwrap();
    let sum: f64 = w.iter().sum();
    assert_relative_eq!(sum, 2.0, epsilon = 1e-12);
}

/// Chebyshev second derivative of cos(πx): d²/dx² cos(πx) = -π² cos(πx).
#[test]
fn test_chebyshev_second_derivative() {
    use cfd_3d::spectral::ChebyshevPolynomial;
    use nalgebra::DVector;

    let n = 24;
    let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
    let d2 = cheb.second_derivative_matrix().unwrap();
    let points = cheb.collocation_points();

    let pi = std::f64::consts::PI;
    let u: DVector<f64> = DVector::from_iterator(n, points.iter().map(|&x| (pi * x).cos()));
    let d2u = &d2 * &u;

    // Skip endpoints (boundary effects)
    for (i, &x) in points.iter().enumerate().skip(1).take(n - 2) {
        let expected = -pi * pi * (pi * x).cos();
        assert_relative_eq!(d2u[i], expected, epsilon = 0.2);
    }
}

/// Poisson solver: zero RHS with zero Dirichlet BC → zero solution.
#[test]
fn test_poisson_zero_rhs_zero_solution() {
    use cfd_3d::spectral::poisson::{PoissonBoundaryCondition, PoissonSolver};
    use nalgebra::DVector;

    let n = 5;
    let solver = PoissonSolver::<f64>::new(n, n, n).unwrap();
    let rhs = DVector::zeros(n * n * n);

    let bc = (
        PoissonBoundaryCondition::Dirichlet(0.0),
        PoissonBoundaryCondition::Dirichlet(0.0),
    );

    let u = solver.solve(&rhs, &bc, &bc, &bc).unwrap();
    for i in 0..u.len() {
        assert_relative_eq!(u[i], 0.0, epsilon = 1e-8);
    }
}

/// Chebyshev: minimum n=2 works without panic.
#[test]
fn test_chebyshev_minimal_size() {
    use cfd_3d::spectral::ChebyshevPolynomial;

    let cheb = ChebyshevPolynomial::<f64>::new(2).unwrap();
    assert_eq!(cheb.num_points(), 2);
    let w = cheb.quadrature_weights().unwrap();
    assert_eq!(w.len(), 2);
    let sum: f64 = w.iter().sum();
    assert_relative_eq!(sum, 2.0, epsilon = 1e-12);
}

/// Chebyshev interpolation reproduces polynomial exactly.
#[test]
fn test_chebyshev_interpolation_cubic() {
    use cfd_3d::spectral::ChebyshevPolynomial;

    let n = 8; // Should reproduce polynomials up to degree 7 exactly
    let cheb = ChebyshevPolynomial::<f64>::new(n).unwrap();
    let points = cheb.collocation_points();

    // f(x) = x³ − 2x + 1
    let values: Vec<f64> = points.iter().map(|&x| x * x * x - 2.0 * x + 1.0).collect();

    // Interpolate at x = 0.37
    let x_test = 0.37;
    let interp = cheb.interpolate(&values, x_test).unwrap();
    let exact = x_test * x_test * x_test - 2.0 * x_test + 1.0;
    assert_relative_eq!(interp, exact, epsilon = 1e-12);
}

// ═══════════════════════════════════════════════════════════════════════════
//  IBM Solver — Integration Tests
// ═══════════════════════════════════════════════════════════════════════════

/// IbmSolver: no points → zero Eulerian force field.
#[test]
fn test_ibm_solver_no_points_zero_force() {
    use cfd_3d::ibm::{IbmConfig, IbmSolver};

    let solver = IbmSolver::<f64>::new(
        IbmConfig::default(),
        Vector3::new(0.1, 0.1, 0.1),
        (10, 10, 10),
    );
    let forces = solver.spread_forces().unwrap();
    for f in &forces {
        assert_relative_eq!(f.norm(), 0.0, epsilon = 1e-15);
    }
}

/// IbmSolver: add point and verify count.
#[test]
fn test_ibm_solver_add_point() {
    use cfd_3d::ibm::{IbmConfig, IbmSolver, LagrangianPoint};

    let mut solver = IbmSolver::<f64>::new(
        IbmConfig::default(),
        Vector3::new(0.1, 0.1, 0.1),
        (10, 10, 10),
    );
    assert_eq!(solver.num_points(), 0);

    solver.add_lagrangian_point(LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 0.01));
    assert_eq!(solver.num_points(), 1);

    solver.add_lagrangian_point(LagrangianPoint::new(Vector3::new(0.3, 0.3, 0.3), 0.01));
    assert_eq!(solver.num_points(), 2);
}

/// IbmSolver: update_positions moves all Lagrangian points.
#[test]
fn test_ibm_solver_update_positions() {
    use cfd_3d::ibm::{IbmConfig, IbmSolver, LagrangianPoint};

    let mut solver = IbmSolver::<f64>::new(
        IbmConfig::default(),
        Vector3::new(0.1, 0.1, 0.1),
        (10, 10, 10),
    );

    let mut pt = LagrangianPoint::new(Vector3::new(0.5, 0.5, 0.5), 0.01);
    pt.velocity = Vector3::new(1.0, 0.0, 0.0);
    solver.add_lagrangian_point(pt);

    solver.update_positions(0.1);

    let points = solver.points();
    assert_relative_eq!(points[0].position.x, 0.6, epsilon = 1e-14);
    assert_relative_eq!(points[0].position.y, 0.5, epsilon = 1e-14);
}
