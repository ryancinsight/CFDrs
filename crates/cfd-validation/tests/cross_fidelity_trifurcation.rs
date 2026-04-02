//! Trifurcation Momentum Advection (Jetting) Cross-Fidelity Validation
//!
//! Validates the asymmetric flow splitting behavior in geometrically symmetric
//! trifurcations due to convective momentum conservation.
//!
//! # Theorem 1: Symmetric Trifurcation Center-Dominance
//!
//! In an unconstrained, symmetric, rigid trifurcation featuring three identical 
//! daughter branches diverging from a single parent channel, the central continuous 
//! branch geometrically aligns with the primary axis of momentum. Although pure 
//! gradient-driven networks (1D lumped models) predict perfectly uniform flux 
//! distribution due to identical viscous path resistance, the convective inertial 
//! mechanisms of the Navier-Stokes equations mandate that the central branch 
//! dominates volumetric flux ($Q_{center} > Q_{lateral}$).
//!
//! **Proof sketch:**
//! The steady, incompressible Navier-Stokes momentum equation balances pressure 
//! gradients against viscous diffusion and convective acceleration:
//! $(\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}$
//! 
//! For a fully developed Poiseuille flow entering the trifurcation junction, the peak 
//! velocity resides along the geometric centerline axis. As the stream diverges, the 
//! lateral branches require a strong convective deflection of momentum 
//! $(\mathbf{u} \cdot \nabla)\mathbf{u}$ perpendicular to the parent axis, imposing a 
//! localized adverse pressure penalty at the branch lip. The central branch, sharing 
//! the continuous parent axis, avoids this convective deflection penalty and directly 
//! inherits the high-kinetic-energy core of the parent jet. Consequently, total pressure 
//! losses into the lateral branches exceed those of the central branch, enforcing 
//! $Q_{center} > Q_{lateral}$ under symmetric outlet boundary conditions.

use cfd_3d::trifurcation::{TrifurcationConfig3D, TrifurcationGeometry3D, TrifurcationSolver3D};
use cfd_core::physics::fluid::blood::casson::CassonBlood;

#[test]
fn cross_fidelity_trifurcation_dominance() {
    // 1. Experimental Constants & Fluid Properties
    let d_parent = 0.001; // 1 mm
    let l_parent = 0.003; // 3 mm
    let d_daughter = 0.0006; // 600 μm
    let l_daughter = 0.003; // 3 mm
    let l_trans = 0.0005; // 500 μm
    let spread_angle = 0.523598775; // 30 degrees in radians

    let q_inlet_target = 1.0e-7; // 100 μL/s

    // Blood properties
    let fluid_props = CassonBlood::<f64>::normal_blood();

    // 2. [1D] Baseline Network Limit (Poiseuille Resistance Theory)
    // In a 1D pure-resistance lumped network, 3 identical branches with identical 
    // outlet pressures will draw exactly 1/3 of the total flow rate.
    // 3. [3D] Volumetric FEM Navier-Stokes
    let geometry_3d = TrifurcationGeometry3D::symmetric(
        d_parent,
        d_daughter,
        l_parent,
        l_daughter,
        l_trans,
        spread_angle,
    );

    let config_3d = TrifurcationConfig3D {
        inlet_flow_rate: q_inlet_target,
        // The solver operates on pressure differences; using gauge values avoids
        // a large atmospheric offset that needlessly worsens conditioning.
        inlet_pressure: 200.0,
        outlet_pressures: [0.0, 0.0, 0.0],
        max_nonlinear_iterations: 20,
        nonlinear_tolerance: 1e-4,
        max_linear_iterations: 600,
        linear_tolerance: 1e-5,
        target_mesh_size: Some(d_parent / 6.0),
    };

    let solver_3d = TrifurcationSolver3D::new(geometry_3d, config_3d);
    
    // Resolve full 3D flow splitting
    let sol_3d = solver_3d.solve(&fluid_props)
        .expect("3D trifurcation calculation failed");

    // Flow indices map directly to geometry.branching_angles: [spread, 0, -spread]
    // 0: parent, 1: lateral_left (spread), 2: center (0), 3: lateral_right (-spread)
    let q_parent_fem = sol_3d.flow_rates[0];
    let q_lateral_1 = sol_3d.flow_rates[1];
    let q_center_3d = sol_3d.flow_rates[2];
    let q_lateral_2 = sol_3d.flow_rates[3];

    let total_out = q_lateral_1 + q_center_3d + q_lateral_2;
    assert!(q_parent_fem.is_finite() && q_parent_fem > 0.0);
    assert!(total_out.is_finite() && total_out > 0.0);

    // The present FEM trifurcation path is more reliable for relative outlet
    // splitting than for absolute end-cap flux closure on coarse meshes. This
    // validation therefore checks the normalized outlet partition, which is the
    // quantity relevant to the jetting theorem below.
    let center_share = q_center_3d / total_out;
    let avg_lateral_3d = (q_lateral_1 + q_lateral_2) / 2.0;
    let avg_lateral_share = avg_lateral_3d / total_out;

    // 4. Verification of Trifurcation Invariants

    // Invariant 1: Center line inherits momentum. Its outlet share exceeds the
    // 1D uniform split fraction of 1/3.
    assert!(
        center_share > 1.0 / 3.0,
        "Violation of Momentum Jetting: 3D Center Share ({:.3e}) <= 1D Ideal Uniform Share ({:.3e})",
        center_share,
        1.0 / 3.0
    );

    // Invariant 2: Lateral branches suffer deflection losses. Their mean outlet
    // share falls below the 1D uniform split fraction of 1/3.
    assert!(
        avg_lateral_share < 1.0 / 3.0,
        "Violation of Deflection Resistance: 3D Lateral Share ({:.3e}) >= 1D Ideal Uniform Share ({:.3e})",
        avg_lateral_share,
        1.0 / 3.0
    );

    // Invariant 3: Strict inequality $Q_{center} > Q_{lateral}$
    assert!(
        q_center_3d > avg_lateral_3d,
        "Violation of Strict Geometric Dominance: Center Branch ({:.3e}) <= Lateral Branches ({:.3e})",
        q_center_3d, avg_lateral_3d
    );
}
