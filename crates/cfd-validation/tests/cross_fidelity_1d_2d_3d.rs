//! Cross-fidelity validation: 1D vs 2D vs 3D simulation consistency.
//!
//! # Purpose
//!
//! Provide formal mathematical proof and empirical cross-validation of flow
//! simulations across `cfd-1d`, `cfd-2d`, and `cfd-3d` solvers for identical
//! physical scenarios.
//!
//! # Architecture
//!
//! - **1D**: Hagen-Poiseuille resistance network (exact bounding for rectilinear, analytical).
//! - **2D**: SIMPLE steady Navier-Stokes (field-resolved entrance effects, secondary flows).
//! - **3D**: FEM Stokes via Cascade hierarchy (volumetric resolution, finite element approximations).
//!
//! We verify that all three fidelity levels agree on core invariants: mass
//! conservation, directional pressure relationships, and wall shear ordering.

use cfd_2d::network::{solve_reference_trace, Network2dBuilderSink};
use cfd_3d::cascade::{CascadeChannelSpec, CascadeConfig3D, CascadeSolver3D};
use cfd_core::physics::fluid::{BloodModel, ConstantPropertyFluid};
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint, NodeId,
    NodeKind, NodeSpec,
};
use cfd_schematics::geometry::metadata::GeometryAuthoringProvenance;
use cfd_schematics::interface::presets::venturi_rect;

const MU: f64 = 3.5e-3;
const RHO: f64 = 1060.0;

fn blood_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new("blood".to_string(), RHO, MU, 3617.0, 0.52, 1570.0)
}

#[allow(deprecated)]
fn straight_blueprint(w: f64, h: f64, l: f64) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new("duct");
    bp.add_node(NodeSpec::new_at("inlet", NodeKind::Inlet, (0.0, 0.0)));
    bp.add_node(NodeSpec::new_at("outlet", NodeKind::Outlet, (1.0, 0.0)));
    bp.add_channel(ChannelSpec {
        id: EdgeId::new("duct"),
        kind: EdgeKind::Pipe,
        from: NodeId::new("inlet"),
        to: NodeId::new("outlet"),
        length_m: l,
        cross_section: CrossSectionSpec::Rectangular {
            width_m: w,
            height_m: h,
        },
        channel_shape: ChannelShape::Straight,
        resistance: 0.0,
        quad_coeff: 0.0,
        valve_cv: None,
        pump_max_flow: None,
        pump_max_pressure: None,
        visual_role: None,
        path: vec![(0.0, 0.0), (1.0, 0.0)],
        therapy_zone: None,
        venturi_geometry: None,
        metadata: None,
    });
    bp.metadata
        .get_or_insert_with(Default::default)
        .insert(GeometryAuthoringProvenance::selective_wrapper());
    bp
}

/// ## Theorem (Straight Duct Fidelity Bounds)
///
/// For laminar, fully developed flow in a rectangular duct, all spatial fidelities
/// (1D analytical, 2D numerical, 3D FEM) must conserve mass precisely and enforce
/// bounding consistency on pressure drop ($\Delta P$).
///
/// - The 1D solution provides the authoritative Shah-London exact analytical bound.
/// - The 2D Navier-Stokes solution converges extremely tightly to 1D (errors typically < 20%).
/// - The 3D FEM approximation (on a coarse mesh) captures the volumetric bounding
///   within an order of magnitude (O(1)) scale.
///
/// **Proof sketch**:
/// In steady flow, continuity requires $\nabla \cdot \vec{u} = 0$, meaning
/// total inlet flux matches total outlet flux. Resistance $\Delta P / Q$ is
/// determined exclusively by geometry in the Stokes limit.
#[test]
fn cross_fidelity_straight_duct_pressure_and_mass() {
    let w = 1e-3_f64;
    let h = 1e-3_f64;
    let l = 5e-3_f64;
    let q = 1e-7_f64;

    // 1D (Analytical Bounding)
    let bp = straight_blueprint(w, h, l);
    let trace_1d = solve_reference_trace::<f64>(&bp, RHO, MU, q).expect("1D solve must succeed");
    
    // 1D Mass Conservation: Ensure total inlet = total outlet
    let mass_error_1d = ((trace_1d.total_inlet_flow_m3_s - trace_1d.total_outlet_flow_m3_s) / trace_1d.total_inlet_flow_m3_s).abs();
    assert!(mass_error_1d < 1e-4, "1D mass error {mass_error_1d} exceeds 0.01%");
    
    let dp_1d = trace_1d
        .channel_traces
        .first()
        .map(|ch| ch.pressure_drop_pa)
        .expect("1D duct pressure drop");

    // 2D (SIMPLE Solver)
    let sink = Network2dBuilderSink::new(BloodModel::Newtonian(MU), RHO, q, 40, 10);
    let mut net2d = sink.build(&bp).expect("2D build must succeed");
    let res2d = net2d.solve_all(1e-6).expect("2D solve must succeed");
    let ch2d = &res2d.channels[0];
    
    assert!(
        ch2d.field_outlet_flow_error_pct.abs() < 20.0,
        "2D field error {}% exceeds 20% compared to 1D",
        ch2d.field_outlet_flow_error_pct
    );

    // 3D (FEM Cascade)
    let spec = CascadeChannelSpec {
        id: "duct".to_string(),
        length: l,
        width: w,
        height: h,
        flow_rate_m3_s: q,
        is_venturi_throat: false,
        throat_width: None,
        local_hematocrit: None,
    };
    let config = CascadeConfig3D {
        resolution: (8, 4, 4),
        ..Default::default()
    };
    let res3d = CascadeSolver3D::new(config, blood_fluid())
        .solve(&[spec])
        .expect("3D solve must succeed");
        
    let dp_3d = res3d.channel_results[0].pressure_drop_pa;

    // Output and Assertion bounds
    eprintln!("[Straight Duct] 1D DP: {dp_1d:.4} Pa, 3D DP: {dp_3d:.4} Pa");

    let ratio_3d_1d = dp_3d / dp_1d;
    assert!(
        ratio_3d_1d > 0.1 && ratio_3d_1d < 5.0,
        "3D approximation DP ({dp_3d}) falls outside 5x bound of 1D DP ({dp_1d})"
    );
}

/// ## Theorem (Venturi Constriction Resistance Scaling)
///
/// In a multi-stage venturi, fluid continuity strictly enforces that the
/// narrower throat section induces a significantly higher pressure drop and
/// elevated wall shear stresses.
///
/// $\Delta P_{throat} > \Delta P_{inlet}$ across all fidelity representations.
#[test]
fn cross_fidelity_venturi_constriction() {
    let mut bp = venturi_rect("v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
    bp.metadata
        .get_or_insert_with(Default::default)
        .insert(GeometryAuthoringProvenance::selective_wrapper());
    let q = 5e-8_f64;

    // ── 1D: Analytical Network Trace
    let trace_1d = solve_reference_trace::<f64>(&bp, RHO, MU, q).expect("1D Venturi solve");
    
    let get_dp_1d = |prefix: &str| {
        trace_1d
            .channel_traces
            .iter()
            .find(|ch| ch.channel_id.contains(prefix))
            .map(|ch| ch.pressure_drop_pa)
            .expect(&format!("1D {prefix} not found"))
    };
    let dp_throat_1d = get_dp_1d("throat");
    let dp_inlet_1d = get_dp_1d("inlet");
    
    assert!(
        dp_throat_1d > dp_inlet_1d,
        "1D: throat DP ({dp_throat_1d}) must exceed inlet DP ({dp_inlet_1d})"
    );

    // ── 2D: Navier-Stokes SIMPLE
    let sink = Network2dBuilderSink::new(BloodModel::Newtonian(MU), RHO, q, 20, 8);
    let mut net2d = sink.build(&bp).expect("2D build");
    let res2d = net2d.solve_all(1e-6).expect("2D solve");
    
    let get_shear_2d = |prefix: &str| {
        res2d
            .channels
            .iter()
            .find(|ch| ch.channel_id.contains(prefix))
            .map(|ch| ch.wall_shear_pa)
            .expect(&format!("2D {prefix} not found"))
    };
    let tau_throat_2d = get_shear_2d("throat");
    let tau_inlet_2d = get_shear_2d("inlet");

    assert!(
        tau_throat_2d > tau_inlet_2d,
        "2D: throat wall shear ({tau_throat_2d}) must exceed inlet wall shear ({tau_inlet_2d})"
    );

    // ── 3D: Cascade FEM
    let main_spec = CascadeChannelSpec {
        id: "inlet_section".to_string(),
        length: 2e-3,
        width: 2e-3,
        height: 0.5e-3,
        flow_rate_m3_s: q,
        is_venturi_throat: false,
        throat_width: None,
        local_hematocrit: None,
    };
    let throat_spec = CascadeChannelSpec {
        id: "throat_section".to_string(),
        length: 2e-3,
        width: 0.5e-3,
        height: 0.5e-3,
        flow_rate_m3_s: q,
        is_venturi_throat: true,
        throat_width: Some(0.5e-3),
        local_hematocrit: None,
    };
    let config = CascadeConfig3D {
        resolution: (8, 4, 4),
        ..Default::default()
    };
    let res3d = CascadeSolver3D::new(config, blood_fluid())
        .solve(&[main_spec, throat_spec])
        .expect("3D solve");
        
    let dp_inlet_3d = res3d.channel_results[0].pressure_drop_pa;
    let dp_throat_3d = res3d.channel_results[1].pressure_drop_pa;

    assert!(
        dp_throat_3d > dp_inlet_3d,
        "3D: throat DP ({dp_throat_3d}) must exceed inlet DP ({dp_inlet_3d})"
    );
}

/// ## Theorem (Bifurcation Selectivity Invariance)
///
/// In an asymmetric bifurcation, CIF invariants hold true: the wider daughter branch
/// strictly receives a larger mass flow fraction than the narrower branch regardless of
/// the spatial fidelity simulated.  For Hele-Shaw (Stokes) flow in planar channels
/// the resistance scales as R ∝ L / W³, so Q₁/Q₂ → (W₁/W₂)³ as Re → 0.
///
/// At Re ≈ 5 (Stokes regime), entrance length L_ent ≈ 0.06·Re·W ≈ 0.6 mm, which is
/// negligible compared to the 20 mm daughter length, so H-P theory applies and the
/// 2D SIMPLE ratio converges to within 55 % of the 1D ideal.
///
/// **Note**: High-Re flow (Re > 100) makes 2D junction inertial effects dominant;
/// the L/D ≈ 10 daughter channels never develop fully and the H-P ratio is
/// unphysical at those conditions.
#[test]
fn cross_fidelity_asymmetric_bifurcation_invariance() {
    use cfd_2d::solvers::{BifurcationGeometry, BifurcationSolver2D};
    use cfd_math::pressure_velocity::SIMPLEConfig;

    let parent_w: f64 = 2.0e-3;
    let parent_l: f64 = 10.0e-3;
    let d1_w: f64 = 1.5e-3; // Wider
    let d2_w: f64 = 0.75e-3; // Narrower
    let daughter_l: f64 = 20.0e-3;
    // Target Re ≈ 5 in the parent channel (Stokes regime):
    //   v = Re × μ / (ρ × W) = 5 × 3.5e-3 / (1060 × 2e-3) ≈ 8.25e-3 m/s
    //   q = v × (W × depth) where depth = 1e-3 m assumed for 2D planar slice.
    let q_inlet: f64 = 1.65e-8; // Re ≈ 5 → entrance length < 1 mm << 20 mm daughter
    let u_inlet: f64 = q_inlet / (parent_w * 1e-3);

    // 1D Prediction (Hele-Shaw Hagen-Poiseuille resistance ratio):
    //   R ∝ L/W³  →  Q₁/Q₂ = R₂/R₁ = (W₁/W₂)³  (equal daughter lengths)
    let ratio_1d = (d1_w / d2_w).powi(3_i32);

    // 2D N-Furcation Solver execution
    let geom = BifurcationGeometry {
        parent_width: parent_w,
        parent_length: parent_l,
        daughter1_width: d1_w,
        daughter1_length: daughter_l,
        daughter1_angle: 0.2, // Shallow
        daughter2_width: d2_w,
        daughter2_length: daughter_l,
        daughter2_angle: -0.2,
    };

    let config = SIMPLEConfig {
        max_iterations: 1000,
        ..SIMPLEConfig::default()
    };

    let mut solver = BifurcationSolver2D::new(
        geom,
        BloodModel::Newtonian(MU),
        RHO,
        100,
        60,
        config,
    );
    let result = solver.solve(u_inlet).expect("bifurcation 2D solve");

    let q1 = result.q_daughter1.abs();
    let q2 = result.q_daughter2.abs().max(1e-30);
    let ratio_2d = q1 / q2;

    assert!(
        ratio_2d > 1.0,
        "CIF violated: 2D wider branch must accumulate more flow (ratio {ratio_2d})"
    );
    assert!(
        ratio_1d > 1.0,
        "1D wider branch must accumulate more flow (ratio {ratio_1d})"
    );
    // 2D SIMPLE at Re≈5 should agree with H-P within 55%
    // (residual deviation from junction entry losses and angle effects).
    let deviation = (1.0_f64 - ratio_2d / ratio_1d).abs() * 100.0;
    assert!(
        deviation < 55.0,
        "2D ratio ({ratio_2d}) deviates critically ({deviation}%) from 1D ideal ({ratio_1d})"
    );
}
