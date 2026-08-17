//! Cross-fidelity validation: venturi calibration across 1D, 2D, and 3D.
//!
//! This complements the existing directional venturi test with a quantitative
//! pressure-loss calibration around a shared total loss coefficient.

use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity, Velocity,
};
use cfd_2d::network::solve_reference_trace;
use cfd_2d::solvers::venturi_flow::{VenturiGeometry, VenturiSolver2D};
use cfd_3d::venturi::{VenturiConfig3D, VenturiSolver3D};
use cfd_core::physics::fluid::{BloodModel, ConstantPropertyFluid};
use cfd_mesh::VenturiMeshBuilder;
use cfd_schematics::domain::model::{CrossSectionSpec, NetworkBlueprint, NodeKind};
use cfd_schematics::geometry::metadata::GeometryAuthoringProvenance;
use cfd_schematics::interface::presets::venturi_rect;

const MU: f64 = 3.5e-3;
const RHO: f64 = 1060.0;
const FLOW_RATE_M3_S: f64 = 5.0e-8;

#[derive(Clone, Copy)]
struct CalibrationGeometry {
    width_inlet_m: f64,
    width_throat_m: f64,
    height_m: f64,
    inlet_length_m: f64,
    throat_length_m: f64,
    diffuser_length_m: f64,
}

fn blood_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new(
        "blood".to_string(),
        MassDensity::from_base(RHO),
        DynamicViscosity::from_base(MU),
        SpecificHeatCapacity::from_base(3617.0),
        ThermalConductivity::from_base(0.52),
        Velocity::from_base(1570.0),
    )
}

fn rectangular_dims(cross_section: &CrossSectionSpec) -> (f64, f64) {
    match cross_section {
        CrossSectionSpec::Rectangular { width_m, height_m } => {
            (width_m.into_base(), height_m.into_base())
        }
        other => panic!("Expected rectangular venturi section, got {other:?}"),
    }
}

fn channel<'a>(
    blueprint: &'a NetworkBlueprint,
    channel_id: &str,
) -> &'a cfd_schematics::domain::model::ChannelSpec {
    blueprint
        .channels
        .iter()
        .find(|channel| channel.id.as_str() == channel_id)
        .unwrap_or_else(|| panic!("Missing blueprint channel '{channel_id}'"))
}

fn inlet_outlet_pressures(trace: &cfd_2d::network::NetworkReferenceTrace<f64>) -> (f64, f64) {
    let inlet_pressure = trace
        .node_traces
        .iter()
        .filter(|node| matches!(node.node_kind, NodeKind::Inlet))
        .map(|node| node.pressure_pa)
        .fold(f64::NEG_INFINITY, f64::max);
    let outlet_pressure = trace
        .node_traces
        .iter()
        .filter(|node| matches!(node.node_kind, NodeKind::Outlet))
        .map(|node| node.pressure_pa)
        .fold(f64::INFINITY, f64::min);

    assert!(
        inlet_pressure.is_finite(),
        "Missing inlet pressure in 1D venturi trace"
    );
    assert!(
        outlet_pressure.is_finite(),
        "Missing outlet pressure in 1D venturi trace"
    );
    (inlet_pressure, outlet_pressure)
}

fn calibration_case() -> (NetworkBlueprint, CalibrationGeometry) {
    let mut blueprint = venturi_rect("venturi-calibration", 2.0e-3, 0.5e-3, 0.5e-3, 2.0e-3);
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(GeometryAuthoringProvenance::selective_wrapper());

    let inlet_section = channel(&blueprint, "inlet_section");
    let throat_section = channel(&blueprint, "throat_section");
    let diffuser_section = channel(&blueprint, "diffuser_section");
    let (width_inlet_m, height_m) = rectangular_dims(&inlet_section.cross_section);
    let (width_throat_m, _) = rectangular_dims(&throat_section.cross_section);
    let inlet_length_m = inlet_section.length_m.into_base();
    let throat_length_m = throat_section.length_m.into_base();
    let diffuser_length_m = diffuser_section.length_m.into_base();

    (
        blueprint,
        CalibrationGeometry {
            width_inlet_m,
            width_throat_m,
            height_m,
            inlet_length_m,
            throat_length_m,
            diffuser_length_m,
        },
    )
}

fn one_dimensional_loss_coefficient(
    blueprint: &NetworkBlueprint,
    geometry: CalibrationGeometry,
) -> (f64, f64) {
    let inlet_area_m2 = geometry.width_inlet_m * geometry.height_m;
    let inlet_mean_velocity_m_s = FLOW_RATE_M3_S / inlet_area_m2;
    let dynamic_pressure_pa = 0.5 * RHO * inlet_mean_velocity_m_s * inlet_mean_velocity_m_s;
    let trace_1d = solve_reference_trace::<f64>(blueprint, RHO, MU, FLOW_RATE_M3_S)
        .expect("1D venturi reference trace");
    let (p_inlet_1d, p_outlet_1d) = inlet_outlet_pressures(&trace_1d);
    (
        (p_inlet_1d - p_outlet_1d) / dynamic_pressure_pa,
        inlet_mean_velocity_m_s,
    )
}

/// ## Theorem (Venturi Loss-Coefficient Cross-Fidelity Consistency)
///
/// For a steady laminar rectangular venturi with fixed contraction ratio and
/// throat length, the 1D reference trace, 2D venturi solver, and 3D venturi
/// solver must all preserve:
///
/// - a positive total inlet-to-outlet pressure loss coefficient,
/// - throat acceleration relative to the inlet section,
/// - coefficient agreement to the same order of magnitude across fidelities.
///
/// **Proof sketch**:
/// The total loss coefficient $K = \Delta P / (\tfrac{1}{2}\rho u_{in}^2)$ is a
/// dimensionless pressure-loss invariant for a fixed contraction geometry. The
/// reduced 1D model aggregates viscous losses, while the 2D/3D solvers resolve
/// the contraction and recovery fields. Although absolute values differ due to
/// resolved geometry and coarse-mesh error, all fidelities should preserve the
/// same positive-loss regime and remain within an $O(1)$ envelope.
#[test]
fn cross_fidelity_venturi_two_dimensional_loss_coefficient() {
    let (blueprint, geometry) = calibration_case();
    let (total_loss_coeff_1d, inlet_mean_velocity_m_s) =
        one_dimensional_loss_coefficient(&blueprint, geometry);

    let l_inlet_2d = geometry.inlet_length_m * 0.5;
    let l_converge_2d = geometry.inlet_length_m * 0.5;
    let geometry_2d = VenturiGeometry::new(
        geometry.width_inlet_m,
        geometry.width_throat_m,
        l_inlet_2d,
        l_converge_2d,
        geometry.throat_length_m,
        geometry.diffuser_length_m,
        geometry.height_m,
    );
    let mut solver_2d = VenturiSolver2D::new(geometry_2d, BloodModel::Newtonian(MU), RHO, 120, 48);
    let solution_2d = solver_2d
        .solve(inlet_mean_velocity_m_s)
        .expect("2D venturi solve");
    let total_loss_coeff_2d = -solution_2d.cp_recovery;

    assert!(total_loss_coeff_1d.is_finite() && total_loss_coeff_1d > 0.0);
    assert!(total_loss_coeff_2d.is_finite() && total_loss_coeff_2d > 0.0);

    let throat_area_ratio = geometry.width_inlet_m / geometry.width_throat_m;
    let ratio_2d_to_1d = total_loss_coeff_2d / total_loss_coeff_1d;
    println!(
        "CALIBRATION 1D/2D: 1D={}, 2D={}, ratio={}",
        total_loss_coeff_1d, total_loss_coeff_2d, ratio_2d_to_1d
    );

    assert!(
        solution_2d.u_throat_mean > inlet_mean_velocity_m_s,
        "2D venturi throat mean velocity must exceed inlet mean velocity"
    );
    let throat_mean_velocity_ratio_2d = solution_2d.u_throat_mean / inlet_mean_velocity_m_s;
    let throat_mean_velocity_ratio_error_2d =
        ((throat_mean_velocity_ratio_2d - throat_area_ratio) / throat_area_ratio).abs();
    assert!(
        throat_mean_velocity_ratio_error_2d < 0.10,
        "2D venturi throat mean velocity ratio {} deviates by more than 10% from the area ratio {}",
        throat_mean_velocity_ratio_2d,
        throat_area_ratio
    );
    assert!(
        ratio_2d_to_1d > 0.25 && ratio_2d_to_1d < 0.75,
        "2D/1D venturi loss-coefficient ratio {ratio_2d_to_1d} falls outside the calibrated envelope"
    );
}

#[test]
fn cross_fidelity_venturi_three_dimensional_loss_coefficient() {
    let (blueprint, geometry) = calibration_case();
    let (total_loss_coeff_1d, _) = one_dimensional_loss_coefficient(&blueprint, geometry);

    let l_inlet_3d = geometry.inlet_length_m * 0.5;
    let l_converge_3d = geometry.inlet_length_m * 0.5;
    let l_diverge_3d = geometry.diffuser_length_m * 0.5;
    let l_outlet_3d = geometry.diffuser_length_m * 0.5;
    let builder_3d = VenturiMeshBuilder::new(
        geometry.width_inlet_m,
        geometry.width_throat_m,
        l_inlet_3d,
        l_converge_3d,
        geometry.throat_length_m,
        l_diverge_3d,
        l_outlet_3d,
    );
    let config_3d = VenturiConfig3D {
        inlet_flow_rate: FLOW_RATE_M3_S,
        inlet_pressure: 100.0,
        outlet_pressure: 0.0,
        max_nonlinear_iterations: 15,
        nonlinear_tolerance: 1e-4,
        resolution: (30, 5),
        circular: false,
        rect_height: Some(geometry.height_m),
    };
    let solution_3d = VenturiSolver3D::new(builder_3d, config_3d)
        .solve(blood_fluid())
        .expect("3D venturi solve");
    let throat_area_ratio = geometry.width_inlet_m / geometry.width_throat_m;
    let total_loss_coeff_3d = -solution_3d.cp_recovery * throat_area_ratio.powi(2);

    assert!(total_loss_coeff_1d.is_finite() && total_loss_coeff_1d > 0.0);
    assert!(total_loss_coeff_3d.is_finite() && total_loss_coeff_3d > 0.0);

    let ratio_3d_to_1d = total_loss_coeff_3d / total_loss_coeff_1d;
    println!(
        "CALIBRATION 1D/3D: 1D={}, 3D={}, ratio={}",
        total_loss_coeff_1d, total_loss_coeff_3d, ratio_3d_to_1d
    );

    assert!(
        solution_3d.u_throat / solution_3d.u_inlet.max(1e-30) > 1.5,
        "3D venturi throat velocity must exceed inlet velocity by a meaningful margin"
    );
    assert!(
        solution_3d.mass_error.abs() < 0.50,
        "3D venturi mass error {} exceeds 50%",
        solution_3d.mass_error
    );

    // Envelope derivation (recalibrated 2026-07-31): the 1D coefficient comes
    // from literature venturi loss correlations; the exact-reduction FEM with
    // strict Picard tolerances (63e49604) and the no-slip inlet/wall rim
    // agrees with it to 3.3% (319.39 vs 308.70), where the previous
    // inexact-Picard solver (100x-loosened early tolerances) under-predicted
    // viscous loss by 56% (3D = 139.82) — an under-converged velocity field
    // is smoother and shears less. Cross-fidelity agreement near 1 is the
    // correct expectation for two faithful models; +/-30% bounds the coarse
    // (30, 5)-resolution discretization error while still failing on a
    // relapse to the under-converged regime (0.44) or a blowup.
    assert!(
        ratio_3d_to_1d > 0.70 && ratio_3d_to_1d < 1.30,
        "3D/1D venturi loss-coefficient ratio {ratio_3d_to_1d} falls outside the calibrated envelope"
    );
}
