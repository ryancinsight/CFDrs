//! Cross-fidelity validation: venturi calibration across 1D, 2D, and 3D.
//!
//! This complements the existing directional venturi test with a quantitative
//! pressure-loss calibration around a shared total loss coefficient.

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

fn blood_fluid() -> ConstantPropertyFluid<f64> {
    ConstantPropertyFluid::new("blood".to_string(), RHO, MU, 3617.0, 0.52, 1570.0)
}

fn rectangular_dims(cross_section: &CrossSectionSpec) -> (f64, f64) {
    match cross_section {
        CrossSectionSpec::Rectangular { width_m, height_m } => (*width_m, *height_m),
        other => panic!("Expected rectangular venturi section, got {other:?}"),
    }
}

fn channel<'a>(blueprint: &'a NetworkBlueprint, channel_id: &str) -> &'a cfd_schematics::domain::model::ChannelSpec {
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

    assert!(inlet_pressure.is_finite(), "Missing inlet pressure in 1D venturi trace");
    assert!(outlet_pressure.is_finite(), "Missing outlet pressure in 1D venturi trace");
    (inlet_pressure, outlet_pressure)
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
fn cross_fidelity_venturi_total_loss_coefficient() {
    let mut blueprint = venturi_rect("venturi-calibration", 2.0e-3, 0.5e-3, 0.5e-3, 2.0e-3);
    blueprint
        .metadata
        .get_or_insert_with(Default::default)
        .insert(GeometryAuthoringProvenance::selective_wrapper());

    let inlet_section = channel(&blueprint, "inlet_section");
    let throat_section = channel(&blueprint, "throat_section");
    let diffuser_section = channel(&blueprint, "diffuser_section");
    let (w_inlet, height_m) = rectangular_dims(&inlet_section.cross_section);
    let (w_throat, _) = rectangular_dims(&throat_section.cross_section);

    let flow_rate_m3_s = 5.0e-8_f64;
    let inlet_area_m2 = w_inlet * height_m;
    let inlet_mean_velocity_m_s = flow_rate_m3_s / inlet_area_m2;
    let dynamic_pressure_pa = 0.5 * RHO * inlet_mean_velocity_m_s * inlet_mean_velocity_m_s;

    let trace_1d = solve_reference_trace::<f64>(&blueprint, RHO, MU, flow_rate_m3_s)
        .expect("1D venturi reference trace");
    let (p_inlet_1d, p_outlet_1d) = inlet_outlet_pressures(&trace_1d);
    let total_loss_coeff_1d = (p_inlet_1d - p_outlet_1d) / dynamic_pressure_pa;

    let l_inlet_2d = inlet_section.length_m * 0.5;
    let l_converge_2d = inlet_section.length_m * 0.5;
    let l_throat = throat_section.length_m;
    let l_diverge = diffuser_section.length_m;
    let geometry_2d = VenturiGeometry::new(
        w_inlet,
        w_throat,
        l_inlet_2d,
        l_converge_2d,
        l_throat,
        l_diverge,
        height_m,
    );
    let mut solver_2d = VenturiSolver2D::new(
        geometry_2d,
        BloodModel::Newtonian(MU),
        RHO,
        120,
        48,
    );
    let solution_2d = solver_2d
        .solve(inlet_mean_velocity_m_s)
        .expect("2D venturi solve");
    let total_loss_coeff_2d = -solution_2d.cp_recovery;

    let l_inlet_3d = inlet_section.length_m * 0.5;
    let l_converge_3d = inlet_section.length_m * 0.5;
    let l_diverge_3d = diffuser_section.length_m * 0.5;
    let l_outlet_3d = diffuser_section.length_m * 0.5;
    let builder_3d = VenturiMeshBuilder::new(
        w_inlet,
        w_throat,
        l_inlet_3d,
        l_converge_3d,
        l_throat,
        l_diverge_3d,
        l_outlet_3d,
    );
    let config_3d = VenturiConfig3D {
        inlet_flow_rate: flow_rate_m3_s,
        inlet_pressure: 100.0,
        outlet_pressure: 0.0,
        max_nonlinear_iterations: 15,
        nonlinear_tolerance: 1e-4,
        resolution: (30, 5),
        circular: false,
        rect_height: Some(height_m),
    };
    let solution_3d = VenturiSolver3D::new(builder_3d, config_3d)
        .solve(blood_fluid())
        .expect("3D venturi solve");
    let total_loss_coeff_3d = -solution_3d.cp_recovery;
    let throat_area_ratio = w_inlet / w_throat;

    assert!(total_loss_coeff_1d.is_finite() && total_loss_coeff_1d > 0.0);
    assert!(total_loss_coeff_2d.is_finite() && total_loss_coeff_2d > 0.0);
    assert!(total_loss_coeff_3d.is_finite() && total_loss_coeff_3d > 0.0);

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
        solution_3d.u_throat / solution_3d.u_inlet.max(1e-30) > 1.5,
        "3D venturi throat velocity must exceed inlet velocity by a meaningful margin"
    );
    assert!(
        solution_3d.mass_error.abs() < 0.50,
        "3D venturi mass error {} exceeds 50%",
        solution_3d.mass_error
    );

    let ratio_2d_to_1d = total_loss_coeff_2d / total_loss_coeff_1d;
    let ratio_3d_to_1d = total_loss_coeff_3d / total_loss_coeff_1d;
    let ratio_3d_to_2d = total_loss_coeff_3d / total_loss_coeff_2d;

    assert!(
        ratio_2d_to_1d > 0.25 && ratio_2d_to_1d < 0.75,
        "2D/1D venturi loss-coefficient ratio {ratio_2d_to_1d} falls outside the calibrated envelope"
    );
    assert!(
        ratio_3d_to_1d > 0.15 && ratio_3d_to_1d < 0.60,
        "3D/1D venturi loss-coefficient ratio {ratio_3d_to_1d} falls outside the calibrated envelope"
    );
    assert!(
        ratio_3d_to_2d > 0.50 && ratio_3d_to_2d < 1.25,
        "3D/2D venturi loss-coefficient ratio {ratio_3d_to_2d} falls outside the calibrated envelope"
    );
}