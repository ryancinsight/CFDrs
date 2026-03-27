//! Tests for channel solver decomposed modules.

use crate::domain::channel::flow::{Channel, FlowRegime};
use crate::domain::channel::geometry::ChannelGeometry;
use approx::assert_relative_eq;

#[test]
fn test_circular_shape_factor_is_64() {
    let geom = ChannelGeometry::circular(0.01_f64, 0.001, 0.0);
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
    let mut ch = Channel::new(geom);
    ch.flow_state.reynolds_number = Some(100.0);
    ch.flow_state.flow_regime = FlowRegime::Laminar;
    let r = ch.laminar_resistance(&fluid).unwrap();

    // R = Po * mu * L / (2 * A * Dh^2)
    // For circular: Po = 64, A = π D²/4, Dh = D
    // R = 128 μ L / (π D⁴)
    let d = 0.001_f64;
    let l = 0.01_f64;
    let mu = cfd_core::physics::fluid::ConstantFluid::dynamic_viscosity(&fluid);
    let expected = 128.0 * mu * l / (std::f64::consts::PI * d.powi(4));
    assert_relative_eq!(r, expected, epsilon = 1e-6);
}

#[test]
fn test_flow_regime_stokes_below_1() {
    let regime = FlowRegime::from_reynolds_number(0.5_f64);
    assert_eq!(regime, FlowRegime::Stokes);
}

#[test]
fn test_flow_regime_laminar_below_2300() {
    let regime = FlowRegime::from_reynolds_number(500.0_f64);
    assert_eq!(regime, FlowRegime::Laminar);
}

#[test]
fn test_flow_regime_transitional() {
    let regime = FlowRegime::from_reynolds_number(3000.0_f64);
    assert_eq!(regime, FlowRegime::Transitional);
}

#[test]
fn test_flow_regime_turbulent_above_4000() {
    let regime = FlowRegime::from_reynolds_number(5000.0_f64);
    assert_eq!(regime, FlowRegime::Turbulent);
}

#[test]
fn test_square_channel_shape_factor_approx_56_9() {
    // Square channel (AR=1): Po should be ~56.9 per Shah-London
    let geom = ChannelGeometry::rectangular(0.01_f64, 0.001, 0.001, 0.0);
    let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().unwrap();
    let mut ch = Channel::new(geom);
    ch.flow_state.reynolds_number = Some(100.0);
    ch.flow_state.flow_regime = FlowRegime::Laminar;
    let r = ch.laminar_resistance(&fluid).unwrap();

    // Extract Po from measured resistance: R = Po · μ · L / (2·A·D_h²)
    let w = 0.001_f64;
    let h = 0.001_f64;
    let a = w * h;
    let dh = 4.0 * a / (2.0 * (w + h));
    let mu = cfd_core::physics::fluid::ConstantFluid::dynamic_viscosity(&fluid);
    let l = 0.01_f64;
    let po = r * 2.0 * a * dh * dh / (mu * l);
    assert_relative_eq!(po, 56.908, epsilon = 0.1);
}
