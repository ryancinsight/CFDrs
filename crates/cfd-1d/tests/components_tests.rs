//! Comprehensive tests for all microfluidic component types.
//!
//! Sprint 8 updates: Micromixer and FlowSensor use first-principles
//! constructors (validated, returns Result) with Idelchik physics model.

use approx::assert_relative_eq;
use cfd_1d::components::{
    CircularChannel, Micropump, Microvalve, OrganCompartment,
    PorousMembrane, RectangularChannel,
};
use cfd_1d::components::mixers::{Micromixer, MixerType};
use cfd_1d::components::sensors::FlowSensor;
use cfd_1d::components::Component;
use cfd_1d::solver::ConvergenceChecker;
use cfd_core::physics::fluid::database::water_20c;

fn water() -> cfd_core::physics::fluid::ConstantPropertyFluid<f64> {
    water_20c::<f64>().unwrap()
}

// ========================  RectangularChannel  ============

#[test]
fn test_rect_channel_resistance_positive_finite() {
    let fluid = water();
    let chan = RectangularChannel::new(0.1, 1e-3, 1e-3, 0.0);
    let r = chan.resistance(&fluid);
    assert!(r > 0.0 && r.is_finite());
}

#[test]
fn test_rect_channel_linear_in_length() {
    let fluid = water();
    let c1 = RectangularChannel::new(0.1, 1e-3, 1e-3, 0.0);
    let c2 = RectangularChannel::new(0.2, 1e-3, 1e-3, 0.0);
    let r1 = c1.resistance(&fluid);
    let r2 = c2.resistance(&fluid);
    assert_relative_eq!(r2 / r1, 2.0, max_relative = 0.001);
}

#[test]
fn test_rect_channel_coefficients_pure_linear() {
    let fluid = water();
    let chan = RectangularChannel::new(0.1, 1e-3, 1e-3, 0.0);
    let (r, k) = chan.coefficients(&fluid);
    assert!(r > 0.0);
    assert_relative_eq!(k, 0.0, epsilon = 1e-30);
}

#[test]
fn test_rect_channel_pressure_drop() {
    let fluid = water();
    let chan = RectangularChannel::new(0.1, 1e-3, 1e-3, 0.0);
    let r = chan.resistance(&fluid);
    let q = 1e-9;
    assert_relative_eq!(chan.pressure_drop(q, &fluid), r * q, epsilon = 1e-30);
}

// ========================  CircularChannel  ===============

#[test]
fn test_circ_channel_resistance_positive_finite() {
    let fluid = water();
    let chan = CircularChannel::new(0.1, 1e-3, 0.0);
    let r = chan.resistance(&fluid);
    assert!(r > 0.0 && r.is_finite());
}

#[test]
fn test_circ_channel_inverse_d4_scaling() {
    let fluid = water();
    let r1 = CircularChannel::new(0.1, 1e-3, 0.0).resistance(&fluid);
    let r2 = CircularChannel::new(0.1, 2e-3, 0.0).resistance(&fluid);
    // R ∝ D^(-4): R(D) / R(2D) = 16
    assert_relative_eq!(r1 / r2, 16.0, max_relative = 1e-6);
}

// ========================  Microvalve  ====================

#[test]
fn test_valve_resistance_open() {
    let fluid = water();
    let valve = Microvalve::new(0.01_f64);
    assert_relative_eq!(valve.resistance(&fluid), 0.0, epsilon = 1e-30);
}

#[test]
fn test_valve_resistance_closed() {
    let fluid = water();
    let mut valve = Microvalve::new(0.01_f64);
    valve.set_parameter("opening", 0.0).unwrap();
    assert!(valve.resistance(&fluid) > 1e10);
}

#[test]
fn test_valve_pressure_drop_at_full_open() {
    let fluid = water();
    let cv = 1e-6_f64; // 1 ul/s/√Pa
    let valve = Microvalve::new(cv);
    let q = 1e-9; // 1 nL/s
    let dp = valve.pressure_drop(q, &fluid);
    // ΔP = k * Q^2 = Q^2 / Cv^2
    let expected = q * q / (cv * cv);
    assert_relative_eq!(dp, expected, epsilon = 1e-30);
}

// ========================  Micropump  =====================

#[test]
fn test_pump_resistance_is_zero() {
    let fluid = water();
    let pump = Micropump::new(1e-6_f64, 1000.0);
    // Pumps must not contribute negative (or any non-zero) passive resistance
    assert_relative_eq!(pump.resistance(&fluid), 0.0, epsilon = 1e-30);
}

#[test]
fn test_pump_is_active() {
    let pump = Micropump::new(1e-6_f64, 1000.0);
    assert!(pump.is_active());
}

#[test]
fn test_pump_max_pressure_and_flow_rate_parameters() {
    let pump = Micropump::new(1e-6_f64, 1000.0);
    assert_relative_eq!(pump.max_flow_rate, 1e-6, epsilon = 1e-30);
    assert_relative_eq!(pump.max_pressure, 1000.0, epsilon = 1e-30);
}

// ========================  FlowSensor  ====================
// (G5: validated constructor, theorem docs)

/// Negative insertion resistance must be rejected.
#[test]
fn test_flow_sensor_rejects_negative_resistance() {
    assert!(FlowSensor::<f64>::new(-1.0, 1e-6).is_err());
}

/// Zero range must be rejected.
#[test]
fn test_flow_sensor_rejects_zero_range() {
    assert!(FlowSensor::<f64>::new(0.0, 0.0).is_err());
}

/// Ideal (zero-insertion) sensor: R = 0.
#[test]
fn test_flow_sensor_ideal_zero_resistance() {
    let fluid = water();
    let sensor = FlowSensor::<f64>::new(0.0, 1e-3).unwrap();
    assert_eq!(sensor.resistance(&fluid), 0.0);
}

/// Insertion resistance returned exactly.
#[test]
fn test_sensor_resistance_pass_through() {
    let fluid = water();
    let sensor = FlowSensor::<f64>::new(1e6, 1e-6).unwrap();
    assert_relative_eq!(sensor.resistance(&fluid), 1e6, epsilon = 1e-30);
}

/// set_parameter rejects negative resistance.
#[test]
fn test_sensor_parameter_update_validates() {
    let mut sensor = FlowSensor::<f64>::new(1e6, 1e-6).unwrap();
    assert!(sensor.set_parameter("resistance", -1.0).is_err());
    sensor.set_parameter("resistance", 2e6).unwrap();
    assert_relative_eq!(sensor.resistance, 2e6, epsilon = 1e-30);
}

/// is_overrange detects flow above measurement range.
#[test]
fn test_flow_sensor_overrange_detection() {
    let sensor = FlowSensor::<f64>::new(0.0, 1e-6).unwrap();
    assert!(!sensor.is_overrange(0.5e-6));
    assert!(sensor.is_overrange(2e-6));
}

// ========================  Micromixer  ====================
// (G1/G6: Idelchik physics model, theorem docs)

/// T-junction resistance must be strictly positive.
///
/// **Theorem**: R_mixer = R_straight + R_minor > 0 for μ > 0, L > 0, D > 0.
#[test]
fn test_mixer_t_junction_resistance_positive() {
    let fluid = water();
    let mixer = Micromixer::<f64>::new(MixerType::TJunction, 200e-6, 5e-3, 1).unwrap();
    let r = mixer.resistance(&fluid);
    assert!(r > 0.0 && r.is_finite(), "R = {r}");
}

/// Serpentine mixer: more bends → strictly higher resistance.
///
/// **Invariant**: ∂R/∂n_bends > 0.
#[test]
fn test_mixer_serpentine_increases_with_bends() {
    let fluid = water();
    let r2 = Micromixer::<f64>::new(MixerType::Serpentine, 200e-6, 5e-3, 2).unwrap().resistance(&fluid);
    let r8 = Micromixer::<f64>::new(MixerType::Serpentine, 200e-6, 5e-3, 8).unwrap().resistance(&fluid);
    assert!(r8 > r2, "R(8 bends)={r8} must exceed R(2 bends)={r2}");
}

/// Resistance increases monotonically with length.
#[test]
fn test_mixer_resistance_increases_with_length() {
    let fluid = water();
    let r1 = Micromixer::<f64>::new(MixerType::TJunction, 200e-6, 1e-3, 1).unwrap().resistance(&fluid);
    let r10 = Micromixer::<f64>::new(MixerType::TJunction, 200e-6, 10e-3, 1).unwrap().resistance(&fluid);
    assert!(r10 > r1, "R(10mm)={r10} must exceed R(1mm)={r1}");
}

/// Larger hydraulic diameter → lower resistance.
#[test]
fn test_mixer_resistance_decreases_with_diameter() {
    let fluid = water();
    let r_small = Micromixer::<f64>::new(MixerType::TJunction, 50e-6, 5e-3, 1).unwrap().resistance(&fluid);
    let r_large = Micromixer::<f64>::new(MixerType::TJunction, 500e-6, 5e-3, 1).unwrap().resistance(&fluid);
    assert!(r_large < r_small, "R(500μm)={r_large} must be < R(50μm)={r_small}");
}

/// Invalid geometry (D=0) must be rejected.
#[test]
fn test_mixer_rejects_zero_diameter() {
    assert!(Micromixer::<f64>::new(MixerType::TJunction, 0.0, 5e-3, 1).is_err());
}

/// Invalid geometry (L=0) must be rejected.
#[test]
fn test_mixer_rejects_zero_length() {
    assert!(Micromixer::<f64>::new(MixerType::TJunction, 200e-6, 0.0, 1).is_err());
}

/// T-junction K_loss=1.5 > Y-junction K_loss=0.9 → R_T > R_Y for same geometry.
#[test]
fn test_mixer_t_junction_higher_than_y_junction() {
    let fluid = water();
    let rt = Micromixer::<f64>::new(MixerType::TJunction, 200e-6, 5e-3, 1).unwrap().resistance(&fluid);
    let ry = Micromixer::<f64>::new(MixerType::YJunction, 200e-6, 5e-3, 1).unwrap().resistance(&fluid);
    assert!(rt > ry, "R_T={rt} must exceed R_Y={ry} (K_T=1.5 > K_Y=0.9)");
}

/// Efficiency clamping above 1.
#[test]
fn test_mixer_efficiency_clamped_above_one() {
    let mut mixer = Micromixer::<f64>::new(MixerType::Serpentine, 200e-6, 5e-3, 4).unwrap();
    mixer.set_parameter("efficiency", 1.5).unwrap();
    assert_relative_eq!(mixer.efficiency, 1.0, epsilon = 1e-15);
}

/// Efficiency clamping below zero.
#[test]
fn test_mixer_efficiency_clamped_below_zero() {
    let mut mixer = Micromixer::<f64>::new(MixerType::Herringbone, 200e-6, 5e-3, 4).unwrap();
    mixer.set_parameter("efficiency", -0.1).unwrap();
    assert_relative_eq!(mixer.efficiency, 0.0, epsilon = 1e-15);
}

// ========================  PorousMembrane  ================

#[test]
fn test_porous_membrane_resistance_positive() {
    let fluid = water();
    let membrane = PorousMembrane::new(1e-4, 1e-3, 1e-3, 5e-7, 0.3);
    let r = membrane.resistance(&fluid);
    assert!(r > 0.0 && r.is_finite());
}

#[test]
fn test_porous_membrane_larger_pore_less_resistance() {
    let fluid = water();
    let m_small = PorousMembrane::new(1e-4, 1e-3, 1e-3, 1e-7, 0.3);
    let m_large = PorousMembrane::new(1e-4, 1e-3, 1e-3, 5e-7, 0.3);
    // Larger pore radius → lower R (HP pore: R ∝ 1/r^2 per unit porosity)
    assert!(m_large.resistance(&fluid) < m_small.resistance(&fluid));
}

// ========================  OrganCompartment  ==============

#[test]
fn test_organ_compartment_resistance_exact() {
    let fluid = water();
    let compartment = OrganCompartment::new(0.01, 2e-3, 1e-3, 1.5e8);
    assert_relative_eq!(compartment.resistance(&fluid), 1.5e8, epsilon = 1e-30);
}

#[test]
fn test_organ_compartment_volume() {
    let compartment = OrganCompartment::new(0.01, 2e-3, 1e-3, 1.5e8_f64);
    let expected = 0.01 * 2e-3 * 1e-3;
    assert_relative_eq!(compartment.volume().unwrap(), expected, epsilon = 1e-25);
}

// ========================  ConvergenceChecker  ============
// (G7)

/// max_iterations_reached returns true at and beyond the threshold.
///
/// **Theorem**: `current_iteration >= max_iterations` ⟹ `true`.
#[test]
fn test_convergence_checker_max_iterations_boundary() {
    let checker = ConvergenceChecker::<f64>::new(1e-6);
    // Default max_iterations = 1000
    assert!(!checker.max_iterations_reached(999));
    assert!(checker.max_iterations_reached(1000));
    assert!(checker.max_iterations_reached(1001));
}

/// Not reached at iteration 0 for any positive max.
#[test]
fn test_convergence_checker_not_reached_at_zero() {
    let checker = ConvergenceChecker::<f64>::new(1e-12);
    assert!(!checker.max_iterations_reached(0));
}
