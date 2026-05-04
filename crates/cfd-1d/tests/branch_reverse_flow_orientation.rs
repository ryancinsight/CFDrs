use cfd_1d::domain::channel::{Channel, ChannelGeometry};
use cfd_1d::domain::junctions::branching::{ThreeWayBranchJunction, TwoWayBranchJunction};
use cfd_core::physics::fluid::blood::CassonBlood;

fn circular_channel(length: f64, diameter: f64) -> Channel<f64> {
    Channel::new(ChannelGeometry::<f64>::circular(length, diameter, 1.0e-7))
}

fn assert_relative_eq(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    let bound = 64.0 * f64::EPSILON * scale;
    assert!(
        (actual - expected).abs() <= bound,
        "actual={actual}, expected={expected}, bound={bound}"
    );
}

#[test]
fn pressure_balanced_two_way_branch_preserves_reverse_flow_orientation() {
    let junction = TwoWayBranchJunction::new(
        circular_channel(1.0e-3, 120.0e-6),
        circular_channel(1.0e-3, 80.0e-6),
        circular_channel(1.0e-3, 100.0e-6),
        0.5,
    );
    let blood = CassonBlood::<f64>::normal_blood();

    let forward = junction
        .solve(blood, 4.0e-9, 100.0, 310.0, 100.0)
        .expect("forward bifurcation solve");
    let reverse = junction
        .solve(blood, -4.0e-9, 100.0, 310.0, 100.0)
        .expect("reverse bifurcation solve");

    assert_relative_eq(reverse.q_1, -forward.q_1);
    assert_relative_eq(reverse.q_2, -forward.q_2);
    assert_relative_eq(reverse.dp_1, -forward.dp_1);
    assert_relative_eq(reverse.dp_2, -forward.dp_2);
    assert_relative_eq(reverse.gamma_1, forward.gamma_1);
    assert_relative_eq(reverse.gamma_2, forward.gamma_2);
    assert_relative_eq(reverse.mu_1, forward.mu_1);
    assert_relative_eq(reverse.mu_2, forward.mu_2);
    assert!(reverse.p_junction > reverse.p_parent);
    assert!(reverse.mass_conservation_error < 1.0e-12);
}

#[test]
fn prescribed_two_way_branch_preserves_reverse_flow_orientation() {
    let junction = TwoWayBranchJunction::new(
        circular_channel(1.0e-3, 120.0e-6),
        circular_channel(1.0e-3, 80.0e-6),
        circular_channel(1.0e-3, 100.0e-6),
        0.35,
    );
    let blood = CassonBlood::<f64>::normal_blood();

    let solution = junction
        .solve_with_prescribed_split(blood, -4.0e-9, 100.0, 310.0, 100.0)
        .expect("reverse prescribed bifurcation solve");

    assert_relative_eq(solution.q_1, -1.4e-9);
    assert_relative_eq(solution.q_2, -2.6e-9);
    assert!(solution.gamma_1 > 0.0);
    assert!(solution.gamma_2 > 0.0);
    assert!(solution.mass_conservation_error < 1.0e-12);
}

#[test]
fn pressure_balanced_three_way_branch_preserves_reverse_flow_orientation() {
    let junction = ThreeWayBranchJunction::new(
        circular_channel(1.0e-3, 140.0e-6),
        circular_channel(1.0e-3, 90.0e-6),
        circular_channel(1.0e-3, 80.0e-6),
        circular_channel(1.0e-3, 70.0e-6),
        (0.4, 0.35, 0.25),
    );
    let blood = CassonBlood::<f64>::normal_blood();

    let forward = junction
        .solve(blood, 6.0e-9, 100.0, 310.0, 100.0)
        .expect("forward trifurcation solve");
    let reverse = junction
        .solve(blood, -6.0e-9, 100.0, 310.0, 100.0)
        .expect("reverse trifurcation solve");

    assert_relative_eq(reverse.q_1, -forward.q_1);
    assert_relative_eq(reverse.q_2, -forward.q_2);
    assert_relative_eq(reverse.q_3, -forward.q_3);
    assert_relative_eq(reverse.dp_1, -forward.dp_1);
    assert_relative_eq(reverse.dp_2, -forward.dp_2);
    assert_relative_eq(reverse.dp_3, -forward.dp_3);
    assert_relative_eq(reverse.gamma_1, forward.gamma_1);
    assert_relative_eq(reverse.gamma_2, forward.gamma_2);
    assert_relative_eq(reverse.gamma_3, forward.gamma_3);
    assert_relative_eq(reverse.mu_1, forward.mu_1);
    assert_relative_eq(reverse.mu_2, forward.mu_2);
    assert_relative_eq(reverse.mu_3, forward.mu_3);
    assert!(reverse.mass_conservation_error < 1.0e-12);
}

#[test]
fn prescribed_three_way_branch_preserves_reverse_flow_orientation() {
    let junction = ThreeWayBranchJunction::new(
        circular_channel(1.0e-3, 140.0e-6),
        circular_channel(1.0e-3, 90.0e-6),
        circular_channel(1.0e-3, 80.0e-6),
        circular_channel(1.0e-3, 70.0e-6),
        (0.5, 0.3, 0.2),
    );
    let blood = CassonBlood::<f64>::normal_blood();

    let solution = junction
        .solve_with_prescribed_split(blood, -6.0e-9, 100.0, 310.0, 100.0)
        .expect("reverse prescribed trifurcation solve");

    assert_relative_eq(solution.q_1, -3.0e-9);
    assert_relative_eq(solution.q_2, -1.8e-9);
    assert_relative_eq(solution.q_3, -1.2e-9);
    assert!(solution.gamma_1 > 0.0);
    assert!(solution.gamma_2 > 0.0);
    assert!(solution.gamma_3 > 0.0);
    assert!(solution.mass_conservation_error < 1.0e-12);
}
