//! Two-way and three-way branch junction models with full conservation equations.
//!
//! This module implements the junction flow models governing pressure and flow distribution
//! at branching points in vascular networks. All equations are derived from first principles
//! with references to literature validation.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`two_way_junction`] | `TwoWayBranchJunction<T>` — bifurcation solver |
//! | [`two_way_solution`] | `TwoWayBranchSolution<T>` — solution type |
//! | [`three_way_junction`] | `ThreeWayBranchJunction<T>` + `ThreeWayBranchSolution<T>` |

mod pressure_balance;
pub mod three_way_junction;
pub mod two_way_junction;
pub mod two_way_solution;

pub use three_way_junction::{ThreeWayBranchJunction, ThreeWayBranchSolution};
pub use two_way_junction::TwoWayBranchJunction;
pub use two_way_solution::TwoWayBranchSolution;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::channel::{Channel, ChannelGeometry};
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::blood::CassonBlood;

    #[test]
    fn test_two_way_branch_mass_conservation() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6));

        let branch = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();
        let solution = branch
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert_relative_eq!(solution.q_1 + solution.q_2, 1.0e-6, epsilon = 1e-10);
    }

    #[test]
    fn test_two_way_branch_blood_flow() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 100.0e-6, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6));

        let branch = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();
        let solution = branch
            .solve(blood, 1.0e-8, 100.0, 310.15, 101325.0)
            .unwrap();

        assert!(solution.mu_1 > 0.0);
        assert!(solution.mu_2 > 0.0);
        assert!(solution.gamma_1 > 0.0);
        assert!(solution.gamma_2 > 0.0);
    }

    #[test]
    fn test_murrary_law() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6));

        let branch = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let deviation = branch.murray_law_deviation();
        assert!(deviation < 0.2);
    }

    #[test]
    fn test_three_way_branch_mass_conservation() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 120.0e-6, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 90.0e-6, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 70.0e-6, 1e-6));

        let branch = ThreeWayBranchJunction::new(parent, d1, d2, d3, (0.4, 0.35, 0.25));
        let blood = CassonBlood::<f64>::normal_blood();
        let solution = branch
            .solve(blood, 9.0e-9, 120.0, 310.15, 101325.0)
            .unwrap();

        assert_relative_eq!(
            solution.q_1 + solution.q_2 + solution.q_3,
            solution.q_parent,
            epsilon = 1e-10
        );
        assert!(solution.p_1 <= solution.p_parent);
        assert!(solution.p_2 <= solution.p_parent);
        assert!(solution.p_3 <= solution.p_parent);
    }

    #[test]
    fn test_three_way_murray_law_extension() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let daughter_diameter = 2.0e-3 / 3.0_f64.cbrt();
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(
            1.0e-2,
            daughter_diameter,
            1e-6,
        ));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(
            1.0e-2,
            daughter_diameter,
            1e-6,
        ));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(
            1.0e-2,
            daughter_diameter,
            1e-6,
        ));

        let branch =
            ThreeWayBranchJunction::new(parent, d1, d2, d3, (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
        assert!(branch.murray_law_deviation() < 1e-12);
    }

    #[test]
    fn test_two_way_symmetric_pressure_continuity() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6);
        let d1 = Channel::new(d1_geom.clone());
        let d2 = Channel::new(d1_geom);

        let branch = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();
        let solution = branch
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert_relative_eq!(solution.dp_1, solution.dp_2, epsilon = 1e-10);
        assert_relative_eq!(solution.p_1, solution.p_2, epsilon = 1e-10);
        assert!(solution.junction_pressure_error < 1e-12);
    }

    #[test]
    fn test_two_way_murray_optimal_flow_split() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let r1 = 1.6e-3;
        let r2 = 1.2e-3;
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, r1, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, r2, 1e-6));

        let optimal_split = r1.powi(3) / (r1.powi(3) + r2.powi(3));
        let branch = TwoWayBranchJunction::new(parent, d1, d2, optimal_split);
        let blood = CassonBlood::<f64>::normal_blood();
        let solution = branch
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert!(
            solution.junction_pressure_error < 0.15,
            "Error: {}",
            solution.junction_pressure_error
        );
    }

    #[test]
    fn test_two_way_pressure_balanced_solution_is_seed_independent() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.7e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.1e-3, 1e-6));

        let branch_low_seed =
            TwoWayBranchJunction::new(parent.clone(), d1.clone(), d2.clone(), 0.2);
        let branch_high_seed = TwoWayBranchJunction::new(parent, d1, d2, 0.8);
        let blood = CassonBlood::<f64>::normal_blood();

        let low_seed_solution = branch_low_seed
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();
        let high_seed_solution = branch_high_seed
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert!(low_seed_solution.junction_pressure_error < 1e-8);
        assert!(high_seed_solution.junction_pressure_error < 1e-8);
        assert_relative_eq!(
            low_seed_solution.q_1,
            high_seed_solution.q_1,
            epsilon = 1e-16
        );
        assert_relative_eq!(
            low_seed_solution.p_1,
            high_seed_solution.p_1,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_two_way_prescribed_split_keeps_requested_ratio() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.7e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.1e-3, 1e-6));
        let branch = TwoWayBranchJunction::new(parent, d1, d2, 0.2);
        let blood = CassonBlood::<f64>::normal_blood();

        let prescribed = branch
            .solve_with_prescribed_split(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();
        let balanced = branch
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert_relative_eq!(prescribed.q_1, 2.0e-7, epsilon = 1e-18);
        assert!(balanced.junction_pressure_error < prescribed.junction_pressure_error);
    }

    #[test]
    fn test_three_way_flow_split_sums_to_one() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.0e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.0e-3, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.0e-3, 1e-6));

        let splits = [(0.2, 0.3, 0.5), (0.1, 0.1, 0.8), (0.3333, 0.3333, 0.3334)];
        let blood = CassonBlood::<f64>::normal_blood();
        let q_parent = 1.0e-6;

        for (s1, s2, s3) in splits {
            let branch = ThreeWayBranchJunction::new(
                parent.clone(),
                d1.clone(),
                d2.clone(),
                d3.clone(),
                (s1, s2, s3),
            );
            let solution = branch
                .solve_with_prescribed_split(blood, q_parent, 1000.0, 310.15, 101325.0)
                .unwrap();
            let q_sum = solution.q_1 + solution.q_2 + solution.q_3;
            assert_relative_eq!(q_sum, q_parent, epsilon = 1e-10);
            assert!(solution.mass_conservation_error < 1e-10);
        }
    }

    #[test]
    fn test_three_way_pressure_balanced_solution_is_seed_independent() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.2e-3, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 0.9e-3, 1e-6));

        let low_seed = ThreeWayBranchJunction::new(
            parent.clone(),
            d1.clone(),
            d2.clone(),
            d3.clone(),
            (0.6, 0.25, 0.15),
        );
        let high_seed = ThreeWayBranchJunction::new(parent, d1, d2, d3, (0.2, 0.3, 0.5));
        let blood = CassonBlood::<f64>::normal_blood();

        let low_seed_solution = low_seed
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();
        let high_seed_solution = high_seed
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert!(low_seed_solution.junction_pressure_error < 1e-8);
        assert!(high_seed_solution.junction_pressure_error < 1e-8);
        assert_relative_eq!(
            low_seed_solution.q_1,
            high_seed_solution.q_1,
            epsilon = 1e-16
        );
        assert_relative_eq!(
            low_seed_solution.q_2,
            high_seed_solution.q_2,
            epsilon = 1e-16
        );
        assert_relative_eq!(
            low_seed_solution.p_3,
            high_seed_solution.p_3,
            epsilon = 1e-10
        );
    }

    #[test]
    fn test_three_way_prescribed_split_keeps_requested_ratios() {
        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.5e-3, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 1.2e-3, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-2, 0.9e-3, 1e-6));
        let branch = ThreeWayBranchJunction::new(parent, d1, d2, d3, (0.2, 0.3, 0.5));
        let blood = CassonBlood::<f64>::normal_blood();

        let prescribed = branch
            .solve_with_prescribed_split(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();
        let balanced = branch
            .solve(blood, 1.0e-6, 1000.0, 310.15, 101325.0)
            .unwrap();

        assert_relative_eq!(prescribed.q_1, 2.0e-7, epsilon = 1e-18);
        assert_relative_eq!(prescribed.q_2, 3.0e-7, epsilon = 1e-18);
        assert_relative_eq!(prescribed.q_3, 5.0e-7, epsilon = 1e-18);
        assert!(balanced.junction_pressure_error < prescribed.junction_pressure_error);
    }
}
