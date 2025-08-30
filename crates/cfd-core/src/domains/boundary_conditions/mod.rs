//! Boundary conditions domain - Physical constraints and boundary condition management.
//!
//! This module encapsulates boundary condition knowledge following DDD principles.
//! It provides abstractions for different boundary condition types and their application.

pub mod applicator;
pub mod error;
pub mod geometry;
pub mod ghost_cells;
pub mod manager;
pub mod specification;
pub mod time_dependent;
pub mod types;

pub use applicator::BoundaryConditionApplicator;
pub use geometry::{BoundaryGeometry, BoundaryRegion};
pub use manager::BoundaryConditionManager;
pub use specification::BoundaryConditionSpec;
pub use time_dependent::{TimeDependentSpec, TimeFunctionType};
pub use types::{DirichletApplicator, NeumannApplicator, RobinApplicator};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::boundary::BoundaryCondition;

    #[test]
    fn test_dirichlet_boundary_application() {
        let mut field = vec![0.0, 0.0, 0.0, 0.0, 0.0];
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Dirichlet { value: 5.0 },
            "west".to_string(),
        );

        let applicator = DirichletApplicator::new();
        applicator.apply(&mut field, &spec, 0.0).unwrap();

        // Check that the boundary value was applied
        assert_eq!(field[0], 5.0, "West boundary should be set to 5.0");
        assert_eq!(field[1], 0.0, "Interior should remain unchanged");
    }

    #[test]
    fn test_neumann_boundary_application() {
        let mut field = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Neumann { gradient: 1.0 },
            "east".to_string(),
        );

        let applicator = NeumannApplicator::new();
        applicator.apply(&mut field, &spec, 0.0).unwrap();

        // Check that gradient is applied: u[n-1] = u[n-2] + gradient
        assert_eq!(field[4], 5.0, "East boundary should be u[3] + 1.0");
    }

    #[test]
    fn test_robin_boundary_application() {
        let mut field = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let spec = BoundaryConditionSpec::new(
            BoundaryCondition::Robin {
                alpha: 1.0,
                beta: 1.0,
                gamma: 2.0,
            },
            "west".to_string(),
        );

        let applicator = RobinApplicator::new();
        applicator.apply(&mut field, &spec, 0.0).unwrap();

        // Check Robin condition: u[0] = (gamma + beta*u[1]/dx) / (alpha + beta/dx)
        // With dx=1: u[0] = (2.0 + 1.0*1.0) / (1.0 + 1.0) = 3.0/2.0 = 1.5
        assert_eq!(
            field[0], 1.5,
            "West boundary should satisfy Robin condition"
        );
    }
}
