//! Boundary condition applicator trait

use super::specification::BoundaryConditionSpec;
use crate::physics::boundary::BoundaryCondition;
use eunomia::FloatElement;
use eunomia::RealField;

/// Boundary condition applicator abstraction
pub trait BoundaryConditionApplicator<T: RealField + FloatElement + Copy>: Send + Sync {
    /// Apply boundary condition to field
    ///
    /// # Errors
    /// Returns error if boundary condition application fails or field is incompatible
    fn apply(
        &self,
        field: &mut [T],
        boundary_spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String>;

    /// Get applicator name
    fn name(&self) -> &str;

    /// Check if this applicator supports the given boundary condition
    fn supports(&self, condition: &BoundaryCondition<T>) -> bool;

    /// Get the order of accuracy of the boundary condition implementation
    fn order_of_accuracy(&self) -> usize {
        2 // Default second-order accuracy
    }

    /// Check if the applicator preserves conservation properties
    fn is_conservative(&self) -> bool {
        true // Most boundary conditions should be conservative
    }

    /// Apply boundary condition to a ghost cell (for finite volume methods)
    fn apply_to_ghost_cell(
        &self,
        interior_value: T,
        boundary_value: T,
        distance_to_boundary: T,
    ) -> T {
        // Linear extrapolation accounting for distance
        // Ghost value = boundary_value + (boundary_value - interior_value) * (distance_ghost / distance_interior)
        // For uniform grid where distance_to_boundary is normalized to cell size:
        let two = <T as FloatElement>::from_f64(2.0);
        if distance_to_boundary > T::EPSILON {
            boundary_value + (boundary_value - interior_value) * distance_to_boundary
        } else {
            // Fallback to standard ghost cell when distance is zero (cell-centered)
            two * boundary_value - interior_value
        }
    }

    /// Compute boundary flux (for finite volume methods)
    fn compute_flux(
        &self,
        interior_value: T,
        boundary_spec: &BoundaryConditionSpec<T>,
        normal_gradient: T,
        time: T,
    ) -> T {
        let condition = boundary_spec.evaluate_at_time(time);

        match condition.as_ref() {
            BoundaryCondition::Dirichlet { value, .. } => {
                // Flux based on gradient to boundary value
                normal_gradient * (*value - interior_value)
            }
            BoundaryCondition::Neumann { gradient } => *gradient,
            BoundaryCondition::Robin { alpha, beta, gamma } => {
                // Robin condition: alpha*u + beta*du/dn = gamma
                // Solving for du/dn: du/dn = (gamma - alpha*u)/beta
                if beta.abs() > T::EPSILON {
                    (*gamma - *alpha * interior_value) / *beta
                } else if alpha.abs() > T::EPSILON {
                    // If beta is zero, we have a Dirichlet-like condition: u = gamma/alpha
                    normal_gradient * (*gamma / *alpha - interior_value)
                } else {
                    T::ZERO
                }
            }
            _ => T::ZERO,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::BoundaryConditionApplicator;
    use crate::physics::boundary::{BoundaryCondition, BoundaryConditionSpec};

    struct TestApplicator;

    impl BoundaryConditionApplicator<f64> for TestApplicator {
        fn apply(
            &self,
            _field: &mut [f64],
            _boundary_spec: &BoundaryConditionSpec<f64>,
            _time: f64,
        ) -> Result<(), String> {
            Ok(())
        }

        fn name(&self) -> &'static str {
            "test"
        }

        fn supports(&self, _condition: &BoundaryCondition<f64>) -> bool {
            true
        }
    }

    #[test]
    fn ghost_cell_uses_standard_reflection_at_zero_distance() {
        let applicator = TestApplicator;

        let ghost = applicator.apply_to_ghost_cell(3.0, 10.0, 0.0);

        assert_eq!(ghost, 17.0);
    }

    #[test]
    fn ghost_cell_uses_distance_weighted_extrapolation() {
        let applicator = TestApplicator;

        let ghost = applicator.apply_to_ghost_cell(3.0, 10.0, 0.25);

        assert_eq!(ghost, 11.75);
    }
}
