//! Boundary condition applicator trait

use super::specification::BoundaryConditionSpec;
use crate::boundary::BoundaryCondition;
use nalgebra::RealField;

/// Boundary condition applicator abstraction
pub trait BoundaryConditionApplicator<T: RealField + Copy>: Send + Sync {
    /// Apply boundary condition to field
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
        // Default linear extrapolation
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        two * boundary_value - interior_value
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

        match condition {
            BoundaryCondition::Dirichlet { value } => {
                // Flux based on gradient to boundary value
                normal_gradient * (value - interior_value)
            }
            BoundaryCondition::Neumann { gradient } => gradient,
            BoundaryCondition::Robin { alpha, beta, gamma } => {
                // Robin condition: alpha*u + beta*du/dn = g
                // Flux = -beta/alpha * (g/beta - u)
                if alpha.abs() > T::default_epsilon() {
                    -beta / alpha * interior_value
                } else {
                    T::zero()
                }
            }
            _ => T::zero(),
        }
    }
}
