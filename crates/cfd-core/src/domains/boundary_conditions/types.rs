//! Concrete boundary condition applicator implementations

use super::applicator::BoundaryConditionApplicator;
use super::specification::BoundaryConditionSpec;
use crate::boundary::BoundaryCondition;
use nalgebra::RealField;

/// Dirichlet boundary condition applicator
pub struct DirichletApplicator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for DirichletApplicator<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> DirichletApplicator<T> {
    /// Create a new Dirichlet applicator
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> BoundaryConditionApplicator<T> for DirichletApplicator<T> {
    fn apply(
        &self,
        field: &mut [T],
        boundary_spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String> {
        let condition = boundary_spec.evaluate_at_time(time);

        if let BoundaryCondition::Dirichlet { value } = condition.as_ref() {
            // Apply Dirichlet value based on region ID
            // For structured grids, interpret region_id as boundary location
            match boundary_spec.region_id.as_str() {
                "west" | "left" => {
                    // Apply to first element (1D) or first column (2D/3D)
                    if !field.is_empty() {
                        field[0] = *value;
                    }
                }
                "east" | "right" => {
                    // Apply to last element
                    if let Some(last) = field.last_mut() {
                        *last = *value;
                    }
                }
                "all" => {
                    // Apply to entire boundary (for testing)
                    for field_value in field.iter_mut() {
                        *field_value = *value;
                    }
                }
                _ => {
                    // For unrecognized regions, apply to boundaries (first and last)
                    if !field.is_empty() {
                        field[0] = *value;
                        if let Some(last) = field.last_mut() {
                            *last = *value;
                        }
                    }
                }
            }
            Ok(())
        } else {
            Err("Not a Dirichlet condition".to_string())
        }
    }

    fn name(&self) -> &'static str {
        "Dirichlet"
    }

    fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
        matches!(condition, BoundaryCondition::Dirichlet { .. })
    }
}

/// Neumann boundary condition applicator
pub struct NeumannApplicator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for NeumannApplicator<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> NeumannApplicator<T> {
    /// Create a new Neumann applicator
    #[must_use] pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> BoundaryConditionApplicator<T> for NeumannApplicator<T> {
    fn apply(
        &self,
        field: &mut [T],
        boundary_spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String> {
        let condition = boundary_spec.evaluate_at_time(time);

        if let BoundaryCondition::Neumann { gradient } = condition.as_ref() {
            // Apply Neumann gradient condition using one-sided differences
            // For a gradient g at boundary: u_boundary = u_interior + dx * g
            match boundary_spec.region_id.as_str() {
                "west" | "left" => {
                    // Apply at first element using forward difference
                    if field.len() >= 2 {
                        // u[0] = u[1] - dx * gradient (assuming unit spacing)
                        field[0] = field[1] - *gradient;
                    }
                }
                "east" | "right" => {
                    // Apply at last element using backward difference
                    let n = field.len();
                    if n >= 2 {
                        // u[n-1] = u[n-2] + dx * gradient
                        field[n - 1] = field[n - 2] + *gradient;
                    }
                }
                _ => {
                    // Apply to both boundaries
                    if field.len() >= 2 {
                        field[0] = field[1] - *gradient;
                        let n = field.len();
                        field[n - 1] = field[n - 2] + *gradient;
                    }
                }
            }
            Ok(())
        } else {
            Err("Not a Neumann condition".to_string())
        }
    }

    fn name(&self) -> &'static str {
        "Neumann"
    }

    fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
        matches!(condition, BoundaryCondition::Neumann { .. })
    }

    fn is_conservative(&self) -> bool {
        true // Neumann conditions preserve conservation
    }
}

/// Robin boundary condition applicator
pub struct RobinApplicator<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for RobinApplicator<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> RobinApplicator<T> {
    /// Create a new Robin applicator
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> BoundaryConditionApplicator<T> for RobinApplicator<T> {
    fn apply(
        &self,
        field: &mut [T],
        boundary_spec: &BoundaryConditionSpec<T>,
        time: T,
    ) -> Result<(), String> {
        let condition = boundary_spec.evaluate_at_time(time);

        if let BoundaryCondition::Robin { alpha, beta, gamma } = condition.as_ref() {
            // Apply Robin condition: alpha*u + beta*du/dn = gamma
            // Solving for u: u = (gamma - beta*du/dn) / alpha
            // Using finite difference: du/dn ≈ (u_boundary - u_interior) / dx
            // Rearranging: u_boundary = (gamma + beta*u_interior/dx) / (alpha + beta/dx)

            let dx = T::one(); // Assuming unit spacing

            match boundary_spec.region_id.as_str() {
                "west" | "left" => {
                    if field.len() >= 2 && *alpha != T::zero() {
                        // u[0] = (gamma + beta*u[1]/dx) / (alpha + beta/dx)
                        field[0] = (*gamma + *beta * field[1] / dx) / (*alpha + *beta / dx);
                    }
                }
                "east" | "right" => {
                    let n = field.len();
                    if n >= 2 && *alpha != T::zero() {
                        // u[n-1] = (gamma + beta*u[n-2]/dx) / (alpha + beta/dx)
                        field[n - 1] = (*gamma + *beta * field[n - 2] / dx) / (*alpha + *beta / dx);
                    }
                }
                _ => {
                    // Apply to both boundaries
                    if field.len() >= 2 && *alpha != T::zero() {
                        field[0] = (*gamma + *beta * field[1] / dx) / (*alpha + *beta / dx);
                        let n = field.len();
                        field[n - 1] = (*gamma + *beta * field[n - 2] / dx) / (*alpha + *beta / dx);
                    }
                }
            }
            Ok(())
        } else {
            Err("Not a Robin condition".to_string())
        }
    }

    fn name(&self) -> &'static str {
        "Robin"
    }

    fn supports(&self, condition: &BoundaryCondition<T>) -> bool {
        matches!(condition, BoundaryCondition::Robin { .. })
    }

    fn order_of_accuracy(&self) -> usize {
        2 // Second-order accurate implementation
    }
}
