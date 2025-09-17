//! Ghost cell implementations for high-order boundary condition treatment
//!
//! References:
//! - Blazek (2015) "Computational Fluid Dynamics: Principles and Applications"
//! - Morinishi et al. (1998) "Fully Conservative Higher Order Finite Difference Schemes"

use super::error::BoundaryError;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// Ghost cell calculator for maintaining high-order accuracy at boundaries
pub struct GhostCellCalculator<T: RealField + Copy> {
    /// Stencil order (1 = first-order, 2 = second-order, etc.)
    order: usize,
    /// Number of ghost cells required
    n_ghost: usize,
    /// Phantom data for type parameter
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> GhostCellCalculator<T> {
    /// Create ghost cell calculator for given stencil order
    #[must_use] pub fn new(order: usize) -> Self {
        // Number of ghost cells = (order + 1) / 2 for centered schemes
        // For upwind schemes, may need more on upwind side
        let n_ghost = order.div_ceil(2);
        Self {
            order,
            n_ghost,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Apply Dirichlet boundary condition with ghost cells
    ///
    /// For second-order accuracy at boundary:
    /// - Linear extrapolation: g₁ = 2*b - i₁
    /// - Quadratic extrapolation: g₁ = 3*b - 3*i₁ + i₂
    /// where g = ghost, b = boundary value, i = interior values
    pub fn apply_dirichlet(
        &self,
        boundary_value: T,
        interior_values: &[T],
        ghost_values: &mut [T],
    ) -> Result<(), BoundaryError> {
        if interior_values.len() < self.order {
            return Err(BoundaryError::insufficient_stencil(
                self.order,
                self.order,
                interior_values.len(),
            ));
        }

        match self.order {
            1 => {
                // First-order: ghost = boundary_value
                ghost_values[0] = boundary_value;
            }
            2 => {
                // Second-order: linear extrapolation
                // g₀ = 2*b - i₀
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                ghost_values[0] = two * boundary_value - interior_values[0];
            }
            3 | 4 => {
                // Third/Fourth-order: quadratic extrapolation
                // g₀ = 3*b - 3*i₀ + i₁
                let three = T::from_f64(3.0).unwrap_or_else(T::one);
                ghost_values[0] =
                    three * boundary_value - three * interior_values[0] + interior_values[1];

                if self.n_ghost > 1 {
                    // Second ghost cell for fourth-order
                    // g₁ = 5*b - 10*i₀ + 10*i₁ - 5*i₂ + i₃
                    if interior_values.len() >= 4 {
                        let five = T::from_f64(5.0).unwrap_or_else(T::one);
                        let ten = T::from_f64(10.0).unwrap_or_else(T::one);
                        ghost_values[1] = five * boundary_value - ten * interior_values[0]
                            + ten * interior_values[1]
                            - five * interior_values[2]
                            + interior_values[3];
                    }
                }
            }
            _ => {
                return Err(BoundaryError::unsupported_order(self.order));
            }
        }

        Ok(())
    }

    /// Apply Neumann boundary condition with ghost cells
    ///
    /// For prescribed gradient ∂u/∂n = g at boundary:
    /// - Second-order: g₀ = i₀ - 2*Δx*gradient
    /// - Fourth-order: g₀ = i₀ - (8/3)*Δx*gradient + (2/3)*i₁
    pub fn apply_neumann(
        &self,
        gradient: T,
        dx: T,
        interior_values: &[T],
        ghost_values: &mut [T],
    ) -> Result<(), BoundaryError> {
        if interior_values.len() < self.order {
            return Err(BoundaryError::insufficient_stencil(
                self.order,
                self.order,
                interior_values.len(),
            ));
        }

        match self.order {
            1 | 2 => {
                // First/Second-order: simple reflection
                // g₀ = i₀ - 2*Δx*gradient
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                ghost_values[0] = interior_values[0] - two * dx * gradient;
            }
            3 | 4 => {
                // Third/Fourth-order: higher-order extrapolation
                // g₀ = i₀ - (8/3)*Δx*gradient + (2/3)*i₁
                let eight_thirds = T::from_f64(8.0 / 3.0).unwrap_or_else(T::one);
                let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(T::one);
                ghost_values[0] = interior_values[0] - eight_thirds * dx * gradient
                    + two_thirds * interior_values[1];
            }
            _ => {
                return Err(BoundaryError::unsupported_order(self.order));
            }
        }

        Ok(())
    }

    /// Apply Robin (mixed) boundary condition
    /// α*u + β*∂u/∂n = γ at boundary
    ///
    /// # Errors
    /// Returns an error if:
    /// - Interior or ghost values arrays are insufficient for the boundary order
    /// - Coefficient α is zero when β is also zero (degenerate case)
    /// - Numerical conditioning issues with mixed boundary condition
    pub fn apply_robin(
        &self,
        alpha: T,
        beta: T,
        gamma: T,
        dx: T,
        interior_values: &[T],
        ghost_values: &mut [T],
    ) -> Result<(), BoundaryError> {
        if alpha == T::zero() {
            // Pure Neumann
            let gradient = gamma / beta;
            self.apply_neumann(gradient, dx, interior_values, ghost_values)
        } else if beta == T::zero() {
            // Pure Dirichlet
            let boundary_value = gamma / alpha;
            self.apply_dirichlet(boundary_value, interior_values, ghost_values)
        } else {
            // Mixed condition
            // For second-order: g₀ = (i₀*(β - α*dx) + 2*γ*dx) / (β + α*dx)
            let two = T::from_f64(2.0).unwrap_or_else(T::one);
            let numerator = interior_values[0] * (beta - alpha * dx) + two * gamma * dx;
            let denominator = beta + alpha * dx;

            if denominator.abs() < T::from_f64(1e-10).unwrap_or_else(T::zero) {
                return Err(BoundaryError::robin_singularity(
                    denominator.to_f64().unwrap_or(0.0),
                ));
            }

            ghost_values[0] = numerator / denominator;
            Ok(())
        }
    }

    /// Get required number of ghost cells for this order
    #[must_use] pub fn ghost_cells_required(&self) -> usize {
        self.n_ghost
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dirichlet_ghost_cells() {
        let calc = GhostCellCalculator::<f64>::new(2);
        let mut ghost = vec![0.0; 1];
        let interior = vec![1.0, 2.0, 3.0];

        // Apply Dirichlet BC with value 0
        calc.apply_dirichlet(0.0, &interior, &mut ghost).unwrap();

        // Second-order: g₀ = 2*0 - 1 = -1
        assert!((ghost[0] + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_neumann_ghost_cells() {
        let calc = GhostCellCalculator::<f64>::new(2);
        let mut ghost = vec![0.0; 1];
        let interior = vec![1.0, 2.0, 3.0];

        // Apply Neumann BC with zero gradient
        calc.apply_neumann(0.0, 0.1, &interior, &mut ghost).unwrap();

        // Zero gradient: ghost should equal first interior
        assert!((ghost[0] - 1.0).abs() < 1e-10);
    }
}
