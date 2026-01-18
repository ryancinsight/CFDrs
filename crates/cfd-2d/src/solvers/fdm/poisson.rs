//! Poisson equation solver using Finite Difference Method.
//!
//! Solves the Poisson equation: ∇²φ = f

use cfd_core::error::{Error, Result};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;

use super::config::FdmConfig;
use super::linear_solver::solve_gauss_seidel;
use crate::grid::{Grid2D, StructuredGrid2D};

/// Poisson equation solver
pub struct PoissonSolver<T: RealField + Copy> {
    config: FdmConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> PoissonSolver<T> {
    /// Create new Poisson solver
    pub fn new(config: FdmConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self {
        Self::new(FdmConfig::<T>::default())
    }

    /// Solve Poisson equation on structured grid (Dirichlet only)
    pub fn solve(
        &self,
        grid: &StructuredGrid2D<T>,
        source: &HashMap<(usize, usize), T>,
        boundary_values: &HashMap<(usize, usize), T>,
    ) -> Result<HashMap<(usize, usize), T>> {
        self.solve_with_neumann(grid, source, boundary_values, &HashMap::new())
    }

    /// Solve Poisson equation with mixed Dirichlet and Neumann boundary conditions
    pub fn solve_with_neumann(
        &self,
        grid: &StructuredGrid2D<T>,
        source: &HashMap<(usize, usize), T>,
        dirichlet_boundaries: &HashMap<(usize, usize), T>,
        neumann_boundaries: &HashMap<(usize, usize), T>,
    ) -> Result<HashMap<(usize, usize), T>> {
        let n = grid.nx() * grid.ny();
        let mut matrix_builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::from_element(n, T::zero());

        // Build system matrix and RHS vector
        for (linear_idx, (i, j)) in grid.iter().enumerate() {
            if let Some(boundary_value) = dirichlet_boundaries.get(&(i, j)).copied() {
                // Dirichlet boundary condition: φ = boundary_value
                matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
                rhs[linear_idx] = boundary_value;
            } else if let Some(gradient) = neumann_boundaries.get(&(i, j)).copied() {
                // Neumann boundary condition: ∂φ/∂n = gradient
                self.add_neumann_bc(
                    &mut matrix_builder,
                    &mut rhs,
                    grid,
                    i,
                    j,
                    linear_idx,
                    gradient,
                )?;
            } else {
                // Interior point: discretize Laplacian
                self.add_laplacian_stencil(
                    &mut matrix_builder,
                    &mut rhs,
                    grid,
                    i,
                    j,
                    linear_idx,
                    source,
                )?;
            }
        }

        // Solve the linear system
        let matrix = matrix_builder.build()?;
        let solution = solve_gauss_seidel(&matrix, &rhs, &self.config, "Poisson")?;

        // Convert solution back to grid coordinates
        let mut result = HashMap::new();
        for (linear_idx, (i, j)) in grid.iter().enumerate() {
            result.insert((i, j), solution[linear_idx]);
        }

        Ok(result)
    }

    /// Add Neumann boundary condition
    fn add_neumann_bc(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        grid: &StructuredGrid2D<T>,
        i: usize,
        j: usize,
        linear_idx: usize,
        gradient: T,
    ) -> Result<()> {
        let (dx, dy) = grid.spacing();

        if i == 0 {
            // Left wall, normal -x: T_{0,j} - T_{1,j} = g * dx
            let neighbor_idx = Self::linear_index(grid, i + 1, j);
            matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -T::one())?;
            rhs[linear_idx] = gradient * dx;
        } else if i == grid.nx() - 1 {
            // Right wall, normal +x: T_{nx-1,j} - T_{nx-2,j} = g * dx
            let neighbor_idx = Self::linear_index(grid, i - 1, j);
            matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -T::one())?;
            rhs[linear_idx] = gradient * dx;
        } else if j == 0 {
            // Bottom wall, normal -y: T_{i,0} - T_{i,1} = g * dy
            let neighbor_idx = Self::linear_index(grid, i, j + 1);
            matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -T::one())?;
            rhs[linear_idx] = gradient * dy;
        } else if j == grid.ny() - 1 {
            // Top wall, normal +y: T_{i,ny-1} - T_{i,ny-2} = g * dy
            let neighbor_idx = Self::linear_index(grid, i, j - 1);
            matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -T::one())?;
            rhs[linear_idx] = gradient * dy;
        } else {
            return Err(Error::InvalidConfiguration(
                "Neumann boundary condition applied to interior point".into(),
            ));
        }

        Ok(())
    }

    /// Add 5-point Laplacian stencil for interior points
    fn add_laplacian_stencil(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        grid: &StructuredGrid2D<T>,
        i: usize,
        j: usize,
        linear_idx: usize,
        source: &HashMap<(usize, usize), T>,
    ) -> Result<()> {
        let (dx, dy) = grid.spacing();
        let dx2 = dx * dx;
        let dy2 = dy * dy;

        // Laplacian can include non-uniform spacing
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Correct 5-point stencil for ∇²φ = f:
        // Center: φ_i,j coefficient = -(2/dx² + 2/dy²)
        let center_coeff = -two / dx2 - two / dy2;
        matrix_builder.add_entry(linear_idx, linear_idx, center_coeff)?;

        // Neighbor contributions: +1/dx² and +1/dy²
        // Left neighbor (i-1)
        if i > 0 {
            let neighbor_idx = Self::linear_index(grid, i - 1, j);
            matrix_builder.add_entry(linear_idx, neighbor_idx, T::one() / dx2)?;
        }

        // Right neighbor (i+1)
        if i < grid.nx() - 1 {
            let neighbor_idx = Self::linear_index(grid, i + 1, j);
            matrix_builder.add_entry(linear_idx, neighbor_idx, T::one() / dx2)?;
        }

        // Bottom neighbor (j-1)
        if j > 0 {
            let neighbor_idx = Self::linear_index(grid, i, j - 1);
            matrix_builder.add_entry(linear_idx, neighbor_idx, T::one() / dy2)?;
        }

        // Top neighbor (j+1)
        if j < grid.ny() - 1 {
            let neighbor_idx = Self::linear_index(grid, i, j + 1);
            matrix_builder.add_entry(linear_idx, neighbor_idx, T::one() / dy2)?;
        }

        // Set RHS: source term f
        rhs[linear_idx] = source.get(&(i, j)).copied().unwrap_or_else(T::zero);

        Ok(())
    }

    /// Convert 2D grid indices to linear index
    fn linear_index(grid: &StructuredGrid2D<T>, i: usize, j: usize) -> usize {
        j * grid.nx() + i
    }
}
