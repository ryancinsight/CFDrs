//! Advection-diffusion equation solver using Finite Difference Method.
//!
//! Solves: ∂φ/∂t + u·∇φ = α∇²φ + S

use cfd_core::error::Result;
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::collections::HashMap;

use super::config::FdmConfig;
use super::linear_solver::solve_gauss_seidel;
use crate::grid::{Grid2D, StructuredGrid2D};

/// Advection-diffusion equation solver
pub struct AdvectionDiffusionSolver<T: RealField + Copy> {
    config: FdmConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> AdvectionDiffusionSolver<T> {
    /// Create a new advection-diffusion solver
    pub fn new(config: FdmConfig<T>) -> Self {
        Self { config }
    }

    /// Solve steady-state advection-diffusion equation
    pub fn solve_steady(
        &self,
        grid: &StructuredGrid2D<T>,
        velocity_x: &HashMap<(usize, usize), T>,
        velocity_y: &HashMap<(usize, usize), T>,
        diffusivity: T,
        source: &HashMap<(usize, usize), T>,
        boundary_values: &HashMap<(usize, usize), T>,
    ) -> Result<HashMap<(usize, usize), T>> {
        let n = grid.nx() * grid.ny();
        let mut matrix_builder = SparseMatrixBuilder::new(n, n);
        let mut rhs = DVector::from_element(n, T::zero());

        // Build system matrix and RHS vector
        for (linear_idx, (i, j)) in grid.iter().enumerate() {
            if let Some(boundary_value) = boundary_values.get(&(i, j)).copied() {
                // Dirichlet boundary condition
                matrix_builder.add_entry(linear_idx, linear_idx, T::one())?;
                rhs[linear_idx] = boundary_value;
            } else {
                // Interior point: discretize advection-diffusion operator
                self.add_advection_diffusion_stencil(
                    &mut matrix_builder,
                    &mut rhs,
                    grid,
                    i,
                    j,
                    linear_idx,
                    velocity_x,
                    velocity_y,
                    diffusivity,
                    source,
                )?;
            }
        }

        // Solve the linear system
        let matrix = matrix_builder.build()?;
        let solution = solve_gauss_seidel(&matrix, &rhs, &self.config, "AdvectionDiffusion")?;

        // Convert solution back to grid coordinates
        let mut result = HashMap::new();
        for (linear_idx, (i, j)) in grid.iter().enumerate() {
            result.insert((i, j), solution[linear_idx]);
        }

        Ok(result)
    }

    /// Add advection-diffusion stencil for interior points
    fn add_advection_diffusion_stencil(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        grid: &StructuredGrid2D<T>,
        i: usize,
        j: usize,
        linear_idx: usize,
        velocity_x: &HashMap<(usize, usize), T>,
        velocity_y: &HashMap<(usize, usize), T>,
        diffusivity: T,
        source: &HashMap<(usize, usize), T>,
    ) -> Result<()> {
        let (dx, dy) = grid.spacing();
        let dx2 = dx * dx;
        let dy2 = dy * dy;

        let u = velocity_x.get(&(i, j)).copied().unwrap_or(T::zero());
        let v = velocity_y.get(&(i, j)).copied().unwrap_or(T::zero());

        // Central coefficient (diffusion part): -2α/dx² - 2α/dy²
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let mut center_coeff = -two * diffusivity / dx2 - two * diffusivity / dy2;

        // Add neighbor contributions with upwind scheme for advection
        let neighbors = grid.neighbors(i, j);
        for &(ni, nj) in &neighbors {
            let neighbor_idx = Self::linear_index(grid, ni, nj);
            let mut coeff = T::zero();

            if ni == i + 1 {
                // Right neighbor: diffusion + upwind advection
                coeff += diffusivity / dx2;
                if u < T::zero() {
                    // Negative velocity: flow from right to left
                    coeff += u / dx;
                    center_coeff -= u / dx;
                }
            } else if ni + 1 == i {
                // Left neighbor: diffusion + upwind advection
                coeff += diffusivity / dx2;
                if u > T::zero() {
                    // Positive velocity: flow from left to right
                    coeff += u / dx;
                    center_coeff -= u / dx;
                }
            } else if nj == j + 1 {
                // Top neighbor: diffusion + upwind advection
                coeff += diffusivity / dy2;
                if v < T::zero() {
                    // Negative velocity: flow from top to bottom
                    coeff += v / dy;
                    center_coeff -= v / dy;
                }
            } else if nj + 1 == j {
                // Bottom neighbor: diffusion + upwind advection
                coeff += diffusivity / dy2;
                if v > T::zero() {
                    // Positive velocity: flow from bottom to top
                    coeff += v / dy;
                    center_coeff -= v / dy;
                }
            }

            matrix_builder.add_entry(linear_idx, neighbor_idx, coeff)?;
        }

        matrix_builder.add_entry(linear_idx, linear_idx, center_coeff)?;

        // Set RHS from source term
        if let Some(source_value) = source.get(&(i, j)).copied() {
            rhs[linear_idx] = source_value;
        }

        Ok(())
    }

    /// Convert 2D grid indices to linear index
    fn linear_index(grid: &StructuredGrid2D<T>, i: usize, j: usize) -> usize {
        j * grid.nx() + i
    }
}
