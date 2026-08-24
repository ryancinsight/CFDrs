//! Advection-diffusion equation solver using Finite Difference Method.
//!
//! Solves: ∂φ/∂t + u·∇φ = α∇²φ + S
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

use cfd_core::error::{Error, Result};
use cfd_core::CfdScalar;
use cfd_math::sparse::SparseMatrixBuilder;
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};
use leto::Array1;
use std::collections::HashMap;

use super::config::FdmConfig;
use super::linear_solver::solve_gauss_seidel;
use crate::grid::StructuredGrid2D;
use crate::scalar;

/// Advection-diffusion equation solver
pub struct AdvectionDiffusionSolver<T: CfdScalar + EunomiaRealField + Copy> {
    config: FdmConfig<T>,
    matrix_builder: core::cell::RefCell<Option<SparseMatrixBuilder<T>>>,
}

impl<T: CfdScalar + EunomiaRealField + Copy + FloatElement> AdvectionDiffusionSolver<T> {
    /// Create a new advection-diffusion solver.
    ///
    /// # Panics
    /// Panics if `config` violates any invariant (see [`Self::try_new`]).
    #[must_use]
    pub fn new(config: FdmConfig<T>) -> Self {
        Self::try_new(config).unwrap_or_else(|error| {
            panic!("AdvectionDiffusionSolver::new called with invalid config: {error}");
        })
    }

    /// Create a new advection-diffusion solver with invariant validation.
    ///
    /// Validates the SSOT [`cfd_core::compute::solver::SolverConfig`]
    /// wrapped inside `FdmConfig`: tolerance must be finite and positive,
    /// `max_iterations` must be ≥ 1, and the relaxation factor must be
    /// finite and in `(0, 2]` (the same convergence range as Gauss–Seidel
    /// in [`super::PoissonSolver`]).
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if any invariant is violated.
    pub fn try_new(config: FdmConfig<T>) -> Result<Self> {
        if !<T as NumericElement>::is_finite(config.tolerance())
            || config.tolerance() <= <T as NumericElement>::ZERO
        {
            return Err(Error::InvalidConfiguration(format!(
                "AdvectionDiffusionSolver::try_new: tolerance must be finite and positive, got {:?}",
                config.tolerance()
            )));
        }
        if config.max_iterations() == 0 {
            return Err(Error::InvalidConfiguration(
                "AdvectionDiffusionSolver::try_new: max_iterations must be at least 1".to_string(),
            ));
        }
        let relaxation = config.relaxation_factor();
        let two = <T as FloatElement>::from_f64(2.0);
        if !<T as NumericElement>::is_finite(relaxation)
            || relaxation <= <T as NumericElement>::ZERO
            || relaxation > two
        {
            return Err(Error::InvalidConfiguration(format!(
                "AdvectionDiffusionSolver::try_new: relaxation_factor must be finite and in (0, 2] for Gauss-Seidel convergence, got {relaxation:?}"
            )));
        }
        Ok(Self {
            config,
            matrix_builder: core::cell::RefCell::new(None),
        })
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
        let mut matrix_builder = self
            .matrix_builder
            .borrow_mut()
            .take()
            .unwrap_or_else(|| SparseMatrixBuilder::new(n, n));
        let zero: T = scalar::zero();
        let one: T = scalar::one();
        let mut rhs = Array1::from_elem([n], zero);

        // Build system matrix and RHS vector
        for (linear_idx, (i, j)) in grid.iter().enumerate() {
            if let Some(boundary_value) = boundary_values.get(&(i, j)).copied() {
                // Dirichlet boundary condition
                matrix_builder.add_entry(linear_idx, linear_idx, one)?;
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
        *self.matrix_builder.borrow_mut() = Some(SparseMatrixBuilder::new(n, n));
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
        rhs: &mut Array1<T>,
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

        let zero: T = scalar::zero();
        let u = velocity_x.get(&(i, j)).copied().unwrap_or(zero);
        let v = velocity_y.get(&(i, j)).copied().unwrap_or(zero);

        // Central coefficient for u·∇φ - α∇²φ = S.
        let two = <T as FloatElement>::from_f64(2.0);
        let mut center_coeff = two * diffusivity / dx2 + two * diffusivity / dy2;

        // Add neighbor contributions with upwind scheme for advection
        let neighbors = grid.neighbors(i, j);
        for &(ni, nj) in &neighbors {
            let neighbor_idx = Self::linear_index(grid, ni, nj);
            let mut coeff = zero;

            if ni == i + 1 {
                // Right neighbor: diffusion + upwind advection
                coeff -= diffusivity / dx2;
                if u < zero {
                    // Negative velocity: flow from right to left
                    coeff += u / dx;
                    center_coeff -= u / dx;
                }
            } else if ni + 1 == i {
                // Left neighbor: diffusion + upwind advection
                coeff -= diffusivity / dx2;
                if u > zero {
                    // Positive velocity: flow from left to right
                    coeff -= u / dx;
                    center_coeff += u / dx;
                }
            } else if nj == j + 1 {
                // Top neighbor: diffusion + upwind advection
                coeff -= diffusivity / dy2;
                if v < zero {
                    // Negative velocity: flow from top to bottom
                    coeff += v / dy;
                    center_coeff -= v / dy;
                }
            } else if nj + 1 == j {
                // Bottom neighbor: diffusion + upwind advection
                coeff -= diffusivity / dy2;
                if v > zero {
                    // Positive velocity: flow from bottom to top
                    coeff -= v / dy;
                    center_coeff += v / dy;
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

#[cfg(test)]
mod tests {
    use super::*;

    /// **Positive**: `try_new` accepts a default configuration.
    #[test]
    fn advection_diffusion_try_new_accepts_default() {
        let _solver = AdvectionDiffusionSolver::<f64>::try_new(FdmConfig::<f64>::default())
            .expect("default must succeed");
    }

    /// **Adversarial**: zero `tolerance` is rejected.
    #[test]
    fn advection_diffusion_try_new_rejects_zero_tolerance() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.tolerance = 0.0;
        match AdvectionDiffusionSolver::try_new(config) {
            Err(e) => assert!(
                e.to_string().contains("tolerance"),
                "error must mention tolerance: {e}"
            ),
            Ok(_) => panic!("zero tolerance must be rejected"),
        }
    }

    /// **Adversarial**: zero `max_iterations` is rejected.
    #[test]
    fn advection_diffusion_try_new_rejects_zero_max_iterations() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.max_iterations = 0;
        match AdvectionDiffusionSolver::try_new(config) {
            Err(e) => assert!(
                e.to_string().contains("max_iterations"),
                "error must mention max_iterations: {e}"
            ),
            Ok(_) => panic!("zero max_iterations must be rejected"),
        }
    }

    /// **Adversarial**: relaxation factor outside (0, 2] is rejected.
    #[test]
    fn advection_diffusion_try_new_rejects_invalid_relaxation() {
        for &bad in &[0.0_f64, -1.0, 3.0, f64::NAN] {
            let mut config = FdmConfig::<f64>::default();
            config.base.numerical.relaxation = bad;
            match AdvectionDiffusionSolver::try_new(config) {
                Err(e) => assert!(
                    e.to_string().contains("relaxation_factor"),
                    "error must mention relaxation_factor for input {bad}: {e}"
                ),
                Ok(_) => panic!("relaxation {bad} must be rejected"),
            }
        }
    }

    /// **Boundary**: `new` panics on invalid `tolerance` (thin wrapper contract).
    #[test]
    #[should_panic(expected = "tolerance")]
    fn advection_diffusion_new_panics_on_invalid_tolerance() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.tolerance = f64::NAN;
        let _ = AdvectionDiffusionSolver::new(config);
    }
}
