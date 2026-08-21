//! Poisson equation solver using Finite Difference Method.
//!
//! Solves the Poisson equation: ∇²φ = f
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

use cfd_core::CfdScalar;
use cfd_core::error::{Error, Result};
use cfd_math::sparse::SparseMatrixBuilder;
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};
use leto::Array1;
use std::collections::HashMap;

use super::config::FdmConfig;
use super::linear_solver::solve_gauss_seidel;
use crate::grid::StructuredGrid2D;
use crate::scalar;

/// Poisson equation solver
pub struct PoissonSolver<T: CfdScalar + EunomiaRealField + Copy> {
    config: FdmConfig<T>,
    matrix_builder: core::cell::RefCell<Option<SparseMatrixBuilder<T>>>,
}

impl<T: CfdScalar + EunomiaRealField + Copy + FloatElement> PoissonSolver<T> {
    /// Create new Poisson solver.
    ///
    /// # Panics
    /// Panics if `config` violates any invariant (see [`Self::try_new`]).
    #[must_use]
    pub fn new(config: FdmConfig<T>) -> Self {
        Self::try_new(config).unwrap_or_else(|error| {
            panic!("PoissonSolver::new called with invalid config: {error}");
        })
    }

    /// Create new Poisson solver with invariant validation.
    ///
    /// Validates the SSOT [`cfd_core::compute::solver::SolverConfig`]
    /// wrapped inside `FdmConfig`: tolerance must be finite and positive,
    /// `max_iterations` must be ≥ 1, and the relaxation factor must be
    /// finite and in `(0, 2]` (Gauss–Seidel converges iff the relaxation
    /// factor lies in that interval).
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if any invariant is violated.
    pub fn try_new(config: FdmConfig<T>) -> Result<Self> {
        if !<T as NumericElement>::is_finite(config.tolerance())
            || config.tolerance() <= <T as NumericElement>::ZERO
        {
            return Err(Error::InvalidConfiguration(format!(
                "PoissonSolver::try_new: tolerance must be finite and positive, got {:?}",
                config.tolerance()
            )));
        }
        if config.max_iterations() == 0 {
            return Err(Error::InvalidConfiguration(
                "PoissonSolver::try_new: max_iterations must be at least 1".to_string(),
            ));
        }
        let relaxation = config.relaxation_factor();
        let one = <T as FloatElement>::from_f64(1.0);
        let two = <T as FloatElement>::from_f64(2.0);
        if !<T as NumericElement>::is_finite(relaxation)
            || relaxation <= <T as NumericElement>::ZERO
            || relaxation > two
        {
            return Err(Error::InvalidConfiguration(format!(
                "PoissonSolver::try_new: relaxation_factor must be finite and in (0, 2] for Gauss-Seidel convergence, got {relaxation:?}"
            )));
        }
        // Ensure one is used (warning suppressor for cfg(conditionals)).
        let _ = one;
        Ok(Self {
            config,
            matrix_builder: core::cell::RefCell::new(None),
        })
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
            if let Some(boundary_value) = dirichlet_boundaries.get(&(i, j)).copied() {
                // Dirichlet boundary condition: φ = boundary_value
                matrix_builder.add_entry(linear_idx, linear_idx, one)?;
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
        *self.matrix_builder.borrow_mut() = Some(SparseMatrixBuilder::new(n, n));
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
        rhs: &mut Array1<T>,
        grid: &StructuredGrid2D<T>,
        i: usize,
        j: usize,
        linear_idx: usize,
        gradient: T,
    ) -> Result<()> {
        let (dx, dy) = grid.spacing();

        if i == 0 {
            let one: T = scalar::one();
            // Left wall, normal -x: T_{0,j} - T_{1,j} = g * dx
            let neighbor_idx = Self::linear_index(grid, i + 1, j);
            matrix_builder.add_entry(linear_idx, linear_idx, one)?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -one)?;
            rhs[linear_idx] = gradient * dx;
        } else if i == grid.nx() - 1 {
            let one: T = scalar::one();
            // Right wall, normal +x: T_{nx-1,j} - T_{nx-2,j} = g * dx
            let neighbor_idx = Self::linear_index(grid, i - 1, j);
            matrix_builder.add_entry(linear_idx, linear_idx, one)?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -one)?;
            rhs[linear_idx] = gradient * dx;
        } else if j == 0 {
            let one: T = scalar::one();
            // Bottom wall, normal -y: T_{i,0} - T_{i,1} = g * dy
            let neighbor_idx = Self::linear_index(grid, i, j + 1);
            matrix_builder.add_entry(linear_idx, linear_idx, one)?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -one)?;
            rhs[linear_idx] = gradient * dy;
        } else if j == grid.ny() - 1 {
            let one: T = scalar::one();
            // Top wall, normal +y: T_{i,ny-1} - T_{i,ny-2} = g * dy
            let neighbor_idx = Self::linear_index(grid, i, j - 1);
            matrix_builder.add_entry(linear_idx, linear_idx, one)?;
            matrix_builder.add_entry(linear_idx, neighbor_idx, -one)?;
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
        rhs: &mut Array1<T>,
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
        let one: T = scalar::one();
        let two = <T as FloatElement>::from_f64(2.0);

        // Correct 5-point stencil for ∇²φ = f:
        // Center: φ_i,j coefficient = -(2/dx² + 2/dy²)
        let center_coeff = -two / dx2 - two / dy2;
        matrix_builder.add_entry(linear_idx, linear_idx, center_coeff)?;

        // Neighbor contributions: +1/dx² and +1/dy²
        // Left neighbor (i-1)
        if i > 0 {
            let neighbor_idx = Self::linear_index(grid, i - 1, j);
            matrix_builder.add_entry(linear_idx, neighbor_idx, one / dx2)?;
        }

        // Right neighbor (i+1)
        if i < grid.nx() - 1 {
            let neighbor_idx = Self::linear_index(grid, i + 1, j);
            matrix_builder.add_entry(linear_idx, neighbor_idx, one / dx2)?;
        }

        // Bottom neighbor (j-1)
        if j > 0 {
            let neighbor_idx = Self::linear_index(grid, i, j - 1);
            matrix_builder.add_entry(linear_idx, neighbor_idx, one / dy2)?;
        }

        // Top neighbor (j+1)
        if j < grid.ny() - 1 {
            let neighbor_idx = Self::linear_index(grid, i, j + 1);
            matrix_builder.add_entry(linear_idx, neighbor_idx, one / dy2)?;
        }

        // Set RHS: source term f
        rhs[linear_idx] = source
            .get(&(i, j))
            .copied()
            .unwrap_or_else(scalar::zero::<T>);

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
    fn poisson_try_new_accepts_default() {
        let _solver = PoissonSolver::<f64>::try_new(FdmConfig::<f64>::default())
            .expect("default must succeed");
    }

    /// **Adversarial**: zero `tolerance` is rejected.
    #[test]
    fn poisson_try_new_rejects_zero_tolerance() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.tolerance = 0.0;
        match PoissonSolver::try_new(config) {
            Err(e) => assert!(
                e.to_string().contains("tolerance"),
                "error must mention tolerance: {e}"
            ),
            Ok(_) => panic!("zero tolerance must be rejected"),
        }
    }

    /// **Adversarial**: zero `max_iterations` is rejected.
    #[test]
    fn poisson_try_new_rejects_zero_max_iterations() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.max_iterations = 0;
        match PoissonSolver::try_new(config) {
            Err(e) => assert!(
                e.to_string().contains("max_iterations"),
                "error must mention max_iterations: {e}"
            ),
            Ok(_) => panic!("zero max_iterations must be rejected"),
        }
    }

    /// **Adversarial**: relaxation factor outside (0, 2] is rejected
    /// (Gauss–Seidel convergence range).
    #[test]
    fn poisson_try_new_rejects_invalid_relaxation() {
        for &bad in &[0.0_f64, -1.0, 3.0, f64::NAN] {
            let mut config = FdmConfig::<f64>::default();
            config.base.numerical.relaxation = bad;
            match PoissonSolver::try_new(config) {
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
    fn poisson_new_panics_on_invalid_tolerance() {
        let mut config = FdmConfig::<f64>::default();
        config.base.convergence.tolerance = f64::NAN;
        let _ = PoissonSolver::new(config);
    }
}
