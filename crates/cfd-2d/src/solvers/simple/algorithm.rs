use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::{validate_boundary_consistency, MomentumSolver};
use crate::solvers::fdm::PoissonSolver;
use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::sparse::SparseMatrix;
use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use std::collections::HashMap;

/// SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm
pub struct SimpleAlgorithm<T: RealField + Copy + FromPrimitive + std::fmt::Debug> {
    pub(crate) pressure_relaxation: T,
    pub(crate) velocity_relaxation: T,
    pub(crate) max_iterations: usize,
    pub(crate) tolerance: T,
    pub(crate) pressure_matrix: Option<SparseMatrix<T>>,
    pub(crate) matrix_builder: Option<cfd_math::sparse::SparseMatrixBuilder<T>>,
    pub(crate) rhs: Option<DVector<T>>,
    pub(crate) p_prime: Option<DVector<T>>,
    pub(crate) d_u: Option<Field2D<T>>,
    pub(crate) d_v: Option<Field2D<T>>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
    /// Construct with Patankar-recommended defaults: α_u = 0.7, α_p = 0.3.
    pub fn new() -> Self {
        Self {
            pressure_relaxation: T::from_f64(0.3).expect("Exact mathematically representable f64"),
            velocity_relaxation: T::from_f64(0.7).expect("Exact mathematically representable f64"),
            max_iterations: 50,
            tolerance: T::from_f64(1e-6).expect("Exact mathematically representable f64"),
            pressure_matrix: None,
            matrix_builder: None,
            rhs: None,
            p_prime: None,
            d_u: None,
            d_v: None,
        }
    }

    /// Set pressure under-relaxation factor α_p ∈ (0, 1).
    pub fn with_pressure_relaxation(mut self, alpha_p: T) -> Self {
        self.pressure_relaxation = alpha_p;
        self
    }

    /// Set velocity under-relaxation factor α_u ∈ (0, 1).
    pub fn with_velocity_relaxation(mut self, alpha_u: T) -> Self {
        self.velocity_relaxation = alpha_u;
        self
    }

    /// Set the maximum number of outer SIMPLE iterations.
    pub fn with_max_iterations(mut self, max_iter: usize) -> Self {
        self.max_iterations = max_iter;
        self
    }

    /// Set convergence tolerance for ‖∇·u‖_∞.
    pub fn with_tolerance(mut self, tol: T) -> Self {
        self.tolerance = tol;
        self
    }

    pub(crate) fn ensure_buffers(&mut self, nx: usize, ny: usize) {
        let n = nx * ny;
        if self.pressure_matrix.as_ref().is_none_or(|b| b.nrows() != n) {
            self.pressure_matrix = None; // Reset if matrix geometry violates new bounds
        }
        if self.rhs.as_ref().is_none_or(|v| v.len() != n) {
            self.rhs = Some(DVector::zeros(n));
        }
        if self.p_prime.as_ref().is_none_or(|v| v.len() != n) {
            self.p_prime = Some(DVector::zeros(n));
        }
        if self
            .d_u
            .as_ref()
            .is_none_or(|f| f.nx() != nx || f.ny() != ny)
        {
            self.d_u = Some(Field2D::new(nx, ny, T::zero()));
        }
        if self
            .d_v
            .as_ref()
            .is_none_or(|f| f.nx() != nx || f.ny() != ny)
        {
            self.d_v = Some(Field2D::new(nx, ny, T::zero()));
        }
    }

    /// Execute one SIMPLE outer iteration.
    ///
    /// ## Steps
    ///
    /// 1. **Momentum predictor** — solve u*, v* from linearised momentum equation.
    /// 2. **D-coefficient computation** — D = Δy/A_P for u-faces, Δx/A_P for v-faces.
    /// 3. **Rhie-Chow face velocities** — apply Theorem 2 correction to suppress
    ///    checker-board modes (correction factor = 1.0, mathematically exact).
    /// 4. **Pressure-correction system** — assemble M-matrix Poisson system for p'.
    /// 5. **Correction step** — update p and u,v using p'.
    ///
    /// Returns `(max_continuity_residual, converged)`.
    pub fn simple_iteration(
        &mut self,
        momentum_solver: &mut MomentumSolver<T>,
        _poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<(T, bool)> {
        validate_boundary_consistency(boundary_conditions, grid)
            .map_err(|error| cfd_core::error::Error::InvalidConfiguration(error.to_string()))?;
        self.ensure_buffers(grid.nx, grid.ny);
        self.predict_momentum(momentum_solver, fields, dt)?;
        self.compute_d_coefficients(momentum_solver, grid);
        let max_residual = self.assemble_pressure_correction(fields, grid, boundary_conditions)?;
        self.solve_pressure_correction()?;
        self.apply_corrections(fields, grid);
        let converged = max_residual < self.tolerance;
        Ok((max_residual, converged))
    }

    /// Run the full SIMPLE outer loop until convergence or `max_iterations`.
    ///
    /// Returns the number of iterations taken on success, or
    /// `Err(ConvergenceError::MaxIterationsExceeded)` if the residual does not
    /// fall below `tolerance`.
    pub fn solve_simple(
        &mut self,
        momentum_solver: &mut MomentumSolver<T>,
        poisson_solver: &mut PoissonSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
    ) -> Result<usize> {
        let mut iteration = 0;
        let mut converged = false;
        while iteration < self.max_iterations && !converged {
            let (_, conv) = self.simple_iteration(
                momentum_solver,
                poisson_solver,
                fields,
                dt,
                grid,
                boundary_conditions,
            )?;
            converged = conv;
            iteration += 1;
        }
        if !converged {
            return Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.max_iterations,
                },
            ));
        }
        Ok(iteration)
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> Default
    for SimpleAlgorithm<T>
{
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::StructuredGrid2D;
    use crate::physics::momentum::MomentumSolver;
    use crate::solvers::fdm::PoissonSolver;
    use cfd_core::physics::boundary::BoundaryCondition;
    use std::collections::HashMap;

    #[test]
    fn test_simple_algorithm_creation() {
        let simple = SimpleAlgorithm::<f64>::new();
        assert_eq!(simple.max_iterations, 50);
        assert!(simple.pressure_relaxation > 0.0 && simple.pressure_relaxation < 1.0);
        assert!(simple.velocity_relaxation > 0.0 && simple.velocity_relaxation < 1.0);
        // Invariant: Patankar's recommendation α_u + α_p ≈ 1
        assert!(
            (simple.velocity_relaxation + simple.pressure_relaxation - 1.0).abs() < 0.01,
            "Expected α_u + α_p ≈ 1 (Patankar 1980 recommendation)"
        );
    }

    #[test]
    fn test_rhie_chow_uniform_pressure_zero_correction() {
        // Theorem 2, Uniform Field Invariant:
        // For a uniform pressure field p = const, the RC correction term
        //   d_e * (avg_grad_pe - face_grad_pe)
        // must be identically zero, since both gradients equal zero.
        let dx = 0.1_f64;
        let _dy = 0.1_f64;
        let two = 2.0_f64;

        // All pressures equal → all gradients zero
        let (p_p, p_e, p_w, p_ee) = (1.0, 1.0, 1.0, 1.0);
        let grad_p_p_x = (p_e - p_w) / (two * dx);
        let grad_p_e_x = (p_ee - p_p) / (two * dx);
        let avg_grad_pe = (grad_p_p_x + grad_p_e_x) * 0.5;
        let face_grad_pe = (p_e - p_p) / dx;
        let rc_correction = avg_grad_pe - face_grad_pe;

        assert!(
            rc_correction.abs() < 1e-14,
            "RC correction for uniform pressure must be zero; got {rc_correction}"
        );
    }

    #[test]
    fn test_rhie_chow_checkerboard_nonzero_correction() {
        // Theorem 2, Checker-Board Mode:
        // For alternating p = 1, -1, 1 the RC correction must be non-zero.
        let dx = 0.1_f64;
        let dy = 0.1_f64;
        let two = 2.0_f64;

        // Checker-board pattern: p_P = 1, p_E = -1, p_W = -1, p_EE = 1
        let (p_p, p_e, p_w, p_ee) = (1.0_f64, -1.0, -1.0, 1.0);
        let _ = dy; // used conceptually
        let grad_p_p_x = (p_e - p_w) / (two * dx); // (-1 - -1)/(0.2) = 0
        let grad_p_e_x = (p_ee - p_p) / (two * dx); // (1 - 1)/(0.2)   = 0
        let avg_grad_pe = (grad_p_p_x + grad_p_e_x) * 0.5; // = 0
        let face_grad_pe = (p_e - p_p) / dx; // (-1 - 1)/0.1    = -20
        let rc_correction = avg_grad_pe - face_grad_pe; // 0 - (-20) = 20 ≠ 0

        assert!(
            rc_correction.abs() > 1e-10,
            "RC correction for checker-board pressure must be non-zero; got {rc_correction}"
        );
    }

    #[test]
    fn test_simple_iteration_rejects_missing_required_boundary() {
        let grid = StructuredGrid2D::<f64>::new(3, 3, 0.0, 1.0, 0.0, 1.0).unwrap();
        let mut simple = SimpleAlgorithm::<f64>::new();
        let mut momentum_solver = MomentumSolver::<f64>::new(&grid);
        let mut poisson_solver = PoissonSolver::<f64>::default();
        let mut fields = SimulationFields::<f64>::new(3, 3);
        let boundary_conditions = HashMap::from([
            ("west".to_string(), BoundaryCondition::wall_no_slip()),
            ("north".to_string(), BoundaryCondition::Outflow),
            ("south".to_string(), BoundaryCondition::Outflow),
        ]);

        let result = simple.simple_iteration(
            &mut momentum_solver,
            &mut poisson_solver,
            &mut fields,
            0.01,
            &grid,
            &boundary_conditions,
        );

        assert!(result.is_err());
    }
}
