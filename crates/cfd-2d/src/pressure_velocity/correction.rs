//! Pressure correction equation assembly and solution
//!
//! Assembles the discrete Laplacian system from the momentum-predicted velocity
//! field and dispatches to the configured linear solver.
//!
//! # Theorem (Pressure Poisson Equation — Patankar 1980)
//!
//! For incompressible flows, substituting the velocity correction into the
//! continuity constraint `∇·u = 0` yields the pressure correction equation:
//!
//! ```text
//! ∇·(d ∇p') = ∇·u*
//! ```
//!
//! where `d = Vol / a_P` is the pressure-gradient coefficient and `u*` is the
//! momentum-predicted velocity. This is a symmetric positive-definite Poisson
//! equation on the fluid domain, guaranteeing existence and uniqueness of `p'`
//! up to an additive constant (fixed by a reference-pressure cell).
//!
//! **Proof sketch**: Diagonal dominance of the momentum matrix ensures `a_P > 0`,
//! so `d > 0` everywhere. The discrete Laplacian `∇·(d ∇·)` with homogeneous
//! Neumann BCs is then SPD on the complement of the null-space spanned by the
//! constant vector. Pinning one cell removes this null-space, yielding a unique
//! solution.

use super::config::PressureLinearSolver;
use super::pressure::PressureCorrectionSolver;
use crate::grid::array2d::Array2D;
use crate::physics::momentum::validate_boundary_consistency;
use cfd_math::linear_solver::preconditioners::{AlgebraicMultigrid, IdentityPreconditioner};
use cfd_math::linear_solver::{DirectSparseSolver, IterativeLinearSolver};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::fmt::Debug;

impl<T: RealField + Copy + FromPrimitive + Debug + num_traits::ToPrimitive>
    PressureCorrectionSolver<T>
{
    /// Dispatch a linear solve to the configured solver backend
    pub(super) fn dispatch_solve(
        &self,
        matrix: &cfd_math::sparse::SparseMatrix<T>,
        rhs: &DVector<T>,
        solution: &mut DVector<T>,
    ) -> cfd_core::error::Result<()> {
        // Phase 8: Deep Optimization & AMG Caching
        //
        // Theorem — Sparse Galerkin Caching (Ruge & Stüben 1987)
        // Re-computing algebraic multigrid transfer operators (R, P) takes O(N) but with a large constant.
        // For pressure Poisson equations on fixed grid topologies, the non-zero pattern remains invariant
        // while coefficients change slowly. We exactly recompute A_c = R A_f P without graph re-coarsening.
        let mut amg_cache = self._amg_preconditioner.borrow_mut();
        if amg_cache.is_none() {
            let config = cfd_math::linear_solver::AMGConfig::default();
            // Initialize the AMG preconditioner and cache it
            if let Ok(amg) = AlgebraicMultigrid::new(matrix, config) {
                *amg_cache = Some(amg);
            }
        } else {
            // Hot-path optimization: $A_c = R A_f P$
            // Recalculates coarse matrices exactly mapping changing continuity coefficients
            // while bypassing O(N) Ruge-Stüben strong-connection heuristics.
            let _ = amg_cache
                .as_mut()
                .expect("guarded by is_none check")
                .recompute(matrix);
        }

        let amg_preconditioner: Option<&AlgebraicMultigrid<T>> = amg_cache.as_ref();

        let solve_result = match self.solver_type {
            PressureLinearSolver::ConjugateGradient => {
                let result = if let Some(amg) = amg_preconditioner {
                    self.cg_solver.solve(matrix, rhs, solution, Some(amg))
                } else {
                    self.cg_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                };
                match result {
                    Ok(_) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => self
                        .cg_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ()),
                    Err(e) => Err(e),
                }
            }
            PressureLinearSolver::BiCGSTAB => {
                let result = if let Some(amg) = amg_preconditioner {
                    self.bicgstab_solver.solve(matrix, rhs, solution, Some(amg))
                } else {
                    self.bicgstab_solver.solve(
                        matrix,
                        rhs,
                        solution,
                        None::<&IdentityPreconditioner>,
                    )
                };
                match result {
                    Ok(_) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => self
                        .bicgstab_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ()),
                    Err(e) => Err(e),
                }
            }
            PressureLinearSolver::GMRES { .. } => {
                let Some(ref solver) = self.gmres_solver else {
                    return Err(cfd_core::error::Error::InvalidConfiguration(
                        "GMRES solver not initialized".to_string(),
                    ));
                };
                let result = if let Some(amg) = amg_preconditioner {
                    solver.solve(matrix, rhs, solution, Some(amg))
                } else {
                    solver.solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                };
                match result {
                    Ok(_) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ()),
                    Err(e) => Err(e),
                }
            }
        };

        match solve_result {
            Ok(()) => Ok(()),
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { max },
            )) => {
                tracing::warn!(
                    solver = ?self.solver_type,
                    max_iterations = max,
                    "Pressure solve stalled; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                *solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Pressure direct sparse solve failed for {:?}: {error}",
                        self.solver_type
                    ))
                })?;
                Ok(())
            }
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Breakdown,
            )) => {
                tracing::warn!(
                    solver = ?self.solver_type,
                    "Pressure solve breakdown; falling back to direct sparse solve"
                );
                let direct_solver = DirectSparseSolver::default();
                *solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                    cfd_core::error::Error::Solver(format!(
                        "Pressure direct sparse solve failed for {:?}: {error}",
                        self.solver_type
                    ))
                })?;
                Ok(())
            }
            Err(error) => Err(error),
        }
    }

    /// Solve pressure correction equation from cell-centred velocities
    ///
    /// The pressure correction equation is derived from the continuity equation:
    /// `∇²p' = ρ/Δt · ∇·u*`
    pub fn solve_pressure_correction(
        &self,
        fields: &crate::fields::SimulationFields<T>,
        dt: T,
        rho: T,
        _rebuild_matrix: bool,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        let n = (nx - 2) * (ny - 2);

        if n <= 1 {
            for i in 0..nx {
                for j in 0..ny {
                    output_correction[(i, j)] = T::zero();
                }
            }
            return Ok(());
        }

        let system_size = n - 1;
        let reference_idx = (1..nx - 1)
            .flat_map(|i| (1..ny - 1).map(move |j| (i, j)))
            .enumerate()
            .find(|(_, (i, j))| fields.mask.at(*i, *j))
            .map_or(0usize, |(idx, _)| idx);

        let map_index = |idx: usize| -> Option<usize> {
            match idx.cmp(&reference_idx) {
                std::cmp::Ordering::Equal => None,
                std::cmp::Ordering::Less => Some(idx),
                std::cmp::Ordering::Greater => Some(idx - 1),
            }
        };

        let mut builder = self.take_matrix_builder(system_size, system_size);
        let mut rhs = self
            ._rhs_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.len() == system_size)
            .unwrap_or_else(|| DVector::zeros(system_size));
        rhs.fill(T::zero());

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);
        let coeff = rho / dt;

        let two = T::from_f64(cfd_core::physics::constants::mathematical::numeric::TWO)
            .unwrap_or_else(|| T::one() + T::one());

        let mut n_fluid = 0;
        let mut n_solid = 0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                if idx == reference_idx {
                    continue;
                }
                let row_idx = map_index(idx).expect("row index must exist");

                if !fields.mask.at(i, j) {
                    n_solid += 1;
                    builder.add_entry(row_idx, row_idx, T::one())?;
                    rhs[row_idx] = T::zero();
                    continue;
                }
                n_fluid += 1;

                let mut ap = T::zero();

                // West neighbour
                if i > 1 && fields.mask.at(i - 1, j) {
                    ap += dx2_inv;
                    let ni = idx - (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dx2_inv)?;
                    }
                } else {
                    ap += dx2_inv;
                }
                // East neighbour
                if i < nx - 2 && fields.mask.at(i + 1, j) {
                    ap += dx2_inv;
                    let ni = idx + (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dx2_inv)?;
                    }
                } else {
                    ap += dx2_inv;
                }
                // South neighbour
                if j > 1 && fields.mask.at(i, j - 1) {
                    ap += dy2_inv;
                    let ni = idx - 1;
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dy2_inv)?;
                    }
                } else {
                    ap += dy2_inv;
                }
                // North neighbour
                if j < ny - 2 && fields.mask.at(i, j + 1) {
                    ap += dy2_inv;
                    let ni = idx + 1;
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dy2_inv)?;
                    }
                } else {
                    ap += dy2_inv;
                }

                builder.add_entry(row_idx, row_idx, ap)?;

                let div_u = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (two * dx)
                    + (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (two * dy);
                rhs[row_idx] = -coeff * div_u;
            }
        }

        let rhs_norm = rhs.norm();
        tracing::debug!(
            "Pressure Solve: n_fluid={n_fluid}, n_solid={n_solid}, \
             system_size={system_size}, rhs_norm={rhs_norm:?}"
        );

        let matrix = builder.build()?;
        self.reset_matrix_builder_cache(system_size, system_size);
        let mut p_correction_vec = self
            ._solution_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.len() == matrix.nrows())
            .unwrap_or_else(|| DVector::zeros(matrix.nrows()));
        p_correction_vec.fill(T::zero());
        self.dispatch_solve(&matrix, &rhs, &mut p_correction_vec)?;

        self.scatter_correction(
            &p_correction_vec,
            reference_idx,
            &map_index,
            output_correction,
        )?;

        *self._rhs_cache.borrow_mut() = Some(rhs);
        *self._solution_cache.borrow_mut() = Some(p_correction_vec);

        Ok(())
    }

    /// Solve pressure correction equation using face velocities (Rhie-Chow)
    pub fn solve_pressure_correction_from_faces(
        &self,
        u_face: &Array2D<T>,
        v_face: &Array2D<T>,
        d_x: &Array2D<T>,
        d_y: &Array2D<T>,
        rho: T,
        fields: &crate::fields::SimulationFields<T>,
        boundary_conditions: &std::collections::HashMap<
            String,
            cfd_core::physics::boundary::BoundaryCondition<T>,
        >,
        _rebuild_matrix: bool,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        validate_boundary_consistency(boundary_conditions, &self.grid)
            .map_err(|error| cfd_core::error::Error::InvalidConfiguration(error.to_string()))?;

        let n = (nx - 2) * (ny - 2);
        if n <= 1 {
            for i in 0..nx {
                for j in 0..ny {
                    output_correction[(i, j)] = T::zero();
                }
            }
            return Ok(());
        }

        let system_size = n - 1;
        let reference_idx = (1..nx - 1)
            .flat_map(|i| (1..ny - 1).map(move |j| (i, j)))
            .enumerate()
            .find(|(_, (i, j))| fields.mask.at(*i, *j))
            .map_or(0usize, |(idx, _)| idx);

        let map_index = |idx: usize| -> Option<usize> {
            match idx.cmp(&reference_idx) {
                std::cmp::Ordering::Equal => None,
                std::cmp::Ordering::Less => Some(idx),
                std::cmp::Ordering::Greater => Some(idx - 1),
            }
        };

        let mut builder = self.take_matrix_builder(system_size, system_size);
        let mut rhs = self
            ._rhs_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.len() == system_size)
            .unwrap_or_else(|| DVector::zeros(system_size));
        rhs.fill(T::zero());

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);

        let mut max_residual = T::zero();

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                if idx == reference_idx {
                    continue;
                }
                let row_idx = map_index(idx).expect("row index must exist");

                if !fields.mask.at(i, j) {
                    builder.add_entry(row_idx, row_idx, T::one())?;
                    rhs[row_idx] = T::zero();
                    continue;
                }

                let mut ap = T::zero();

                if i > 1 && fields.mask.at(i - 1, j) {
                    let aw = d_x[(i - 1, j)] * dx2_inv;
                    ap += aw;
                    if let Some(ci) = map_index(idx - (ny - 2)) {
                        builder.add_entry(row_idx, ci, -aw)?;
                    }
                } else {
                    ap += d_x[(i - 1, j)] * dx2_inv;
                }
                if i < nx - 2 && fields.mask.at(i + 1, j) {
                    let ae = d_x[(i, j)] * dx2_inv;
                    ap += ae;
                    if let Some(ci) = map_index(idx + (ny - 2)) {
                        builder.add_entry(row_idx, ci, -ae)?;
                    }
                } else {
                    ap += d_x[(i, j)] * dx2_inv;
                }
                if j > 1 && fields.mask.at(i, j - 1) {
                    let as_ = d_y[(i, j - 1)] * dy2_inv;
                    ap += as_;
                    if let Some(ci) = map_index(idx - 1) {
                        builder.add_entry(row_idx, ci, -as_)?;
                    }
                } else {
                    ap += d_y[(i, j - 1)] * dy2_inv;
                }
                if j < ny - 2 && fields.mask.at(i, j + 1) {
                    let an = d_y[(i, j)] * dy2_inv;
                    ap += an;
                    if let Some(ci) = map_index(idx + 1) {
                        builder.add_entry(row_idx, ci, -an)?;
                    }
                } else {
                    ap += d_y[(i, j)] * dy2_inv;
                }

                builder.add_entry(row_idx, row_idx, ap)?;

                let div_u = (u_face[(i, j)] - u_face[(i - 1, j)]) / dx
                    + (v_face[(i, j)] - v_face[(i, j - 1)]) / dy;
                rhs[row_idx] = -rho * div_u;

                if rhs[row_idx].abs() > max_residual {
                    max_residual = rhs[row_idx].abs();
                }
            }
        }

        tracing::debug!(
            "Pressure Solve (faces): n={n}, rhs_norm={:?}, max_residual={max_residual:?}",
            rhs.norm()
        );

        let matrix = builder.build()?;
        self.reset_matrix_builder_cache(system_size, system_size);
        let mut p_correction_vec = self
            ._solution_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.len() == matrix.nrows())
            .unwrap_or_else(|| DVector::zeros(matrix.nrows()));
        p_correction_vec.fill(T::zero());
        self.dispatch_solve(&matrix, &rhs, &mut p_correction_vec)?;

        self.scatter_correction(
            &p_correction_vec,
            reference_idx,
            &map_index,
            output_correction,
        )?;

        *self._rhs_cache.borrow_mut() = Some(rhs);
        *self._solution_cache.borrow_mut() = Some(p_correction_vec);

        Ok(())
    }

    /// Scatter solution vector back to 2D grid with Neumann BCs
    fn scatter_correction(
        &self,
        solution: &DVector<T>,
        reference_idx: usize,
        map_index: &dyn Fn(usize) -> Option<usize>,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                let value = if idx == reference_idx {
                    T::zero()
                } else if let Some(col_idx) = map_index(idx) {
                    solution[col_idx]
                } else {
                    T::zero()
                };
                output_correction[(i, j)] = value;
            }
        }

        // Neumann BCs (zero gradient)
        for i in 0..nx {
            output_correction[(i, 0)] = output_correction[(i, 1)];
            output_correction[(i, ny - 1)] = output_correction[(i, ny - 2)];
        }
        for j in 0..ny {
            output_correction[(0, j)] = output_correction[(1, j)];
            output_correction[(nx - 1, j)] = output_correction[(nx - 2, j)];
        }

        Ok(())
    }

    fn take_matrix_builder(&self, rows: usize, cols: usize) -> SparseMatrixBuilder<T> {
        self._matrix_builder_cache
            .borrow_mut()
            .take()
            .filter(|builder| builder.num_rows() == rows && builder.num_cols() == cols)
            .unwrap_or_else(|| SparseMatrixBuilder::new(rows, cols))
    }

    fn reset_matrix_builder_cache(&self, rows: usize, cols: usize) {
        *self._matrix_builder_cache.borrow_mut() = Some(SparseMatrixBuilder::new(rows, cols));
    }
}

#[cfg(test)]
mod tests {
    use crate::fields::SimulationFields;
    use crate::grid::array2d::Array2D;
    use crate::grid::StructuredGrid2D;
    use crate::pressure_velocity::config::PressureLinearSolver;
    use crate::pressure_velocity::pressure::PressureCorrectionSolver;
    use cfd_core::physics::boundary::BoundaryCondition;
    use std::collections::HashMap;

    /// Helper: create a small grid and pressure correction solver
    fn make_solver(nx: usize, ny: usize) -> (PressureCorrectionSolver<f64>, StructuredGrid2D<f64>) {
        let grid = StructuredGrid2D::new(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
        let solver =
            PressureCorrectionSolver::new(grid.clone(), PressureLinearSolver::ConjugateGradient)
                .unwrap();
        (solver, grid)
    }

    #[test]
    fn zero_divergence_produces_near_zero_correction() {
        let (solver, _grid) = make_solver(8, 8);
        // All-zero velocity => zero divergence
        let fields: SimulationFields<f64> = SimulationFields::new(8, 8);
        let dt = 0.01;
        let rho = 1.0;
        let mut p_corr = Array2D::new(8, 8, 0.0);
        solver
            .solve_pressure_correction(&fields, dt, rho, true, &mut p_corr)
            .unwrap();

        for &val in p_corr.iter() {
            assert!(
                val.abs() < 1e-10,
                "Expected near-zero pressure correction, got {val}"
            );
        }
    }

    #[test]
    fn pressure_correction_field_is_finite() {
        let (solver, _grid) = make_solver(8, 8);
        let mut fields: SimulationFields<f64> = SimulationFields::new(8, 8);

        // Set a non-trivial velocity field
        for i in 0..8 {
            for j in 0..8 {
                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = 0.1 * (i as f64);
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = -0.05 * (j as f64);
                }
            }
        }

        let dt = 0.01;
        let rho = 1000.0;
        let mut p_corr = Array2D::new(8, 8, 0.0);
        solver
            .solve_pressure_correction(&fields, dt, rho, true, &mut p_corr)
            .unwrap();

        for &val in p_corr.iter() {
            assert!(val.is_finite(), "Pressure correction contains NaN or Inf");
        }
    }

    #[test]
    fn face_based_pressure_correction_with_walls_is_well_posed() {
        let grid = StructuredGrid2D::new(6, 6, 0.0, 1.0, 0.0, 1.0).unwrap();
        let solver = PressureCorrectionSolver::new(
            grid.clone(),
            PressureLinearSolver::GMRES { restart_dim: 10 },
        )
        .unwrap();
        let fields: SimulationFields<f64> = SimulationFields::new(6, 6);

        let boundary_conditions = HashMap::from([
            ("west".to_string(), BoundaryCondition::wall_no_slip()),
            ("east".to_string(), BoundaryCondition::wall_no_slip()),
            ("north".to_string(), BoundaryCondition::wall_no_slip()),
            ("south".to_string(), BoundaryCondition::wall_no_slip()),
        ]);

        let u_face = Array2D::new(grid.nx - 1, grid.ny, 0.0);
        let v_face = Array2D::new(grid.nx, grid.ny - 1, 0.0);
        let d_x = Array2D::new(grid.nx - 1, grid.ny, 1.0);
        let d_y = Array2D::new(grid.nx, grid.ny - 1, 1.0);
        let mut p_corr = Array2D::new(grid.nx, grid.ny, 1.0);

        solver
            .solve_pressure_correction_from_faces(
                &u_face,
                &v_face,
                &d_x,
                &d_y,
                1.0,
                &fields,
                &boundary_conditions,
                true,
                &mut p_corr,
            )
            .unwrap();

        for i in 0..grid.nx {
            for j in 0..grid.ny {
                let value = p_corr[(i, j)];
                assert!(value.is_finite(), "pressure correction must be finite");
                assert!(value.abs() < 1e-10, "expected zero correction, got {value}");
            }
        }
    }

    #[test]
    fn basic_creation_and_initialization() {
        for solver_type in [
            PressureLinearSolver::ConjugateGradient,
            PressureLinearSolver::BiCGSTAB,
        ] {
            let grid = StructuredGrid2D::new(6, 6, 0.0, 1.0, 0.0, 1.0).unwrap();
            let solver = PressureCorrectionSolver::new(grid, solver_type);
            assert!(
                solver.is_ok(),
                "Failed to create solver with {solver_type:?}"
            );
        }
    }
}
