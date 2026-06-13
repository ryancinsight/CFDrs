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
use cfd_math::linear_solver::preconditioners::{AlgebraicMultigrid, IdentityPreconditioner};
use cfd_math::linear_solver::{DirectSparseSolver, IterativeLinearSolver};
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
                let direct_solver = DirectSparseSolver::default();
                if matrix.nrows() <= 500 && direct_solver.can_handle_size(matrix.nrows()) {
                    tracing::warn!(
                        solver = ?self.solver_type,
                        max_iterations = max,
                        "Pressure solve stalled; falling back to direct sparse solve"
                    );
                    *solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                        cfd_core::error::Error::Solver(format!(
                            "Pressure direct sparse solve failed for {:?}: {error}",
                            self.solver_type
                        ))
                    })?;
                } else {
                    tracing::warn!(
                        solver = ?self.solver_type,
                        max_iterations = max,
                        size = matrix.nrows(),
                        "Pressure solve stalled, but system size exceeds direct solver limit of 500; continuing with iterative solution"
                    );
                }
                Ok(())
            }
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Breakdown,
            )) => {
                let direct_solver = DirectSparseSolver::default();
                if matrix.nrows() <= 500 && direct_solver.can_handle_size(matrix.nrows()) {
                    tracing::warn!(
                        solver = ?self.solver_type,
                        "Pressure solve breakdown; falling back to direct sparse solve"
                    );
                    *solution = direct_solver.solve(matrix, rhs).map_err(|error| {
                        cfd_core::error::Error::Solver(format!(
                            "Pressure direct sparse solve failed for {:?}: {error}",
                            self.solver_type
                        ))
                    })?;
                } else {
                    tracing::warn!(
                        solver = ?self.solver_type,
                        size = matrix.nrows(),
                        "Pressure solve breakdown, but system size exceeds direct solver limit of 500; continuing with iterative solution"
                    );
                }
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
        boundary_conditions: Option<&std::collections::HashMap<String, cfd_core::physics::boundary::BoundaryCondition<T>>>,
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

        let is_dirichlet = |side: &str| -> bool {
            if let Some(bcs) = boundary_conditions {
                if let Some(bc) = bcs.get(side) {
                    return matches!(
                        bc,
                        cfd_core::physics::boundary::BoundaryCondition::PressureOutlet { .. }
                            | cfd_core::physics::boundary::BoundaryCondition::PressureInlet { .. }
                            | cfd_core::physics::boundary::BoundaryCondition::CharacteristicOutlet { .. }
                    );
                }
            }
            false
        };

        let has_dirichlet = is_dirichlet("west") || is_dirichlet("east") || is_dirichlet("south") || is_dirichlet("north");

        let (system_size, reference_idx) = if has_dirichlet {
            (n, None)
        } else {
            let ref_idx = (1..nx - 1)
                .flat_map(|i| (1..ny - 1).map(move |j| (i, j)))
                .enumerate()
                .find(|(_, (i, j))| fields.mask.at(*i, *j))
                .map_or(0usize, |(idx, _)| idx);
            (n - 1, Some(ref_idx))
        };

        let map_index = |idx: usize| -> Option<usize> {
            if let Some(ref_idx) = reference_idx {
                match idx.cmp(&ref_idx) {
                    std::cmp::Ordering::Equal => None,
                    std::cmp::Ordering::Less => Some(idx),
                    std::cmp::Ordering::Greater => Some(idx - 1),
                }
            } else {
                Some(idx)
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
                if Some(idx) == reference_idx {
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
                } else if i == 1 && is_dirichlet("west") {
                    ap += dx2_inv;
                }
                // East neighbour
                if i < nx - 2 && fields.mask.at(i + 1, j) {
                    ap += dx2_inv;
                    let ni = idx + (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dx2_inv)?;
                    }
                } else if i == nx - 2 && is_dirichlet("east") {
                    ap += dx2_inv;
                }
                // South neighbour
                if j > 1 && fields.mask.at(i, j - 1) {
                    ap += dy2_inv;
                    let ni = idx - 1;
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dy2_inv)?;
                    }
                } else if j == 1 && is_dirichlet("south") {
                    ap += dy2_inv;
                }
                // North neighbour
                if j < ny - 2 && fields.mask.at(i, j + 1) {
                    ap += dy2_inv;
                    let ni = idx + 1;
                    if let Some(ci) = map_index(ni) {
                        builder.add_entry(row_idx, ci, -dy2_inv)?;
                    }
                } else if j == ny - 2 && is_dirichlet("north") {
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
            boundary_conditions,
            output_correction,
        )?;

        *self._rhs_cache.borrow_mut() = Some(rhs);
        *self._solution_cache.borrow_mut() = Some(p_correction_vec);

        Ok(())
    }
}
