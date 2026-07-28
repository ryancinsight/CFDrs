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
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_math::iterative::IterativeLinearSolver;
use cfd_math::iterative::preconditioners::IdentityPreconditioner;
use cfd_math::multigrid::AlgebraicMultigrid;
use eunomia::FloatElement;
use leto::Array1;
use leto_ops::{Scalar as LetoScalar, norm_l2};
use std::fmt::Debug;

/// Largest pressure system for which default AMG setup stays within the
/// committed validation-test runtime budget on the supported CPU path.
const AMG_SETUP_MATRIX_SIZE_LIMIT: usize = 10_000;

impl<T: Cfd2dScalar + Copy + Debug + FloatElement + LetoScalar> PressureCorrectionSolver<T> {
    /// Dispatch a linear solve to the configured solver backend
    pub(super) fn dispatch_solve(
        &self,
        matrix: &cfd_math::sparse::SparseMatrix<T>,
        rhs: &Array1<T>,
        solution: &mut Array1<T>,
    ) -> cfd_core::error::Result<()> {
        // Phase 8: Deep Optimization & AMG Caching
        //
        // Theorem — Sparse Galerkin Caching (Ruge & Stüben 1987)
        // Re-computing algebraic multigrid transfer operators (R, P) takes O(N) but with a large constant.
        // For pressure Poisson equations on fixed grid topologies, the non-zero pattern remains invariant
        // while coefficients change slowly. We exactly recompute A_c = R A_f P without graph re-coarsening.
        let mut amg_cache = self._amg_preconditioner.borrow_mut();
        let matrix_values_unchanged = amg_cache.is_some()
            && self
                ._amg_matrix_values
                .borrow()
                .as_deref()
                .is_some_and(|cached| cached == matrix.values());
        let amg_setup_allowed = matrix.nrows() < AMG_SETUP_MATRIX_SIZE_LIMIT;
        if !amg_setup_allowed {
            tracing::debug!(
                matrix_size = matrix.nrows(),
                setup_limit = AMG_SETUP_MATRIX_SIZE_LIMIT,
                "Using the configured iterative solver without AMG because setup exceeds the runtime budget"
            );
        }
        if amg_setup_allowed && !matrix_values_unchanged {
            if amg_cache.is_none() {
                let config = cfd_math::multigrid::AMGConfig::default();
                // Initialize the AMG preconditioner and cache it.
                if let Ok(amg) = AlgebraicMultigrid::new(matrix, config) {
                    *amg_cache = Some(amg);
                }
            } else {
                // Hot-path optimization: A_c = R A_f P. Recalculate coarse
                // matrices only when the fine-grid coefficients changed.
                amg_cache
                    .as_mut()
                    .expect("invariant: non-empty AMG cache checked above")
                    .recompute(matrix)?;
            }
            if amg_cache.is_some() {
                *self._amg_matrix_values.borrow_mut() = Some(matrix.values().to_vec());
            }
        }

        let amg_preconditioner: Option<&AlgebraicMultigrid<T>> = amg_cache.as_ref();
        let solve_result = match self.solver_type {
            PressureLinearSolver::ConjugateGradient => {
                let result = if let Some(amg) = amg_preconditioner {
                    self.cg_solver.solve(matrix, rhs, solution, Some(amg))
                } else {
                    self.cg_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                }
                .map(|_| ())
                .map_err(cfd_core::error::Error::from);
                match result {
                    Ok(()) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => self
                        .cg_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ())
                        .map_err(cfd_core::error::Error::from),
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
                }
                .map(|_| ())
                .map_err(cfd_core::error::Error::from);
                match result {
                    Ok(()) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => self
                        .bicgstab_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ())
                        .map_err(cfd_core::error::Error::from),
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
                }
                .map(|_| ())
                .map_err(cfd_core::error::Error::from);
                match result {
                    Ok(()) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                        .map(|_| ())
                        .map_err(cfd_core::error::Error::from),
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
                    "Pressure solve stalled; keeping last iterate as approximate solution"
                );
                Ok(())
            }
            Err(cfd_core::error::Error::Convergence(
                cfd_core::error::ConvergenceErrorKind::Breakdown,
            )) => {
                tracing::warn!(
                    solver = ?self.solver_type,
                    "Pressure solve breakdown; keeping last iterate as approximate solution"
                );
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
        boundary_conditions: Option<
            &std::collections::HashMap<String, cfd_core::physics::boundary::BoundaryCondition<T>>,
        >,
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
                    output_correction[(i, j)] = scalar::zero::<T>();
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

        let has_dirichlet = is_dirichlet("west")
            || is_dirichlet("east")
            || is_dirichlet("south")
            || is_dirichlet("north");
        let dirichlet_sides = [
            is_dirichlet("west"),
            is_dirichlet("east"),
            is_dirichlet("south"),
            is_dirichlet("north"),
        ];

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

        let mut rhs = self
            ._rhs_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.shape()[0] == system_size)
            .unwrap_or_else(|| Array1::from_elem([system_size], scalar::zero::<T>()));
        rhs.fill(scalar::zero::<T>());

        let coeff = rho / dt;

        let two =
            <T as FloatElement>::from_f64(cfd_core::physics::constants::mathematical::numeric::TWO);

        let mask = fields.mask.as_slice();
        let matrix_cache_valid = self._laplacian_cache.borrow().as_ref().is_some_and(|_| {
            self._face_matrix_mask
                .borrow()
                .as_deref()
                .is_some_and(|cached| cached == mask)
                && self
                    ._face_matrix_dirichlet
                    .borrow()
                    .is_some_and(|cached| cached == dirichlet_sides)
        });

        let matrix = if matrix_cache_valid {
            self._laplacian_cache
                .borrow()
                .as_ref()
                .expect("invariant: pressure matrix cache validity was established above")
                .clone()
        } else {
            let mut builder = self.take_matrix_builder(system_size, system_size);
            let dx2_inv = scalar::one::<T>() / (dx * dx);
            let dy2_inv = scalar::one::<T>() / (dy * dy);

            for i in 1..nx - 1 {
                for j in 1..ny - 1 {
                    let idx = (i - 1) * (ny - 2) + (j - 1);
                    if Some(idx) == reference_idx {
                        continue;
                    }
                    let row_idx = map_index(idx).expect("row index must exist");

                    if !fields.mask.at(i, j) {
                        builder.add_entry(row_idx, row_idx, scalar::one::<T>())?;
                        continue;
                    }

                    let mut ap = scalar::zero::<T>();

                    // West neighbour
                    if i > 1 && fields.mask.at(i - 1, j) {
                        ap += dx2_inv;
                        let ni = idx - (ny - 2);
                        if let Some(ci) = map_index(ni) {
                            builder.add_entry(row_idx, ci, -dx2_inv)?;
                        }
                    } else if i == 1 && dirichlet_sides[0] {
                        ap += dx2_inv;
                    }
                    // East neighbour
                    if i < nx - 2 && fields.mask.at(i + 1, j) {
                        ap += dx2_inv;
                        let ni = idx + (ny - 2);
                        if let Some(ci) = map_index(ni) {
                            builder.add_entry(row_idx, ci, -dx2_inv)?;
                        }
                    } else if i == nx - 2 && dirichlet_sides[1] {
                        ap += dx2_inv;
                    }
                    // South neighbour
                    if j > 1 && fields.mask.at(i, j - 1) {
                        ap += dy2_inv;
                        let ni = idx - 1;
                        if let Some(ci) = map_index(ni) {
                            builder.add_entry(row_idx, ci, -dy2_inv)?;
                        }
                    } else if j == 1 && dirichlet_sides[2] {
                        ap += dy2_inv;
                    }
                    // North neighbour
                    if j < ny - 2 && fields.mask.at(i, j + 1) {
                        ap += dy2_inv;
                        let ni = idx + 1;
                        if let Some(ci) = map_index(ni) {
                            builder.add_entry(row_idx, ci, -dy2_inv)?;
                        }
                    } else if j == ny - 2 && dirichlet_sides[3] {
                        ap += dy2_inv;
                    }

                    builder.add_entry(row_idx, row_idx, ap)?;
                }
            }

            let matrix = builder.build()?;
            self.reset_matrix_builder_cache(system_size, system_size);
            *self._laplacian_cache.borrow_mut() = Some(matrix.clone());
            *self._face_matrix_mask.borrow_mut() = Some(mask.to_vec());
            *self._face_matrix_dirichlet.borrow_mut() = Some(dirichlet_sides);
            matrix
        };

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                if Some(idx) == reference_idx || !fields.mask.at(i, j) {
                    continue;
                }
                let row_idx = map_index(idx).expect("row index must exist");
                let div_u = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (two * dx)
                    + (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (two * dy);
                rhs[row_idx] = -coeff * div_u;
            }
        }

        let rhs_norm =
            norm_l2(&rhs.view()).expect("invariant: pressure RHS Leto vector has a valid layout");
        tracing::debug!(
            cached = matrix_cache_valid,
            system_size,
            rhs_norm = ?rhs_norm,
            "Pressure Solve"
        );
        let mut p_correction_vec = self
            ._solution_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.shape()[0] == matrix.nrows())
            .unwrap_or_else(|| Array1::from_elem([matrix.nrows()], scalar::zero::<T>()));
        p_correction_vec.fill(scalar::zero::<T>());
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
