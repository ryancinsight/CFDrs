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
use cfd_math::linear_solver::preconditioners::{AlgebraicMultigrid, IdentityPreconditioner};
use cfd_math::linear_solver::IterativeLinearSolver;
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::fmt::Debug;

impl<T: RealField + Copy + FromPrimitive + Debug> PressureCorrectionSolver<T> {
    /// Dispatch a linear solve to the configured solver backend
    pub(super) fn dispatch_solve(
        &self,
        matrix: &cfd_math::sparse::SparseMatrix<T>,
        rhs: &DVector<T>,
        solution: &mut DVector<T>,
    ) -> cfd_core::error::Result<()> {
        let amg_preconditioner: Option<AlgebraicMultigrid<T>> = None;

        match self.solver_type {
            PressureLinearSolver::ConjugateGradient => {
                let result = if let Some(ref amg) = amg_preconditioner {
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
                let result = if let Some(ref amg) = amg_preconditioner {
                    self.bicgstab_solver
                        .solve(matrix, rhs, solution, Some(amg))
                } else {
                    self.bicgstab_solver
                        .solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                };
                match result {
                    Ok(_) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => self.bicgstab_solver.solve(
                        matrix,
                        rhs,
                        solution,
                        None::<&IdentityPreconditioner>,
                    ).map(|_| ()),
                    Err(e) => Err(e),
                }
            }
            PressureLinearSolver::GMRES { .. } => {
                let Some(ref solver) = self.gmres_solver else {
                    return Err(cfd_core::error::Error::InvalidConfiguration(
                        "GMRES solver not initialized".to_string(),
                    ));
                };
                let result = if let Some(ref amg) = amg_preconditioner {
                    solver.solve(matrix, rhs, solution, Some(amg))
                } else {
                    solver.solve(matrix, rhs, solution, None::<&IdentityPreconditioner>)
                };
                match result {
                    Ok(_) => Ok(()),
                    Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::Breakdown,
                    )) if amg_preconditioner.is_some() => solver.solve(
                        matrix,
                        rhs,
                        solution,
                        None::<&IdentityPreconditioner>,
                    ).map(|_| ()),
                    Err(e) => Err(e),
                }
            }
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
    ) -> cfd_core::error::Result<Vec<Vec<T>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        let n = (nx - 2) * (ny - 2);

        if n <= 1 {
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }

        let reference_idx = 0usize;
        let system_size = n - 1;

        let map_index = |idx: usize| -> Option<usize> {
            match idx.cmp(&reference_idx) {
                std::cmp::Ordering::Equal => None,
                std::cmp::Ordering::Less => Some(idx),
                std::cmp::Ordering::Greater => Some(idx - 1),
            }
        };

        let mut builder = SparseMatrixBuilder::new(system_size, system_size);
        let mut rhs = DVector::zeros(system_size);

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

                let ap = two * (dx2_inv + dy2_inv);
                builder.add_entry(row_idx, row_idx, ap)?;

                // West neighbour
                if i > 1 {
                    let ni = idx - (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i - 1, j) {
                            builder.add_entry(row_idx, ci, -dx2_inv)?;
                        }
                    }
                }
                // East neighbour
                if i < nx - 2 {
                    let ni = idx + (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i + 1, j) {
                            builder.add_entry(row_idx, ci, -dx2_inv)?;
                        }
                    }
                }
                // South neighbour
                if j > 1 {
                    let ni = idx - 1;
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i, j - 1) {
                            builder.add_entry(row_idx, ci, -dy2_inv)?;
                        }
                    }
                }
                // North neighbour
                if j < ny - 2 {
                    let ni = idx + 1;
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i, j + 1) {
                            builder.add_entry(row_idx, ci, -dy2_inv)?;
                        }
                    }
                }

                let div_u = (fields.u.at(i + 1, j) - fields.u.at(i - 1, j)) / (two * dx)
                    + (fields.v.at(i, j + 1) - fields.v.at(i, j - 1)) / (two * dy);
                rhs[row_idx] = -coeff * div_u;
            }
        }

        let rhs_norm = rhs.norm();
        eprintln!(
            "Pressure Solve: n_fluid={n_fluid}, n_solid={n_solid}, \
             system_size={system_size}, rhs_norm={rhs_norm:?}"
        );

        let matrix = builder.build()?;
        let mut p_correction_vec = DVector::zeros(matrix.nrows());
        self.dispatch_solve(&matrix, &rhs, &mut p_correction_vec)?;

        self.scatter_correction(&p_correction_vec, reference_idx, &map_index)
    }

    /// Solve pressure correction equation using face velocities (Rhie-Chow)
    pub fn solve_pressure_correction_from_faces(
        &self,
        u_face: &[Vec<T>],
        v_face: &[Vec<T>],
        d_x: &[Vec<T>],
        d_y: &[Vec<T>],
        rho: T,
        fields: &crate::fields::SimulationFields<T>,
    ) -> cfd_core::error::Result<Vec<Vec<T>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        let n = (nx - 2) * (ny - 2);
        if n <= 1 {
            return Ok(vec![vec![T::zero(); ny]; nx]);
        }

        let reference_idx = 0usize;
        let system_size = n - 1;

        let map_index = |idx: usize| -> Option<usize> {
            match idx.cmp(&reference_idx) {
                std::cmp::Ordering::Equal => None,
                std::cmp::Ordering::Less => Some(idx),
                std::cmp::Ordering::Greater => Some(idx - 1),
            }
        };

        let mut builder = SparseMatrixBuilder::new(system_size, system_size);
        let mut rhs = DVector::zeros(system_size);

        let dx2_inv = T::one() / (dx * dx);
        let dy2_inv = T::one() / (dy * dy);

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

                let de = d_x[i][j];
                let dw = d_x[i - 1][j];
                let dn = d_y[i][j];
                let ds = d_y[i][j - 1];

                let ae = de * dx2_inv;
                let aw = dw * dx2_inv;
                let an = dn * dy2_inv;
                let as_ = ds * dy2_inv;
                let ap = ae + aw + an + as_;

                builder.add_entry(row_idx, row_idx, ap)?;

                if i > 1 {
                    let ni = idx - (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i - 1, j) {
                            builder.add_entry(row_idx, ci, -aw)?;
                        }
                    }
                }
                if i < nx - 2 {
                    let ni = idx + (ny - 2);
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i + 1, j) {
                            builder.add_entry(row_idx, ci, -ae)?;
                        }
                    }
                }
                if j > 1 {
                    let ni = idx - 1;
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i, j - 1) {
                            builder.add_entry(row_idx, ci, -as_)?;
                        }
                    }
                }
                if j < ny - 2 {
                    let ni = idx + 1;
                    if let Some(ci) = map_index(ni) {
                        if fields.mask.at(i, j + 1) {
                            builder.add_entry(row_idx, ci, -an)?;
                        }
                    }
                }

                let div_u = (u_face[i][j] - u_face[i - 1][j]) / dx
                    + (v_face[i][j] - v_face[i][j - 1]) / dy;
                rhs[row_idx] = -rho * div_u;
            }
        }

        let matrix = builder.build()?;
        let mut p_correction_vec = DVector::zeros(matrix.nrows());
        self.dispatch_solve(&matrix, &rhs, &mut p_correction_vec)?;

        self.scatter_correction(&p_correction_vec, reference_idx, &map_index)
    }

    /// Scatter solution vector back to 2D grid with Neumann BCs
    fn scatter_correction(
        &self,
        solution: &DVector<T>,
        reference_idx: usize,
        map_index: &dyn Fn(usize) -> Option<usize>,
    ) -> cfd_core::error::Result<Vec<Vec<T>>> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;

        let mut p_correction = vec![vec![T::zero(); ny]; nx];
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
                p_correction[i][j] = value;
            }
        }

        // Neumann BCs (zero gradient)
        for i in 0..nx {
            p_correction[i][0] = p_correction[i][1];
            p_correction[i][ny - 1] = p_correction[i][ny - 2];
        }
        for j in 0..ny {
            p_correction[0][j] = p_correction[1][j];
            p_correction[nx - 1][j] = p_correction[nx - 2][j];
        }

        Ok(p_correction)
    }
}
