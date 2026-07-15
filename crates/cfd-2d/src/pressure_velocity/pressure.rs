//! Pressure correction solver — struct definition and field corrections
//!
//! # Theorem (Pressure-Velocity Coupling — Patankar & Spalding 1972)
//!
//! For a fixed linearization of the momentum equations, solving the anchored
//! pressure-correction Poisson problem and applying the corresponding velocity
//! update enforces discrete mass conservation on each control volume to the
//! linear-solver tolerance:
//!
//! ```text
//! u = u* − (Vol / a_P) ∇p'   (velocity correction)
//! p = p  + α_p · p'           (pressure correction)
//! ```
//!
//! See [`correction`](super::correction) for the matrix assembly and solve details.

use super::config::PressureLinearSolver;
use crate::fields::Field2D;
use crate::grid::array2d::Array2D;
use crate::grid::StructuredGrid2D;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_math::linear_solver::preconditioners::AlgebraicMultigrid;
use cfd_math::linear_solver::{BiCGSTAB, ConjugateGradient, GMRES};
use cfd_math::sparse::SparseMatrixBuilder;
use eunomia::FloatElement;
use leto::{geometry::Vector2, Array1};
use std::fmt::Debug;

/// Pressure correction solver supporting multiple linear solver backends
pub struct PressureCorrectionSolver<T: Cfd2dScalar + Copy> {
    pub(super) grid: StructuredGrid2D<T>,
    pub(super) solver_type: PressureLinearSolver,
    pub(super) cg_solver: ConjugateGradient<T>,
    pub(super) bicgstab_solver: BiCGSTAB<T>,
    pub(super) gmres_solver: Option<GMRES<T>>,
    pub(super) _amg_preconditioner: std::cell::RefCell<Option<AlgebraicMultigrid<T>>>,
    pub(super) _laplacian_cache: std::cell::RefCell<Option<cfd_math::sparse::SparseMatrix<T>>>,
    pub(super) _rhs_cache: std::cell::RefCell<Option<Array1<T>>>,
    pub(super) _solution_cache: std::cell::RefCell<Option<Array1<T>>>,
    pub(super) _p_correction_cache: std::cell::RefCell<Option<Array2D<T>>>,
    pub(super) _matrix_builder_cache: std::cell::RefCell<Option<SparseMatrixBuilder<T>>>,
}

impl<T: Cfd2dScalar + Copy + Debug + FloatElement> PressureCorrectionSolver<T> {
    /// Create new pressure correction solver with specified linear solver
    ///
    /// ## Recommended Solver Selection
    ///
    /// - **GMRES**: Industry standard for SIMPLE/PISO (default).
    ///   Reference: Saad (2003), Patankar (1980).
    /// - **BiCGSTAB**: Good alternative for non-symmetric systems.
    /// - **CG**: Only for symmetric systems.
    pub fn new(
        grid: StructuredGrid2D<T>,
        solver_type: PressureLinearSolver,
    ) -> cfd_core::error::Result<Self> {
        let config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: 200,
            tolerance: <T as FloatElement>::from_f64(1e-3),
            use_preconditioner: false,
        };

        let gmres_solver = match solver_type {
            PressureLinearSolver::GMRES { restart_dim } => Some(GMRES::new(config, restart_dim)),
            _ => None,
        };

        Ok(Self {
            grid,
            solver_type,
            cg_solver: ConjugateGradient::new(config),
            bicgstab_solver: BiCGSTAB::new(config),
            gmres_solver,
            _amg_preconditioner: std::cell::RefCell::new(None),
            _laplacian_cache: std::cell::RefCell::new(None),
            _rhs_cache: std::cell::RefCell::new(None),
            _solution_cache: std::cell::RefCell::new(None),
            _p_correction_cache: std::cell::RefCell::new(None),
            _matrix_builder_cache: std::cell::RefCell::new(None),
        })
    }

    /// Correct velocity field using pressure correction gradient
    pub fn correct_velocity(
        &self,
        u_star: &mut Array2D<Vector2<T>>,
        p_correction: &Array2D<T>,
        ap_u: &Field2D<T>,
        ap_v: &Field2D<T>,
        _rho: T,
        alpha: T,
        fields: &crate::fields::SimulationFields<T>,
    ) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;
        let volume = dx * dy;
        let two =
            <T as FloatElement>::from_f64(cfd_core::physics::constants::mathematical::numeric::TWO);

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                if !fields.mask.at(i, j) {
                    u_star[(i, j)][0] = scalar::zero::<T>();
                    u_star[(i, j)][1] = scalar::zero::<T>();
                    continue;
                }

                let left_fluid = fields.mask.at(i - 1, j);
                let right_fluid = fields.mask.at(i + 1, j);
                let dp_dx = if left_fluid && right_fluid {
                    (p_correction[(i + 1, j)] - p_correction[(i - 1, j)]) / (two * dx)
                } else if right_fluid {
                    (p_correction[(i + 1, j)] - p_correction[(i, j)]) / (two * dx)
                } else if left_fluid {
                    (p_correction[(i, j)] - p_correction[(i - 1, j)]) / (two * dx)
                } else {
                    scalar::zero::<T>()
                };

                let south_fluid = fields.mask.at(i, j - 1);
                let north_fluid = fields.mask.at(i, j + 1);
                let dp_dy = if south_fluid && north_fluid {
                    (p_correction[(i, j + 1)] - p_correction[(i, j - 1)]) / (two * dy)
                } else if north_fluid {
                    (p_correction[(i, j + 1)] - p_correction[(i, j)]) / (two * dy)
                } else if south_fluid {
                    (p_correction[(i, j)] - p_correction[(i, j - 1)]) / (two * dy)
                } else {
                    scalar::zero::<T>()
                };

                let factor_u = volume / ap_u.at(i, j);
                let factor_v = volume / ap_v.at(i, j);

                u_star[(i, j)][0] -= alpha * factor_u * dp_dx;
                u_star[(i, j)][1] -= alpha * factor_v * dp_dy;
            }
        }
    }

    /// Correct pressure field with under-relaxation
    pub fn correct_pressure(&self, p: &mut Array2D<T>, p_correction: &Array2D<T>, alpha: T) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for i in 0..nx {
            for j in 0..ny {
                p[(i, j)] += alpha * p_correction[(i, j)];
            }
        }
    }
}
