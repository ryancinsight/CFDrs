//! Tiered linear solver fallback chain for saddle-point systems.
//!
//! # Algorithm — Tiered Linear Solver Fallback (Robust Sparse Solve)
//!
//! Attempts solvers in priority order; returns the first successful solution.
//! Rationale: saddle-point systems from mixed FEM (Taylor-Hood P2-P1) are
//! non-SPD; the best preconditioner depends on problem size and conditioning.
//!
//! ## Algorithm (Solve Priority Order)
//!
//! ```text
//! 1. DirectSparseSolver (LU)             — exact, O(n^1.5), used when n < threshold
//! 2. GMRES + BlockDiagonalPreconditioner — best for large saddle-point systems
//! 3. GMRES (unpreconditioned)            — fallback when block preconditioner fails
//! 4. GMRES + ILU preconditioner          — fallback with incomplete LU
//! 5. BiCGSTAB (unpreconditioned)         — last resort for extreme non-symmetry
//! ```
//!
//! **Rationale.** Direct LU is exact and preferred for small systems; for large
//! saddle-point systems from mixed FEM the block diagonal preconditioner
//! (Elman, Silvester & Wathen 2014, §3.4) exploits the 2×2 block structure
//! to give mesh-independent GMRES convergence.  The ILU fallback addresses
//! highly anisotropic or nearly singular systems where block preconditioning
//! fails.  BiCGSTAB provides a last resort for strongly non-normal operators.
//!
//! ## References
//!
//! - Elman, H., Silvester, D. & Wathen, A. (2014). *Finite Elements and Fast
//!   Iterative Solvers.* Oxford University Press, §3.4.
//! - Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.).
//!   SIAM, Chapter 6.
//! - Benzi, M., Golub, G.H. & Liesen, J. (2005). "Numerical solution of
//!   saddle point problems." *Acta Numerica* 14:1–137.

use crate::linear_solver::{
    BiCGSTAB, BlockDiagonalPreconditioner, DirectSparseSolver, GMRES, IncompleteLU,
    IterativeSolverConfig,
};
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use num_traits::{Float, FromPrimitive};
use std::fmt::Debug;

use crate::sparse::SparseMatrix;

/// Tiered linear solver fallback chain for saddle-point systems arising from
/// mixed FEM discretizations of incompressible flow.
///
/// The chain attempts solvers in priority order (direct LU → GMRES/block →
/// GMRES/unpreconditioned → GMRES/ILU → BiCGSTAB), returning the first
/// successful solution.  This eliminates duplicated solver fallback logic
/// across domain-specific solvers.
///
/// # Usage
///
/// ```rust,ignore
/// let chain = LinearSolverChain::new(solver_config);
/// let x = chain.solve(&matrix, &rhs, n_velocity_dof)?;
/// ```
pub struct LinearSolverChain<T: RealField + Copy + Float + FromPrimitive + Debug> {
    /// Configuration (tolerance, max iterations) for all iterative solvers.
    config: IterativeSolverConfig<T>,
    /// DOF count below which the direct LU solver is preferred.
    ///
    /// Default: 100,000.  Direct LU is exact but has O(n^1.5) cost;
    /// for large systems iterative solvers are more efficient.
    direct_threshold: usize,
    /// GMRES restart parameter m (maximum Krylov subspace dimension).
    ///
    /// Default: 100.  Larger values reduce restart overhead but increase
    /// memory O(n·m); for saddle-point systems m ∈ [50, 200] is typical.
    krylov_restart: usize,
}

impl<T: RealField + Copy + Float + FromPrimitive + Debug> LinearSolverChain<T> {
    /// Create a new solver chain with the given iterative solver configuration.
    ///
    /// Defaults: direct_threshold = 100,000;  krylov_restart = 100.
    #[must_use]
    pub fn new(config: IterativeSolverConfig<T>) -> Self {
        Self {
            config,
            direct_threshold: 100_000,
            krylov_restart: 100,
        }
    }

    /// Override the DOF count threshold below which direct LU is used.
    #[must_use]
    pub fn with_direct_threshold(mut self, threshold: usize) -> Self {
        self.direct_threshold = threshold;
        self
    }

    /// Override the GMRES restart parameter.
    #[must_use]
    pub fn with_krylov_restart(mut self, restart: usize) -> Self {
        self.krylov_restart = restart;
        self
    }

    /// Solve the sparse linear system `A·x = b` using a tiered fallback strategy.
    ///
    /// # Algorithm
    ///
    /// Attempts solvers in order:
    /// 1. **Direct LU** — used when `rhs.len() < direct_threshold` (exact, fast for small n).
    /// 2. **GMRES + BlockDiagonal** — exploits 2×2 saddle-point block structure.
    /// 3. **GMRES (unpreconditioned)** — fallback if block preconditioner fails.
    /// 4. **GMRES + ILU** — fallback for highly anisotropic systems.
    /// 5. **BiCGSTAB (unpreconditioned)** — last resort for strongly non-normal operators.
    ///
    /// # Arguments
    /// * `matrix` — Sparse coefficient matrix A (CSR format)
    /// * `rhs` — Right-hand side vector b
    /// * `n_velocity_dof` — Number of velocity DOFs; determines the 2×2 block split
    ///   for the BlockDiagonal preconditioner (velocity block = `n_velocity_dof`,
    ///   pressure block = `rhs.len() − n_velocity_dof`)
    ///
    /// # Errors
    /// Returns `Error::Solver` only if all five solver tiers fail.
    pub fn solve(
        &self,
        matrix: &SparseMatrix<T>,
        rhs: &DVector<T>,
        n_velocity_dof: usize,
    ) -> Result<DVector<T>> {
        let n_total_dof = rhs.len();
        let n_pressure_dof = n_total_dof.saturating_sub(n_velocity_dof);
        let mut x = DVector::zeros(n_total_dof);

        // ── Tier 1: Direct sparse LU (exact, preferred for small systems) ─────
        if n_total_dof < self.direct_threshold {
            let direct = DirectSparseSolver::default();
            match direct.solve(matrix, rhs) {
                Ok(x_direct) => {
                    tracing::debug!("LinearSolverChain: direct LU succeeded (n={n_total_dof})");
                    return Ok(x_direct);
                }
                Err(e) => {
                    tracing::warn!("LinearSolverChain: direct LU failed ({e}); trying iterative");
                }
            }
        }

        let restart = std::cmp::min(self.krylov_restart, n_total_dof.max(1));
        let solver = GMRES::new(self.config, restart);

        // ── Tier 2: GMRES + BlockDiagonal preconditioner (saddle-point) ───────
        match BlockDiagonalPreconditioner::new(matrix, n_velocity_dof, n_pressure_dof) {
            Ok(block_precond) => {
                match solver.solve_preconditioned(matrix, rhs, &block_precond, &mut x) {
                    Ok(monitor) => {
                        tracing::debug!(
                            "LinearSolverChain: GMRES+BlockDiag converged in {} iters",
                            monitor.iteration
                        );
                        return Ok(x);
                    }
                    Err(e) => {
                        tracing::warn!("LinearSolverChain: GMRES+BlockDiag failed ({e})");
                    }
                }
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain: BlockDiag preconditioner failed ({e})");
            }
        }

        // ── Tier 3: GMRES unpreconditioned ────────────────────────────────────
        x.fill(T::zero());
        match solver.solve_unpreconditioned(matrix, rhs, &mut x) {
            Ok(monitor) => {
                tracing::debug!(
                    "LinearSolverChain: GMRES (unpreconditioned) converged in {} iters",
                    monitor.iteration
                );
                return Ok(x);
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain: GMRES unpreconditioned failed ({e})");
            }
        }

        // ── Tier 4: GMRES + ILU preconditioner ───────────────────────────────
        x.fill(T::zero());
        match IncompleteLU::new(matrix) {
            Ok(ilu) => match solver.solve_preconditioned(matrix, rhs, &ilu, &mut x) {
                Ok(monitor) => {
                    tracing::debug!(
                        "LinearSolverChain: GMRES+ILU converged in {} iters",
                        monitor.iteration
                    );
                    return Ok(x);
                }
                Err(e) => {
                    tracing::warn!("LinearSolverChain: GMRES+ILU failed ({e})");
                }
            },
            Err(e) => {
                tracing::warn!("LinearSolverChain: ILU construction failed ({e})");
            }
        }

        // ── Tier 5: BiCGSTAB (last resort) ────────────────────────────────────
        x.fill(T::zero());
        let bicg = BiCGSTAB::new(self.config);
        bicg.solve_unpreconditioned(matrix, rhs, &mut x)
            .map_err(|e| {
                Error::Solver(
                    format!(
                        "LinearSolverChain: all solver tiers failed. \
                         Final BiCGSTAB error: {e}"
                    )
                )
            })?;

        tracing::debug!("LinearSolverChain: BiCGSTAB (last resort) converged");
        Ok(x)
    }

    /// Solve with optional warm-start initial guess for Picard/continuation methods.
    ///
    /// When `initial_guess` is `Some`, iterative solvers begin from that vector
    /// rather than zero, dramatically reducing iteration counts for successive
    /// solves on slowly-varying systems (e.g., Picard viscosity updates).
    ///
    /// Additional optimizations over [`Self::solve`]:
    /// - On tier failure, resets to `initial_guess` (not zero) before next tier.
    /// - Skips unpreconditioned GMRES (Tier 3) when block preconditioning was
    ///   constructed but GMRES stagnated — for saddle-point systems unpreconditioned
    ///   is strictly worse.
    pub fn solve_with_guess(
        &self,
        matrix: &SparseMatrix<T>,
        rhs: &DVector<T>,
        n_velocity_dof: usize,
        initial_guess: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n_total_dof = rhs.len();
        let n_pressure_dof = n_total_dof.saturating_sub(n_velocity_dof);
        let mut x = initial_guess.cloned().unwrap_or_else(|| DVector::zeros(n_total_dof));

        let reset = |x: &mut DVector<T>| {
            if let Some(guess) = initial_guess {
                x.copy_from(guess);
            } else {
                x.fill(T::zero());
            }
        };

        // ── Tier 1: Direct sparse LU ──────────────────────────────────────────
        if n_total_dof < self.direct_threshold {
            let direct = DirectSparseSolver::default();
            match direct.solve(matrix, rhs) {
                Ok(x_direct) => {
                    tracing::debug!("LinearSolverChain: direct LU succeeded (n={n_total_dof})");
                    return Ok(x_direct);
                }
                Err(e) => {
                    tracing::warn!("LinearSolverChain: direct LU failed ({e}); trying iterative");
                }
            }
        }

        let restart = std::cmp::min(self.krylov_restart, n_total_dof.max(1));
        let solver = GMRES::new(self.config, restart);
        let mut block_precond_constructed = false;

        // ── Tier 2: GMRES + BlockDiagonal preconditioner ──────────────────────
        match BlockDiagonalPreconditioner::new(matrix, n_velocity_dof, n_pressure_dof) {
            Ok(block_precond) => {
                block_precond_constructed = true;
                match solver.solve_preconditioned(matrix, rhs, &block_precond, &mut x) {
                    Ok(monitor) => {
                        tracing::debug!(
                            "LinearSolverChain(warm): GMRES+BlockDiag converged in {} iters",
                            monitor.iteration
                        );
                        return Ok(x);
                    }
                    Err(e) => {
                        tracing::warn!("LinearSolverChain(warm): GMRES+BlockDiag failed ({e})");
                        reset(&mut x);
                    }
                }
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain(warm): BlockDiag construction failed ({e})");
            }
        }

        // ── Tier 3: GMRES unpreconditioned ────────────────────────────────────
        // Skip if block preconditioner was built but GMRES stagnated — for
        // saddle-point systems, unpreconditioned will be strictly worse.
        if !block_precond_constructed {
            match solver.solve_unpreconditioned(matrix, rhs, &mut x) {
                Ok(monitor) => {
                    tracing::debug!(
                        "LinearSolverChain(warm): GMRES unpreconditioned converged in {} iters",
                        monitor.iteration
                    );
                    return Ok(x);
                }
                Err(e) => {
                    tracing::warn!("LinearSolverChain(warm): GMRES unpreconditioned failed ({e})");
                    reset(&mut x);
                }
            }
        }

        // ── Tier 4: GMRES + ILU preconditioner ───────────────────────────────
        reset(&mut x);
        match IncompleteLU::new(matrix) {
            Ok(ilu) => match solver.solve_preconditioned(matrix, rhs, &ilu, &mut x) {
                Ok(monitor) => {
                    tracing::debug!(
                        "LinearSolverChain(warm): GMRES+ILU converged in {} iters",
                        monitor.iteration
                    );
                    return Ok(x);
                }
                Err(e) => {
                    tracing::warn!("LinearSolverChain(warm): GMRES+ILU failed ({e})");
                }
            },
            Err(e) => {
                tracing::warn!("LinearSolverChain(warm): ILU construction failed ({e})");
            }
        }

        // ── Tier 5: BiCGSTAB (last resort) ────────────────────────────────────
        reset(&mut x);
        let bicg = BiCGSTAB::new(self.config);
        bicg.solve_unpreconditioned(matrix, rhs, &mut x)
            .map_err(|e| {
                Error::Solver(
                    format!(
                        "LinearSolverChain: all solver tiers failed. \
                         Final BiCGSTAB error: {e}"
                    )
                )
            })?;

        tracing::debug!("LinearSolverChain(warm): BiCGSTAB (last resort) converged");
        Ok(x)
    }
}
