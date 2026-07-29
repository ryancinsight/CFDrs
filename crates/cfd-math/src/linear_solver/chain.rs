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

use crate::linear_solver::krylov::{self, SolveOutcome};
use crate::linear_solver::preconditioners::multigrid::AMGHierarchy;
use crate::linear_solver::{
    AMGConfig, AlgebraicMultigrid, BlockDiagonalPreconditioner, DirectSparseSolver,
    IterativeSolverConfig,
};
use athena_leto::IncompleteLu;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use leto_ops::RealScalar as LetoRealScalar;
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
pub struct LinearSolverChain<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> {
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

/// Reusable state for successive solves with the same sparse matrix pattern.
///
/// The state caches only AMG transfer operators. Matrix values are rebuilt for
/// every solve, while exact CSR row and column arrays gate reuse across mesh or
/// boundary-topology changes.
pub struct LinearSolverState<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> {
    amg: Option<CachedAmg<T>>,
}

struct CachedAmg<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> {
    hierarchy: AMGHierarchy<T>,
    row_ptr: Vec<usize>,
    col_indices: Vec<usize>,
}

impl<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> Default for LinearSolverState<T> {
    fn default() -> Self {
        Self { amg: None }
    }
}

impl<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> CachedAmg<T> {
    fn matches(&self, matrix: &SparseMatrix<T>) -> bool {
        self.row_ptr == matrix.row_ptr() && self.col_indices == matrix.col_indices()
    }
}

impl<T: RealField + Copy + FloatElement + LetoRealScalar + Debug> LinearSolverChain<T> {
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
        rhs: &Array1<T>,
        n_velocity_dof: usize,
    ) -> Result<Array1<T>> {
        let n_total_dof = rhs.shape()[0];
        let n_pressure_dof = n_total_dof.saturating_sub(n_velocity_dof);
        let mut x = Array1::zeros([n_total_dof]);

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

        // ── Tier 2: GMRES + Algebraic Multigrid (AMG) ─────────────────────────
        // AMG is exceptionally fast for Poisson and SPD systems (like pressure correction
        // in fractional step algorithms). For saddle-point systems it may fail setup or diverge,
        // safely falling back to block preconditioning.
        match AlgebraicMultigrid::new(matrix, AMGConfig::default()) {
            Ok(amg) => {
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChain: GMRES+AMG",
                    krylov::gmres_preconditioned(matrix, rhs, &amg, &mut x, &self.config, restart),
                ) {
                    tracing::debug!(
                        "LinearSolverChain: GMRES+AMG converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
            }
            Err(e) => {
                tracing::debug!("LinearSolverChain: AMG setup skipped/failed ({e})");
            }
        }

        let restart = std::cmp::min(self.krylov_restart, n_total_dof.max(1));

        // ── Tier 3: GMRES + BlockDiagonal preconditioner (saddle-point) ───────
        match BlockDiagonalPreconditioner::new(matrix, n_velocity_dof, n_pressure_dof) {
            Ok(block_precond) => {
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChain: GMRES+BlockDiag",
                    krylov::gmres_preconditioned(
                        matrix,
                        rhs,
                        &block_precond,
                        &mut x,
                        &self.config,
                        restart,
                    ),
                ) {
                    tracing::debug!(
                        "LinearSolverChain: GMRES+BlockDiag converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain: BlockDiag preconditioner failed ({e})");
            }
        }

        // ── Tier 4: GMRES unpreconditioned ────────────────────────────────────
        x.fill(<T as NumericElement>::ZERO);
        if let Some(report) = krylov::converged_or_none(
            "LinearSolverChain: GMRES unpreconditioned",
            krylov::gmres(matrix, rhs, &mut x, &self.config, restart),
        ) {
            tracing::debug!(
                "LinearSolverChain: GMRES (unpreconditioned) converged in {} iters",
                report.iterations
            );
            return Ok(x);
        }

        // ── Tier 5: GMRES + ILU preconditioner ───────────────────────────────
        x.fill(<T as NumericElement>::ZERO);
        match IncompleteLu::from_csr(matrix) {
            Ok(ilu) => {
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChain: GMRES+ILU",
                    krylov::gmres_preconditioned(matrix, rhs, &ilu, &mut x, &self.config, restart),
                ) {
                    tracing::debug!(
                        "LinearSolverChain: GMRES+ILU converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain: ILU construction failed ({e})");
            }
        }

        // ── Tier 6: BiCGSTAB (last resort) ────────────────────────────────────
        x.fill(<T as NumericElement>::ZERO);
        // Last resort: unlike the tiers above this has nowhere to fall through
        // to, so a non-convergence fails the chain. The specific outcome is
        // logged by `converged_or_none`.
        if krylov::converged_or_none(
            "LinearSolverChain: BiCGSTAB (last resort)",
            krylov::bicgstab(matrix, rhs, &mut x, &self.config),
        )
        .is_none()
        {
            return Err(Error::Solver(
                "LinearSolverChain: all solver tiers failed".to_string(),
            ));
        }

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
        rhs: &Array1<T>,
        n_velocity_dof: usize,
        initial_guess: Option<&Array1<T>>,
    ) -> Result<Array1<T>> {
        let mut state = LinearSolverState::default();
        self.solve_with_guess_state(matrix, rhs, n_velocity_dof, initial_guess, &mut state)
    }

    /// Solve with a warm start and reusable AMG state.
    ///
    /// `state` may be retained across Picard or continuation solves. AMG
    /// transfer operators are reused only when the complete CSR sparsity
    /// pattern matches; changing mesh connectivity or constrained topology
    /// rebuilds the hierarchy before solving.
    pub fn solve_with_guess_state(
        &self,
        matrix: &SparseMatrix<T>,
        rhs: &Array1<T>,
        n_velocity_dof: usize,
        initial_guess: Option<&Array1<T>>,
        state: &mut LinearSolverState<T>,
    ) -> Result<Array1<T>> {
        let n_total_dof = rhs.shape()[0];
        let n_pressure_dof = n_total_dof.saturating_sub(n_velocity_dof);
        let initial_guess_vector = initial_guess.cloned();
        let mut x = initial_guess
            .cloned()
            .unwrap_or_else(|| Array1::zeros([n_total_dof]));

        let reset = |x: &mut Array1<T>| {
            if let Some(guess) = initial_guess_vector.as_ref() {
                for row in 0..guess.shape()[0] {
                    x[row] = guess[row];
                }
            } else {
                x.fill(<T as NumericElement>::ZERO);
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

        // ── Tier 2: GMRES + Algebraic Multigrid (AMG) ─────────────────────────
        let cached_hierarchy = state
            .amg
            .take()
            .filter(|cached| cached.matches(matrix))
            .map(|cached| cached.hierarchy);
        let amg = match cached_hierarchy {
            Some(hierarchy) => {
                match AlgebraicMultigrid::with_hierarchy(matrix, AMGConfig::default(), hierarchy) {
                    Ok(amg) => {
                        tracing::debug!("LinearSolverChain(warm): reused AMG hierarchy");
                        Ok(amg)
                    }
                    Err(e) => {
                        tracing::warn!(
                        "LinearSolverChain(warm): cached AMG hierarchy failed ({e}); rebuilding"
                    );
                        AlgebraicMultigrid::new(matrix, AMGConfig::default())
                    }
                }
            }
            None => AlgebraicMultigrid::new(matrix, AMGConfig::default()),
        };
        match amg {
            Ok(amg) => {
                state.amg = Some(CachedAmg {
                    hierarchy: amg.get_hierarchy(),
                    row_ptr: matrix.row_ptr().to_vec(),
                    col_indices: matrix.col_indices().to_vec(),
                });
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChain(warm): GMRES+AMG",
                    krylov::gmres_preconditioned(matrix, rhs, &amg, &mut x, &self.config, restart),
                ) {
                    tracing::debug!(
                        "LinearSolverChain(warm): GMRES+AMG converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
                reset(&mut x);
            }
            Err(e) => {
                tracing::debug!("LinearSolverChain(warm): AMG setup skipped/failed ({e})");
            }
        }

        let restart = std::cmp::min(self.krylov_restart, n_total_dof.max(1));
        let mut block_precond_constructed = false;

        // ── Tier 3: GMRES + BlockDiagonal preconditioner ──────────────────────
        match BlockDiagonalPreconditioner::new(matrix, n_velocity_dof, n_pressure_dof) {
            Ok(block_precond) => {
                block_precond_constructed = true;
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChain(warm): GMRES+BlockDiag",
                    krylov::gmres_preconditioned(
                        matrix,
                        rhs,
                        &block_precond,
                        &mut x,
                        &self.config,
                        restart,
                    ),
                ) {
                    tracing::debug!(
                        "LinearSolverChain(warm): GMRES+BlockDiag converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
                reset(&mut x);
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain(warm): BlockDiag construction failed ({e})");
            }
        }

        // ── Tier 4: GMRES unpreconditioned ────────────────────────────────────
        // Skip if block preconditioner was built but GMRES stagnated — for
        // saddle-point systems, unpreconditioned will be strictly worse.
        if !block_precond_constructed {
            if let Some(report) = krylov::converged_or_none(
                "LinearSolverChain(warm): GMRES unpreconditioned",
                krylov::gmres(matrix, rhs, &mut x, &self.config, restart),
            ) {
                tracing::debug!(
                    "LinearSolverChain(warm): GMRES unpreconditioned converged in {} iters",
                    report.iterations
                );
                return Ok(x);
            }
            reset(&mut x);
        }

        // ── Tier 5: GMRES + ILU preconditioner ───────────────────────────────
        reset(&mut x);
        match IncompleteLu::from_csr(matrix) {
            Ok(ilu) => {
                if let Some(report) = krylov::converged_or_none(
                    "LinearSolverChainwarm): GMRES+ILU",
                    krylov::gmres_preconditioned(matrix, rhs, &ilu, &mut x, &self.config, restart),
                ) {
                    tracing::debug!(
                        "LinearSolverChain(warm): GMRES+ILU converged in {} iters",
                        report.iterations
                    );
                    return Ok(x);
                }
            }
            Err(e) => {
                tracing::warn!("LinearSolverChain(warm): ILU construction failed ({e})");
            }
        }

        // ── Tier 6: BiCGSTAB (last resort) ────────────────────────────────────
        reset(&mut x);
        // Last resort: unlike the tiers above this has nowhere to fall through
        // to, so a non-convergence fails the chain. The specific outcome is
        // logged by `converged_or_none`.
        if krylov::converged_or_none(
            "LinearSolverChain: BiCGSTAB (last resort)",
            krylov::bicgstab(matrix, rhs, &mut x, &self.config),
        )
        .is_none()
        {
            return Err(Error::Solver(
                "LinearSolverChain: all solver tiers failed".to_string(),
            ));
        }

        tracing::debug!("LinearSolverChain(warm): BiCGSTAB (last resort) converged");
        Ok(x)
    }
}
