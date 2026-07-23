//! # Finite Element Method (FEM) Solver for 3D Incompressible Flow
//!
//! This module implements a high-performance FEM solver for the incompressible
//! Navier-Stokes equations using mixed Taylor-Hood (P2-P1) formulation.
//!
//! # Theorem — Taylor–Hood Inf-Sup Stability (Brezzi 1974, Taylor & Hood 1973)
//!
//! The P2-P1 mixed element pair (quadratic velocity, linear pressure) satisfies
//! the Ladyzhenskaya–Babuška–Brezzi (LBB) inf-sup condition:
//!
//! ```text
//! inf_{q ∈ Q_h} sup_{v ∈ V_h} [(∇·v, q) / (‖v‖_1 ‖q‖_0)] ≥ β > 0
//! ```
//!
//! with $\beta$ independent of $h$. This guarantees unique solvability of the
//! discrete saddle-point system and optimal-order convergence:
//!
//! ```text
//! ‖u − u_h‖_1 + ‖p − p_h‖_0 ≤ C h² (‖u‖_3 + ‖p‖_2)
//! ```
//!
//! **Reference:** Brezzi, F., "On the existence, uniqueness and discretization of
//! saddle-point problems arising from Lagrangian multipliers", RAIRO Anal. Numér. 8(R-2), 1974.

use cfd_core::error::Result;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::{LinearSolverChain, GMRES};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use eunomia::{FloatElement, NumericElement};
use leto::{Array1, Vector3};
use tracing;

use crate::fem::leto_bridge::build_with_vector_rhs;
use crate::fem::mid_node_cache::MidNodeCache;
use crate::fem::quadrature::TetrahedronQuadrature;
use crate::fem::shape_functions::LagrangeTet10;
use crate::fem::{scalar, FemConfig, StokesFlowProblem, StokesFlowSolution};
use crate::linalg::{
    array1_copy, array1_l2_norm, array1_len, array1_subarray, matrix3_determinant,
    matrix3_from_columns, matrix3_try_inverse, reference_tet_gradients, vector3_from_indexed,
};
use crate::scalar::Cfd3dScalar;
use moirai::{fold_reduce_with, Adaptive};
use std::collections::HashMap;

// Re-export mesh utility functions that were previously defined here.
pub(crate) use super::mesh_utils::compute_mesh_scale;
pub use super::mesh_utils::{extract_vertex_indices, extract_vertex_indices_cached};

/// Finite Element Method solver for 3D incompressible flow
pub struct FemSolver<T: Cfd3dScalar> {
    config: FemConfig<T>,
    _linear_solver: GMRES<T>,
    /// Reusable matrix builder to avoid O(N) allocations per iteration
    matrix_builder: Option<SparseMatrixBuilder<T>>,
    /// Reusable RHS vector
    rhs: Option<Array1<T>>,
    /// Reusable edge→mid-node lookup cache (GAP-PERF-001).
    ///
    /// Built once per mesh (keyed by `n_corner_nodes` + vertex count) and reused
    /// across Picard iterations. The mid-edge scan inside `extract_vertex_indices`
    /// is O(n_mid) per cell; the cached variant `extract_vertex_indices_cached`
    /// shrinks that to O(1) amortised. Without this hoist the per-iteration cost
    /// is O(n_elements × n_mid) which dominates large FEM assemblies.
    mid_node_cache: Option<MidNodeCache>,
    /// Cache key for `mid_node_cache`: `(n_corner_nodes, vertex_count)`.
    mid_cache_key: Option<(usize, usize)>,
    /// Reusable vertex-position array (one entry per mesh node, in mesh order).
    /// Avoids rebuilding the `Vec<Vector3<T>>` per Picard iteration since the
    /// mesh is immutable across Picard iterations of the same problem.
    vertex_positions: Option<Vec<Vector3<T>>>,
}

impl<T: Cfd3dScalar> FemSolver<T> {
    /// Create a new FEM solver with the given configuration
    pub fn new(config: FemConfig<T>) -> Self {
        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: config.base.convergence.max_iterations,
            tolerance: config.base.convergence.tolerance,
            ..cfd_math::linear_solver::IterativeSolverConfig::default()
        };
        let linear_solver = GMRES::new(solver_config, 100);
        Self {
            _linear_solver: linear_solver,
            config,
            matrix_builder: None,
            rhs: None,
            mid_node_cache: None,
            mid_cache_key: None,
            vertex_positions: None,
        }
    }

    /// Stokes Flow Problem ($-\mu \nabla^2 \mathbf{u} + \nabla p = \mathbf{f}$, $\nabla \cdot \mathbf{u} = 0$)
    ///
    /// ### Mathematical Invariants
    /// 1. **LBB (Babuška-Brezzi) Stability**:
    ///    The velocity-pressure space pair $(\mathcal{V}_h, \mathcal{Q}_h)$ must satisfy
    ///    $\inf_{q \in \mathcal{Q}_h} \sup_{v \in \mathcal{V}_h} \frac{\int q \nabla \cdot v}{\|v\|_{\mathcal{V}}\|q\|_{\mathcal{Q}}} \geq \beta > 0$.
    ///    Ensured here by Taylor-Hood $P_k / P_{k-1}$ elements (Quad/Linear).
    /// 2. **Mass Conservation**:
    ///    $\oint_{\partial \Omega} \mathbf{u} \cdot \mathbf{n} dA = \int_{\Omega} \nabla \cdot \mathbf{u} d\Omega = 0$.
    ///    Monitored via `continuity_residual` during assembly.
    /// 3. **Force Balance**:
    ///    Local momentum residual must vanish in the sense of distributions: $\langle \mathcal{R}, \phi \rangle = 0$.
    #[allow(clippy::too_many_lines)]
    pub fn solve(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<StokesFlowSolution<T>> {
        tracing::info!("Starting Taylor-Hood Stokes solver");

        // Strictly enforce mathematical well-posedness (Lax-Milgram requirement)
        // If this fails, the geometric mesh builder missed assigning boundary conditions
        // to some topological boundary nodes (e.g., P2 mid-edge nodes), making the BVP ill-posed.
        problem.validate()?;

        let (matrix, rhs) = self.assemble_system(problem, previous_solution)?;

        let n_total_dof = array1_len(&rhs);
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;

        // ── Linear solve: tiered fallback chain (Direct LU → GMRES/block → BiCGSTAB) ──
        //
        // # Algorithm — Tiered Linear Solver Fallback (see cfd_math::linear_solver::chain)
        //
        // The LinearSolverChain encapsulates the 5-tier fallback strategy that
        // was previously duplicated across fem/solver.rs and fem/projection_solver.rs.
        // Tier order: Direct LU → GMRES+BlockDiag → GMRES (unprec) → GMRES+ILU → BiCGSTAB.
        let rel_tol = self.config.base.convergence.tolerance;
        let abs_tol =
            (rel_tol * array1_l2_norm(&rhs)).max_scalar(<T as FloatElement>::from_f64(1e-14));
        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: self.config.base.convergence.max_iterations,
            tolerance: abs_tol,
            ..cfd_math::linear_solver::IterativeSolverConfig::default()
        };
        tracing::info!(
            n_total_dof,
            n_velocity_dof,
            n_corner_nodes,
            "FemSolver: starting linear solve"
        );
        let chain = LinearSolverChain::new(solver_config)
            // Direct LU threshold lowered from 100_000 to 512: the upstream
            // `leto_ops::SparseLuSolver` is a misnamed dense partial-pivoting LU
            // (see crates/cfd-math/src/linear_solver/direct_solver.rs:3-7), which
            // is O(n^3) and dominates per-Picard-iteration cost for any meaningfully
            // sized FEM system (1700-DOF saddle-point: ~3s/iter). For n <= 512 the
            // dense cost is <~0.05s, so direct LU stays the right call; above that,
            // GMRES+AMG / GMRES+BlockDiag (Tier 2/3) is faster for the saddle-point
            // structure. The strategic fix is a real sparse LU upstream in leto-ops
            // (arch board item); the threshold lowers is the tactical routing.
            .with_direct_threshold(512)
            .with_krylov_restart(std::cmp::min(200, n_total_dof.max(1)));
        let x_array = chain.solve(&matrix, &rhs, n_velocity_dof)?;

        // Guard: detect NaN/Inf from ill-conditioned or singular systems before
        // they silently propagate into the velocity/pressure solution fields.
        if x_array.iter().any(|v| !NumericElement::is_finite(*v)) {
            return Err(cfd_core::error::Error::Solver(
                "FEM linear solve produced non-finite values (NaN or Inf)".into(),
            ));
        }

        let velocity = array1_subarray(&x_array, 0, n_velocity_dof);
        let pressure = array1_subarray(&x_array, n_velocity_dof, n_corner_nodes);
        let solution =
            StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes);
        self.print_continuity_residual_stats(problem, &solution)?;

        Ok(solution)
    }

    /// Picard-optimized solve with warm-start and adaptive linear tolerance.
    ///
    /// # Optimizations over [`Self::solve`]
    ///
    /// 1. **Warm-start**: Uses the previous Picard solution as GMRES initial
    ///    guess, dramatically reducing iteration counts on subsequent iterations.
    /// 2. **Adaptive tolerance**: Early Picard iterations use 100× looser
    ///    linear solve tolerance (inexact Picard, cf. Eisenstat–Walker 1996),
    ///    since the viscosity field is still changing significantly.
    /// 3. **Reduced iteration budget**: 10K max linear iterations (vs 50K in
    ///    [`Self::solve`]) — with warm-starting, convergence should occur in
    ///    hundreds of iterations, not tens of thousands.
    /// 4. **Smart tier progression**: Skips unpreconditioned GMRES tier for
    ///    saddle-point systems when block-preconditioned GMRES stagnates.
    /// 5. **Timing diagnostics**: Logs assembly and linear solve wall-clock
    ///    time for performance monitoring.
    pub fn solve_picard(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
        picard_iteration: usize,
        max_picard_iterations: usize,
    ) -> Result<StokesFlowSolution<T>> {
        tracing::info!(
            picard_iteration,
            "Starting Taylor-Hood Stokes solver (Picard mode)"
        );

        problem.validate()?;

        let assembly_start = std::time::Instant::now();
        let (matrix, rhs) = self.assemble_system(problem, previous_solution)?;
        let assembly_elapsed = assembly_start.elapsed();

        let n_total_dof = array1_len(&rhs);
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;

        // Adaptive tolerance (inexact Picard / Eisenstat–Walker strategy):
        // Early iterations use 100× looser tolerance since viscosity is still changing.
        let base_tol = NumericElement::to_f64(self.config.base.convergence.tolerance);
        let adaptive_factor = if picard_iteration < max_picard_iterations / 2 {
            100.0_f64
        } else {
            1.0_f64
        };
        let rel_tol = <T as FloatElement>::from_f64(base_tol * adaptive_factor);
        let abs_tol =
            (rel_tol * array1_l2_norm(&rhs)).max_scalar(<T as FloatElement>::from_f64(1e-14));

        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: self.config.base.convergence.max_iterations,
            tolerance: abs_tol,
            ..cfd_math::linear_solver::IterativeSolverConfig::default()
        };

        let chain = LinearSolverChain::new(solver_config)
            // See `solve` above for the threshold rationale: dense-LU-as-sparse-LU
            // defeats O(n^3) cost above ~512 DOF, so route saddle-point Picard
            // iters to GMRES+AMG / GMRES+BlockDiag (Tier 2/3) instead.
            .with_direct_threshold(512)
            .with_krylov_restart(std::cmp::min(200, n_total_dof.max(1)));

        // Warm-start: reconstruct DOF vector from previous Picard solution.
        let initial_guess = previous_solution.map(|prev| {
            let mut x0 = Array1::zeros([n_total_dof]);
            let vel_len = n_velocity_dof.min(prev.velocity.len());
            let pres_len = n_corner_nodes.min(prev.pressure.len());
            let x0_slice = x0
                .as_slice_mut()
                .expect("invariant: warm-start vectors are dense one-dimensional Leto arrays");
            x0_slice[..vel_len].copy_from_slice(&prev.velocity.as_slice()[..vel_len]);
            x0_slice[n_velocity_dof..n_velocity_dof + pres_len]
                .copy_from_slice(&prev.pressure.as_slice()[..pres_len]);
            x0
        });
        let initial_guess_array = initial_guess.as_ref();

        tracing::info!(
            n_total_dof,
            n_velocity_dof,
            picard_iteration,
            ?abs_tol,
            max_iter = self.config.base.convergence.max_iterations,
            assembly_secs = assembly_elapsed.as_secs_f64(),
            "FemSolver: starting linear solve (Picard mode)"
        );

        tracing::info!(
            assembly_secs = format!("{:.2}", assembly_elapsed.as_secs_f64()).as_str(),
            ?abs_tol,
            max_iter = self.config.base.convergence.max_iterations,
            restart = std::cmp::min(200, n_total_dof.max(1)),
            "Assembly and linear solve config"
        );

        let solve_start = std::time::Instant::now();
        let x_array = chain.solve_with_guess(&matrix, &rhs, n_velocity_dof, initial_guess_array)?;
        let solve_elapsed = solve_start.elapsed();

        tracing::info!(
            solve_secs = format!("{:.2}", solve_elapsed.as_secs_f64()).as_str(),
            n_total_dof,
            "Linear solve complete"
        );

        // Guard: detect NaN/Inf from ill-conditioned or singular systems before
        // they silently propagate into the velocity/pressure solution fields.
        if x_array.iter().any(|v| !NumericElement::is_finite(*v)) {
            return Err(cfd_core::error::Error::Solver(
                "FEM linear solve (Picard) produced non-finite values (NaN or Inf)".into(),
            ));
        }

        let velocity = array1_subarray(&x_array, 0, n_velocity_dof);
        let pressure = array1_subarray(&x_array, n_velocity_dof, n_corner_nodes);
        let solution =
            StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes);
        self.print_continuity_residual_stats(problem, &solution)?;

        Ok(solution)
    }

    /// Compute and log continuity residual (∇·u) statistics.
    ///
    /// # Parallelization (GAP-PERF-008)
    ///
    /// Cells are processed in parallel via Moirai `fold_reduce_with`, accumulating
    /// per-worker partial accumulators `(residual_vec, max_abs, sum, l2, net)` that
    /// are merged hierarchically. This matches the pattern used in `assemble_system`.
    fn print_continuity_residual_stats(
        &self,
        problem: &StokesFlowProblem<T>,
        solution: &StokesFlowSolution<T>,
    ) -> Result<()> {
        let n_corner = problem.n_corner_nodes;

        // Per-cell contribution is purely read-only on mesh and solution;
        // accumulate per-thread local residual vectors, then merge.
        //
        // If `assemble_system` already ran on this `FemSolver` for the same mesh,
        // the `MidNodeCache` is cached on `self` and we use the cache-accelerated
        // index extraction; otherwise we fall back to the uncached variant to
        // preserve independent callability and identical element-wise output.
        let mid_cache = self.mid_node_cache.as_ref();

        let (residual, max_abs, sum_abs, l2, net) = fold_reduce_with::<Adaptive, _, _, _, _>(
            problem.mesh.cells.len(),
            || {
                (
                    vec![scalar::zero::<T>(); n_corner],
                    scalar::zero::<T>(),
                    scalar::zero::<T>(),
                    scalar::zero::<T>(),
                    scalar::zero::<T>(),
                )
            },
            |(mut res, mut mx, mut sm, mut l2acc, mut nt), cell_idx| {
                let cell = &problem.mesh.cells[cell_idx];
                let idxs = match mid_cache {
                    Some(cache) => {
                        extract_vertex_indices_cached(cell, &problem.mesh, n_corner, cache)
                    }
                    None => extract_vertex_indices(cell, &problem.mesh, n_corner),
                };
                let idxs = match idxs {
                    Ok(v) => v,
                    Err(_) => return (res, mx, sm, l2acc, nt),
                };
                if idxs.len() < 4 {
                    return (res, mx, sm, l2acc, nt);
                }

                let verts: Vec<Vector3<T>> = idxs
                    .iter()
                    .map(|&idx| {
                        vector3_from_indexed(
                            &problem
                                .mesh
                                .vertices
                                .position(cfd_mesh::domain::core::index::VertexId::from_usize(idx))
                                .coords,
                        )
                    })
                    .collect();

                let j_mat = matrix3_from_columns(
                    verts[1] - verts[0],
                    verts[2] - verts[0],
                    verts[3] - verts[0],
                );
                let det_j = matrix3_determinant(&j_mat);
                let abs_det = NumericElement::abs(det_j);
                if abs_det <= <T as FloatElement>::from_f64(1e-24) {
                    return (res, mx, sm, l2acc, nt);
                }
                let j_inv_t = match matrix3_try_inverse(&j_mat) {
                    Some(ji) => ji.transpose(),
                    None => return (res, mx, sm, l2acc, nt),
                };

                let grad_ref_p1 = reference_tet_gradients();
                let p1_gradients_phys = j_inv_t * grad_ref_p1;
                let shape = LagrangeTet10::new(p1_gradients_phys);

                let quad = TetrahedronQuadrature::keast_degree_3();
                for (qp, &qw) in quad.points().iter().zip(quad.weights().iter()) {
                    let weight = qw * abs_det;
                    let l = [scalar::one::<T>() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
                    let grad_p2 = shape.gradients(&l);

                    let mut div_u = scalar::zero::<T>();
                    for i in 0..idxs.len().min(10) {
                        let vel = solution.get_velocity(idxs[i]);
                        let grad_i = if idxs.len() == 4 {
                            Vector3::new(
                                p1_gradients_phys[(0, i)],
                                p1_gradients_phys[(1, i)],
                                p1_gradients_phys[(2, i)],
                            )
                        } else {
                            Vector3::new(grad_p2[[0, i]], grad_p2[[1, i]], grad_p2[[2, i]])
                        };
                        div_u += grad_i.x * vel.x + grad_i.y * vel.y + grad_i.z * vel.z;
                    }

                    for j in 0..4 {
                        let p_idx = idxs[j];
                        if p_idx < n_corner {
                            res[p_idx] += l[j] * div_u * weight;
                        }
                    }
                }

                // Update running stats (local thread, no synchronization needed)
                for &r in &res {
                    let a = NumericElement::abs(r);
                    if a > mx {
                        mx = a;
                    }
                    sm += a;
                    l2acc += r * r;
                    nt += r;
                }

                (res, mx, sm, l2acc, nt)
            },
            |(mut r1, mx1, sm1, l2a1, nt1), (r2, mx2, sm2, l2a2, nt2)| {
                for i in 0..n_corner {
                    r1[i] += r2[i];
                }
                (
                    r1,
                    if mx2 > mx1 { mx2 } else { mx1 },
                    sm1 + sm2,
                    l2a1 + l2a2,
                    nt1 + nt2,
                )
            },
        );

        let n = residual.len();
        if n > 0 {
            let n_as_u64 = u64::try_from(n).expect("residual slice length is representable as u64");
            let mean_abs =
                sum_abs / <T as FloatElement>::from_f64(NumericElement::to_f64(n_as_u64));
            let l2_norm = NumericElement::sqrt(l2);
            tracing::debug!(
                max_abs = ?max_abs, mean_abs = ?mean_abs,
                l2 = ?l2_norm, net = ?net, n,
                "Continuity Residual (Bu)"
            );
        }

        Ok(())
    }

    /// Assemble the global saddle-point system from element-level contributions.
    ///
    /// # Per-Element Viscosity (Generalised-Newtonian Extension)
    ///
    /// When `problem.element_viscosities` is `Some(vec)`, the assembly uses
    /// per-element viscosity `μ_e = vec[e]` instead of the global constant
    /// `problem.fluid.viscosity`. This enables generalised-Newtonian
    /// constitutive models (Carreau-Yasuda, Casson) via outer Picard
    /// iteration without modifying the assembly kernel.
    ///
    /// The mathematical equivalence: the global stiffness matrix becomes
    ///
    /// ```text
    /// K = Σ_e  μ_e · K̃_e
    /// ```
    ///
    /// where $\tilde{K}_e$ is the unit-viscosity element stiffness matrix.
    /// This is standard globalisation with element-varying material
    /// properties (Zienkiewicz & Taylor, 2000, §10.3).
    ///
    /// # Parallelization
    ///
    /// Element-level assembly is data-parallel via Moirai `fold_reduce_with`.
    /// Each worker accumulates into a local `HashMap<(row, col), T>` and
    /// `DVector<T>`. Worker-local maps are merged, then inserted into the
    /// global `SparseMatrixBuilder`.
    fn assemble_system(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<(SparseMatrix<T>, Array1<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;
        let n_total_dof = n_velocity_dof + n_corner_nodes;

        if self
            .matrix_builder
            .as_ref()
            .is_none_or(|b| b.num_rows() != n_total_dof)
        {
            self.matrix_builder = Some(SparseMatrixBuilder::new(n_total_dof, n_total_dof));
        } else {
            self.matrix_builder
                .as_mut()
                .expect("checked Some above")
                .clear();
        }

        if self
            .rhs
            .as_ref()
            .is_none_or(|r| array1_len(r) != n_total_dof)
        {
            self.rhs = Some(Array1::zeros([n_total_dof]));
        } else {
            self.rhs
                .as_mut()
                .expect("checked Some above")
                .fill(scalar::zero::<T>());
        }

        let mut matrix_builder = self
            .matrix_builder
            .take()
            .expect("matrix_builder initialized above");
        let mut rhs_store = self.rhs.take().expect("rhs initialized above");

        // GAP-PERF-001 (audited Session 12): `MidNodeCache::build` is O(n_mid × ~6)
        // and was being recomputed per Picard iteration only to be discarded — the
        // worker closure was calling the un-cached `extract_vertex_indices`, paying
        // O(n_mid) per cell. Hoist the cache to `FemSolver` and rebuild it only
        // when the mesh shape changes. The cached variant `extract_vertex_indices_cached`
        // is element-wise identical by the cache invariant (see mesh_utils docs).
        let mid_cache_key_now = (n_corner_nodes, n_nodes);
        if self
            .mid_node_cache
            .as_ref()
            .is_none_or(|_| self.mid_cache_key != Some(mid_cache_key_now))
        {
            self.mid_node_cache = Some(MidNodeCache::build(&problem.mesh, n_corner_nodes));
            self.mid_cache_key = Some(mid_cache_key_now);
        }
        let mid_cache = self
            .mid_node_cache
            .as_ref()
            .expect("mid_node_cache built above");

        // Hoist `vertex_positions` likewise: it depends only on the mesh, so an
        // immutable mesh across Picard iterations yields an immutable mapping.
        // Invalidate when the vertex count changes (the only free dimension of the
        // `vertex_positions[idx]` indexing used by the worker closure below).
        if self
            .vertex_positions
            .as_ref()
            .is_none_or(|v| v.len() != n_nodes)
        {
            self.vertex_positions = Some(
                problem
                    .mesh
                    .vertices
                    .iter()
                    .map(|v| vector3_from_indexed(&v.1.position.coords))
                    .collect(),
            );
        }
        let vertex_positions = self
            .vertex_positions
            .as_ref()
            .expect("vertex_positions built above");

        tracing::info!(
            n_elements = problem.mesh.cells.len(),
            "Assembling elements in parallel"
        );

        let (entry_map, mut rhs) = fold_reduce_with::<Adaptive, _, _, _, _>(
            problem.mesh.cells.len(),
            || (HashMap::with_capacity(512), Array1::zeros([n_total_dof])),
            |(mut local_map, local_rhs), i| {
                let cell = &problem.mesh.cells[i];
                let viscosity = problem
                    .element_viscosities
                    .as_ref()
                    .map_or(problem.fluid.viscosity, |v| v[i]);
                // Cache-accelerated index extraction: O(1) per edge vs O(N_mid)
                let idxs = match extract_vertex_indices_cached(
                    cell,
                    &problem.mesh,
                    problem.n_corner_nodes,
                    mid_cache,
                ) {
                    Ok(v) => v,
                    Err(_) => return (local_map, local_rhs), // skip degenerate cell
                };
                let local_verts: Vec<Vector3<T>> =
                    idxs.iter().map(|&idx| vertex_positions[idx]).collect();

                if local_verts.len() < 4_usize {
                    // Only tetrahedral cells are supported by this assembler.
                    // Mesh extraction can occasionally surface malformed cells
                    // with insufficient recovered corner vertices; they carry no
                    // valid tetra contribution and must be skipped.
                    return (local_map, local_rhs);
                }

                let v0 = local_verts[0];
                let v1 = local_verts[1];
                let v2 = local_verts[2];
                let v3 = local_verts[3];

                let six = <T as FloatElement>::from_f64(6.0);
                let vol = ((v1 - v0).cross(v2 - v0)).dot(v3 - v0) / six;
                let vol_tol = <T as FloatElement>::from_f64(1e-22);
                if NumericElement::abs(vol) < vol_tol {
                    // Skip degenerate element — zero volume means zero
                    // contribution to the global stiffness matrix
                    return (local_map, local_rhs);
                }

                let u_avg = self.calculate_u_avg(&idxs, previous_solution);

                self.assemble_element_local(
                    &mut local_map,
                    &idxs,
                    &local_verts,
                    viscosity,
                    problem.fluid.density,
                    u_avg,
                    n_nodes,
                );

                (local_map, local_rhs)
            },
            |(mut map1, mut rhs1), (map2, rhs2)| {
                for (k, v) in map2 {
                    *map1.entry(k).or_insert(scalar::zero::<T>()) += v;
                }
                for index in 0..array1_len(&rhs1) {
                    rhs1[index] += rhs2[index];
                }
                (map1, rhs1)
            },
        );

        tracing::debug!("Assembly map-reduce complete. Applying boundary conditions");
        // Populate builder with accumulated map entries
        for ((row, col), val) in entry_map {
            matrix_builder.add_entry(row, col, val)?;
        }

        self.apply_boundary_conditions_block(&mut matrix_builder, &mut rhs, problem, n_nodes)?;

        let velocity_dofs_constrained = problem.boundary_conditions.len() * 3;
        tracing::debug!(
            constrained = velocity_dofs_constrained,
            total = n_velocity_dof,
            "Velocity DOFs constrained"
        );

        if velocity_dofs_constrained == n_velocity_dof {
            tracing::warn!("All velocity DOFs constrained — may cause incompressibility conflict");
        }

        let diag_eps = problem.fluid.viscosity * <T as FloatElement>::from_f64(1e-12);
        for i in n_velocity_dof..n_total_dof {
            let _ = matrix_builder.add_entry(i, i, diag_eps);
        }

        let (matrix, rhs_after_assembly) =
            build_with_vector_rhs(matrix_builder, rhs, "FEM saddle-point RHS")?;
        rhs = rhs_after_assembly;
        array1_copy(&rhs, &mut rhs_store);
        self.rhs = Some(rhs_store);
        Ok((matrix, rhs))
    }

    fn assemble_element_local(
        &self,
        local_map: &mut HashMap<(usize, usize), T>,
        idxs: &[usize],
        verts: &[Vector3<T>],
        viscosity: T,
        density: T,
        u_avg: Vector3<T>,
        n_nodes: usize,
    ) {
        let quad = TetrahedronQuadrature::keast_degree_3();
        let points = quad.points();
        let weights = quad.weights();
        let grad_div_penalty = self.config.grad_div_penalty;

        let j_mat = matrix3_from_columns(
            verts[1] - verts[0],
            verts[2] - verts[0],
            verts[3] - verts[0],
        );
        let det_j = matrix3_determinant(&j_mat);
        let abs_det = NumericElement::abs(det_j);
        let j_inv_t = match matrix3_try_inverse(&j_mat) {
            Some(inv) => inv.transpose(),
            None => return, // Degenerate element — skip assembly
        };

        let grad_ref_p1 = reference_tet_gradients();
        let p1_gradients_phys = j_inv_t * grad_ref_p1;
        let shape = LagrangeTet10::new(p1_gradients_phys);

        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;

        for (qp, &qw) in points.iter().zip(weights.iter()) {
            let weight = qw * abs_det;
            let l = [scalar::one::<T>() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
            let n_p2 = shape.values(&l);
            let grad_p2_mat = shape.gradients(&l);
            let n_p1 = l;

            for i in 0..idxs.len().min(10) {
                let gi = idxs[i];
                let grad_i = if idxs.len() == 4 {
                    Vector3::new(
                        p1_gradients_phys[(0, i)],
                        p1_gradients_phys[(1, i)],
                        p1_gradients_phys[(2, i)],
                    )
                } else {
                    Vector3::new(
                        grad_p2_mat[[0, i]],
                        grad_p2_mat[[1, i]],
                        grad_p2_mat[[2, i]],
                    )
                };
                let n_i = if idxs.len() == 4 { n_p1[i] } else { n_p2[i] };

                for d in 0..3 {
                    let gv_i = gi + d * v_offset;
                    for j in 0..idxs.len().min(10) {
                        let gj = idxs[j];
                        let gv_j = gj + d * v_offset;
                        let grad_j = if idxs.len() == 4 {
                            Vector3::new(
                                p1_gradients_phys[(0, j)],
                                p1_gradients_phys[(1, j)],
                                p1_gradients_phys[(2, j)],
                            )
                        } else {
                            Vector3::new(
                                grad_p2_mat[[0, j]],
                                grad_p2_mat[[1, j]],
                                grad_p2_mat[[2, j]],
                            )
                        };
                        let _n_j = if idxs.len() == 4 { n_p1[j] } else { n_p2[j] };

                        // GAP-PERF-002: Fuse viscous + advection into a single HashMap probe.
                        // Before: two separate `entry().or_insert() += visc` and `+= adv`
                        // After: one probe computes visc+adv sum, stored in one HashMap op.
                        //
                        // # Proof of numerical equivalence:
                        // visc + adv = μ(∇Nᵢ · ∇Nⱼ)w + ρ Nᵢ (u_avg · ∇Nⱼ)w
                        // Accumulated in same float accumulator — identical to two separate adds.
                        let visc = viscosity * grad_i.dot(grad_j) * weight;
                        let adv = density * n_i * u_avg.dot(grad_j) * weight;
                        *local_map.entry((gv_i, gv_j)).or_insert(scalar::zero::<T>()) += visc + adv;
                        if grad_div_penalty > scalar::zero::<T>() {
                            for e in 0..3 {
                                let gv_j_e = gj + e * v_offset;
                                let grad_div = grad_div_penalty * grad_i[d] * grad_j[e] * weight;
                                *local_map
                                    .entry((gv_i, gv_j_e))
                                    .or_insert(scalar::zero::<T>()) += grad_div;
                            }
                        }
                    }
                    for j in 0..4 {
                        let gj = idxs[j];
                        let gp_j = p_offset + gj;
                        let b_val = n_p1[j] * grad_i[d] * weight;
                        *local_map.entry((gv_i, gp_j)).or_insert(scalar::zero::<T>()) -= b_val;
                        *local_map.entry((gp_j, gv_i)).or_insert(scalar::zero::<T>()) += b_val;
                    }
                }
            }

            // PSPG pressure stabilization (Brezzi-Pitkäranta)
            // Adds τ_BP * ∫ ∇q_i · ∇q_j dΩ to the pressure-pressure block.
            // τ_BP = h_e² / (12 * μ), where h_e = (6V)^(1/3).
            // Circumvents the LBB inf-sup condition for equal-order P1-P1 elements.
            if viscosity > scalar::zero::<T>() {
                let one_third = <T as FloatElement>::from_f64(1.0 / 3.0);
                let h_e = FloatElement::powf(abs_det, one_third); // (6V)^(1/3) ≈ element diameter
                let twelve = <T as FloatElement>::from_f64(12.0);
                let tau_bp = h_e * h_e / (twelve * viscosity);
                let vol_e = abs_det / <T as FloatElement>::from_f64(6.0);

                for i in 0..4 {
                    let gi = idxs[i];
                    let gp_i = p_offset + gi;
                    let grad_p_i = Vector3::new(
                        p1_gradients_phys[(0, i)],
                        p1_gradients_phys[(1, i)],
                        p1_gradients_phys[(2, i)],
                    );
                    for j in 0..4 {
                        let gj = idxs[j];
                        let gp_j = p_offset + gj;
                        let grad_p_j = Vector3::new(
                            p1_gradients_phys[(0, j)],
                            p1_gradients_phys[(1, j)],
                            p1_gradients_phys[(2, j)],
                        );
                        let pspg = tau_bp * grad_p_i.dot(grad_p_j) * vol_e;
                        *local_map.entry((gp_i, gp_j)).or_insert(scalar::zero::<T>()) += pspg;
                    }
                }
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    fn apply_boundary_conditions_block(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut Array1<T>,
        problem: &StokesFlowProblem<T>,
        n_nodes: usize,
    ) -> Result<()> {
        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;
        let mesh_scale = compute_mesh_scale(&problem.mesh);
        let diag_scale = problem.fluid.viscosity * mesh_scale;

        let mut vel_dofs = std::collections::HashSet::new();
        let mut p_dofs = std::collections::HashSet::new();
        let mut inlet_nodes = 0usize;
        let mut wall_nodes = 0usize;
        let mut outlet_nodes = 0usize;
        let mut dirichlet_nodes = 0usize;
        let mut unconstrained_boundary_nodes = Vec::new();

        let boundary_nodes = problem.get_boundary_nodes();
        for &node_idx in &boundary_nodes {
            if !problem.boundary_conditions.contains_key(&node_idx) {
                unconstrained_boundary_nodes.push(node_idx);
            }
        }

        if !unconstrained_boundary_nodes.is_empty() {
            tracing::warn!(
                count = unconstrained_boundary_nodes.len(),
                first_5 = ?&unconstrained_boundary_nodes[..unconstrained_boundary_nodes.len().min(5)],
                "Boundary Leak: boundary nodes have no BCs"
            );
        }

        // Track if any pressure boundary condition is applied
        let mut has_pressure_bc = false;

        // Extract and uniquely sort boundary conditions to guarantee exact deterministic linear matrix
        // assembly across identical geometries on multi-threaded parallel executors with randomized `HashMap`s.
        let mut sorted_bcs: Vec<_> = problem.boundary_conditions.iter().collect();
        sorted_bcs.sort_unstable_by_key(|(&k, _)| k);

        for (&node_idx, bc) in sorted_bcs {
            match bc {
                BoundaryCondition::VelocityInlet { velocity } => {
                    inlet_nodes += 1;
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        builder.set_dirichlet_row(dof, diag_scale, velocity[d]);
                        rhs[dof] = velocity[d] * diag_scale;
                        vel_dofs.insert(dof);
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    wall_nodes += 1;
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        builder.set_dirichlet_row(dof, diag_scale, scalar::zero::<T>());
                        rhs[dof] = scalar::zero::<T>();
                        vel_dofs.insert(dof);
                    }
                }
                BoundaryCondition::PressureOutlet { pressure }
                | BoundaryCondition::PressureInlet { pressure, .. } => {
                    outlet_nodes += 1;
                    if node_idx < problem.n_corner_nodes {
                        has_pressure_bc = true;
                        let dof = p_offset + node_idx;
                        builder.set_dirichlet_row(dof, diag_scale, *pressure);
                        rhs[dof] = *pressure * diag_scale;
                        p_dofs.insert(dof);
                    }
                }
                BoundaryCondition::Dirichlet {
                    value,
                    component_values,
                } => {
                    dirichlet_nodes += 1;
                    if let Some(comps) = component_values {
                        for d in 0..3 {
                            if let Some(Some(val)) = comps.get(d) {
                                let dof = node_idx + d * v_offset;
                                builder.set_dirichlet_row(dof, diag_scale, *val);
                                rhs[dof] = *val * diag_scale;
                                vel_dofs.insert(dof);
                            }
                        }
                        if let Some(Some(p_val)) = comps.get(3) {
                            if node_idx < problem.n_corner_nodes {
                                has_pressure_bc = true;
                                let dof = p_offset + node_idx;
                                builder.set_dirichlet_row(dof, diag_scale, *p_val);
                                rhs[dof] = *p_val * diag_scale;
                                p_dofs.insert(dof);
                            }
                        }
                    } else {
                        // Scalar Dirichlet: apply to all velocity components (standard wall/inlet)
                        // This usually isn't what's desired for pressure, but we follow the old logic.
                        for d in 0..3 {
                            let dof = node_idx + d * v_offset;
                            builder.set_dirichlet_row(dof, diag_scale, *value);
                            rhs[dof] = *value * diag_scale;
                            vel_dofs.insert(dof);
                        }
                    }
                }
                _ => {}
            }
        }

        // CRITICAL FIX: For incompressible flow without pressure BC, pressure is
        // undetermined up to a constant. Pin first corner node pressure to zero
        // as reference to remove this null space and make the system non-singular.
        if !has_pressure_bc && problem.n_corner_nodes > 0 {
            let reference_pressure_dof = p_offset; // First corner node
            builder.set_dirichlet_row(reference_pressure_dof, diag_scale, scalar::zero::<T>());
            rhs[reference_pressure_dof] = scalar::zero::<T>();
            tracing::info!(
                reference_pressure_dof,
                "Pinned pressure DOF to zero (no pressure BC specified)"
            );
        }

        tracing::debug!(
            inlet_nodes,
            wall_nodes,
            outlet_nodes,
            dirichlet_nodes,
            "BC Diagnostics: node counts"
        );
        tracing::debug!(
            velocity_dofs_set = vel_dofs.len(),
            pressure_dofs_set = p_dofs.len(),
            n_velocity_dof = n_nodes * 3,
            n_pressure_dof = problem.n_corner_nodes,
            "BC Diagnostics: DOF counts"
        );

        Ok(())
    }

    fn calculate_u_avg(
        &self,
        nodes: &[usize],
        solution: Option<&StokesFlowSolution<T>>,
    ) -> Vector3<T> {
        if let Some(sol) = solution {
            if nodes.is_empty() {
                return Vector3::zeros();
            }
            let mut sum = Vector3::zeros();
            for &n in nodes {
                sum += sol.get_velocity(n);
            }
            let node_count =
                u64::try_from(nodes.len()).expect("node count is representable as u64");
            sum / <T as FloatElement>::from_f64(NumericElement::to_f64(node_count))
        } else {
            Vector3::zeros()
        }
    }
}
