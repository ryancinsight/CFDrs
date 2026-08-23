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
use cfd_math::linear_solver::{LinearSolverChain, LinearSolverState};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use eunomia::{FloatElement, NumericElement};
use leto::{Array1, Vector3};
use tracing;

use crate::fem::leto_bridge::build_with_vector_rhs;
use crate::fem::mid_node_cache::MidNodeCache;
use crate::fem::quadrature::TetrahedronQuadrature;
use crate::fem::shape_functions::LagrangeTet10;
use crate::fem::{FemConfig, StokesFlowProblem, StokesFlowSolution};
use crate::linalg::{
    array1_l2_norm, array1_len, array1_subarray, matrix3_determinant, matrix3_from_columns,
    matrix3_try_inverse, reference_tet_gradients, vector3_from_indexed, Matrix3x4,
};
use crate::scalar;
use cfd_core::CfdScalar;
use moirai::{fold_reduce_with, Adaptive};
use std::collections::HashMap;

// Re-export mesh utility functions that were previously defined here.
pub(crate) use super::mesh_utils::compute_mesh_scale;
pub use super::mesh_utils::{extract_vertex_indices, extract_vertex_indices_cached};

/// Finite Element Method solver for 3D incompressible flow
pub struct FemSolver<T: CfdScalar + cfd_mesh::domain::core::Scalar> {
    config: FemConfig<T>,
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
    /// Mesh geometry terms reused across Picard assemblies.
    element_geometry: Option<Vec<Option<ElementGeometry<T>>>>,
    /// Cache key for `element_geometry`: `(corner nodes, vertices, cells)`.
    element_geometry_key: Option<(usize, usize, usize)>,
    /// Immutable Keast degree-3 quadrature rule shared by all elements.
    quadrature: TetrahedronQuadrature<T>,
    /// AMG transfer state reused across Picard solves with the same CSR pattern.
    linear_solver_state: LinearSolverState<T>,
}

struct ElementGeometry<T: CfdScalar + cfd_mesh::domain::core::Scalar> {
    indices: Vec<usize>,
    abs_det: T,
    p1_gradients_phys: Matrix3x4<T>,
}

struct AssembledSystem<T> {
    matrix: SparseMatrix<T>,
    rhs: Array1<T>,
    constrained_dofs: Vec<(usize, T)>,
}

impl<T: CfdScalar + cfd_mesh::domain::core::Scalar> FemSolver<T> {
    /// Create a new FEM solver with the given configuration
    pub fn new(config: FemConfig<T>) -> Self {
        Self {
            config,
            matrix_builder: None,
            rhs: None,
            mid_node_cache: None,
            mid_cache_key: None,
            vertex_positions: None,
            element_geometry: None,
            element_geometry_key: None,
            quadrature: TetrahedronQuadrature::keast_degree_3(),
            linear_solver_state: LinearSolverState::default(),
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

        let AssembledSystem {
            matrix,
            rhs,
            constrained_dofs,
        } = self.assemble_system(problem, previous_solution)?;

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
            relative_tolerance: self.config.base.convergence.relative_tolerance,
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
        let x_array = self.solve_reduced_system(
            &chain,
            &matrix,
            &rhs,
            &constrained_dofs,
            n_velocity_dof,
            None,
        )?;

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
        if tracing::enabled!(tracing::Level::DEBUG) {
            self.print_continuity_residual_stats(problem, &solution)?;
        }

        Ok(solution)
    }

    /// Picard solve with warm-start and configured linear tolerance.
    ///
    /// # Optimizations over [`Self::solve`]
    ///
    /// 1. **Warm-start**: Uses the previous Picard solution as GMRES initial
    ///    guess, dramatically reducing iteration counts on subsequent iterations.
    /// 2. **Configured tolerance**: Every Picard linear system uses the
    ///    configured relative tolerance so the returned velocity satisfies the
    ///    continuity contract before viscosity is updated.
    /// 3. **Smart tier progression**: Skips unpreconditioned GMRES tier for
    ///    saddle-point systems when block-preconditioned GMRES stagnates.
    /// 4. **Timing diagnostics**: Logs assembly and linear solve wall-clock
    ///    time for performance monitoring.
    pub fn solve_picard(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
        picard_iteration: usize,
        _max_picard_iterations: usize,
    ) -> Result<StokesFlowSolution<T>> {
        tracing::info!(
            picard_iteration,
            "Starting Taylor-Hood Stokes solver (Picard mode)"
        );

        problem.validate()?;

        let assembly_start = std::time::Instant::now();
        let AssembledSystem {
            matrix,
            rhs,
            constrained_dofs,
        } = self.assemble_system(problem, previous_solution)?;
        let assembly_elapsed = assembly_start.elapsed();

        let n_total_dof = array1_len(&rhs);
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;

        // The linear solve must satisfy the configured relative tolerance at
        // every Picard step. An inexact early solve can preserve the previous
        // velocity trace while the viscosity field has already changed,
        // violating the inlet/outlet continuity contract.
        let base_tol = NumericElement::to_f64(self.config.base.convergence.relative_tolerance);
        let rel_tol = <T as FloatElement>::from_f64(base_tol);
        let abs_tol =
            (rel_tol * array1_l2_norm(&rhs)).max_scalar(<T as FloatElement>::from_f64(1e-14));

        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: self.config.base.convergence.max_iterations,
            tolerance: abs_tol,
            relative_tolerance: rel_tol,
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
            ?rel_tol,
            max_iter = self.config.base.convergence.max_iterations,
            restart = std::cmp::min(200, n_total_dof.max(1)),
            "Assembly and linear solve config"
        );

        let solve_start = std::time::Instant::now();
        let x_array = self.solve_reduced_system(
            &chain,
            &matrix,
            &rhs,
            &constrained_dofs,
            n_velocity_dof,
            initial_guess_array,
        )?;
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
        if tracing::enabled!(tracing::Level::DEBUG) {
            self.print_continuity_residual_stats(problem, &solution)?;
        }

        Ok(solution)
    }

    /// Solve only the unconstrained Schur system.
    ///
    /// Dirichlet enforcement has already eliminated constrained columns and
    /// replaced constrained rows with prescribed-value identities. Removing
    /// those identities before Krylov iteration is algebraically exact: the
    /// free rows see the same RHS, while constrained values are restored after
    /// the solve. The reduced ordering preserves the velocity/pressure block
    /// split required by the saddle-point preconditioners.
    fn solve_reduced_system(
        &mut self,
        chain: &LinearSolverChain<T>,
        matrix: &SparseMatrix<T>,
        rhs: &Array1<T>,
        constrained_dofs: &[(usize, T)],
        n_velocity_dof: usize,
        initial_guess: Option<&Array1<T>>,
    ) -> Result<Array1<T>> {
        let n_total_dof = rhs.shape()[0];
        let mut constrained = vec![false; n_total_dof];
        let mut prescribed = vec![scalar::zero::<T>(); n_total_dof];
        for &(dof, value) in constrained_dofs {
            if let (Some(is_constrained), Some(prescribed_value)) =
                (constrained.get_mut(dof), prescribed.get_mut(dof))
            {
                *is_constrained = true;
                *prescribed_value = value;
            }
        }

        let free_indices: Vec<usize> = (0..n_total_dof).filter(|&dof| !constrained[dof]).collect();
        let mut reduced_index = vec![None; n_total_dof];
        for (reduced, &full) in free_indices.iter().enumerate() {
            reduced_index[full] = Some(reduced);
        }
        let reduced_velocity_dof = free_indices
            .iter()
            .take_while(|&&dof| dof < n_velocity_dof)
            .count();
        let reduced_size = free_indices.len();
        let mut reduced_scales = vec![scalar::one::<T>(); reduced_size];
        for (reduced_row, &full_row) in free_indices.iter().enumerate() {
            let row = matrix.row(full_row);
            let diagonal = row
                .col_indices()
                .iter()
                .zip(row.values())
                .find_map(|(&column, &value)| (column == full_row).then_some(value));
            if let Some(diagonal) = diagonal {
                let magnitude = NumericElement::abs(diagonal);
                if magnitude > scalar::zero::<T>() {
                    reduced_scales[reduced_row] =
                        <T as NumericElement>::ONE / NumericElement::sqrt(magnitude);
                }
            }
        }
        let mut reduced_builder = SparseMatrixBuilder::new(reduced_size, reduced_size);
        let mut reduced_rhs = Array1::zeros([reduced_size]);

        for (reduced_row, &full_row) in free_indices.iter().enumerate() {
            let row_scale = reduced_scales[reduced_row];
            reduced_rhs[reduced_row] = rhs[full_row] * row_scale;
            let row = matrix.row(full_row);
            for (&full_col, &value) in row.col_indices().iter().zip(row.values()) {
                if let Some(reduced_col) = reduced_index.get(full_col).and_then(|index| *index) {
                    reduced_builder.add_entry(
                        reduced_row,
                        reduced_col,
                        row_scale * value * reduced_scales[reduced_col],
                    )?;
                }
            }
        }
        let (reduced_matrix, reduced_rhs) =
            build_with_vector_rhs(reduced_builder, reduced_rhs, "reduced FEM saddle-point RHS")?;

        let reduced_guess = initial_guess.map(|guess| {
            let mut reduced = Array1::zeros([reduced_size]);
            for (reduced_index, &full_index) in free_indices.iter().enumerate() {
                // The Krylov system is D A D y = D b and the reconstructed
                // solution is x = D y. Convert the previous full-space x into
                // the scaled-space y before using it as a warm start.
                reduced[reduced_index] = guess[full_index] / reduced_scales[reduced_index];
            }
            reduced
        });
        let reduced_solution = chain.solve_with_guess_state(
            &reduced_matrix,
            &reduced_rhs,
            reduced_velocity_dof,
            reduced_guess.as_ref(),
            &mut self.linear_solver_state,
        )?;

        let mut solution = Array1::zeros([n_total_dof]);
        for (dof, value) in solution.iter_mut().zip(prescribed) {
            *dof = value;
        }
        for (reduced_index, &full_index) in free_indices.iter().enumerate() {
            solution[full_index] = reduced_scales[reduced_index] * reduced_solution[reduced_index];
        }
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

                for (qp, &qw) in self
                    .quadrature
                    .points()
                    .iter()
                    .zip(self.quadrature.weights().iter())
                {
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
    ) -> Result<AssembledSystem<T>> {
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
        let rhs_store = self.rhs.take().expect("rhs initialized above");

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

        let geometry_key = (n_corner_nodes, n_nodes, problem.mesh.cells.len());
        if self
            .element_geometry
            .as_ref()
            .is_none_or(|_| self.element_geometry_key != Some(geometry_key))
        {
            self.element_geometry = Some(
                problem
                    .mesh
                    .cells
                    .iter()
                    .map(|cell| {
                        let indices = extract_vertex_indices_cached(
                            cell,
                            &problem.mesh,
                            n_corner_nodes,
                            mid_cache,
                        )
                        .ok()?;
                        if indices.len() < 4 {
                            return None;
                        }

                        let vertices = indices
                            .iter()
                            .map(|&idx| vertex_positions[idx])
                            .collect::<Vec<_>>();
                        let v0 = vertices[0];
                        let v1 = vertices[1];
                        let v2 = vertices[2];
                        let v3 = vertices[3];
                        let six = <T as FloatElement>::from_f64(6.0);
                        let volume = ((v1 - v0).cross(v2 - v0)).dot(v3 - v0) / six;
                        let vol_tol = <T as FloatElement>::from_f64(1e-22);
                        if NumericElement::abs(volume) < vol_tol {
                            return None;
                        }

                        let jacobian = matrix3_from_columns(v1 - v0, v2 - v0, v3 - v0);
                        let abs_det = NumericElement::abs(matrix3_determinant(&jacobian));
                        let inverse_transpose = matrix3_try_inverse(&jacobian)?.transpose();
                        let p1_gradients_phys = inverse_transpose * reference_tet_gradients();

                        Some(ElementGeometry {
                            indices,
                            abs_det,
                            p1_gradients_phys,
                        })
                    })
                    .collect(),
            );
            self.element_geometry_key = Some(geometry_key);
        }
        let element_geometry = self
            .element_geometry
            .as_ref()
            .expect("element geometry cache built above");

        tracing::info!(
            n_elements = problem.mesh.cells.len(),
            "Assembling elements in parallel"
        );

        let entry_map = fold_reduce_with::<Adaptive, _, _, _, _>(
            problem.mesh.cells.len(),
            || HashMap::with_capacity(512),
            |mut local_map, i| {
                let viscosity = problem
                    .element_viscosities
                    .as_ref()
                    .map_or(problem.fluid.viscosity.into_base(), |v| v[i]);
                let Some(geometry) = element_geometry[i].as_ref() else {
                    return local_map;
                };

                let u_avg = self.calculate_u_avg(&geometry.indices, previous_solution);

                self.assemble_element_local(
                    &mut local_map,
                    geometry,
                    viscosity,
                    problem.fluid.density.into_base(),
                    u_avg,
                    n_nodes,
                );

                local_map
            },
            |mut map1, map2| {
                for (k, v) in map2 {
                    *map1.entry(k).or_insert(scalar::zero::<T>()) += v;
                }
                map1
            },
        );

        tracing::debug!("Assembly map-reduce complete. Applying boundary conditions");
        // Populate builder with accumulated map entries
        for ((row, col), val) in entry_map {
            matrix_builder.add_entry(row, col, val)?;
        }

        let mut rhs = rhs_store;
        let constrained_dofs =
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

        let diag_eps = problem.fluid.viscosity.into_base() * <T as FloatElement>::from_f64(1e-12);
        for i in n_velocity_dof..n_total_dof {
            let _ = matrix_builder.add_entry(i, i, diag_eps);
        }

        let (matrix, rhs_after_assembly) =
            build_with_vector_rhs(matrix_builder, rhs, "FEM saddle-point RHS")?;
        rhs = rhs_after_assembly;
        self.rhs = Some(rhs.clone());
        Ok(AssembledSystem {
            matrix,
            rhs,
            constrained_dofs,
        })
    }

    fn assemble_element_local(
        &self,
        local_map: &mut HashMap<(usize, usize), T>,
        geometry: &ElementGeometry<T>,
        viscosity: T,
        density: T,
        u_avg: Vector3<T>,
        n_nodes: usize,
    ) {
        let points = self.quadrature.points();
        let weights = self.quadrature.weights();
        let grad_div_penalty = self.config.grad_div_penalty;

        let idxs = geometry.indices.as_slice();
        let abs_det = geometry.abs_det;
        let p1_gradients_phys = geometry.p1_gradients_phys;
        if idxs.len() == 4 {
            self.assemble_element_linear(
                local_map,
                idxs,
                p1_gradients_phys,
                abs_det,
                viscosity,
                density,
                u_avg,
                n_nodes,
            );
            return;
        }

        let p2_shape = (idxs.len() != 4).then(|| LagrangeTet10::new(p1_gradients_phys));

        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;

        for (qp, &qw) in points.iter().zip(weights.iter()) {
            let weight = qw * abs_det;
            let l = [scalar::one::<T>() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
            let zero = scalar::zero::<T>();
            let mut shape_values = [zero; 10];
            let mut shape_gradients = [Vector3::new(zero, zero, zero); 10];
            if let Some(shape) = p2_shape.as_ref() {
                let n_p2 = shape.values(&l);
                let grad_p2_mat = shape.gradients(&l);
                for i in 0..10 {
                    shape_values[i] = n_p2[i];
                    shape_gradients[i] = Vector3::new(
                        grad_p2_mat[[0, i]],
                        grad_p2_mat[[1, i]],
                        grad_p2_mat[[2, i]],
                    );
                }
            } else {
                shape_values[..4].copy_from_slice(&l);
                for i in 0..4 {
                    shape_gradients[i] = Vector3::new(
                        p1_gradients_phys[(0, i)],
                        p1_gradients_phys[(1, i)],
                        p1_gradients_phys[(2, i)],
                    );
                }
            }

            for i in 0..idxs.len().min(10) {
                let gi = idxs[i];
                let grad_i = shape_gradients[i];
                let n_i = shape_values[i];

                for d in 0..3 {
                    let gv_i = gi + d * v_offset;
                    for j in 0..idxs.len().min(10) {
                        let gj = idxs[j];
                        let gv_j = gj + d * v_offset;
                        let grad_j = shape_gradients[j];

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
                        let b_val = shape_values[j] * grad_i[d] * weight;
                        *local_map.entry((gv_i, gp_j)).or_insert(scalar::zero::<T>()) -= b_val;
                        *local_map.entry((gp_j, gv_i)).or_insert(scalar::zero::<T>()) += b_val;
                    }
                }
            }

            // PSPG pressure stabilization (Brezzi-Pitkäranta)
            // Adds τ_BP * ∫ ∇q_i · ∇q_j dΩ to the pressure-pressure block.
            // τ_BP = h_e² / (12 * μ), where h_e = (6V)^(1/3).
            // Circumvents the LBB inf-sup condition for equal-order P1-P1 elements.
            // Taylor-Hood P2-P1 is inf-sup stable and does not use this term.
            if idxs.len() == 4 && viscosity > scalar::zero::<T>() {
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

    fn assemble_element_linear(
        &self,
        local_map: &mut HashMap<(usize, usize), T>,
        idxs: &[usize],
        p1_gradients_phys: Matrix3x4<T>,
        abs_det: T,
        viscosity: T,
        density: T,
        u_avg: Vector3<T>,
        n_nodes: usize,
    ) {
        let points = self.quadrature.points();
        let weights = self.quadrature.weights();
        let zero = scalar::zero::<T>();
        let mut weight_sum = zero;
        let mut basis_integrals = [zero; 4];
        for (qp, &qw) in points.iter().zip(weights.iter()) {
            weight_sum += qw;
            let l = [scalar::one::<T>() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
            for i in 0..4 {
                basis_integrals[i] += qw * l[i];
            }
        }

        let gradients = [
            Vector3::new(
                p1_gradients_phys[(0, 0)],
                p1_gradients_phys[(1, 0)],
                p1_gradients_phys[(2, 0)],
            ),
            Vector3::new(
                p1_gradients_phys[(0, 1)],
                p1_gradients_phys[(1, 1)],
                p1_gradients_phys[(2, 1)],
            ),
            Vector3::new(
                p1_gradients_phys[(0, 2)],
                p1_gradients_phys[(1, 2)],
                p1_gradients_phys[(2, 2)],
            ),
            Vector3::new(
                p1_gradients_phys[(0, 3)],
                p1_gradients_phys[(1, 3)],
                p1_gradients_phys[(2, 3)],
            ),
        ];
        let integrated_weight = abs_det * weight_sum;
        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;

        for i in 0..4 {
            let gi = idxs[i];
            let grad_i = gradients[i];
            for d in 0..3 {
                let gv_i = gi + d * v_offset;
                for j in 0..4 {
                    let gj = idxs[j];
                    let grad_j = gradients[j];
                    let gv_j = gj + d * v_offset;
                    let visc = viscosity * grad_i.dot(grad_j) * integrated_weight;
                    let adv = density * basis_integrals[i] * u_avg.dot(grad_j) * abs_det;
                    *local_map.entry((gv_i, gv_j)).or_insert(zero) += visc + adv;
                    if self.config.grad_div_penalty > zero {
                        for e in 0..3 {
                            let gv_j_e = gj + e * v_offset;
                            let grad_div = self.config.grad_div_penalty
                                * grad_i[d]
                                * grad_j[e]
                                * integrated_weight;
                            *local_map.entry((gv_i, gv_j_e)).or_insert(zero) += grad_div;
                        }
                    }
                }
                for j in 0..4 {
                    let gj = idxs[j];
                    let gp_j = p_offset + gj;
                    let b_val = basis_integrals[j] * grad_i[d] * abs_det;
                    *local_map.entry((gv_i, gp_j)).or_insert(zero) -= b_val;
                    *local_map.entry((gp_j, gv_i)).or_insert(zero) += b_val;
                }
            }
        }

        // Add the pressure-stabilization bilinear form once per element. The
        // element volume is already represented by
        // `vol_e`; repeating this contribution for every quadrature point
        // would scale the incompressibility defect with the quadrature rule
        // instead of the physical element.
        if viscosity > zero {
            let one_third = <T as FloatElement>::from_f64(1.0 / 3.0);
            let h_e = FloatElement::powf(abs_det, one_third);
            let twelve = <T as FloatElement>::from_f64(12.0);
            let tau_bp = h_e * h_e / (twelve * viscosity);
            let vol_e = abs_det / <T as FloatElement>::from_f64(6.0);
            for i in 0..4 {
                let gi = idxs[i];
                let gp_i = p_offset + gi;
                for j in 0..4 {
                    let gj = idxs[j];
                    let gp_j = p_offset + gj;
                    let pspg = tau_bp * gradients[i].dot(gradients[j]) * vol_e;
                    *local_map.entry((gp_i, gp_j)).or_insert(zero) += pspg;
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
    ) -> Result<Vec<(usize, T)>> {
        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;
        let mesh_scale = compute_mesh_scale(&problem.mesh);
        let diag_scale = problem.fluid.viscosity.into_base() * mesh_scale;

        let mut vel_dofs = std::collections::HashSet::new();
        let mut p_dofs = std::collections::HashSet::new();
        let mut has_pressure_bc = false;
        let mut inlet_nodes = 0usize;
        let mut wall_nodes = 0usize;
        let mut outlet_nodes = 0usize;
        let mut dirichlet_nodes = 0usize;
        let mut unconstrained_boundary_nodes = Vec::new();
        let mut constrained_dofs = Vec::new();

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
                        constrained_dofs.push((dof, velocity[d]));
                        vel_dofs.insert(dof);
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    wall_nodes += 1;
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        builder.set_dirichlet_row(dof, diag_scale, scalar::zero::<T>());
                        rhs[dof] = scalar::zero::<T>();
                        constrained_dofs.push((dof, scalar::zero::<T>()));
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
                        constrained_dofs.push((dof, *pressure));
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
                                constrained_dofs.push((dof, *val));
                                vel_dofs.insert(dof);
                            }
                        }
                        if let Some(Some(p_val)) = comps.get(3) {
                            if node_idx < problem.n_corner_nodes {
                                has_pressure_bc = true;
                                let dof = p_offset + node_idx;
                                builder.set_dirichlet_row(dof, diag_scale, *p_val);
                                rhs[dof] = *p_val * diag_scale;
                                constrained_dofs.push((dof, *p_val));
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
                            constrained_dofs.push((dof, *value));
                            vel_dofs.insert(dof);
                        }
                    }
                }
                _ => {}
            }
        }

        // In the natural-traction case pressure is defined only up to a
        // constant. Pin one pressure DOF for the gauge; the pressure
        // stabilization retains the remaining continuity equations.
        if !has_pressure_bc && problem.n_corner_nodes > 0 {
            let reference_pressure_dof = p_offset;
            builder.set_dirichlet_row(reference_pressure_dof, diag_scale, scalar::zero::<T>());
            rhs[reference_pressure_dof] = scalar::zero::<T>();
            constrained_dofs.push((reference_pressure_dof, scalar::zero::<T>()));
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

        Ok(constrained_dofs)
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
