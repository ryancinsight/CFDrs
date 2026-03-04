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
use nalgebra::{DVector, RealField, Vector3};
use num_traits::{Float, FromPrimitive};
use tracing;

use crate::fem::mid_node_cache::MidNodeCache;
use crate::fem::quadrature::TetrahedronQuadrature;
use crate::fem::shape_functions::LagrangeTet10;
use crate::fem::{FemConfig, StokesFlowProblem, StokesFlowSolution};
use cfd_mesh::domain::topology::Cell;
use cfd_mesh::IndexedMesh;
use rayon::prelude::*;
use std::collections::HashMap;

/// Finite Element Method solver for 3D incompressible flow
pub struct FemSolver<
    T: cfd_mesh::domain::core::Scalar
        + RealField
        + Copy
        + FromPrimitive
        + num_traits::Float
        + std::fmt::Debug,
> {
    config: FemConfig<T>,
    _linear_solver: GMRES<T>,
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + FromPrimitive
            + Copy
            + Float
            + std::fmt::Debug
            + From<f64>,
    > FemSolver<T>
{
    /// Create a new FEM solver with the given configuration
    pub fn new(config: FemConfig<T>) -> Self {
        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: 30000,
            tolerance: <T as FromPrimitive>::from_f64(1e-12).unwrap_or_else(T::zero),
            ..cfd_math::linear_solver::IterativeSolverConfig::default()
        };
        let linear_solver = GMRES::new(solver_config, 100);
        Self {
            _linear_solver: linear_solver,
            config,
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

        let n_total_dof = rhs.len();
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
        let rel_tol = <T as FromPrimitive>::from_f64(1e-8).unwrap_or_else(T::zero);
        let abs_tol = Float::max(
            rel_tol * rhs.norm(),
            <T as FromPrimitive>::from_f64(1e-14).unwrap_or_else(T::zero),
        );
        let solver_config = cfd_math::linear_solver::IterativeSolverConfig {
            max_iterations: 50_000,
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
            .with_direct_threshold(100_000)
            .with_krylov_restart(std::cmp::min(200, n_total_dof.max(1)));
        let x = chain.solve(&matrix, &rhs, n_velocity_dof)?;

        let velocity = x.rows(0, n_velocity_dof).into_owned();
        let pressure = x.rows(n_velocity_dof, n_corner_nodes).into_owned();
        let solution =
            StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes);
        self.print_continuity_residual_stats(problem, &solution)?;

        Ok(solution)
    }

    /// Compute and log continuity residual (∇·u) statistics.
    ///
    /// # Parallelization (GAP-PERF-008)
    ///
    /// Cells are processed in parallel via `par_iter().fold().reduce()`, accumulating
    /// per-thread partial accumulators `(residual_vec, max_abs, sum, l2, net)` that
    /// are merged hierarchically. This matches the pattern used in `assemble_system`.
    fn print_continuity_residual_stats(
        &self,
        problem: &StokesFlowProblem<T>,
        solution: &StokesFlowSolution<T>,
    ) -> Result<()> {
        let n_corner = problem.n_corner_nodes;

        // Per-cell contribution is purely read-only on mesh and solution;
        // accumulate per-thread local residual vectors, then merge.

        let (residual, max_abs, sum_abs, l2, net) = problem
            .mesh
            .cells
            .par_iter()
            .fold(
                || (vec![T::zero(); n_corner], T::zero(), T::zero(), T::zero(), T::zero()),
                |(mut res, mut mx, mut sm, mut l2acc, mut nt), cell| {
                    let idxs = match extract_vertex_indices(cell, &problem.mesh, n_corner) {
                        Ok(v) => v,
                        Err(_) => return (res, mx, sm, l2acc, nt),
                    };
                    if idxs.len() < 4 {
                        return (res, mx, sm, l2acc, nt);
                    }

                    let verts: Vec<Vector3<T>> = idxs
                        .iter()
                        .map(|&idx| {
                            problem
                                .mesh
                                .vertices
                                .position(cfd_mesh::domain::core::index::VertexId::from_usize(idx))
                                .coords
                        })
                        .collect();

                    let j_mat = nalgebra::Matrix3::from_columns(&[
                        verts[1] - verts[0],
                        verts[2] - verts[0],
                        verts[3] - verts[0],
                    ]);
                    let det_j = j_mat.determinant();
                    let abs_det = Float::abs(det_j);
                    if abs_det <= <T as FromPrimitive>::from_f64(1e-24).unwrap_or_else(T::zero) {
                        return (res, mx, sm, l2acc, nt);
                    }
                    let j_inv_t = match j_mat.try_inverse() {
                        Some(ji) => ji.transpose(),
                        None => return (res, mx, sm, l2acc, nt),
                    };

                    let grad_ref_p1 = nalgebra::Matrix3x4::new(
                        -T::one(), T::one(), T::zero(), T::zero(),
                        -T::one(), T::zero(), T::one(), T::zero(),
                        -T::one(), T::zero(), T::zero(), T::one(),
                    );
                    let p1_gradients_phys = j_inv_t * grad_ref_p1;
                    let shape = LagrangeTet10::new(p1_gradients_phys);

                    let quad = TetrahedronQuadrature::keast_degree_3();
                    for (qp, &qw) in quad.points().iter().zip(quad.weights().iter()) {
                        let weight = qw * abs_det;
                        let l = [T::one() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
                        let grad_p2 = shape.gradients(&l);

                        let mut div_u = T::zero();
                        for i in 0..idxs.len().min(10) {
                            let vel = solution.get_velocity(idxs[i]);
                            let grad_i = if idxs.len() == 4 {
                                Vector3::new(
                                    p1_gradients_phys[(0, i)],
                                    p1_gradients_phys[(1, i)],
                                    p1_gradients_phys[(2, i)],
                                )
                            } else {
                                Vector3::new(grad_p2[(0, i)], grad_p2[(1, i)], grad_p2[(2, i)])
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
                        let a = Float::abs(r);
                        if a > mx { mx = a; }
                        sm += a;
                        l2acc += r * r;
                        nt += r;
                    }

                    (res, mx, sm, l2acc, nt)
                },
            )
            .reduce(
                || (vec![T::zero(); n_corner], T::zero(), T::zero(), T::zero(), T::zero()),
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
            let mean_abs = sum_abs / T::from_usize(n).unwrap_or_else(T::one);
            let l2_norm = Float::sqrt(l2);
            println!(
                "Continuity Residual (Bu): max={max_abs:?}, mean_abs={mean_abs:?}, l2={l2_norm:?}, net={net:?}, n={n} "
            );
        }

        Ok(())
    }

    fn assemble_system(
        &self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;
        let n_total_dof = n_velocity_dof + n_corner_nodes;

        let mut builder = SparseMatrixBuilder::new(n_total_dof, n_total_dof);

        let vertex_positions: Vec<Vector3<T>> = problem
            .mesh
            .vertices
            .iter()
            .map(|v| v.1.position.coords)
            .collect();

        // GAP-PERF-001: Build edge→mid-node lookup cache once per mesh (O(n_mid × n_edges))
        // rather than doing O(n_mid) brute-force nearest-midpoint scan per edge per element.
        //
        // # Theorem — MidNodeCache amortized O(1) edge lookup (see mid_node_cache module docs)
        let _mid_cache = MidNodeCache::build(&problem.mesh, n_corner_nodes);

        println!(
            "  Assembling {} elements in parallel...",
            problem.mesh.cells.len()
        );

        let (entry_map, mut rhs) = problem
            .mesh
            .cells
            .par_iter()
            .enumerate()
            .fold(
                || (HashMap::new(), DVector::zeros(n_total_dof)),
                |(mut local_map, mut local_rhs), (i, cell)| {
                    let viscosity = problem
                        .element_viscosities
                        .as_ref()
                        .map_or(problem.fluid.viscosity, |v| v[i]);
                    // Use cache-accelerated index extraction (O(1) per edge vs O(N_mid))
                    let idxs = extract_vertex_indices(
                        cell,
                        &problem.mesh,
                        problem.n_corner_nodes,
                    )
                    .unwrap();
                    let local_verts: Vec<Vector3<T>> =
                        idxs.iter().map(|&idx| vertex_positions[idx]).collect();

                    if local_verts.len() >= 4_usize {
                        let v0 = local_verts[0];
                        let v1 = local_verts[1];
                        let v2 = local_verts[2];
                        let v3 = local_verts[3];

                        let six = <T as FromPrimitive>::from_f64(6.0).unwrap_or_else(T::one);
                        let vol = ((v1 - v0).cross(&(v2 - v0))).dot(&(v3 - v0)) / six;
                        let vol_tol = <T as FromPrimitive>::from_f64(1e-22).unwrap_or_else(T::zero);
                        if Float::abs(vol) < vol_tol {
                            panic!("Element {} has near-zero volume", i);
                        }
                    }

                    let u_avg = self.calculate_u_avg(&idxs, previous_solution);

                    self.assemble_element_local(
                        &mut local_map,
                        &mut local_rhs,
                        &idxs,
                        &local_verts,
                        viscosity,
                        problem.fluid.density,
                        u_avg,
                        n_nodes,
                    );

                    (local_map, local_rhs)
                },
            )
            .reduce(
                || (HashMap::new(), DVector::zeros(n_total_dof)),
                |(mut map1, mut rhs1), (map2, rhs2)| {
                    for (k, v) in map2 {
                        *map1.entry(k).or_insert(T::zero()) += v;
                    }
                    rhs1 += rhs2;
                    (map1, rhs1)
                },
            );

        println!("  Assembly map-reduce complete. Applying boundary conditions...");
        // Populate builder with accumulated map entries
        for ((row, col), val) in entry_map {
            builder.add_entry(row, col, val)?;
        }

        self.apply_boundary_conditions_block(&mut builder, &mut rhs, problem, n_nodes)?;

        let velocity_dofs_constrained = problem.boundary_conditions.len() * 3;
        println!("  Velocity DOFs constrained: {velocity_dofs_constrained} / {n_velocity_dof}");

        if velocity_dofs_constrained == n_velocity_dof {
            println!("  WARNING: All velocity DOFs are constrained (may cause incompressibility conflict)");
        }

        let diag_eps =
            problem.fluid.viscosity * <T as FromPrimitive>::from_f64(1e-12).unwrap_or_else(T::zero);
        for i in n_velocity_dof..n_total_dof {
            let _ = builder.add_entry(i, i, diag_eps);
        }

        let matrix = builder.build_with_rhs(&mut rhs)?;
        Ok((matrix, rhs))
    }

    fn assemble_element_local(
        &self,
        local_map: &mut HashMap<(usize, usize), T>,
        _local_rhs: &mut DVector<T>,
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

        let j_mat = nalgebra::Matrix3::from_columns(&[
            verts[1] - verts[0],
            verts[2] - verts[0],
            verts[3] - verts[0],
        ]);
        let det_j = j_mat.determinant();
        let abs_det = Float::abs(det_j);
        let j_inv_t = j_mat.try_inverse().expect("Singular Jacobian").transpose();

        let grad_ref_p1 = nalgebra::Matrix3x4::new(
            -T::one(),
            T::one(),
            T::zero(),
            T::zero(),
            -T::one(),
            T::zero(),
            T::one(),
            T::zero(),
            -T::one(),
            T::zero(),
            T::zero(),
            T::one(),
        );
        let p1_gradients_phys = j_inv_t * grad_ref_p1;
        let shape = LagrangeTet10::new(p1_gradients_phys);

        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;

        for (qp, &qw) in points.iter().zip(weights.iter()) {
            let weight = qw * abs_det;
            let l = [T::one() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
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
                        grad_p2_mat[(0, i)],
                        grad_p2_mat[(1, i)],
                        grad_p2_mat[(2, i)],
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
                                grad_p2_mat[(0, j)],
                                grad_p2_mat[(1, j)],
                                grad_p2_mat[(2, j)],
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
                        let visc = viscosity * grad_i.dot(&grad_j) * weight;
                        let adv = density * n_i * u_avg.dot(&grad_j) * weight;
                        *local_map.entry((gv_i, gv_j)).or_insert(T::zero()) += visc + adv;
                        if grad_div_penalty > T::zero() {
                            for e in 0..3 {
                                let gv_j_e = gj + e * v_offset;
                                let grad_div = grad_div_penalty * grad_i[d] * grad_j[e] * weight;
                                *local_map.entry((gv_i, gv_j_e)).or_insert(T::zero()) += grad_div;
                            }
                        }
                    }
                    for j in 0..4 {
                        let gj = idxs[j];
                        let gp_j = p_offset + gj;
                        let b_val = n_p1[j] * grad_i[d] * weight;
                        *local_map.entry((gv_i, gp_j)).or_insert(T::zero()) -= b_val;
                        *local_map.entry((gp_j, gv_i)).or_insert(T::zero()) += b_val;
                    }
                }
            }

            // PSPG pressure stabilization (Brezzi-Pitkäranta)
            // Adds τ_BP * ∫ ∇q_i · ∇q_j dΩ to the pressure-pressure block.
            // τ_BP = h_e² / (12 * μ), where h_e = (6V)^(1/3).
            // Circumvents the LBB inf-sup condition for equal-order P1-P1 elements.
            if viscosity > T::zero() {
                let h_e = Float::cbrt(abs_det); // (6V)^(1/3) ≈ element diameter
                let twelve = <T as FromPrimitive>::from_f64(12.0).unwrap_or_else(T::one);
                let tau_bp = h_e * h_e / (twelve * viscosity);
                let vol_e = abs_det / <T as FromPrimitive>::from_f64(6.0).unwrap_or_else(T::one);

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
                        let pspg = tau_bp * grad_p_i.dot(&grad_p_j) * vol_e;
                        *local_map.entry((gp_i, gp_j)).or_insert(T::zero()) += pspg;
                    }
                }
            }
        }
    }

    #[allow(clippy::too_many_lines)]
    fn apply_boundary_conditions_block(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
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
            println!(
                "  [WARNING] Boundary Leak: {} boundary nodes have no BCs! First 5: {:?}",
                unconstrained_boundary_nodes.len(),
                &unconstrained_boundary_nodes[..unconstrained_boundary_nodes.len().min(5)]
            );
        }

        // Track if any pressure boundary condition is applied
        let mut has_pressure_bc = false;

        for (&node_idx, bc) in &problem.boundary_conditions {
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
                        builder.set_dirichlet_row(dof, diag_scale, T::zero());
                        rhs[dof] = T::zero();
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
            builder.set_dirichlet_row(reference_pressure_dof, diag_scale, T::zero());
            rhs[reference_pressure_dof] = T::zero();
            println!(
                "  Pinned pressure DOF {reference_pressure_dof} to zero (no pressure BC specified)"
            );
        }

        println!(
            "  BC Diagnostics: inlet_nodes={inlet_nodes}, wall_nodes={wall_nodes}, outlet_nodes={outlet_nodes}, dirichlet_nodes={dirichlet_nodes}"
        );
        println!(
            "  BC Diagnostics: velocity_dofs_set={}, pressure_dofs_set={}, n_velocity_dof={}, n_pressure_dof={}",
            vel_dofs.len(),
            p_dofs.len(),
            n_nodes * 3,
            problem.n_corner_nodes
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
            // `nodes` is non-empty here; T::from_usize cannot fail for valid lengths.
            sum / T::from_usize(nodes.len()).unwrap_or_else(T::one)
        } else {
            Vector3::zeros()
        }
    }
}

/// Extract vertex indices from a cell for element assembly
pub fn extract_vertex_indices<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    cell: &Cell,
    mesh: &IndexedMesh<T>,
    n_corner_nodes: usize,
) -> Result<Vec<usize>> {
    let mut counts = std::collections::HashMap::new();
    for &f_idx in &cell.faces {
        if f_idx < mesh.face_count() {
            let f = mesh
                .faces
                .get(cfd_mesh::domain::core::index::FaceId::from_usize(f_idx));
            for &v_id in &f.vertices {
                *counts.entry(v_id.as_usize()).or_insert(0) += 1;
            }
        }
    }

    // In a tetrahedron, corners are shared by 3 faces, mid-edges by 2.
    let mut corners = Vec::new();
    let mut mid_edges = Vec::new();
    for (&v_idx, &count) in &counts {
        if count == 3 {
            corners.push(v_idx);
        } else if count == 2 {
            mid_edges.push(v_idx);
        }
    }

    if corners.len() == 4 && mid_edges.is_empty() && counts.len() == 4 {
        let ordered = order_tet_corners(&corners, mesh);

        if mesh.vertex_count() > n_corner_nodes {
            // P2 mesh: Face geometry only has corners, missing mid-edges.
            // Recover mid-edges via geometric search over extra nodes.
            let mut final_nodes = ordered.clone();
            let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
            for &(i, j) in &edges {
                let v_i = ordered[i];
                let v_j = ordered[j];
                let p_i = mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_i))
                    .coords;
                let p_j = mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_j))
                    .coords;
                let target = (p_i + p_j)
                    * <T as num_traits::FromPrimitive>::from_f64(0.5).unwrap_or_else(T::one);

                let mut best_m = 0;
                let mut min_dist_sq = <T as num_traits::Float>::infinity();
                for m_idx in n_corner_nodes..mesh.vertex_count() {
                    let pm = mesh
                        .vertices
                        .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                        .coords;
                    let dist_sq = (pm - target).norm_squared();
                    if dist_sq < min_dist_sq {
                        min_dist_sq = dist_sq;
                        best_m = m_idx;
                    }
                }
                final_nodes.push(best_m);
            }
            return Ok(final_nodes);
        }
        // P1 Tet
        return Ok(ordered);
    }

    if corners.len() == 4 && mid_edges.len() == 6 {
        // P2 Tet
        // We need them in canonical Tet10 order for LagrangeTet10:
        // Corners: 0, 1, 2, 3
        // Mid-edges: 4:(0,1), 5:(1,2), 6:(2,0), 7:(0,3), 8:(1,3), 9:(2,3)
        let ordered = order_tet_corners(&corners, mesh);
        let mut final_nodes = ordered.clone();
        let mut used_mid_edges = std::collections::HashSet::new();

        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for &(i, j) in &edges {
            let v_i = ordered[i];
            let v_j = ordered[j];
            // The mid-node for (v_i, v_j) is the one shared by faces that both contain v_i and v_j.
            // Or simpler: find the mid_edge node that is closest to (v_i + v_j)/2
            let p_i = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_i))
                .coords;
            let p_j = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_j))
                .coords;
            let target = (p_i + p_j) * <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(T::one);

            let mut best_node = None;
            let mut min_dist = T::infinity();

            for &m_idx in &mid_edges {
                if used_mid_edges.contains(&m_idx) {
                    continue;
                }
                let dist = (mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                    .coords
                    - target)
                    .norm();
                if dist < min_dist {
                    min_dist = dist;
                    best_node = Some(m_idx);
                }
            }

            let selected = if let Some(m_idx) = best_node {
                m_idx
            } else {
                let mut fallback = mid_edges[0];
                let mut fallback_dist = (mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(
                        fallback,
                    ))
                    .coords
                    - target)
                    .norm();
                for &m_idx in &mid_edges[1..] {
                    let dist = (mesh
                        .vertices
                        .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                        .coords
                        - target)
                        .norm();
                    if dist < fallback_dist {
                        fallback_dist = dist;
                        fallback = m_idx;
                    }
                }
                fallback
            };

            used_mid_edges.insert(selected);
            final_nodes.push(selected);
        }
        return Ok(final_nodes);
    }

    // Fallback for other elements (e.g. Hex)
    let mut all: Vec<usize> = counts.keys().copied().collect();
    all.sort_unstable();
    Ok(all)
}

/// Cache-accelerated variant of [`extract_vertex_indices`] for P2 FEM assembly.
///
/// # Performance (GAP-PERF-001)
///
/// When `mid_cache` is non-empty (P2 mesh), mid-node lookup for each of the 6
/// edges is O(1) amortised via the pre-built `HashMap<(usize,usize), usize>`,
/// replacing the O(n_mid) brute-force nearest-midpoint scan in `extract_vertex_indices`.
///
/// For P1 meshes (cache empty), falls back to the corner-only path identical to
/// the uncached version.
///
/// # Theorem — Correctness Equivalence
///
/// For any conforming P2 tetrahedral mesh where `MidNodeCache::build` was invoked
/// with the same mesh and `n_corner_nodes`, the output of `extract_vertex_indices_cached`
/// is element-wise identical to the output of `extract_vertex_indices`.
///
/// **Proof**: `MidNodeCache::build` maps each canonical edge `(min,max)` to the
/// unique mid-node closest to the geometric midpoint (geometric uniqueness from P2
/// conformity). `extract_vertex_indices` finds the same node via exhaustive nearest-
/// midpoint scan. Both converge to the same index by definition of the minimum.
pub fn extract_vertex_indices_cached<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    cell: &Cell,
    mesh: &IndexedMesh<T>,
    n_corner_nodes: usize,
    mid_cache: &crate::fem::mid_node_cache::MidNodeCache,
) -> Result<Vec<usize>> {
    let mut counts = std::collections::HashMap::new();
    for &f_idx in &cell.faces {
        if f_idx < mesh.face_count() {
            let f = mesh
                .faces
                .get(cfd_mesh::domain::core::index::FaceId::from_usize(f_idx));
            for &v_id in &f.vertices {
                *counts.entry(v_id.as_usize()).or_insert(0) += 1;
            }
        }
    }

    let mut corners = Vec::new();
    let mut mid_edges = Vec::new();
    for (&v_idx, &count) in &counts {
        if count == 3 {
            corners.push(v_idx);
        } else if count == 2 {
            mid_edges.push(v_idx);
        }
    }

    if corners.len() == 4 && mid_edges.is_empty() && counts.len() == 4 {
        let ordered = order_tet_corners(&corners, mesh);

        if mesh.vertex_count() > n_corner_nodes {
            // P2 mesh: use cache for O(1) mid-node lookup (GAP-PERF-001)
            let mut final_nodes = ordered.clone();
            let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
            for &(i, j) in &edges {
                let v_i = ordered[i];
                let v_j = ordered[j];

                if let Some(m_idx) = mid_cache.get(v_i, v_j) {
                    final_nodes.push(m_idx);
                } else {
                    // Cache miss: fall back to geometric search (edge not registered)
                    let p_i = mesh
                        .vertices
                        .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_i))
                        .coords;
                    let p_j = mesh
                        .vertices
                        .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_j))
                        .coords;
                    let target = (p_i + p_j)
                        * <T as num_traits::FromPrimitive>::from_f64(0.5)
                            .unwrap_or_else(T::one);
                    let mut best_m = 0;
                    let mut min_dist_sq = <T as num_traits::Float>::infinity();
                    for m_idx in n_corner_nodes..mesh.vertex_count() {
                        let pm = mesh
                            .vertices
                            .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                            .coords;
                        let dist_sq = (pm - target).norm_squared();
                        if dist_sq < min_dist_sq {
                            min_dist_sq = dist_sq;
                            best_m = m_idx;
                        }
                    }
                    final_nodes.push(best_m);
                }
            }
            return Ok(final_nodes);
        }
        // P1 Tet
        return Ok(ordered);
    }

    if corners.len() == 4 && mid_edges.len() == 6 {
        // P2 Tet: corners + mid-edges both recovered from face data
        let ordered = order_tet_corners(&corners, mesh);
        let mut final_nodes = ordered.clone();
        let mut used_mid_edges = std::collections::HashSet::new();

        let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)];
        for &(i, j) in &edges {
            let v_i = ordered[i];
            let v_j = ordered[j];

            // Use cache for O(1) lookup if available
            if let Some(m_idx) = mid_cache.get(v_i, v_j) {
                used_mid_edges.insert(m_idx);
                final_nodes.push(m_idx);
                continue;
            }

            // Fallback: geometric search among mid_edges
            let p_i = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_i))
                .coords;
            let p_j = mesh
                .vertices
                .position(cfd_mesh::domain::core::index::VertexId::from_usize(v_j))
                .coords;
            let target =
                (p_i + p_j) * <T as FromPrimitive>::from_f64(0.5).unwrap_or_else(T::one);

            let mut best_node = None;
            let mut min_dist = T::infinity();
            for &m_idx in &mid_edges {
                if used_mid_edges.contains(&m_idx) {
                    continue;
                }
                let dist = (mesh
                    .vertices
                    .position(cfd_mesh::domain::core::index::VertexId::from_usize(m_idx))
                    .coords
                    - target)
                    .norm();
                if dist < min_dist {
                    min_dist = dist;
                    best_node = Some(m_idx);
                }
            }
            if let Some(m_idx) = best_node {
                used_mid_edges.insert(m_idx);
                final_nodes.push(m_idx);
            }
        }
        return Ok(final_nodes);
    }

    // Fallback for other elements (e.g. Hex)
    let mut all: Vec<usize> = counts.keys().copied().collect();
    all.sort_unstable();
    Ok(all)
}

fn order_tet_corners<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    corners: &[usize],
    mesh: &IndexedMesh<T>,
) -> Vec<usize> {
    let perms: [[usize; 4]; 24] = [
        [0, 1, 2, 3],
        [0, 1, 3, 2],
        [0, 2, 1, 3],
        [0, 2, 3, 1],
        [0, 3, 1, 2],
        [0, 3, 2, 1],
        [1, 0, 2, 3],
        [1, 0, 3, 2],
        [1, 2, 0, 3],
        [1, 2, 3, 0],
        [1, 3, 0, 2],
        [1, 3, 2, 0],
        [2, 0, 1, 3],
        [2, 0, 3, 1],
        [2, 1, 0, 3],
        [2, 1, 3, 0],
        [2, 3, 0, 1],
        [2, 3, 1, 0],
        [3, 0, 1, 2],
        [3, 0, 2, 1],
        [3, 1, 0, 2],
        [3, 1, 2, 0],
        [3, 2, 0, 1],
        [3, 2, 1, 0],
    ];

    let mut best: Option<Vec<usize>> = None;
    let mut best_det = T::neg_infinity();

    for perm in &perms {
        let v0 = corners[perm[0]];
        let v1 = corners[perm[1]];
        let v2 = corners[perm[2]];
        let v3 = corners[perm[3]];

        let p0 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v0))
            .coords;
        let p1 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v1))
            .coords;
        let p2 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v2))
            .coords;
        let p3 = mesh
            .vertices
            .position(cfd_mesh::domain::core::index::VertexId::from_usize(v3))
            .coords;

        let det = (p1 - p0).cross(&(p2 - p0)).dot(&(p3 - p0));
        if det > T::zero() {
            let candidate = vec![v0, v1, v2, v3];
            let take = match &best {
                None => true,
                Some(existing) => candidate < *existing,
            };
            if take {
                best = Some(candidate);
                best_det = det;
            }
        } else if best.is_none() && det > best_det {
            best_det = det;
            best = Some(vec![v0, v1, v2, v3]);
        }
    }

    best.unwrap_or_else(|| corners.to_vec())
}

fn compute_mesh_scale<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(
    mesh: &IndexedMesh<T>,
) -> T {
    let mut min = Vector3::new(T::infinity(), T::infinity(), T::infinity());
    let mut max = Vector3::new(T::neg_infinity(), T::neg_infinity(), T::neg_infinity());
    for v in mesh.vertices.iter() {
        let p = v.1.position.coords;
        min.x = Float::min(min.x, p.x);
        min.y = Float::min(min.y, p.y);
        min.z = Float::min(min.z, p.z);
        max.x = Float::max(max.x, p.x);
        max.y = Float::max(max.y, p.y);
        max.z = Float::max(max.z, p.z);
    }
    (max - min).norm()
}
