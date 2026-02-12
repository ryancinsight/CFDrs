//! # Finite Element Method (FEM) Solver for 3D Incompressible Flow
//!
//! This module implements a high-performance FEM solver for the incompressible
//! Navier-Stokes equations using mixed Taylor-Hood (P2-P1) formulation.

use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::{
    BiCGSTAB, BlockDiagonalPreconditioner, DirectSparseSolver, GMRES, IncompleteLU, LinearSolver,
    Preconditioner,
};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField, Vector3};
use num_traits::{Float, FromPrimitive};
use tracing;

use crate::fem::{FemConfig, StokesFlowProblem, StokesFlowSolution};
use crate::fem::shape_functions::LagrangeTet10;
use crate::fem::quadrature::TetrahedronQuadrature;
use cfd_mesh::mesh::Mesh;
use cfd_mesh::topology::Cell;

/// Finite Element Method solver for 3D incompressible flow
pub struct FemSolver<T: RealField + Copy + FromPrimitive + num_traits::Float + std::fmt::Debug> {
    _config: FemConfig<T>,
    linear_solver: GMRES<T>,
}

impl<T: RealField + FromPrimitive + Copy + Float + std::fmt::Debug + From<f64>> FemSolver<T> {
    pub fn new(config: FemConfig<T>) -> Self {
        let mut solver_config = cfd_math::linear_solver::IterativeSolverConfig::default();
        solver_config.max_iterations = 30000;
        solver_config.tolerance = T::from_f64(1e-12).unwrap_or_else(T::zero);
        let linear_solver = GMRES::new(solver_config, 100);
        Self { linear_solver, _config: config }
    }

    pub fn solve(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<StokesFlowSolution<T>> {
        tracing::info!("Starting Taylor-Hood Stokes solver");
        
        // Validate problem configuration before solving
        problem.validate()?;

        let (matrix, rhs) = self.assemble_system(problem, previous_solution)?;
        
        let n_total_dof = rhs.len();
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;

        let mut x = if let Some(sol) = previous_solution {
            sol.interleave()
        } else {
            DVector::zeros(n_total_dof)
        };

        let direct_solver = DirectSparseSolver::default();
        if direct_solver.can_handle_size(n_total_dof) {
            println!("  Attempting direct sparse LU solve (rsparse)...");
            match direct_solver.solve(&matrix, &rhs) {
                Ok(x_direct) => {
                    x = x_direct;
                    let velocity = x.rows(0, n_velocity_dof).into_owned();
                    let pressure = x.rows(n_velocity_dof, n_corner_nodes).into_owned();
                    println!("  Direct solve completed successfully.");
                    return Ok(StokesFlowSolution::new_with_corners(
                        velocity,
                        pressure,
                        n_nodes,
                        n_corner_nodes,
                    ));
                }
                Err(err) => {
                    println!("  Direct solve failed: {}", err);
                    println!("  Falling back to iterative solvers...");
                }
            }
        } else {
            println!(
                "  Skipping direct solve: system size {} exceeds limit {}",
                n_total_dof, direct_solver.max_size
            );
        }

        let rel_tol = T::from_f64(1e-8).unwrap_or_else(T::zero);
        let abs_tol = Float::max(rel_tol * rhs.norm(), T::from_f64(1e-14).unwrap_or_else(T::zero));
        
        let mut gmres_config = cfd_math::linear_solver::IterativeSolverConfig::default();
        gmres_config.max_iterations = 50000; // Increased for complex 3D problems
        gmres_config.tolerance = abs_tol;

        let restart = 200.min(n_total_dof);
        
        println!("  FEM Solver: System size {}x{}", n_total_dof, n_total_dof);
        println!("  GMRES: restart={}, max_iter={}", restart, gmres_config.max_iterations);
        
        // Try block diagonal preconditioner (specialized for saddle-point systems)
        println!("  Attempting block diagonal preconditioner for saddle-point system...");
        let block_precond = BlockDiagonalPreconditioner::new(&matrix, n_velocity_dof, n_corner_nodes)
            .map_err(|e| Error::Solver(format!("Block preconditioner failed: {}", e)))?;
        
        let solver = GMRES::new(gmres_config.clone(), restart);
        
        // First attempt with block preconditioner
        match solver.solve_preconditioned(&matrix, &rhs, &block_precond, &mut x) {
            Ok(monitor) => {
                println!("  ✓ Block diagonal preconditioner: Converged in {} iterations, resid={:?}", 
                    monitor.iteration, monitor.residual_history.last());
            }
            Err(block_err) => {
                let block_err_msg = format!("{block_err}");
                println!("  ✗ Block diagonal preconditioner failed: {}", block_err_msg);
                println!("  Falling back to unpreconditioned GMRES...");

                // Reset initial guess before each fallback strategy.
                x.fill(T::zero());
                match solver.solve_unpreconditioned(&matrix, &rhs, &mut x) {
                    Ok(monitor) => {
                        println!(
                            "  ✓ Unpreconditioned GMRES: Converged in {} iterations, resid={:?}",
                            monitor.iteration,
                            monitor.residual_history.last()
                        );
                    }
                    Err(gmres_unprec_err) => {
                        let gmres_unprec_err_msg = format!("{gmres_unprec_err}");
                        println!("  ✗ Unpreconditioned GMRES failed: {}", gmres_unprec_err_msg);
                        println!("  Falling back to ILU-preconditioned GMRES...");

                        x.fill(T::zero());
                        let ilu_result = IncompleteLU::new(&matrix);
                        let gmres_ilu_result = if let Ok(ilu_preconditioner) = ilu_result {
                            solver.solve_preconditioned(&matrix, &rhs, &ilu_preconditioner, &mut x)
                        } else {
                            Err(cfd_core::error::Error::Solver(
                                "ILU preconditioner construction failed".to_string(),
                            ))
                        };

                        match gmres_ilu_result {
                            Ok(monitor) => {
                                println!(
                                    "  ✓ ILU preconditioned GMRES: Converged in {} iterations, resid={:?}",
                                    monitor.iteration,
                                    monitor.residual_history.last()
                                );
                            }
                            Err(gmres_ilu_err) => {
                                let gmres_ilu_err_msg = format!("{gmres_ilu_err}");
                                println!("  ✗ ILU-preconditioned GMRES failed: {}", gmres_ilu_err_msg);
                                println!("  Falling back to unpreconditioned BiCGSTAB...");

                                x.fill(T::zero());
                                let bicg = BiCGSTAB::new(gmres_config.clone());
                                let monitor = bicg
                                    .solve_unpreconditioned(&matrix, &rhs, &mut x)
                                    .map_err(|bicg_err| {
                                        Error::Solver(format!(
                                            "All linear solver attempts failed. block={}, gmres_unpreconditioned={}, gmres_ilu={}, bicgstab={}",
                                            block_err_msg, gmres_unprec_err_msg, gmres_ilu_err_msg, bicg_err
                                        ))
                                    })?;

                                println!(
                                    "  ✓ Unpreconditioned BiCGSTAB: Converged in {} iterations, resid={:?}",
                                    monitor.iteration,
                                    monitor.residual_history.last()
                                );
                            }
                        }
                    }
                }
            }
        }

        let velocity = x.rows(0, n_velocity_dof).into_owned();
        let pressure = x.rows(n_velocity_dof, n_corner_nodes).into_owned();

        Ok(StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes))
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
        let mut rhs = DVector::zeros(n_total_dof);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices().iter()
            .map(|v| v.position.coords).collect();

        println!("  Assembling {} elements...", problem.mesh.cells().len());
        
        for (i, cell) in problem.mesh.cells().iter().enumerate() {
            let viscosity = problem.element_viscosities.as_ref().map_or(problem.fluid.viscosity, |v| v[i]);
            let idxs = extract_vertex_indices(cell, &problem.mesh)?;
            let local_verts: Vec<Vector3<T>> = idxs.iter().map(|&idx| vertex_positions[idx]).collect();
            
            // Check for degenerate elements
            if local_verts.len() >= 4 {
                let v0 = local_verts[0];
                let v1 = local_verts[1];
                let v2 = local_verts[2];
                let v3 = local_verts[3];
                
                let vol = ((v1 - v0).cross(&(v2 - v0))).dot(&(v3 - v0)) / T::from_f64(6.0).unwrap();
                if Float::abs(vol) < T::from_f64(1e-15).unwrap() {
                    return Err(Error::Solver(format!(
                        "Element {} has near-zero volume: {:?}. Vertices: {:?}, {:?}, {:?}, {:?}",
                        i, vol, v0, v1, v2, v3
                    )));
                }
            }
            
            let u_avg = self.calculate_u_avg(&idxs, previous_solution);

            self.assemble_element_block(&mut builder, &mut rhs, &idxs, &local_verts, viscosity, problem.fluid.density, u_avg, n_nodes)?;
        }

        println!("  Assembly complete. Applying boundary conditions...");
        self.apply_boundary_conditions_block(&mut builder, &mut rhs, problem, n_nodes)?;

        // Count Dirichlet DOFs
        let velocity_dofs_constrained = problem.boundary_conditions.len() * 3;
        println!("  Velocity DOFs constrained: {} / {}", velocity_dofs_constrained, n_velocity_dof);
        
        if velocity_dofs_constrained == n_velocity_dof {
            println!("  WARNING: All velocity DOFs are constrained (may cause incompressibility conflict)");
        }

        let diag_eps = problem.fluid.viscosity * T::from_f64(1e-12).unwrap();
        for i in n_velocity_dof..n_total_dof {
            let _ = builder.add_entry(i, i, diag_eps);
        }

        let matrix = builder.build_with_rhs(&mut rhs)?;
        Ok((matrix, rhs))
    }

    fn assemble_element_block(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        _rhs: &mut DVector<T>,
        idxs: &[usize],
        verts: &[Vector3<T>],
        viscosity: T,
        density: T,
        u_avg: Vector3<T>,
        n_nodes: usize,
    ) -> Result<()> {
        let quad = TetrahedronQuadrature::keast_degree_3();
        let points = quad.points();
        let weights = quad.weights();

        let j_mat = nalgebra::Matrix3::from_columns(&[
            verts[1]-verts[0],
            verts[2]-verts[0],
            verts[3]-verts[0]
        ]);
        let det_j = j_mat.determinant();
        let abs_det = Float::abs(det_j);
        let j_inv_t = j_mat.try_inverse().ok_or_else(|| Error::Solver("Singular Jacobian".to_string()))?.transpose();

        let grad_ref_p1 = nalgebra::Matrix3x4::new(
            T::one(), T::zero(), T::zero(), -T::one(),
            T::zero(), T::one(), T::zero(), -T::one(),
            T::zero(), T::zero(), T::one(), -T::one(),
        );
        let p1_gradients_phys = j_inv_t * grad_ref_p1;
        let shape = LagrangeTet10::new(p1_gradients_phys);

        let v_offset = n_nodes;
        let p_offset = n_nodes * 3;

        for (qp, &qw) in points.iter().zip(weights.iter()) {
            let weight = qw * abs_det;
            let l = [qp.x, qp.y, qp.z, T::one() - qp.x - qp.y - qp.z];
            let n_p2 = shape.values(&l);
            let grad_p2_mat = shape.gradients(&l);
            let n_p1 = l;

            for i in 0..idxs.len().min(10) {
                let gi = idxs[i];
                let grad_i = grad_p2_mat.column(i);
                for d in 0..3 {
                    let gv_i = gi + d * v_offset;
                    for j in 0..idxs.len().min(10) {
                        let gj = idxs[j];
                        let gv_j = gj + d * v_offset;
                        let grad_j = grad_p2_mat.column(j);
                        let visc = viscosity * grad_i.dot(&grad_j) * weight;
                        builder.add_entry(gv_i, gv_j, visc)?;
                        let adv = density * n_p2[i] * u_avg.dot(&grad_j) * weight;
                        builder.add_entry(gv_i, gv_j, adv)?;
                    }
                    for j in 0..4 {
                        let gj = idxs[j];
                        let gp_j = p_offset + gj;
                        let b_val = -n_p1[j] * grad_i[d] * weight;
                        builder.add_entry(gv_i, gp_j, b_val)?;
                        builder.add_entry(gp_j, gv_i, b_val)?;
                    }
                }
            }
        }
        Ok(())
    }

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
                BoundaryCondition::PressureOutlet { pressure } | BoundaryCondition::PressureInlet { pressure, .. } => {
                    outlet_nodes += 1;
                    if node_idx < problem.n_corner_nodes {
                        has_pressure_bc = true;
                        let dof = p_offset + node_idx;
                        builder.set_dirichlet_row(dof, diag_scale, *pressure);
                        rhs[dof] = *pressure * diag_scale;
                        p_dofs.insert(dof);
                    }
                }
                BoundaryCondition::Dirichlet { value, component_values } => {
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
            println!("  Pinned pressure DOF {} to zero (no pressure BC specified)", reference_pressure_dof);
        }

        println!(
            "  BC Diagnostics: inlet_nodes={}, wall_nodes={}, outlet_nodes={}, dirichlet_nodes={}",
            inlet_nodes, wall_nodes, outlet_nodes, dirichlet_nodes
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

    fn calculate_u_avg(&self, nodes: &[usize], solution: Option<&StokesFlowSolution<T>>) -> Vector3<T> {
        if let Some(sol) = solution {
            let mut sum = Vector3::zeros();
            for &n in nodes { sum += sol.get_velocity(n); }
            sum / T::from_usize(nodes.len()).unwrap()
        } else {
            Vector3::zeros()
        }
    }
}

pub fn extract_vertex_indices<T: RealField + Copy + Float>(cell: &Cell, mesh: &Mesh<T>) -> Result<Vec<usize>> {
    let mut counts = std::collections::HashMap::new();
    for &f_idx in &cell.faces {
        if let Some(f) = mesh.face(f_idx) {
            for &v_idx in &f.vertices {
                *counts.entry(v_idx).or_insert(0) += 1;
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
        // P1 Tet
        let ordered = order_tet_corners(&corners, mesh);
        return Ok(ordered);
    }
    
    if corners.len() == 4 && mid_edges.len() == 6 {
        // P2 Tet
        // We need them in canonical Tet10 order for LagrangeTet10:
        // Corners: 0, 1, 2, 3
        // Mid-edges: 4:(0,1), 5:(1,2), 6:(2,0), 7:(0,3), 8:(1,3), 9:(2,3)
        let ordered = order_tet_corners(&corners, mesh);
        let mut final_nodes = ordered.clone();
        
        let edges = [(0,1), (1,2), (2,0), (0,3), (1,3), (2,3)];
        for &(i, j) in &edges {
            let v_i = ordered[i];
            let v_j = ordered[j];
            // The mid-node for (v_i, v_j) is the one shared by faces that both contain v_i and v_j.
            // Or simpler: find the mid_edge node that is closest to (v_i + v_j)/2
            let p_i = mesh.vertex(v_i).unwrap().position.coords;
            let p_j = mesh.vertex(v_j).unwrap().position.coords;
            let target = (p_i + p_j) * T::from_f64(0.5).unwrap();
            
            let mut best_node = mid_edges[0];
            let mut min_dist = (mesh.vertex(best_node).unwrap().position.coords - target).norm();
            for &m_idx in &mid_edges[1..] {
                let dist = (mesh.vertex(m_idx).unwrap().position.coords - target).norm();
                if dist < min_dist {
                    min_dist = dist;
                    best_node = m_idx;
                }
            }
            final_nodes.push(best_node);
        }
        return Ok(final_nodes);
    }
    
    // Fallback for other elements (e.g. Hex)
    let mut all: Vec<usize> = counts.keys().copied().collect();
    all.sort_unstable();
    Ok(all)
}

fn order_tet_corners<T: RealField + Copy + Float>(corners: &[usize], mesh: &Mesh<T>) -> Vec<usize> {
    let perms: [[usize; 4]; 24] = [
        [0, 1, 2, 3], [0, 1, 3, 2], [0, 2, 1, 3], [0, 2, 3, 1], [0, 3, 1, 2], [0, 3, 2, 1],
        [1, 0, 2, 3], [1, 0, 3, 2], [1, 2, 0, 3], [1, 2, 3, 0], [1, 3, 0, 2], [1, 3, 2, 0],
        [2, 0, 1, 3], [2, 0, 3, 1], [2, 1, 0, 3], [2, 1, 3, 0], [2, 3, 0, 1], [2, 3, 1, 0],
        [3, 0, 1, 2], [3, 0, 2, 1], [3, 1, 0, 2], [3, 1, 2, 0], [3, 2, 0, 1], [3, 2, 1, 0],
    ];

    let mut best: Option<Vec<usize>> = None;
    let mut best_det = T::neg_infinity();

    for perm in perms.iter() {
        let v0 = corners[perm[0]];
        let v1 = corners[perm[1]];
        let v2 = corners[perm[2]];
        let v3 = corners[perm[3]];

        let p0 = mesh.vertex(v0).unwrap().position.coords;
        let p1 = mesh.vertex(v1).unwrap().position.coords;
        let p2 = mesh.vertex(v2).unwrap().position.coords;
        let p3 = mesh.vertex(v3).unwrap().position.coords;

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

fn compute_mesh_scale<T: RealField + Copy + Float>(mesh: &Mesh<T>) -> T {
    let mut min = Vector3::new(T::infinity(), T::infinity(), T::infinity());
    let mut max = Vector3::new(T::neg_infinity(), T::neg_infinity(), T::neg_infinity());
    for v in mesh.vertices() {
        let p = v.position.coords;
        min.x = Float::min(min.x, p.x); 
        min.y = Float::min(min.y, p.y); 
        min.z = Float::min(min.z, p.z);
        max.x = Float::max(max.x, p.x); 
        max.y = Float::max(max.y, p.y); 
        max.z = Float::max(max.z, p.z);
    }
    (max - min).norm()
}
