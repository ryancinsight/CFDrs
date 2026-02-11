//! # Finite Element Method (FEM) Solver for 3D Incompressible Flow
//!
//! This module implements a high-performance FEM solver for the incompressible
//! Navier-Stokes equations using mixed Taylor-Hood (P2-P1) formulation.

use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::linear_solver::{GMRES, LinearSolver, IncompleteLU};
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
        problem.validate()?;

        let (matrix, rhs) = self.assemble_system(problem, previous_solution)?;
        
        let n_total_dof = rhs.len();
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = self.count_corner_nodes(&problem.mesh);
        let n_velocity_dof = n_nodes * 3;

        let mut x = if let Some(sol) = previous_solution {
            sol.interleave()
        } else {
            DVector::zeros(n_total_dof)
        };

        let rel_tol = T::from_f64(1e-8).unwrap_or_else(T::zero);
        let abs_tol = Float::max(rel_tol * rhs.norm(), T::from_f64(1e-14).unwrap_or_else(T::zero));
        
        let mut gmres_config = cfd_math::linear_solver::IterativeSolverConfig::default();
        gmres_config.max_iterations = 10000;
        gmres_config.tolerance = abs_tol;

        let restart = 200.min(n_total_dof);
        
        println!("  FEM Solver: System size {}x{}", n_total_dof, n_total_dof);
        
        let preconditioner = IncompleteLU::new(&matrix)
             .map_err(|e| Error::Solver(format!("ILU failed: {}", e)))?;
             
        let solver = GMRES::new(gmres_config, restart);
        let monitor = solver.solve_preconditioned(&matrix, &rhs, &preconditioner, &mut x)
            .map_err(|e| Error::Solver(format!("GMRES failed: {}", e)))?;
            
        println!("  Taylor-Hood FEM Solver: Converged in {} iterations, resid={:?}", 
            monitor.iteration, monitor.residual_history.last());

        let velocity = x.rows(0, n_velocity_dof).into_owned();
        let pressure = x.rows(n_velocity_dof, n_corner_nodes).into_owned();

        Ok(StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes))
    }

    fn count_corner_nodes(&self, mesh: &Mesh<T>) -> usize {
        mesh.cells().iter()
            .flat_map(|c| extract_vertex_indices(c, mesh).unwrap_or_default().into_iter().take(4))
            .fold(0, |max, v| if v > max { v } else { max }) + 1
    }

    fn assemble_system(
        &self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = self.count_corner_nodes(&problem.mesh);
        let n_velocity_dof = n_nodes * 3;
        let n_total_dof = n_velocity_dof + n_corner_nodes;

        let mut builder = SparseMatrixBuilder::new(n_total_dof, n_total_dof);
        let mut rhs = DVector::zeros(n_total_dof);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices().iter()
            .map(|v| v.position.coords).collect();

        for (i, cell) in problem.mesh.cells().iter().enumerate() {
            let viscosity = problem.element_viscosities.as_ref().map_or(problem.fluid.viscosity, |v| v[i]);
            let idxs = extract_vertex_indices(cell, &problem.mesh)?;
            let local_verts: Vec<Vector3<T>> = idxs.iter().map(|&idx| vertex_positions[idx]).collect();
            let u_avg = self.calculate_u_avg(&idxs, previous_solution);

            self.assemble_element_block(&mut builder, &mut rhs, &idxs, &local_verts, viscosity, problem.fluid.density, u_avg, n_nodes)?;
        }

        self.apply_boundary_conditions_block(&mut builder, &mut rhs, problem, n_nodes)?;

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

        for (&node_idx, bc) in &problem.boundary_conditions {
            match bc {
                BoundaryCondition::VelocityInlet { velocity } => {
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        builder.set_dirichlet_row(dof, diag_scale, velocity[d]);
                        rhs[dof] = velocity[d] * diag_scale;
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        builder.set_dirichlet_row(dof, diag_scale, T::zero());
                        rhs[dof] = T::zero();
                    }
                }
                BoundaryCondition::PressureOutlet { pressure } | BoundaryCondition::PressureInlet { pressure, .. } => {
                    let dof = p_offset + node_idx;
                    builder.set_dirichlet_row(dof, diag_scale, *pressure);
                    rhs[dof] = *pressure * diag_scale;
                }
                BoundaryCondition::Dirichlet { value, component_values } => {
                    if let Some(comps) = component_values {
                        for d in 0..3 {
                            if let Some(Some(val)) = comps.get(d) {
                                let dof = node_idx + d * v_offset;
                                builder.set_dirichlet_row(dof, diag_scale, *val);
                                rhs[dof] = *val * diag_scale;
                            }
                        }
                        if let Some(Some(p_val)) = comps.get(3) {
                            let dof = p_offset + node_idx;
                            builder.set_dirichlet_row(dof, diag_scale, *p_val);
                            rhs[dof] = *p_val * diag_scale;
                        }
                    } else {
                        // Scalar Dirichlet: apply to all velocity components (standard wall/inlet)
                        // This usually isn't what's desired for pressure, but we follow the old logic.
                        for d in 0..3 {
                            let dof = node_idx + d * v_offset;
                            builder.set_dirichlet_row(dof, diag_scale, *value);
                            rhs[dof] = *value * diag_scale;
                        }
                    }
                }
                _ => {}
            }
        }
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

pub fn extract_vertex_indices<T: RealField + Copy>(cell: &Cell, mesh: &Mesh<T>) -> Result<Vec<usize>> {
    let mut indices = Vec::new();
    for &f_idx in &cell.faces {
        if let Some(f) = mesh.face(f_idx) {
            for &v_idx in &f.vertices {
                if !indices.contains(&v_idx) { indices.push(v_idx); }
            }
        }
    }
    Ok(indices)
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
