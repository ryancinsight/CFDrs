//! # Pressure Projection Solver for Incompressible Flow
//!
//! Implements Chorin's projection method (also known as fractional step method)
//! for solving incompressible Navier-Stokes equations. This method decouples
//! the velocity and pressure, avoiding the ill-conditioned saddle-point system.
//!
//! # Algorithm
//!
//! For each time step:
//! 1. **Momentum Prediction**: Solve for intermediate velocity without pressure:
//!    ```text
//!    μ∇²u* - u·∇u* = f  (with BCs)
//!    ```
//!
//! 2. **Pressure Poisson**: Enforce incompressibility by solving for pressure:
//!    ```text
//!    ∇²p = ρ∇·u* / Δt
//!    ```
//!
//! 3. **Velocity Correction**: Update velocity to be divergence-free:
//!    ```text
//!    u^(n+1) = u* - Δt/ρ ∇p
//!    ```
//!
//! # References
//!
//! 1. Chorin, A. J. (1968): "Numerical solution of the Navier-Stokes equations"
//! 2. Temam, R. (1969): "Sur l'approximation de la solution des équations de Navier-Stokes"
//! 3. Guermond, Minev & Shen (2006): "An overview of projection methods"

use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_math::linear_solver::{ConjugateGradient, GMRES, IdentityPreconditioner, IterativeLinearSolver, LinearSolver};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField, Vector3};
use num_traits::{Float, FromPrimitive};
use std::collections::HashSet;

use crate::fem::{FemConfig, StokesFlowProblem, StokesFlowSolution};
use crate::fem::shape_functions::LagrangeTet10;
use crate::fem::quadrature::TetrahedronQuadrature;
use cfd_mesh::IndexedMesh;
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use cfd_mesh::domain::topology::Cell;
use crate::fem::solver::extract_vertex_indices;

/// Pressure projection solver for incompressible Stokes/Navier-Stokes equations
pub struct ProjectionSolver<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + Float + std::fmt::Debug> {
    config: FemConfig<T>,
    /// Time step for transient simulations
    dt: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy + Float + std::fmt::Debug + From<f64>> ProjectionSolver<T> {
    /// Create new projection solver with default time step
    pub fn new(config: FemConfig<T>) -> Self {
        Self { 
            config,
            dt: <T as FromPrimitive>::from_f64(0.001).unwrap_or_else(T::one),
        }
    }

    /// Create projection solver with specified time step
    pub fn with_timestep(config: FemConfig<T>, dt: T) -> Self {
        Self { config, dt }
    }

    /// Solve using projection method
    pub fn solve(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<StokesFlowSolution<T>> {
        println!("  Starting Chorin's projection method (dt={:?})...", self.dt);

        let n_nodes = problem.mesh.vertex_count();
        let n_corner_nodes = problem.n_corner_nodes;
        let n_velocity_dof = n_nodes * 3;

        // Step 1: Momentum prediction (solve for u*)
        println!("  Step 1: Momentum prediction...");
        let (momentum_matrix, momentum_rhs) = self.assemble_momentum_system(problem)?;
        
        let mut u_star = if let Some(sol) = previous_solution {
            sol.velocity.clone()
        } else {
            DVector::zeros(n_velocity_dof)
        };

        // Use GMRES for momentum (non-symmetric due to convection)
        let mut gmres_config = cfd_math::linear_solver::IterativeSolverConfig::default();
        gmres_config.max_iterations = 10000;
        gmres_config.tolerance = <T as FromPrimitive>::from_f64(1e-10).unwrap_or_else(T::zero);
        
        let momentum_solver = GMRES::new(gmres_config.clone(), 100);
        let monitor = momentum_solver.solve(&momentum_matrix, &momentum_rhs, &mut u_star, None::<&IdentityPreconditioner>)
            .map_err(|e| Error::Solver(format!("Momentum solve failed: {}", e)))?;
        
        println!("    ✓ Momentum converged in {} iterations (resid={:?})", 
            monitor.iteration, monitor.residual_history.last());

        // Step 2: Pressure Poisson equation
        println!("  Step 2: Pressure Poisson equation...");
        let (pressure_matrix, pressure_rhs) = self.assemble_pressure_poisson(
            problem,
            &u_star,
        )?;
        
        let mut pressure = if let Some(sol) = previous_solution {
            sol.pressure.clone()
        } else {
            DVector::zeros(n_corner_nodes)
        };

        // Use Conjugate Gradient for pressure (symmetric positive definite after pinning)
        let cg_solver = ConjugateGradient::new(gmres_config);
        let monitor = cg_solver.solve(&pressure_matrix, &pressure_rhs, &mut pressure, None::<&IdentityPreconditioner>)
            .map_err(|e| Error::Solver(format!("Pressure solve failed: {}", e)))?;
        
        println!("    ✓ Pressure converged in {} iterations (resid={:?})", 
            monitor.iteration, monitor.residual_history.last());

        // Step 3: Velocity correction (project onto divergence-free space)
        println!("  Step 3: Velocity correction...");
        let velocity = self.correct_velocity(problem, &u_star, &pressure)?;
        
        // Compute divergence for verification
        let max_div = self.compute_max_divergence(problem, &velocity)?;
        println!("    ✓ Max divergence after correction: {:?}", max_div);
        
        println!("  ✓ Projection method complete");

        Ok(StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes))
    }

    /// Assemble momentum system: (ρ/Δt)u* + μ∇²u - ρ(u·∇)u* = f + (ρ/Δt)u^n (without pressure gradient)
    fn assemble_momentum_system(
        &self,
        problem: &StokesFlowProblem<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_velocity_dof = n_nodes * 3;

        let mut builder = SparseMatrixBuilder::new(n_velocity_dof, n_velocity_dof);
        let mut rhs = DVector::zeros(n_velocity_dof);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices.iter()
            .map(|v| v.1.position.coords).collect();

        println!("    Assembling momentum matrix ({} elements)...", problem.mesh.cells.len());

        // Assemble viscous + transient terms (no pressure gradient)
        for (i, cell) in problem.mesh.cells.iter().enumerate() {
            let viscosity = problem.element_viscosities.as_ref().map_or(problem.fluid.viscosity, |v| v[i]);
            let idxs = extract_vertex_indices(cell, &problem.mesh, problem.n_corner_nodes)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            self.assemble_element_momentum(
                &mut builder, 
                &mut rhs, 
                &idxs, 
                &positions, 
                viscosity,
                problem.fluid.density,
                n_nodes
            )?;
        }

        // Apply boundary conditions
        let matrix = builder.build_with_rhs(&mut rhs)?;
        let (matrix, rhs) = self.apply_velocity_boundary_conditions(matrix, rhs, problem)?;

        Ok((matrix, rhs))
    }

    /// Assemble pressure Poisson equation: ∇²p = (ρ/Δt)∇·u*
    fn assemble_pressure_poisson(
        &self,
        problem: &StokesFlowProblem<T>,
        u_star: &DVector<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_corner_nodes = problem.n_corner_nodes;
        let mut builder = SparseMatrixBuilder::new(n_corner_nodes, n_corner_nodes);
        let mut rhs = DVector::zeros(n_corner_nodes);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices.iter()
            .map(|v| v.1.position.coords).collect();

        println!("    Assembling pressure Poisson matrix...");

        // Assemble Laplacian for pressure using P1 elements
        for cell in problem.mesh.cells.iter() {
            let idxs = extract_vertex_indices(cell, &problem.mesh, problem.n_corner_nodes)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            // Pressure uses only corner nodes (P1 elements)
            let corner_idxs: Vec<usize> = idxs.iter().take(4).copied().collect();
            
            self.assemble_element_pressure_laplacian(
                &mut builder,
                &mut rhs,
                &corner_idxs,
                &positions,
                problem.fluid.density,
                u_star,
            )?;
        }

        let matrix = builder.build_with_rhs(&mut rhs)?;
        
        // Pin one pressure DOF to remove null space
        let (matrix, rhs) = self.pin_pressure_reference(matrix, rhs)?;

        Ok((matrix, rhs))
    }

    /// Correct velocity to enforce incompressibility: u = u* - (Δt/ρ) ∇p
    fn correct_velocity(
        &self,
        problem: &StokesFlowProblem<T>,
        u_star: &DVector<T>,
        pressure: &DVector<T>,
    ) -> Result<DVector<T>> {
        let n_nodes = problem.mesh.vertex_count();
        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices.iter()
            .map(|v| v.1.position.coords).collect();

        let mut velocity = u_star.clone();
        let dt_over_rho = self.dt / problem.fluid.density;

        // For each element, compute pressure gradient and correct velocity
        for cell in problem.mesh.cells.iter() {
            let idxs = extract_vertex_indices(cell, &problem.mesh, problem.n_corner_nodes)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            let corner_idxs: Vec<usize> = idxs.iter().take(4).copied().collect();
            
            // Get pressure gradient at element center (constant for P1 elements)
            let pressure_grad = self.compute_pressure_gradient(&corner_idxs, &positions, pressure)?;
            
            // Correct velocities at all nodes of the element
            for &node_idx in &idxs {
                for d in 0..3 {
                    let vel_idx = node_idx * 3 + d;
                    velocity[vel_idx] = velocity[vel_idx] - dt_over_rho * pressure_grad[d];
                }
            }
        }

        // Re-apply velocity boundary conditions after correction
        self.apply_velocity_correction_bcs(problem, &mut velocity)?;

        Ok(velocity)
    }

    /// Assemble element contribution to momentum matrix (viscous + mass terms)
    /// 
    /// For the momentum prediction step:
    /// (ρ/Δt) ∫ N_i N_j dV + μ ∫ ∇N_i · ∇N_j dV
    fn assemble_element_momentum(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        indices: &[usize],
        positions: &[Vector3<T>],
        viscosity: T,
        density: T,
        n_nodes: usize,
    ) -> Result<()> {
        // Use quadrature for accurate integration
        let quad = TetrahedronQuadrature::keast_degree_3();
        
        // Compute element Jacobian transformation
        let v0 = positions[0];
        let v1 = positions[1];
        let v2 = positions[2];
        let v3 = positions[3];

        let j_mat = nalgebra::Matrix3::from_columns(&[
            v1 - v0,
            v2 - v0,
            v3 - v0,
        ]);
        
        let det_j = j_mat.determinant();
        let abs_det = Float::abs(det_j);
        
        if abs_det < <T as FromPrimitive>::from_f64(1e-20).unwrap() {
            return Err(Error::Solver("Near-zero element volume in momentum assembly".to_string()));
        }

        let j_inv_t = j_mat
            .try_inverse()
            .ok_or_else(|| Error::Solver("Singular Jacobian in momentum assembly".to_string()))?
            .transpose();

        // P1 gradients in reference space
        let grad_ref_p1 = nalgebra::Matrix3x4::new(
            -T::one(), T::one(), T::zero(), T::zero(),
            -T::one(), T::zero(), T::one(), T::zero(),
            -T::one(), T::zero(), T::zero(), T::one(),
        );
        
        // Transform to physical space
        let p1_gradients_phys = j_inv_t * grad_ref_p1;
        let shape = LagrangeTet10::new(p1_gradients_phys);

        let v_offset = n_nodes;
        let mass_coeff = density / self.dt;

        // Integrate using quadrature
        for (qp, &qw) in quad.points().iter().zip(quad.weights().iter()) {
            let weight = qw * abs_det;
            let l = [T::one() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
            
            // P2 shape function values and gradients
            let n_p2 = shape.values(&l);
            let grad_p2 = shape.gradients(&l);

            // Assemble for each velocity component
            for d in 0..3 {
                for i in 0..indices.len().min(10) {
                    let gi = indices[i];
                    let dof_i = gi + d * v_offset;
                    let grad_i = grad_p2.column(i);

                    for j in 0..indices.len().min(10) {
                        let gj = indices[j];
                        let dof_j = gj + d * v_offset;
                        let grad_j = grad_p2.column(j);

                        // Mass matrix term: (ρ/Δt) N_i N_j
                        let mass_term = mass_coeff * n_p2[i] * n_p2[j] * weight;
                        
                        // Viscous term: μ ∇N_i · ∇N_j
                        let visc_term = viscosity * grad_i.dot(&grad_j) * weight;

                        builder.add_entry(dof_i, dof_j, mass_term + visc_term)?;
                    }
                }
            }
        }

        Ok(())
    }

    /// Assemble element Laplacian for pressure Poisson equation
    /// 
    /// ∫ ∇N_i · ∇N_j dV for pressure (P1 elements)
    /// RHS: (ρ/Δt) ∫ N_i ∇·u* dV
    fn assemble_element_pressure_laplacian(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        corner_indices: &[usize],
        positions: &[Vector3<T>],
        density: T,
        u_star: &DVector<T>,
    ) -> Result<()> {
        // Compute element volume
        let v0 = positions[0];
        let v1 = positions[1];
        let v2 = positions[2];
        let v3 = positions[3];

        let j_mat = nalgebra::Matrix3::from_columns(&[
            v1 - v0,
            v2 - v0,
            v3 - v0,
        ]);
        
        let det_j = j_mat.determinant();
        let abs_det = Float::abs(det_j);
        
        if abs_det < <T as FromPrimitive>::from_f64(1e-20).unwrap() {
            return Err(Error::Solver("Near-zero element volume in pressure assembly".to_string()));
        }

        let j_inv_t = j_mat
            .try_inverse()
            .ok_or_else(|| Error::Solver("Singular Jacobian in pressure assembly".to_string()))?
            .transpose();

        // P1 gradients in reference space
        let grad_ref_p1 = nalgebra::Matrix3x4::new(
            -T::one(), T::one(), T::zero(), T::zero(),
            -T::one(), T::zero(), T::one(), T::zero(),
            -T::one(), T::zero(), T::zero(), T::one(),
        );
        
        // Transform to physical space (constant for P1)
        let grad_p1_phys = j_inv_t * grad_ref_p1;

        // For P1 elements, the Laplacian matrix is constant over the element
        // K_ij = ∫ ∇N_i · ∇N_j dV = (∇N_i · ∇N_j) * V/4
        // Using exact integration for linear elements
        let vol = abs_det / <T as FromPrimitive>::from_f64(6.0).unwrap();
        
        for i in 0..4 {
            let grad_i = grad_p1_phys.column(i);
            for j in 0..4 {
                let grad_j = grad_p1_phys.column(j);
                let laplacian_term = grad_i.dot(&grad_j) * vol;
                builder.add_entry(corner_indices[i], corner_indices[j], laplacian_term)?;
            }
        }

        // RHS: (ρ/Δt) ∫ N_i ∇·u* dV
        // For P1, N_i is linear, so we use quadrature
        let quad = TetrahedronQuadrature::keast_degree_3();
        let rho_over_dt = density / self.dt;
        
        for (qp, &qw) in quad.points().iter().zip(quad.weights().iter()) {
            let weight = qw * abs_det;
            let l = [T::one() - qp.x - qp.y - qp.z, qp.x, qp.y, qp.z];
            
            // Compute divergence of u* at this quadrature point
            let div_u = self.compute_divergence_at_quad_point(
                corner_indices, 
                &grad_p1_phys, 
                u_star
            );
            
            for i in 0..4 {
                rhs[corner_indices[i]] += rho_over_dt * l[i] * div_u * weight;
            }
        }

        Ok(())
    }

    /// Compute divergence at a quadrature point using velocity gradient
    fn compute_divergence_at_quad_point(
        &self,
        corner_indices: &[usize],
        grad_p1_phys: &nalgebra::Matrix3x4<T>,
        velocity: &DVector<T>,
    ) -> T {
        let mut div = T::zero();
        
        // ∇·u = Σ (∂u_x/∂x + ∂u_y/∂y + ∂u_z/∂z)
        // For P1 interpolation: u = Σ N_i u_i
        // ∇u = Σ u_i ⊗ ∇N_i
        // ∇·u = Σ u_i · ∇N_i (sum over components)
        
        for (i, &node_idx) in corner_indices.iter().enumerate() {
            let grad_N = grad_p1_phys.column(i);
            
            // Get velocity at this node
            let u_x = velocity[node_idx * 3];
            let u_y = velocity[node_idx * 3 + 1];
            let u_z = velocity[node_idx * 3 + 2];
            
            // Contribution to divergence: u · ∇N = u_x * ∂N/∂x + u_y * ∂N/∂y + u_z * ∂N/∂z
            // But we need ∂u_x/∂x + ∂u_y/∂y + ∂u_z/∂z
            // = Σ u_x_i * ∂N_i/∂x + Σ u_y_i * ∂N_i/∂y + Σ u_z_i * ∂N_i/∂z
            div += u_x * grad_N.x + u_y * grad_N.y + u_z * grad_N.z;
        }
        
        div
    }

    /// Compute pressure gradient from nodal pressures (P1 elements)
    /// 
    /// For P1 elements, pressure gradient is constant within element:
    /// ∇p = Σ p_i ∇N_i
    fn compute_pressure_gradient(
        &self,
        corner_indices: &[usize],
        positions: &[Vector3<T>],
        pressure: &DVector<T>,
    ) -> Result<Vector3<T>> {
        let v0 = positions[0];
        let v1 = positions[1];
        let v2 = positions[2];
        let v3 = positions[3];

        let j_mat = nalgebra::Matrix3::from_columns(&[
            v1 - v0,
            v2 - v0,
            v3 - v0,
        ]);
        
        let j_inv_t = j_mat
            .try_inverse()
            .ok_or_else(|| Error::Solver("Singular Jacobian in pressure gradient".to_string()))?
            .transpose();

        // P1 gradients in reference space
        let grad_ref_p1 = nalgebra::Matrix3x4::new(
            -T::one(), T::one(), T::zero(), T::zero(),
            -T::one(), T::zero(), T::one(), T::zero(),
            -T::one(), T::zero(), T::zero(), T::one(),
        );
        
        let grad_p1_phys = j_inv_t * grad_ref_p1;

        // ∇p = Σ p_i ∇N_i
        let mut grad_p = Vector3::zeros();
        
        for (i, &node_idx) in corner_indices.iter().enumerate() {
            if node_idx < pressure.len() {
                let p_i = pressure[node_idx];
                let grad_N = grad_p1_phys.column(i);
                grad_p += grad_N * p_i;
            }
        }

        Ok(grad_p)
    }

    /// Apply velocity boundary conditions to the momentum system
    fn apply_velocity_boundary_conditions(
        &self,
        matrix: SparseMatrix<T>,
        mut rhs: DVector<T>,
        problem: &StokesFlowProblem<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let v_offset = n_nodes;
        
        // Compute mesh scale for diagonal scaling
        let mesh_scale = compute_mesh_scale(&problem.mesh);
        let diag_scale = problem.fluid.viscosity * mesh_scale;

        let mut builder = csr_to_builder(matrix);
        let mut applied_bcs = HashSet::new();

        for (&node_idx, bc) in &problem.boundary_conditions {
            match bc {
                BoundaryCondition::VelocityInlet { velocity } => {
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        if !applied_bcs.contains(&dof) {
                            builder.set_dirichlet_row(dof, diag_scale, velocity[d]);
                            rhs[dof] = velocity[d] * diag_scale;
                            applied_bcs.insert(dof);
                        }
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    for d in 0..3 {
                        let dof = node_idx + d * v_offset;
                        if !applied_bcs.contains(&dof) {
                            builder.set_dirichlet_row(dof, diag_scale, T::zero());
                            rhs[dof] = T::zero();
                            applied_bcs.insert(dof);
                        }
                    }
                }
                BoundaryCondition::Dirichlet { value, component_values } => {
                    if let Some(comps) = component_values {
                        for d in 0..3 {
                            if let Some(Some(val)) = comps.get(d) {
                                let dof = node_idx + d * v_offset;
                                if !applied_bcs.contains(&dof) {
                                    builder.set_dirichlet_row(dof, diag_scale, *val);
                                    rhs[dof] = *val * diag_scale;
                                    applied_bcs.insert(dof);
                                }
                            }
                        }
                    } else {
                        for d in 0..3 {
                            let dof = node_idx + d * v_offset;
                            if !applied_bcs.contains(&dof) {
                                builder.set_dirichlet_row(dof, diag_scale, *value);
                                rhs[dof] = *value * diag_scale;
                                applied_bcs.insert(dof);
                            }
                        }
                    }
                }
                _ => {
                    // PressureOutlet, PressureInlet, etc. are handled in pressure solve
                }
            }
        }

        let matrix = builder.build_with_rhs(&mut rhs)?;
        Ok((matrix, rhs))
    }

    /// Apply velocity boundary conditions after pressure correction
    fn apply_velocity_correction_bcs(
        &self,
        problem: &StokesFlowProblem<T>,
        velocity: &mut DVector<T>,
    ) -> Result<()> {
        for (&node_idx, bc) in &problem.boundary_conditions {
            match bc {
                BoundaryCondition::VelocityInlet { velocity: inlet_vel } => {
                    for d in 0..3 {
                        velocity[node_idx * 3 + d] = inlet_vel[d];
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    for d in 0..3 {
                        velocity[node_idx * 3 + d] = T::zero();
                    }
                }
                BoundaryCondition::Dirichlet { value, component_values } => {
                    if let Some(comps) = component_values {
                        for d in 0..3 {
                            if let Some(Some(val)) = comps.get(d) {
                                velocity[node_idx * 3 + d] = *val;
                            }
                        }
                    } else {
                        for d in 0..3 {
                            velocity[node_idx * 3 + d] = *value;
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Pin reference pressure to remove null space
    fn pin_pressure_reference(
        &self,
        matrix: SparseMatrix<T>,
        mut rhs: DVector<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_dof = rhs.len();
        if n_dof == 0 {
            return Ok((matrix, rhs));
        }

        // Compute diagonal scale for pressure DOF
        let diag_scale = T::one();
        
        let mut builder = csr_to_builder(matrix);
        
        // Pin first pressure DOF to zero
        builder.set_dirichlet_row(0, diag_scale, T::zero());
        rhs[0] = T::zero();
        
        println!("    Pinned pressure DOF 0 to zero (reference pressure)");

        let matrix = builder.build_with_rhs(&mut rhs)?;
        Ok((matrix, rhs))
    }

    /// Compute maximum divergence over all elements for verification
    fn compute_max_divergence(
        &self,
        problem: &StokesFlowProblem<T>,
        velocity: &DVector<T>,
    ) -> Result<T> {
        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices.iter()
            .map(|v| v.1.position.coords)
            .collect();

        let mut max_div = T::zero();

        for cell in problem.mesh.cells.iter() {
            let idxs = extract_vertex_indices(cell, &problem.mesh, problem.n_corner_nodes)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            let corner_idxs: Vec<usize> = idxs.iter().take(4).copied().collect();

            // Compute element Jacobian
            let v0 = positions[0];
            let v1 = positions[1];
            let v2 = positions[2];
            let v3 = positions[3];

            let j_mat = nalgebra::Matrix3::from_columns(&[
                v1 - v0,
                v2 - v0,
                v3 - v0,
            ]);
            
            let j_inv_t = match j_mat.try_inverse() {
                Some(inv) => inv.transpose(),
                None => continue,
            };

            // P1 gradients
            let grad_ref_p1 = nalgebra::Matrix3x4::new(
                -T::one(), T::one(), T::zero(), T::zero(),
                -T::one(), T::zero(), T::one(), T::zero(),
                -T::one(), T::zero(), T::zero(), T::one(),
            );
            let grad_p1_phys = j_inv_t * grad_ref_p1;

            // Compute divergence at element center
            let div = self.compute_divergence_at_quad_point(&corner_idxs, &grad_p1_phys, velocity);
            
            if Float::abs(div) > max_div {
                max_div = Float::abs(div);
            }
        }

        Ok(max_div)
    }
}

/// Convert a CsrMatrix into a SparseMatrixBuilder so that Dirichlet BCs can be applied.
fn csr_to_builder<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(matrix: SparseMatrix<T>) -> SparseMatrixBuilder<T> {
    let nrows = matrix.nrows();
    let ncols = matrix.ncols();
    let nnz = matrix.nnz();
    let mut builder = SparseMatrixBuilder::with_capacity(nrows, ncols, nnz);
    let offsets = matrix.row_offsets();
    let col_indices = matrix.col_indices();
    let values = matrix.values();
    for row in 0..nrows {
        let start = offsets[row];
        let end = offsets[row + 1];
        for idx in start..end {
            let _ = builder.add_entry(row, col_indices[idx], values[idx]);
        }
    }
    builder
}



/// Compute mesh scale for diagonal scaling
fn compute_mesh_scale<T: cfd_mesh::domain::core::Scalar + RealField + Copy + Float>(mesh: &IndexedMesh<T>) -> T {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_projection_solver_creation() {
        let config = FemConfig::default();
        let _solver = ProjectionSolver::<f64>::new(config);
    }

    #[test]
    fn test_projection_solver_with_timestep() {
        let config = FemConfig::default();
        let solver = ProjectionSolver::with_timestep(config, 0.01);
        assert!(solver.dt > 0.0);
    }
}
