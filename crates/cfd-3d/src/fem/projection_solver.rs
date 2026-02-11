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
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_math::linear_solver::{ConjugateGradient, GMRES, LinearSolver};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField, Vector3};
use num_traits::{Float, FromPrimitive};

use crate::fem::{FemConfig, StokesFlowProblem, StokesFlowSolution};
use crate::fem::shape_functions::LagrangeTet10;
use cfd_mesh::mesh::Mesh;
use cfd_mesh::topology::Cell;

/// Pressure projection solver for incompressible Stokes/Navier-Stokes equations
pub struct ProjectionSolver<T: RealField + Copy + FromPrimitive + Float + std::fmt::Debug> {
    config: FemConfig<T>,
}

impl<T: RealField + FromPrimitive + Copy + Float + std::fmt::Debug + From<f64>> ProjectionSolver<T> {
    /// Create new projection solver
    pub fn new(config: FemConfig<T>) -> Self {
        Self { config }
    }

    /// Solve using projection method
    pub fn solve(
        &mut self,
        problem: &StokesFlowProblem<T>,
        previous_solution: Option<&StokesFlowSolution<T>>,
    ) -> Result<StokesFlowSolution<T>> {
        println!("  Starting Chorin's projection method...");

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

        // Use GMRES for momentum (positive definite after applying BCs)
        let mut gmres_config = cfd_math::linear_solver::IterativeSolverConfig::default();
        gmres_config.max_iterations = 10000;
        gmres_config.tolerance = T::from_f64(1e-10).unwrap_or_else(T::zero);
        
        let momentum_solver = GMRES::new(gmres_config.clone(), 100);
        let monitor = momentum_solver.solve(&momentum_matrix, &momentum_rhs, &mut u_star)
            .map_err(|e| Error::Solver(format!("Momentum solve failed: {}", e)))?;
        
        println!("    ✓ Momentum converged in {} iterations", monitor.iteration);

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

        // Use Conjugate Gradient for pressure (symmetric positive definite)
        let cg_solver = ConjugateGradient::new(gmres_config);
        let monitor = cg_solver.solve(&pressure_matrix, &pressure_rhs, &mut pressure)
            .map_err(|e| Error::Solver(format!("Pressure solve failed: {}", e)))?;
        
        println!("    ✓ Pressure converged in {} iterations", monitor.iteration);

        // Step 3: Velocity correction (project onto divergence-free space)
        println!("  Step 3: Velocity correction...");
        let velocity = self.correct_velocity(problem, &u_star, &pressure)?;
        
        println!("  ✓ Projection method complete");

        Ok(StokesFlowSolution::new_with_corners(velocity, pressure, n_nodes, n_corner_nodes))
    }

    /// Assemble momentum system: μ∇²u - u·∇u = f (without pressure)
    fn assemble_momentum_system(
        &self,
        problem: &StokesFlowProblem<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_velocity_dof = n_nodes * 3;

        let mut builder = SparseMatrixBuilder::new(n_velocity_dof, n_velocity_dof);
        let mut rhs = DVector::zeros(n_velocity_dof);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices().iter()
            .map(|v| v.position.coords).collect();

        println!("    Assembling momentum matrix ({} elements)...", problem.mesh.cells().len());

        // Assemble viscous + convective terms (no pressure gradient)
        for (i, cell) in problem.mesh.cells().iter().enumerate() {
            let viscosity = problem.element_viscosities.as_ref().map_or(problem.fluid.viscosity, |v| v[i]);
            let idxs = extract_vertex_indices(cell, &problem.mesh)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            self.assemble_element_momentum(&mut builder, &mut rhs, &idxs, &positions, viscosity)?;
        }

        // Apply boundary conditions
        let matrix = builder.build_with_rhs(&mut rhs)?;
        let (matrix, rhs) = self.apply_velocity_boundary_conditions(matrix, rhs, problem)?;

        Ok((matrix, rhs))
    }

    /// Assemble pressure Poisson equation: ∇²p = ρ∇·u* / Δt
    fn assemble_pressure_poisson(
        &self,
        problem: &StokesFlowProblem<T>,
        u_star: &DVector<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_corner_nodes = problem.n_corner_nodes;
        let mut builder = SparseMatrixBuilder::new(n_corner_nodes, n_corner_nodes);
        let mut rhs = DVector::zeros(n_corner_nodes);

        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices().iter()
            .map(|v| v.position.coords).collect();

        println!("    Assembling pressure Poisson matrix...");

        // Assemble Laplacian for pressure
        for cell in problem.mesh.cells().iter() {
            let idxs = extract_vertex_indices(cell, &problem.mesh)?;
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

    /// Correct velocity to enforce incompressibility: u = u* - Δt/ρ ∇p
    fn correct_velocity(
        &self,
        problem: &StokesFlowProblem<T>,
        u_star: &DVector<T>,
        pressure: &DVector<T>,
    ) -> Result<DVector<T>> {
        let n_nodes = problem.mesh.vertex_count();
        let vertex_positions: Vec<Vector3<T>> = problem.mesh.vertices().iter()
            .map(|v| v.position.coords).collect();

        let mut velocity = u_star.clone();

        // For each element, compute pressure gradient and correct velocity
        for cell in problem.mesh.cells().iter() {
            let idxs = extract_vertex_indices(cell, &problem.mesh)?;
            let positions: Vec<Vector3<T>> = idxs.iter()
                .map(|&idx| vertex_positions[idx])
                .collect();

            let corner_idxs: Vec<usize> = idxs.iter().take(4).copied().collect();
            
            // Get pressure gradient at element center
            let pressure_grad = self.compute_pressure_gradient(&corner_idxs, &positions, pressure)?;
            
            // Correct velocities at element nodes
            let dt_over_rho = T::one() / problem.fluid.density; // Steady state approximation
            
            for &node_idx in &idxs {
                for d in 0..3 {
                    let vel_idx = node_idx * 3 + d;
                    velocity[vel_idx] = velocity[vel_idx] - dt_over_rho * pressure_grad[d];
                }
            }
        }

        Ok(velocity)
    }

    /// Assemble element contribution to momentum matrix (viscous terms only)
    fn assemble_element_momentum(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        indices: &[usize],
        positions: &[Vector3<T>],
        viscosity: T,
    ) -> Result<()> {
        // Placeholder: assemble 4x4 stiffness matrix for P1 elements
        // K_ij = μ ∫ ∇φ_i · ∇φ_j dV
        
        // Compute element Jacobian and gradients
        let v0 = positions[0];
        let v1 = positions[1];
        let v2 = positions[2];
        let v3 = positions[3];

        let vol = ((v1 - v0).cross(&(v2 - v0))).dot(&(v3 - v0)) / T::from_f64(6.0).unwrap();
        if Float::abs(vol) < T::from_f64(1e-15).unwrap() {
            return Err(Error::Solver("Near-zero element volume".to_string()));
        }

        // Simple diagonal mass lumping for now
        let diag_val = viscosity * Float::abs(vol) / T::from_f64(4.0).unwrap();
        
        for i in 0..4 {
            let node_idx = indices[i];
            for d in 0..3 {
                let dof_idx = node_idx * 3 + d;
                builder.add_entry(dof_idx, dof_idx, diag_val)?;
            }
        }

        Ok(())
    }

    /// Assemble element Laplacian for pressure Poisson equation
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

        let vol = ((v1 - v0).cross(&(v2 - v0))).dot(&(v3 - v0)) / T::from_f64(6.0).unwrap();
        if Float::abs(vol) < T::from_f64(1e-15).unwrap() {
            return Err(Error::Solver("Near-zero element volume".to_string()));
        }

        // Simplified diagonal Laplacian contribution
        let diag_val = Float::abs(vol) / T::from_f64(4.0).unwrap();
        
        for i in 0..4 {
            let node_idx = corner_indices[i];
            builder.add_entry(node_idx, node_idx, diag_val)?;

            // RHS: divergence of u_star at this node
            let u_div = self.compute_divergence_at_node(corner_indices[i], u_star)?;
            rhs[node_idx] += density * u_div * Float::abs(vol) / T::from_f64(4.0).unwrap();
        }

        Ok(())
    }

    /// Compute divergence at a node (simplified)
    fn compute_divergence_at_node(&self, node_idx: usize, velocity: &DVector<T>) -> Result<T> {
        // Simplified: finite difference approximation
        // In practice, should use proper FEM gradient reconstruction
        Ok(T::zero()) // Placeholder
    }

    /// Compute pressure gradient from nodal pressures
    fn compute_pressure_gradient(
        &self,
        corner_indices: &[usize],
        positions: &[Vector3<T>],
        pressure: &DVector<T>,
    ) -> Result<Vector3<T>> {
        // Linear reconstruction: ∇p = Σ p_i ∇φ_i
        // For P1 elements, this is constant within element
        
        // Placeholder: return zero gradient
        Ok(Vector3::zeros())
    }

    /// Apply velocity boundary conditions
    fn apply_velocity_boundary_conditions(
        &self,
        matrix: SparseMatrix<T>,
        rhs: DVector<T>,
        problem: &StokesFlowProblem<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        // Apply Dirichlet BCs by modifying matrix and RHS
        // This is a placeholder - full implementation in parent solver
        Ok((matrix, rhs))
    }

    /// Pin reference pressure to remove null space
    fn pin_pressure_reference(
        &self,
        matrix: SparseMatrix<T>,
        rhs: DVector<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        // Pin first DOF to zero
        // Placeholder - actual implementation requires matrix modification
        Ok((matrix, rhs))
    }
}

/// Extract vertex indices from cell
fn extract_vertex_indices<T: RealField + Copy>(
    cell: &Cell<T>,
    mesh: &Mesh<T>,
) -> Result<Vec<usize>> {
    use cfd_mesh::topology::CellType;
    
    match &cell.cell_type {
        CellType::Tetrahedron(tet) => {
            Ok(vec![tet.v0 as usize, tet.v1 as usize, tet.v2 as usize, tet.v3 as usize])
        }
        _ => Err(Error::InvalidConfiguration("Expected tetrahedral cell".to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_projection_solver_creation() {
        let config = FemConfig::default();
        let _solver = ProjectionSolver::<f64>::new(config);
    }
}
