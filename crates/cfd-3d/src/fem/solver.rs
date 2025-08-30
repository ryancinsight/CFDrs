//! FEM solver implementation

use cfd_core::boundary::BoundaryCondition;
use cfd_core::error::Result;
use cfd_math::linear_solver::{ConjugateGradient, LinearSolver};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use nalgebra::{DVector, RealField, Vector3};
use num_traits::{Float, FromPrimitive};

use crate::fem::constants;
use crate::fem::{ElementMatrices, FemConfig, FluidElement, StokesFlowProblem, StokesFlowSolution};
use cfd_mesh::mesh::Mesh;
use cfd_mesh::topology::Cell;

/// Finite Element Method solver for 3D incompressible flow
pub struct FemSolver<T: RealField + Copy> {
    /// Configuration (stored for future use)
    _config: FemConfig<T>,
    /// Linear solver for the system
    linear_solver: Box<dyn LinearSolver<T>>,
}

/// Extract vertex indices from a cell
fn extract_vertex_indices<T: RealField + Copy>(cell: &Cell, mesh: &Mesh<T>) -> Vec<usize> {
    // For tetrahedral elements, extract 4 unique vertex indices from faces
    let mut indices = Vec::with_capacity(4);
    let mut seen = std::collections::HashSet::new();

    for &face_idx in &cell.faces {
        if let Some(face) = mesh.face(face_idx) {
            for &vertex_idx in &face.vertices {
                if seen.insert(vertex_idx) && indices.len() < 4 {
                    indices.push(vertex_idx);
                }
            }
        }
        if indices.len() >= 4 {
            break;
        }
    }

    // Ensure we have exactly 4 indices for tetrahedral element
    while indices.len() < 4 {
        indices.push(0);
    }

    indices
}

impl<T: RealField + FromPrimitive + Copy + Float> FemSolver<T> {
    /// Create a new FEM solver
    pub fn new(config: FemConfig<T>) -> Self {
        let linear_solver: Box<dyn LinearSolver<T>> = Box::new(ConjugateGradient::new(
            cfd_math::linear_solver::IterativeSolverConfig::default(),
        ));

        Self {
            _config: config,
            linear_solver,
        }
    }

    /// Solve the Stokes flow problem
    pub fn solve(&mut self, problem: &StokesFlowProblem<T>) -> Result<StokesFlowSolution<T>> {
        // Validate problem setup
        problem.validate()?;

        let n_nodes = problem.mesh.vertex_count();
        let n_velocity_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pressure_dof = n_nodes;
        let n_total_dof = n_velocity_dof + n_pressure_dof;

        // Assemble global system
        let (matrix, rhs) = self.assemble_system(problem)?;

        // Validate system dimensions
        debug_assert_eq!(matrix.nrows(), n_total_dof, "Matrix row dimension mismatch");
        debug_assert_eq!(
            matrix.ncols(),
            n_total_dof,
            "Matrix column dimension mismatch"
        );
        debug_assert_eq!(rhs.len(), n_total_dof, "RHS vector dimension mismatch");

        // Solve linear system
        let solution = self.linear_solver.solve(&matrix, &rhs, None)?;

        // Extract velocity and pressure
        let velocity = solution.rows(0, n_velocity_dof).into();
        let pressure = solution.rows(n_velocity_dof, n_pressure_dof).into();

        Ok(StokesFlowSolution::new(velocity, pressure, n_nodes))
    }

    /// Assemble the global system matrix and RHS
    fn assemble_system(
        &self,
        problem: &StokesFlowProblem<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_nodes = problem.mesh.vertex_count();
        let n_velocity_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pressure_dof = n_nodes;
        let n_total_dof = n_velocity_dof + n_pressure_dof;

        let mut builder = SparseMatrixBuilder::new(n_total_dof, n_total_dof);
        let mut rhs = DVector::zeros(n_total_dof);

        // Get fluid properties
        let viscosity = problem.fluid.viscosity;

        // Loop over elements
        for (_elem_idx, cell) in problem.mesh.cells().iter().enumerate() {
            // Get vertex indices for this cell
            let vertex_indices = extract_vertex_indices(cell, &problem.mesh);

            // Create element
            let mut element = FluidElement::new(vertex_indices);

            // Calculate element properties
            // Convert vertices to Vector3 format
            let vertex_positions: Vec<Vector3<T>> = problem
                .mesh
                .vertices()
                .iter()
                .map(|v| v.position.coords)
                .collect();
            element.calculate_volume(&vertex_positions[..4]); // Use first 4 vertices for tetrahedral
            element.calculate_shape_derivatives(&vertex_positions[..4]);

            // Calculate element matrices
            let elem_matrices = self.calculate_element_matrices(&element, viscosity);

            // Assemble into global system
            self.assemble_element(&mut builder, &mut rhs, &element, &elem_matrices)?;
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(&mut builder, &mut rhs, problem)?;

        let matrix = builder.build()?;
        Ok((matrix, rhs))
    }

    /// Calculate element matrices
    fn calculate_element_matrices(
        &self,
        element: &FluidElement<T>,
        viscosity: T,
    ) -> ElementMatrices<T> {
        let n_nodes = element.nodes.len();
        let n_velocity_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pressure_dof = n_nodes;
        let n_dof = n_velocity_dof + n_pressure_dof;

        let mut matrices = ElementMatrices::new(n_dof);

        // Get element volume
        let volume = self.compute_element_volume(element);

        // For Stokes flow: -∇²u + ∇p = f, ∇·u = 0
        // We need K (viscous), B (divergence), and B^T (gradient) matrices

        // Viscous stiffness matrix (K) for velocity DOFs
        // Using simplified Laplacian operator
        let visc_factor =
            viscosity * volume / T::from_usize(n_nodes * n_nodes).unwrap_or_else(T::one);

        for i in 0..n_nodes {
            for j in 0..n_nodes {
                for d in 0..constants::VELOCITY_COMPONENTS {
                    let row = i * constants::VELOCITY_COMPONENTS + d;
                    let col = j * constants::VELOCITY_COMPONENTS + d;

                    if i == j {
                        // Diagonal dominance for stability
                        matrices.k_e[(row, col)] =
                            visc_factor * T::from_f64(2.0).unwrap_or_else(T::one);
                    } else {
                        // Off-diagonal coupling
                        matrices.k_e[(row, col)] =
                            -visc_factor / T::from_usize(n_nodes - 1).unwrap_or_else(T::one);
                    }
                }
            }
        }

        // Divergence matrix (B) and gradient matrix (B^T)
        // These couple velocity and pressure
        let div_factor = volume / T::from_usize(n_nodes).unwrap_or_else(T::one);

        for i in 0..n_nodes {
            for j in 0..n_nodes {
                if i == j {
                    // Gradient operator: pressure to velocity
                    for d in 0..constants::VELOCITY_COMPONENTS {
                        let vel_idx = i * constants::VELOCITY_COMPONENTS + d;
                        let pres_idx = n_velocity_dof + j;

                        // B^T: gradient (pressure affects velocity)
                        matrices.k_e[(vel_idx, pres_idx)] = div_factor;
                        // B: divergence (velocity affects pressure)
                        matrices.k_e[(pres_idx, vel_idx)] = div_factor;
                    }
                }
            }
        }

        // Add small stabilization to pressure block to avoid singular matrix
        let stab_factor = T::from_f64(1e-8).unwrap_or_else(T::zero) * volume;
        for i in 0..n_nodes {
            let pres_idx = n_velocity_dof + i;
            matrices.k_e[(pres_idx, pres_idx)] = matrices.k_e[(pres_idx, pres_idx)] + stab_factor;
        }

        matrices
    }

    fn compute_element_volume(&self, element: &FluidElement<T>) -> T {
        // Use pre-calculated volume if available, otherwise default
        if element.volume > T::zero() {
            element.volume
        } else {
            T::one()
        }
    }

    /// Assemble element contribution into global system
    fn assemble_element(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        _rhs: &mut DVector<T>, // Will be used when body forces are added
        element: &FluidElement<T>,
        matrices: &ElementMatrices<T>,
    ) -> Result<()> {
        // Map local DOFs to global DOFs
        let n_nodes = element.nodes.len();
        let dofs_per_node = constants::VELOCITY_COMPONENTS + 1;

        // Direct assembly - add element stiffness to global matrix
        for i in 0..n_nodes {
            for j in 0..n_nodes {
                for k in 0..dofs_per_node {
                    let local_i = i * dofs_per_node + k;
                    let local_j = j * dofs_per_node + k;
                    let global_i = element.nodes[i] * dofs_per_node + k;
                    let global_j = element.nodes[j] * dofs_per_node + k;

                    if local_i < matrices.k_e.nrows() && local_j < matrices.k_e.ncols() {
                        builder.add_entry(global_i, global_j, matrices.k_e[(local_i, local_j)])?;
                    }
                }
            }
        }
        Ok(())
    }

    /// Apply boundary conditions to the system
    fn apply_boundary_conditions(
        &self,
        builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        problem: &StokesFlowProblem<T>,
    ) -> Result<()> {
        // Apply Dirichlet boundary conditions using penalty method
        let penalty = T::from_f64(1e10).unwrap_or_else(T::one);

        for (node_idx, bc) in &problem.boundary_conditions {
            let dof = *node_idx * (constants::VELOCITY_COMPONENTS + 1);

            match bc {
                BoundaryCondition::VelocityInlet { velocity } => {
                    // Apply penalty method for velocity components
                    for i in 0..constants::VELOCITY_COMPONENTS {
                        let component_dof = dof + i;
                        builder.add_entry(component_dof, component_dof, penalty)?;
                        if component_dof < rhs.len() {
                            rhs[component_dof] = penalty * velocity[i];
                        }
                    }
                }
                BoundaryCondition::Wall { .. } => {
                    // No-slip wall: zero velocity
                    for i in 0..constants::VELOCITY_COMPONENTS {
                        let component_dof = dof + i;
                        builder.add_entry(component_dof, component_dof, penalty)?;
                        if component_dof < rhs.len() {
                            rhs[component_dof] = T::zero();
                        }
                    }
                }
                _ => { /* Other boundary conditions not implemented yet */ }
            }
        }

        Ok(())
    }
}
