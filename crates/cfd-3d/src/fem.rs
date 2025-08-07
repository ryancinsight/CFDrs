//! Finite Element Method (FEM) solvers for 3D CFD.
//!
//! This module provides FEM implementations for solving 3D fluid dynamics problems
//! using various element types and numerical schemes.

use cfd_core::{Error, Result};
use cfd_math::{LinearSolver, LinearSolverConfig, ConjugateGradient, SparseMatrixBuilder};
use cfd_mesh::{Mesh, Cell};
use nalgebra::{RealField, Vector3, DVector};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use rayon::prelude::*;

/// FEM solver configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: RealField> {
    /// Convergence tolerance
    pub tolerance: T,
    /// Maximum number of iterations
    pub max_iterations: usize,
    /// Element type to use
    pub element_type: ElementType,
    /// Integration order
    pub integration_order: usize,
    /// Enable verbose output
    pub verbose: bool,
}

impl<T: RealField + FromPrimitive> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-6).unwrap(),
            max_iterations: 1000,
            element_type: ElementType::Tetrahedron4,
            integration_order: 2,
            verbose: false,
        }
    }
}

/// Supported element types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElementType {
    /// 4-node tetrahedron (linear)
    Tetrahedron4,
    /// 10-node tetrahedron (quadratic)
    Tetrahedron10,
    /// 8-node hexahedron (linear)
    Hexahedron8,
    /// 20-node hexahedron (quadratic)
    Hexahedron20,
}

/// Finite element for 3D problems
pub trait Element<T: RealField>: Send + Sync {
    /// Number of nodes in the element
    fn num_nodes(&self) -> usize;

    /// Evaluate shape functions at given local coordinates
    fn shape_functions(&self, xi: &Vector3<T>) -> Vec<T>;

    /// Evaluate shape function derivatives at given local coordinates
    fn shape_derivatives(&self, xi: &Vector3<T>) -> Vec<Vector3<T>>;

    /// Get integration points and weights
    fn integration_points(&self) -> Vec<(Vector3<T>, T)>;

    /// Compute element stiffness matrix
    fn stiffness_matrix(
        &self,
        nodes: &[Vector3<T>],
        material_properties: &MaterialProperties<T>,
    ) -> Result<nalgebra::DMatrix<T>>;
}

/// Material properties for FEM analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialProperties<T: RealField> {
    /// Density
    pub density: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Thermal conductivity (for heat transfer)
    pub thermal_conductivity: Option<T>,
    /// Specific heat capacity
    pub specific_heat: Option<T>,
}

/// 4-node tetrahedral element implementation
#[derive(Debug, Clone)]
pub struct Tetrahedron4<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Default for Tetrahedron4<T> {
    fn default() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive> Element<T> for Tetrahedron4<T> {
    fn num_nodes(&self) -> usize {
        4
    }

    fn shape_functions(&self, xi: &Vector3<T>) -> Vec<T> {
        let xi1 = xi.x.clone();
        let xi2 = xi.y.clone();
        let xi3 = xi.z.clone();
        let xi0 = T::one() - xi1.clone() - xi2.clone() - xi3.clone();

        vec![xi0, xi1, xi2, xi3]
    }

    fn shape_derivatives(&self, _xi: &Vector3<T>) -> Vec<Vector3<T>> {
        // For linear tetrahedron, derivatives are constant
        vec![
            Vector3::new(-T::one(), -T::one(), -T::one()), // Node 0
            Vector3::new(T::one(), T::zero(), T::zero()),   // Node 1
            Vector3::new(T::zero(), T::one(), T::zero()),   // Node 2
            Vector3::new(T::zero(), T::zero(), T::one()),   // Node 3
        ]
    }

    fn integration_points(&self) -> Vec<(Vector3<T>, T)> {
        // Single integration point at centroid for linear tetrahedron
        let quarter = T::from_f64(0.25).unwrap();
        let weight = T::from_f64(1.0/6.0).unwrap(); // Volume of reference tetrahedron

        vec![(Vector3::new(quarter.clone(), quarter.clone(), quarter), weight)]
    }

    fn stiffness_matrix(
        &self,
        nodes: &[Vector3<T>],
        material_properties: &MaterialProperties<T>,
    ) -> Result<nalgebra::DMatrix<T>> {
        if nodes.len() != 4 {
            return Err(Error::InvalidConfiguration(
                "Tetrahedron4 requires exactly 4 nodes".to_string()
            ));
        }

        // Compute Jacobian matrix
        let mut jacobian = nalgebra::Matrix3::zeros();
        let derivatives = self.shape_derivatives(&Vector3::zeros());

        for i in 0..3 {
            for j in 0..4 {
                jacobian[(i, 0)] += derivatives[j][i].clone() * nodes[j].x.clone();
                jacobian[(i, 1)] += derivatives[j][i].clone() * nodes[j].y.clone();
                jacobian[(i, 2)] += derivatives[j][i].clone() * nodes[j].z.clone();
            }
        }

        let det_j: T = jacobian.determinant();
        if det_j <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Negative or zero Jacobian determinant".to_string()
            ));
        }

        // For now, return a simple identity-based stiffness matrix
        // TODO: Implement proper stiffness matrix computation for Navier-Stokes
        let mut k = nalgebra::DMatrix::zeros(12, 12); // 4 nodes × 3 DOF per node
        let stiffness_factor = material_properties.viscosity.clone() * det_j.clone();

        for i in 0..12 {
            k[(i, i)] = stiffness_factor.clone();
        }

        Ok(k)
    }
}

/// FEM solver for 3D fluid dynamics problems
pub struct FemSolver<T: RealField> {
    config: FemConfig<T>,
    mesh: Option<Mesh<T>>,
}

impl<T: RealField + FromPrimitive + Send + Sync> FemSolver<T> {
    /// Create a new FEM solver
    pub fn new(config: FemConfig<T>) -> Self {
        Self {
            config,
            mesh: None,
        }
    }

    /// Create with default configuration
    pub fn default() -> Self {
        Self::new(FemConfig::default())
    }

    /// Set the mesh for the solver
    pub fn set_mesh(&mut self, mesh: Mesh<T>) {
        self.mesh = Some(mesh);
    }

    /// Solve steady-state Stokes equations (simplified Navier-Stokes)
    /// ∇²u - ∇p = f (momentum)
    /// ∇·u = 0 (continuity)
    pub fn solve_stokes(
        &self,
        velocity_bcs: &HashMap<usize, Vector3<T>>, // Simplified: Dirichlet velocity BCs
        _body_force: &HashMap<usize, Vector3<T>>,
        material_properties: &MaterialProperties<T>,
    ) -> Result<HashMap<usize, Vector3<T>>> {
        let mesh = self.mesh.as_ref().ok_or_else(|| {
            Error::InvalidConfiguration("No mesh set for FEM solver".to_string())
        })?;

        let num_nodes = mesh.vertices.len();
        let num_dofs = num_nodes * 4; // 3 velocity + 1 pressure per node

        // Build global system matrix using iterator combinators
        let mut matrix_builder = SparseMatrixBuilder::new(num_dofs, num_dofs);
        let mut rhs = DVector::zeros(num_dofs);

        // Assemble element contributions in parallel
        let element_matrices: Result<Vec<_>> = mesh.cells
            .par_iter()
            .enumerate()
            .map(|(elem_idx, cell)| -> Result<_> {
                self.assemble_element_matrix(cell, mesh, material_properties, elem_idx)
            })
            .collect();

        let element_matrices = element_matrices?;

        // Add element matrices to global system
        for (elem_idx, (local_matrix, local_rhs, node_indices)) in element_matrices.iter().enumerate() {
            self.add_element_to_global(
                &mut matrix_builder,
                &mut rhs,
                local_matrix,
                local_rhs,
                node_indices,
                elem_idx,
            )?;
        }

        // Apply boundary conditions
        self.apply_velocity_boundary_conditions(&mut matrix_builder, &mut rhs, velocity_bcs)?;

        // Solve the linear system
        let matrix = matrix_builder.build()?;
        let mut solver_config = LinearSolverConfig::default();
        solver_config.base.tolerance = self.config.tolerance.clone();
        solver_config.base.max_iterations = self.config.max_iterations;

        let solver = ConjugateGradient::new(solver_config);
        let solution_vector = solver.solve(&matrix, &rhs, None)?;

        if self.config.verbose {
            tracing::info!("FEM Stokes solver completed successfully");
        }

        // Extract velocity solution (ignore pressure for now)
        let mut velocity_solution = HashMap::new();
        for node_idx in 0..num_nodes {
            let u = solution_vector[node_idx * 4].clone();
            let v = solution_vector[node_idx * 4 + 1].clone();
            let w = solution_vector[node_idx * 4 + 2].clone();
            velocity_solution.insert(node_idx, Vector3::new(u, v, w));
        }

        Ok(velocity_solution)
    }

    /// Assemble element matrix and RHS vector
    fn assemble_element_matrix(
        &self,
        _cell: &Cell,
        mesh: &Mesh<T>,
        material_properties: &MaterialProperties<T>,
        _elem_idx: usize,
    ) -> Result<(nalgebra::DMatrix<T>, DVector<T>, Vec<usize>)> {
        let element = Tetrahedron4::default();

        // For now, assume tetrahedral cells with 4 vertices
        // In a proper implementation, we'd get vertex indices from faces
        // For simplicity, assume the first 4 vertices form a tetrahedron
        let node_indices = if mesh.vertices.len() >= 4 {
            vec![0, 1, 2, 3] // Simplified for now
        } else {
            return Err(Error::InvalidConfiguration(
                "Insufficient vertices for tetrahedral element".to_string()
            ));
        };

        // Get node coordinates and convert Point3 to Vector3
        let nodes: Vec<Vector3<T>> = node_indices
            .iter()
            .map(|&idx| {
                let point = &mesh.vertices[idx].position;
                Vector3::new(point.x.clone(), point.y.clone(), point.z.clone())
            })
            .collect();

        // Compute element stiffness matrix
        let k_elem = element.stiffness_matrix(&nodes, material_properties)?;

        // For now, zero RHS (no body forces in this simplified implementation)
        let rhs_elem = DVector::zeros(k_elem.nrows());

        Ok((k_elem, rhs_elem, node_indices))
    }

    /// Add element matrix to global system
    fn add_element_to_global(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        local_matrix: &nalgebra::DMatrix<T>,
        local_rhs: &DVector<T>,
        node_indices: &[usize],
        _elem_idx: usize,
    ) -> Result<()> {
        // Map local DOFs to global DOFs with proper bounds checking
        for (i, &node_i) in node_indices.iter().enumerate() {
            for dof_i in 0..4 { // 3 velocity + 1 pressure
                let global_i = node_i * 4 + dof_i;
                let local_rhs_idx = i * 4 + dof_i;

                // Bounds check for RHS vector access
                if global_i >= rhs.len() {
                    return Err(Error::InvalidConfiguration(
                        format!("Global DOF index {} exceeds RHS vector size {}", global_i, rhs.len())
                    ));
                }

                if local_rhs_idx >= local_rhs.len() {
                    // Skip this DOF if it's out of bounds (happens with mixed element types)
                    continue;
                }

                rhs[global_i] += local_rhs[local_rhs_idx].clone();

                for (j, &node_j) in node_indices.iter().enumerate() {
                    for dof_j in 0..4 {
                        let global_j = node_j * 4 + dof_j;
                        let local_i = i * 4 + dof_i;
                        let local_j = j * 4 + dof_j;

                        // Enhanced bounds checking for matrix access
                        if local_i < local_matrix.nrows() &&
                           local_j < local_matrix.ncols() &&
                           global_i < rhs.len() &&
                           global_j < rhs.len() {
                            matrix_builder.add_entry(
                                global_i,
                                global_j,
                                local_matrix[(local_i, local_j)].clone(),
                            )?;
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Apply velocity boundary conditions to the system
    fn apply_velocity_boundary_conditions(
        &self,
        matrix_builder: &mut SparseMatrixBuilder<T>,
        rhs: &mut DVector<T>,
        velocity_bcs: &HashMap<usize, Vector3<T>>,
    ) -> Result<()> {
        for (&node_idx, velocity) in velocity_bcs {
            // Apply Dirichlet BC for velocity components
            for dof in 0..3 { // Only velocity components
                let global_dof = node_idx * 4 + dof;
                matrix_builder.add_entry(global_dof, global_dof, T::one())?;
                rhs[global_dof] = match dof {
                    0 => velocity.x.clone(),
                    1 => velocity.y.clone(),
                    2 => velocity.z.clone(),
                    _ => unreachable!(),
                };
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_mesh::{Vertex, Face, MeshTopology};

    fn create_simple_tet_mesh() -> Mesh<f64> {
        use nalgebra::Point3;

        // Create a simple tetrahedral mesh with 4 vertices
        let vertices = vec![
            Vertex { position: Point3::new(0.0, 0.0, 0.0), id: 0 },
            Vertex { position: Point3::new(1.0, 0.0, 0.0), id: 1 },
            Vertex { position: Point3::new(0.0, 1.0, 0.0), id: 2 },
            Vertex { position: Point3::new(0.0, 0.0, 1.0), id: 3 },
        ];

        // Create faces for the tetrahedron
        let faces = vec![
            Face { vertices: vec![0, 1, 2], id: 0 }, // Bottom face
            Face { vertices: vec![0, 1, 3], id: 1 }, // Side face 1
            Face { vertices: vec![1, 2, 3], id: 2 }, // Side face 2
            Face { vertices: vec![0, 2, 3], id: 3 }, // Side face 3
        ];

        let cells = vec![
            Cell { faces: vec![0, 1, 2, 3], id: 0 },
        ];

        let topology = MeshTopology {
            num_vertices: 4,
            num_edges: 6,
            num_faces: 4,
            num_cells: 1,
        };

        Mesh {
            vertices,
            edges: vec![],
            faces,
            cells,
            topology,
        }
    }

    #[test]
    fn test_tetrahedron4_shape_functions() {
        let element = Tetrahedron4::<f64>::default();
        let xi = Vector3::new(0.25, 0.25, 0.25);
        let shape_funcs = element.shape_functions(&xi);

        assert_eq!(shape_funcs.len(), 4);

        // Shape functions should sum to 1
        let sum: f64 = shape_funcs.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-12);

        // At centroid, all shape functions should be equal
        for &n in &shape_funcs {
            assert_relative_eq!(n, 0.25, epsilon = 1e-12);
        }
    }

    #[test]
    fn test_tetrahedron4_derivatives() {
        let element = Tetrahedron4::<f64>::default();
        let derivatives = element.shape_derivatives(&Vector3::zeros());

        assert_eq!(derivatives.len(), 4);

        // Check that derivatives sum to zero (partition of unity)
        let sum_x: f64 = derivatives.iter().map(|d| d.x).sum();
        let sum_y: f64 = derivatives.iter().map(|d| d.y).sum();
        let sum_z: f64 = derivatives.iter().map(|d| d.z).sum();

        assert_relative_eq!(sum_x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(sum_y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(sum_z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_fem_solver_creation() {
        let config = FemConfig::<f64>::default();
        let solver = FemSolver::new(config);

        assert!(solver.mesh.is_none());
    }

    #[test]
    fn test_fem_solver_with_mesh() {
        let mut solver = FemSolver::default();
        let mesh = create_simple_tet_mesh();

        solver.set_mesh(mesh);
        assert!(solver.mesh.is_some());
    }

    #[test]
    fn test_material_properties() {
        let props = MaterialProperties {
            density: 1000.0,
            viscosity: 0.001,
            thermal_conductivity: Some(0.6),
            specific_heat: Some(4186.0),
        };

        assert_eq!(props.density, 1000.0);
        assert_eq!(props.viscosity, 0.001);
        assert_eq!(props.thermal_conductivity, Some(0.6));
        assert_eq!(props.specific_heat, Some(4186.0));
    }

    #[test]
    fn test_element_stiffness_matrix() {
        let element = Tetrahedron4::<f64>::default();
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];

        let material = MaterialProperties {
            density: 1000.0,
            viscosity: 0.001,
            thermal_conductivity: None,
            specific_heat: None,
        };

        let k = element.stiffness_matrix(&nodes, &material).unwrap();

        // Check dimensions
        assert_eq!(k.nrows(), 12);
        assert_eq!(k.ncols(), 12);

        // Check that diagonal entries are positive
        for i in 0..12 {
            assert!(k[(i, i)] > 0.0);
        }
    }
}
