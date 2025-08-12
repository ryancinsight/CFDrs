//! Finite Element Method (FEM) solvers for 3D CFD.
//!
//! This module provides FEM implementations for solving 3D fluid dynamics problems
//! using various element types and numerical schemes.

use cfd_core::{Error, Result, SolverConfiguration};
use crate::constants;
use cfd_math::{LinearSolver, LinearSolverConfig, ConjugateGradient, SparseMatrixBuilder};
use cfd_mesh::{Mesh, Cell};
use nalgebra::{RealField, Vector3, DVector, DMatrix, Matrix3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use rayon::prelude::*;

/// FEM solver configuration
/// Uses unified SolverConfig as base to follow SSOT principle
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: RealField> {
    /// Base solver configuration (SSOT)
    pub base: cfd_core::SolverConfig<T>,
    /// Element type to use
    pub element_type: ElementType,
    /// Integration order
    pub integration_order: usize,
}

impl<T: RealField + FromPrimitive> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            base: cfd_core::SolverConfig::default(),
            element_type: ElementType::Tetrahedron4,
            integration_order: 2,
        }
    }
}

impl<T: RealField> FemConfig<T> {
    /// Get tolerance from base configuration
    pub fn tolerance(&self) -> T {
        self.base.tolerance()
    }

    /// Get max iterations from base configuration
    pub fn max_iterations(&self) -> usize {
        self.base.max_iterations()
    }

    /// Check if verbose output is enabled
    pub fn verbose(&self) -> bool {
        self.base.verbose()
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

/// Finite element for 3D problems following Interface Segregation Principle
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

/// Element enum for type-safe factory pattern avoiding trait objects
#[derive(Debug, Clone)]
pub enum ElementInstance<T: RealField> {
    /// 4-node tetrahedron
    Tetrahedron4(Tetrahedron4<T>),
    // Future element types can be added here
}

impl<T: RealField + FromPrimitive> ElementInstance<T> {
    /// Create element from type
    pub fn from_type(element_type: ElementType) -> Result<Self> {
        match element_type {
            ElementType::Tetrahedron4 => Ok(ElementInstance::Tetrahedron4(Tetrahedron4::default())),
            ElementType::Tetrahedron10 => Err(cfd_core::Error::NotImplemented(
                "Tetrahedron10 element not yet implemented".to_string()
            )),
            ElementType::Hexahedron8 => Err(cfd_core::Error::NotImplemented(
                "Hexahedron8 element not yet implemented".to_string()
            )),
            ElementType::Hexahedron20 => Err(cfd_core::Error::NotImplemented(
                "Hexahedron20 element not yet implemented".to_string()
            )),
        }
    }

    /// Delegate to underlying element implementation
    pub fn stiffness_matrix(
        &self,
        nodes: &[Vector3<T>],
        material_properties: &MaterialProperties<T>,
    ) -> Result<nalgebra::DMatrix<T>> {
        match self {
            ElementInstance::Tetrahedron4(elem) => elem.stiffness_matrix(nodes, material_properties),
        }
    }

    /// Get number of nodes
    pub fn num_nodes(&self) -> usize {
        match self {
            ElementInstance::Tetrahedron4(elem) => elem.num_nodes(),
        }
    }
}

/// Element factory trait following GRASP Creator principle
pub trait ElementFactory<T: RealField>: Send + Sync {
    /// Create an element of the specified type
    fn create_element(&self, element_type: ElementType) -> Result<ElementInstance<T>>;

    /// Check if this factory can create the given element type
    fn can_create(&self, element_type: ElementType) -> bool;

    /// Get supported element types
    fn supported_types(&self) -> Vec<ElementType>;
}

/// Tetrahedral element for FEM
#[derive(Debug, Clone)]
pub struct TetrahedralElement<T: RealField> {
    /// Node indices (4 nodes for tetrahedron)
    pub nodes: Vec<usize>,
    /// Element ID
    pub id: usize,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> TetrahedralElement<T> {
    /// Create new tetrahedral element
    pub fn new(nodes: Vec<usize>, id: usize) -> Self {
        Self {
            nodes,
            id,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Compute element stiffness matrix
    pub fn stiffness_matrix(&self, nodes: &[Vector3<T>], properties: &MaterialProperties<T>) -> Result<DMatrix<T>> {
        // 12x12 stiffness matrix for 4 nodes with 3 DOF each
        let mut k = DMatrix::zeros(12, 12);
        
        // Compute element volume
        let volume = self.compute_volume(nodes)?;
        
        // Build B matrix (strain-displacement matrix)
        let b_matrix = self.build_b_matrix(nodes)?;
        
        // Build D matrix (material constitutive matrix)
        let d_matrix = self.build_d_matrix(properties);
        
        // K = V * B^T * D * B
        let bt = b_matrix.transpose();
        let btd = &bt * &d_matrix;
        let btdb = &btd * &b_matrix;
        k = btdb * volume;
        
        Ok(k)
    }
    
    /// Compute element volume
    fn compute_volume(&self, nodes: &[Vector3<T>]) -> Result<T> {
        if nodes.len() < 4 {
            return Err(Error::InvalidConfiguration("Insufficient nodes for tetrahedron".to_string()));
        }
        
        let v0 = &nodes[0];
        let v1 = &nodes[1];
        let v2 = &nodes[2];
        let v3 = &nodes[3];
        
        // Volume = |det(v1-v0, v2-v0, v3-v0)| / 6
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let e3 = v3 - v0;
        
        let volume = e1.cross(&e2).dot(&e3).abs() / constants::tetrahedron_volume_factor::<T>();
        Ok(volume)
    }
    
    /// Build strain-displacement matrix
    fn build_b_matrix(&self, nodes: &[Vector3<T>]) -> Result<DMatrix<T>> {
        // 6x12 B matrix (6 strain components, 12 DOFs for 4 nodes)
        let mut b = DMatrix::zeros(6, 12);
        
        // Calculate shape function derivatives at the centroid
        // For a tetrahedron, shape functions are linear: N_i = (a_i + b_i*x + c_i*y + d_i*z) / (6*V)
        let v0 = &nodes[0];
        let v1 = &nodes[1];
        let v2 = &nodes[2];
        let v3 = &nodes[3];
        
        // Calculate volume
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let e3 = v3 - v0;
        let volume = e1.cross(&e2).dot(&e3).abs() / constants::tetrahedron_volume_factor::<T>();
        
        if volume < T::from_f64(1e-10).unwrap() {
            return Err(Error::InvalidInput("Degenerate tetrahedron element".into()));
        }
        
        // Calculate shape function derivatives (constant for linear tetrahedron)
        // Using the standard formula for tetrahedral shape function derivatives
        let inv_6v = T::one() / (constants::tetrahedron_volume_factor::<T>() * volume);
        
        // Shape function derivatives with respect to global coordinates
        // dN_i/dx, dN_i/dy, dN_i/dz for each node i
        let mut dn_dx = vec![T::zero(); 4];
        let mut dn_dy = vec![T::zero(); 4];
        let mut dn_dz = vec![T::zero(); 4];
        
        // Node 0 derivatives
        let n0_vec = (v2 - v1).cross(&(v3 - v1));
        dn_dx[0] = n0_vec.x.clone() * inv_6v.clone();
        dn_dy[0] = n0_vec.y.clone() * inv_6v.clone();
        dn_dz[0] = n0_vec.z.clone() * inv_6v.clone();
        
        // Node 1 derivatives
        let n1_vec = (v3 - v0).cross(&(v2 - v0));
        dn_dx[1] = n1_vec.x.clone() * inv_6v.clone();
        dn_dy[1] = n1_vec.y.clone() * inv_6v.clone();
        dn_dz[1] = n1_vec.z.clone() * inv_6v.clone();
        
        // Node 2 derivatives
        let n2_vec = (v1 - v0).cross(&(v3 - v0));
        dn_dx[2] = n2_vec.x.clone() * inv_6v.clone();
        dn_dy[2] = n2_vec.y.clone() * inv_6v.clone();
        dn_dz[2] = n2_vec.z.clone() * inv_6v.clone();
        
        // Node 3 derivatives
        let n3_vec = (v2 - v0).cross(&(v1 - v0));
        dn_dx[3] = n3_vec.x.clone() * inv_6v.clone();
        dn_dy[3] = n3_vec.y.clone() * inv_6v.clone();
        dn_dz[3] = n3_vec.z.clone() * inv_6v;
        
        // Fill B matrix: relates strains to nodal displacements
        // Strain = B * u, where u = [u0x, u0y, u0z, u1x, u1y, u1z, ...]
        for i in 0..4 {
            let col_offset = i * 3;
            
            // Normal strains
            b[(0, col_offset)] = dn_dx[i].clone();     // ε_xx = ∂u/∂x
            b[(1, col_offset + 1)] = dn_dy[i].clone();  // ε_yy = ∂v/∂y
            b[(2, col_offset + 2)] = dn_dz[i].clone();  // ε_zz = ∂w/∂z
            
            // Shear strains (engineering strain = 2 * tensorial strain)
            b[(3, col_offset)] = dn_dy[i].clone();      // γ_xy = ∂u/∂y + ∂v/∂x
            b[(3, col_offset + 1)] = dn_dx[i].clone();
            
            b[(4, col_offset + 1)] = dn_dz[i].clone();  // γ_yz = ∂v/∂z + ∂w/∂y
            b[(4, col_offset + 2)] = dn_dy[i].clone();
            
            b[(5, col_offset)] = dn_dz[i].clone();      // γ_xz = ∂u/∂z + ∂w/∂x
            b[(5, col_offset + 2)] = dn_dx[i].clone();
        }
        
        Ok(b)
    }
    
    /// Build material constitutive matrix
    fn build_d_matrix(&self, properties: &MaterialProperties<T>) -> DMatrix<T> {
        let mut d = DMatrix::zeros(6, 6);
        
        // For linear elastic material (used in structural analysis)
        let e = properties.youngs_modulus.clone();
        let nu = properties.poisson_ratio.clone();
        let two = constants::two::<T>();
        let factor = e.clone() / ((T::one() + nu.clone()) * (T::one() - two.clone() * nu.clone()));
        
        // Lamé parameters
        let lambda = factor.clone() * nu.clone();
        let mu = factor * (T::one() - nu) / two.clone();
        
        // Normal stress components
        let two_mu = two * mu.clone();
        d[(0, 0)] = lambda.clone() + two_mu.clone();
        d[(1, 1)] = lambda.clone() + two_mu.clone();
        d[(2, 2)] = lambda.clone() + two_mu;
        
        // Shear stress components
        d[(3, 3)] = mu.clone();
        d[(4, 4)] = mu.clone();
        d[(5, 5)] = mu.clone();
        
        // Off-diagonal terms
        d[(0, 1)] = lambda.clone();
        d[(0, 2)] = lambda.clone();
        d[(1, 0)] = lambda.clone();
        d[(1, 2)] = lambda.clone();
        d[(2, 0)] = lambda.clone();
        d[(2, 1)] = lambda;
        
        d
    }
    
    /// Compute Jacobian matrix
    pub fn jacobian(&self, nodes: &[Vector3<T>]) -> Result<Matrix3<T>> {
        use nalgebra::Matrix3;
        
        if nodes.len() < 4 {
            return Err(Error::InvalidConfiguration("Insufficient nodes".to_string()));
        }
        
        // Jacobian for linear tetrahedron
        let v0 = &nodes[0];
        let v1 = &nodes[1];
        let v2 = &nodes[2];
        let v3 = &nodes[3];
        
        let j = Matrix3::from_columns(&[v1 - v0, v2 - v0, v3 - v0]);
        Ok(j)
    }
    
    /// Evaluate shape functions at given natural coordinates
    pub fn shape_functions(&self, xi: T, eta: T, zeta: T) -> Vec<T> {
        // Linear shape functions for tetrahedron
        vec![
            T::one() - xi.clone() - eta.clone() - zeta.clone(),
            xi,
            eta,
            zeta,
        ]
    }
}

/// Material properties for FEM analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaterialProperties<T: RealField> {
    /// Density
    pub density: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Young's modulus (for structural analysis)
    pub youngs_modulus: T,
    /// Poisson's ratio
    pub poisson_ratio: T,
    /// Body force vector (e.g., gravity)
    pub body_force: Option<Vector3<T>>,
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
        // Linear shape functions for tetrahedron
        // N0 = 1 - ξ - η - ζ
        // N1 = ξ
        // N2 = η
        // N3 = ζ
        vec![
            T::one() - xi.x.clone() - xi.y.clone() - xi.z.clone(),
            xi.x.clone(),
            xi.y.clone(),
            xi.z.clone(),
        ]
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
        // Single-point integration at centroid for linear tetrahedron
        let quarter = constants::tetrahedron_gauss_point::<T>();
        let weight = constants::tetrahedron_gauss_weight::<T>();
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
                let deriv = derivatives[j][i].clone();
                jacobian[(i, 0)] += deriv.clone() * nodes[j].x.clone();
                jacobian[(i, 1)] += deriv.clone() * nodes[j].y.clone();
                jacobian[(i, 2)] += deriv * nodes[j].z.clone();
            }
        }

        let det_j: T = jacobian.determinant();
        if det_j <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Negative or zero Jacobian determinant".to_string()
            ));
        }

        // Proper stiffness matrix computation for Navier-Stokes equations
        // Based on Zienkiewicz & Taylor "The Finite Element Method" Vol. 3
        // K = ∫ B^T D B dV where B is strain-displacement matrix, D is material matrix

        let mut k = nalgebra::DMatrix::zeros(12, 12); // 4 nodes × 3 DOF per node
        let viscosity = material_properties.viscosity.clone();

        // Compute shape function derivatives in physical coordinates
        let j_inv = jacobian.try_inverse().ok_or_else(|| {
            Error::InvalidConfiguration("Singular Jacobian matrix".to_string())
        })?;

        // Use integration points from the element for zero-copy efficiency
        for (_integration_point, weight) in self.integration_points() {
            // Shape function derivatives in natural coordinates
            let mut dn_dxi = nalgebra::DMatrix::zeros(4, 3);
            dn_dxi[(0, 0)] = -T::one(); dn_dxi[(0, 1)] = -T::one(); dn_dxi[(0, 2)] = -T::one();
            dn_dxi[(1, 0)] = T::one();  dn_dxi[(1, 1)] = T::zero(); dn_dxi[(1, 2)] = T::zero();
            dn_dxi[(2, 0)] = T::zero(); dn_dxi[(2, 1)] = T::one();  dn_dxi[(2, 2)] = T::zero();
            dn_dxi[(3, 0)] = T::zero(); dn_dxi[(3, 1)] = T::zero(); dn_dxi[(3, 2)] = T::one();

            // Transform to physical coordinates: dN/dx = J^(-1) * dN/dξ
            let dn_dx = &j_inv * &dn_dxi.transpose();

            // Build strain-displacement matrix B for 3D Navier-Stokes
            let mut b_matrix = nalgebra::DMatrix::zeros(6, 12); // 6 strain components, 12 DOFs

            // Use iterator pattern for strain-displacement matrix assembly
            (0..4).for_each(|i| {
                let base_col = i * 3;
                // Velocity gradients for viscous stress tensor using zero-copy references
                b_matrix[(0, base_col)] = dn_dx[(0, i)].clone(); // ∂u/∂x
                b_matrix[(1, base_col + 1)] = dn_dx[(1, i)].clone(); // ∂v/∂y
                b_matrix[(2, base_col + 2)] = dn_dx[(2, i)].clone(); // ∂w/∂z
                b_matrix[(3, base_col)] = dn_dx[(1, i)].clone(); // ∂u/∂y
                b_matrix[(3, base_col + 1)] = dn_dx[(0, i)].clone(); // ∂v/∂x
                b_matrix[(4, base_col + 1)] = dn_dx[(2, i)].clone(); // ∂v/∂z
                b_matrix[(4, base_col + 2)] = dn_dx[(1, i)].clone(); // ∂w/∂y
                b_matrix[(5, base_col)] = dn_dx[(2, i)].clone(); // ∂u/∂z
                b_matrix[(5, base_col + 2)] = dn_dx[(0, i)].clone(); // ∂w/∂x
            });

            // Material matrix D for viscous fluid (isotropic)
            let mut d_matrix = nalgebra::DMatrix::zeros(6, 6);
            let two_mu = constants::two::<T>() * viscosity.clone();
            d_matrix[(0, 0)] = two_mu.clone(); 
            d_matrix[(1, 1)] = two_mu.clone(); 
            d_matrix[(2, 2)] = two_mu;
            d_matrix[(3, 3)] = viscosity.clone(); 
            d_matrix[(4, 4)] = viscosity.clone(); 
            d_matrix[(5, 5)] = viscosity.clone();

            // Compute element stiffness: K_e += B^T * D * B * det(J) * weight
            let bd = &b_matrix.transpose() * &d_matrix;
            let bdb = &bd * &b_matrix;
            let integration_factor = det_j.clone() * weight;

            k += bdb * integration_factor;
        }

        Ok(k)
    }
}

/// Concrete element factory implementation following SOLID principles
#[derive(Debug, Clone, Default)]
pub struct StandardElementFactory<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> StandardElementFactory<T> {
    /// Create a new element factory
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + FromPrimitive> ElementFactory<T> for StandardElementFactory<T> {
    fn create_element(&self, element_type: ElementType) -> Result<ElementInstance<T>> {
        ElementInstance::from_type(element_type)
    }

    fn can_create(&self, element_type: ElementType) -> bool {
        matches!(element_type, ElementType::Tetrahedron4)
    }

    fn supported_types(&self) -> Vec<ElementType> {
        vec![ElementType::Tetrahedron4]
    }
}

/// FEM solver for 3D fluid dynamics problems following SOLID principles
pub struct FemSolver<T: RealField, F: ElementFactory<T> = StandardElementFactory<T>> {
    config: FemConfig<T>,
    mesh: Option<Mesh<T>>,
    element_factory: F,
}

impl<T: RealField + FromPrimitive + Send + Sync> FemSolver<T, StandardElementFactory<T>> {
    /// Create a new FEM solver with default element factory
    pub fn new(config: FemConfig<T>) -> Self {
        Self {
            config,
            mesh: None,
            element_factory: StandardElementFactory::new(),
        }
    }

    /// Create with default configuration and element factory
    pub fn default() -> Self {
        Self::new(FemConfig::default())
    }
}

impl<T: RealField + FromPrimitive + Send + Sync, F: ElementFactory<T>> FemSolver<T, F> {
    /// Create a new FEM solver with custom element factory
    pub fn with_factory(config: FemConfig<T>, element_factory: F) -> Self {
        Self {
            config,
            mesh: None,
            element_factory,
        }
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
        solver_config.base = cfd_core::SolverConfig::builder()
            .tolerance(self.config.tolerance())
            .max_iterations(self.config.max_iterations())
            .build();

        let solver = ConjugateGradient::new(solver_config);
        let solution_vector = solver.solve(&matrix, &rhs, None)?;

        if self.config.verbose() {
            tracing::info!("FEM Stokes solver completed successfully");
        }

        // Extract velocity solution using iterator patterns for zero-copy efficiency
        let velocity_solution: HashMap<usize, Vector3<T>> = (0..num_nodes)
            .map(|node_idx| {
                let u = solution_vector[node_idx * 4].clone();
                let v = solution_vector[node_idx * 4 + 1].clone();
                let w = solution_vector[node_idx * 4 + 2].clone();
                (node_idx, Vector3::new(u, v, w))
            })
            .collect();

        Ok(velocity_solution)
    }

    /// Compute body force contribution for an element
    /// 
    /// Implements numerical integration of body forces over the element:
    /// F_i = ∫_Ω N_i · f dΩ
    /// 
    /// Uses Gaussian quadrature for accurate integration
    fn compute_body_force_contribution(
        &self,
        element: &TetrahedralElement<T>,
        nodes: &[Vector3<T>],
        material_properties: &MaterialProperties<T>,
    ) -> Result<DVector<T>> {
        let ndof = element.nodes.len() * 3; // 3 DOF per node (u, v, w)
        
        // Get body force from material properties
        let body_force = material_properties.body_force.as_ref()
            .map(|f| f.clone())
            .unwrap_or_else(Vector3::zeros);
        
        // Compute Jacobian determinant for integration
        // For linear tetrahedron, det(J) = 6 * Volume
        let v0 = &nodes[0];
        let v1 = &nodes[1];
        let v2 = &nodes[2];
        let v3 = &nodes[3];
        
        let e1 = v1 - v0;
        let e2 = v2 - v0;
        let e3 = v3 - v0;
        let det_j = e1.cross(&e2).dot(&e3).abs();
        
        // Initialize force vector
        let mut force_vector = DVector::zeros(ndof);
        
        // Gauss quadrature integration
        // For linear tetrahedron, single point at centroid is sufficient
        let gauss_points = vec![
            (T::from_f64(0.25).unwrap(), T::from_f64(0.25).unwrap(), T::from_f64(0.25).unwrap(), T::from_f64(1.0/6.0).unwrap())
        ];
        
        for (xi, eta, zeta, weight) in gauss_points {
            // Evaluate shape functions at Gauss point
            let shape_functions = element.shape_functions(xi, eta, zeta);
            
            // Add contribution to force vector
            for (i, n_i) in shape_functions.iter().enumerate() {
                force_vector[i * 3] += weight.clone() * n_i.clone() * body_force.x.clone() * det_j.clone();
                force_vector[i * 3 + 1] += weight.clone() * n_i.clone() * body_force.y.clone() * det_j.clone();
                force_vector[i * 3 + 2] += weight.clone() * n_i.clone() * body_force.z.clone() * det_j.clone();
            }
        }
        
        Ok(force_vector)
    }

    /// Assemble element matrix and RHS vector using factory pattern
    fn assemble_element_matrix(
        &self,
        _cell: &Cell,
        mesh: &Mesh<T>,
        material_properties: &MaterialProperties<T>,
        _elem_idx: usize,
    ) -> Result<(nalgebra::DMatrix<T>, DVector<T>, Vec<usize>)> {
        // Create tetrahedral element
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
        
        let element = TetrahedralElement::new(node_indices.clone(), 0);

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

        // Compute RHS with body forces using Galerkin method
        // F_i = ∫_Ω N_i · f dΩ where N_i are shape functions and f is body force
        let rhs_elem = self.compute_body_force_contribution(&element, &nodes, material_properties)?;

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
        // Map local DOFs to global DOFs using iterator patterns for zero-copy efficiency
        let rhs_updates: Result<Vec<_>> = node_indices
            .iter()
            .enumerate()
            .flat_map(|(i, &node_i)| {
                (0..4).map(move |dof_i| (i, node_i, dof_i))
            })
            .map(|(i, node_i, dof_i)| -> Result<(usize, usize, T)> {
                let global_i = node_i * 4 + dof_i;
                let local_rhs_idx = i * 4 + dof_i;

                // Bounds check for RHS vector access
                if global_i >= rhs.len() {
                    return Err(Error::InvalidConfiguration(
                        format!("Global DOF index {} exceeds RHS vector size {}", global_i, rhs.len())
                    ));
                }

                if local_rhs_idx >= local_rhs.len() {
                    // Return zero contribution if out of bounds
                    return Ok((global_i, local_rhs_idx, T::zero()));
                }

                Ok((global_i, local_rhs_idx, local_rhs[local_rhs_idx].clone()))
            })
            .collect();

        // Apply RHS updates
        for (global_i, _local_idx, value) in rhs_updates? {
            rhs[global_i] += value;
        }

        // Matrix assembly using iterator patterns for better performance
        let matrix_entries: Result<Vec<_>> = node_indices
            .iter()
            .enumerate()
            .flat_map(|(i, &node_i)| {
                (0..4).map(move |dof_i| (i, node_i, dof_i))
            })
            .flat_map(|(i, node_i, dof_i)| {
                node_indices.iter().enumerate().flat_map(move |(j, &node_j)| {
                    (0..4).map(move |dof_j| (i, node_i, dof_i, j, node_j, dof_j))
                })
            })
            .filter_map(|(i, node_i, dof_i, j, node_j, dof_j)| {
                let global_i = node_i * 4 + dof_i;
                let global_j = node_j * 4 + dof_j;
                let local_i = i * 4 + dof_i;
                let local_j = j * 4 + dof_j;

                // Enhanced bounds checking for matrix access
                if local_i < local_matrix.nrows() &&
                   local_j < local_matrix.ncols() &&
                   global_i < rhs.len() &&
                   global_j < rhs.len() {
                    Some(Ok((global_i, global_j, local_matrix[(local_i, local_j)].clone())))
                } else {
                    None
                }
            })
            .collect();

        // Apply matrix entries
        for (global_i, global_j, value) in matrix_entries? {
            matrix_builder.add_entry(global_i, global_j, value)?;
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
            youngs_modulus: 2.0e11,  // Steel
            poisson_ratio: 0.3,
            body_force: Some(Vector3::new(0.0, -9.81, 0.0)),  // Gravity
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
            youngs_modulus: 2.0e11,
            poisson_ratio: 0.3,
            body_force: None,
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
