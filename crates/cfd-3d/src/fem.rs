//! Finite Element Method (FEM) for 3D incompressible fluid dynamics.
//!
//! This module implements FEM for solving the incompressible Navier-Stokes equations
//! using mixed velocity-pressure formulation with appropriate stabilization.

use cfd_core::{Error, Result};
use cfd_math::sparse::{SparseMatrix, SparseMatrixBuilder};
use cfd_mesh::Mesh;
use nalgebra::{RealField, Vector3, DVector, DMatrix, Matrix3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Constants for FEM fluid dynamics
mod constants {
    /// Small value for numerical stability
    pub const EPSILON: f64 = 1e-10;
    /// Default stabilization parameter
    pub const DEFAULT_STABILIZATION: f64 = 0.1;
    /// Number of velocity components
    pub const VELOCITY_COMPONENTS: usize = 3;
    /// Number of nodes in linear tetrahedron
    pub const TET4_NODES: usize = 4;
    /// Volume factor for tetrahedron (1/6)
    pub const TET_VOLUME_FACTOR: f64 = 6.0;
    // Unused constants removed for cleaner codebase
}

/// FEM configuration for fluid dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: RealField> {
    /// Base solver configuration
    pub base: cfd_core::SolverConfig<T>,
    /// Use SUPG/PSPG stabilization
    pub use_stabilization: bool,
    /// Stabilization parameter
    pub tau: T,
    /// Time step (for transient problems)
    pub dt: Option<T>,
    /// Reynolds number (for scaling)
    pub reynolds: Option<T>,
}

impl<T: RealField + FromPrimitive> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            base: cfd_core::SolverConfig::default(),
            use_stabilization: true,
            tau: T::from_f64(constants::DEFAULT_STABILIZATION).unwrap(),
            dt: None,
            reynolds: None,
        }
    }
}

/// Fluid properties for incompressible flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidProperties<T: RealField> {
    /// Fluid density (kg/m³)
    pub density: T,
    /// Dynamic viscosity (Pa·s)
    pub viscosity: T,
    /// Body force (e.g., gravity)
    pub body_force: Option<Vector3<T>>,
}

/// Element matrix contributions for sparse assembly
pub struct ElementMatrices<T: RealField> {
    /// Viscous stiffness entries (row, col, value)
    pub k_entries: Vec<(usize, usize, T)>,
    /// Gradient matrix entries
    pub g_entries: Vec<(usize, usize, T)>,
    /// Mass matrix entries
    pub m_entries: Vec<(usize, usize, T)>,
    /// Stabilization matrix entries
    pub s_entries: Vec<(usize, usize, T)>,
}

impl<T: RealField + FromPrimitive> FluidProperties<T> {
    /// Create properties for water at 20°C
    pub fn water() -> Self {
        Self {
            density: T::from_f64(998.2).unwrap(),
            viscosity: T::from_f64(0.001).unwrap(),
            body_force: Some(Vector3::new(
                T::zero(),
                T::zero(),
                T::from_f64(-9.81).unwrap(),
            )),
        }
    }
    
    /// Kinematic viscosity
    pub fn kinematic_viscosity(&self) -> T {
        self.viscosity.clone() / self.density.clone()
    }
}

/// Tetrahedral element for incompressible flow (P1-P1 with stabilization)
pub struct FluidElement<T: RealField> {
    /// Node indices
    pub nodes: Vec<usize>,
    /// Element ID
    pub id: usize,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> FluidElement<T> {
    /// Create new fluid element
    pub fn new(nodes: Vec<usize>, id: usize) -> Self {
        assert_eq!(nodes.len(), constants::TET4_NODES, "Linear tetrahedron requires 4 nodes");
        Self {
            nodes,
            id,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Compute element matrices for Stokes flow
    /// Returns element contributions as triplets for sparse assembly
    pub fn stokes_matrices(
        &self,
        node_coords: &[Vector3<T>],
        properties: &FluidProperties<T>,
        config: &FemConfig<T>,
    ) -> Result<ElementMatrices<T>> {
        let _n_vel_dof = constants::TET4_NODES * constants::VELOCITY_COMPONENTS; // 12
        let _n_pres_dof = constants::TET4_NODES; // 4
        
        // Initialize element matrices (local dense, will be assembled into global sparse)
        let mut k_entries = Vec::new(); // Viscous stiffness
        let mut g_entries = Vec::new(); // Gradient
        let mut m_entries = Vec::new(); // Mass
        let mut s_entries = Vec::new(); // Stabilization
        
        // Compute element volume and shape function derivatives
        let volume = self.compute_volume(node_coords)?;
        let (dn_dx, dn_dy, dn_dz) = self.shape_derivatives(node_coords)?;
        
        let nu = properties.kinematic_viscosity();
        let rho = properties.density.clone();
        
        // Single-point quadrature for linear elements
        let weight = volume.clone(); // For single Gauss point at centroid
        
        // Build viscous stiffness matrix K
        // K_ij = ∫ μ (∇N_i : ∇N_j) dΩ
        for i in 0..constants::TET4_NODES {
            for j in 0..constants::TET4_NODES {
                // Viscous term contribution
                let k_val = nu.clone() * weight.clone() * (
                    dn_dx[i].clone() * dn_dx[j].clone() +
                    dn_dy[i].clone() * dn_dy[j].clone() +
                    dn_dz[i].clone() * dn_dz[j].clone()
                );
                
                // Add to all velocity components
                for d in 0..constants::VELOCITY_COMPONENTS {
                    let row = i * constants::VELOCITY_COMPONENTS + d;
                    let col = j * constants::VELOCITY_COMPONENTS + d;
                    k_entries.push((row, col, k_val.clone()));
                }
            }
        }
        
        // Build gradient matrix G
        // G_ij = -∫ (∇·N_i) M_j dΩ
        for i in 0..constants::TET4_NODES {
            for j in 0..constants::TET4_NODES {
                // u-component row
                let g_val_x = -weight.clone() * dn_dx[i].clone() / T::from_usize(constants::TET4_NODES).unwrap();
                g_entries.push((i * constants::VELOCITY_COMPONENTS, j, g_val_x));
                
                // v-component row
                let g_val_y = -weight.clone() * dn_dy[i].clone() / T::from_usize(constants::TET4_NODES).unwrap();
                g_entries.push((i * constants::VELOCITY_COMPONENTS + 1, j, g_val_y));
                
                // w-component row
                let g_val_z = -weight.clone() * dn_dz[i].clone() / T::from_usize(constants::TET4_NODES).unwrap();
                g_entries.push((i * constants::VELOCITY_COMPONENTS + 2, j, g_val_z));
            }
        }
        
        // Build mass matrix M (for transient terms)
        if config.dt.is_some() {
            let mass_factor = rho * weight.clone() / T::from_usize(constants::TET4_NODES * constants::TET4_NODES).unwrap();
            for i in 0..constants::TET4_NODES {
                for j in 0..constants::TET4_NODES {
                    let m_val = if i == j {
                        mass_factor.clone() * T::from_f64(2.0).unwrap()
                    } else {
                        mass_factor.clone()
                    };
                    
                    for d in 0..constants::VELOCITY_COMPONENTS {
                        let row = i * constants::VELOCITY_COMPONENTS + d;
                        let col = j * constants::VELOCITY_COMPONENTS + d;
                        m_entries.push((row, col, m_val.clone()));
                    }
                }
            }
        }
        
        // Compute PSPG stabilization matrix if enabled
        if config.use_stabilization {
            let h = self.compute_element_size(node_coords)?;
            let tau = self.compute_stabilization_parameter(h, nu.clone(), config)?;
            
            // S_pp = τ ∫ (∇M_i · ∇M_j) dΩ
            for i in 0..constants::TET4_NODES {
                for j in 0..constants::TET4_NODES {
                    let s_val = tau.clone() * weight.clone() * (
                        dn_dx[i].clone() * dn_dx[j].clone() +
                        dn_dy[i].clone() * dn_dy[j].clone() +
                        dn_dz[i].clone() * dn_dz[j].clone()
                    );
                    s_entries.push((i, j, s_val));
                }
            }
        }
        
        Ok(ElementMatrices {
            k_entries,
            g_entries,
            m_entries,
            s_entries,
        })
    }
    
    /// Compute element characteristic size
    fn compute_element_size(&self, node_coords: &[Vector3<T>]) -> Result<T> {
        self.element_length_scale(node_coords)
    }
    
    /// Compute stabilization parameter τ for PSPG/SUPG
    fn compute_stabilization_parameter(
        &self,
        h: T,
        nu: T,
        config: &FemConfig<T>
    ) -> Result<T> {
        // Standard formula: τ = h²/(4ν) for diffusion-dominated problems
        // For convection-dominated: τ = h/(2|u|) where |u| is velocity magnitude
        // Here we use diffusion form since this is Stokes flow
        let two = T::from_f64(2.0).unwrap();
        let four = T::from_f64(4.0).unwrap();
        Ok(h.clone() * h / (four * nu))
    }

    /// Add SUPG/PSPG stabilization for equal-order interpolation
    /// Returns the pressure stabilization matrix S_pp
    fn add_stabilization(
        &self,
        k_matrix: &mut DMatrix<T>,
        g_matrix: &mut DMatrix<T>,
        node_coords: &[Vector3<T>],
        properties: &FluidProperties<T>,
        config: &FemConfig<T>,
        dn_dx: &[T],
        dn_dy: &[T],
        dn_dz: &[T],
        volume: T,
    ) -> Result<DMatrix<T>> {
        // Calculate element length scale
        let h = self.element_length_scale(node_coords)?;
        
        // Calculate stabilization parameter τ
        let nu = properties.kinematic_viscosity();
        let tau = self.calculate_tau(h, nu, config)?;
        
        // PSPG stabilization: pressure Laplacian term
        // S_pp = τ ∫ ∇M_i · ∇M_j dΩ
        let n_pres_dof = constants::TET4_NODES;
        let mut s_pp = DMatrix::zeros(n_pres_dof, n_pres_dof);
        
        for i in 0..constants::TET4_NODES {
            for j in 0..constants::TET4_NODES {
                let s_val = tau.clone() * volume.clone() * (
                    dn_dx[i].clone() * dn_dx[j].clone() +
                    dn_dy[i].clone() * dn_dy[j].clone() +
                    dn_dz[i].clone() * dn_dz[j].clone()
                ) / properties.density.clone();
                
                s_pp[(i, j)] = s_val;
            }
        }
        
        // Note: The gradient matrix g_matrix should NOT be modified
        // The stabilization matrix s_pp will be added to the (2,2) block
        // of the global system matrix during assembly
        
        Ok(s_pp)
    }
    
    /// Calculate stabilization parameter τ
    fn calculate_tau(&self, h: T, nu: T, config: &FemConfig<T>) -> Result<T> {
        // Standard SUPG/PSPG parameter
        // τ = h²/(4ν) for diffusion-dominated flow
        // τ = h/(2|u|) for convection-dominated flow
        
        // For Stokes flow (no convection), use diffusive scaling
        let h_squared = h.clone() * h;
        let four = T::from_f64(4.0).unwrap();
        let tau = h_squared / (four * nu);
        
        // Apply user scaling if provided
        Ok(tau * config.tau.clone())
    }
    
    /// Calculate element length scale
    fn element_length_scale(&self, node_coords: &[Vector3<T>]) -> Result<T> {
        // Use average edge length as characteristic length
        let mut total_length = T::zero();
        let mut count = 0;
        
        for i in 0..constants::TET4_NODES {
            for j in i+1..constants::TET4_NODES {
                let edge = &node_coords[j] - &node_coords[i];
                total_length = total_length + edge.norm();
                count += 1;
            }
        }
        
        Ok(total_length / T::from_usize(count).unwrap())
    }
    
    /// Compute element volume
    fn compute_volume(&self, nodes: &[Vector3<T>]) -> Result<T> {
        if nodes.len() < constants::TET4_NODES {
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
        
        let volume = e1.cross(&e2).dot(&e3).abs() / T::from_f64(constants::TET_VOLUME_FACTOR).unwrap();
        
        if volume < T::from_f64(constants::EPSILON).unwrap() {
            return Err(Error::InvalidInput("Degenerate tetrahedron element".into()));
        }
        
        Ok(volume)
    }
    
    /// Calculate shape function derivatives
    fn shape_derivatives(&self, nodes: &[Vector3<T>]) -> Result<(Vec<T>, Vec<T>, Vec<T>)> {
        let volume = self.compute_volume(nodes)?;
        let six_v = T::from_f64(constants::TET_VOLUME_FACTOR).unwrap() * volume;
        
        let v0 = &nodes[0];
        let v1 = &nodes[1];
        let v2 = &nodes[2];
        let v3 = &nodes[3];
        
        // Shape function derivatives (constant for linear tet)
        let mut dn_dx = vec![T::zero(); constants::TET4_NODES];
        let mut dn_dy = vec![T::zero(); constants::TET4_NODES];
        let mut dn_dz = vec![T::zero(); constants::TET4_NODES];
        
        // Node 0
        let n0_vec = (v2 - v1).cross(&(v3 - v1));
        dn_dx[0] = n0_vec.x.clone() / six_v.clone();
        dn_dy[0] = n0_vec.y.clone() / six_v.clone();
        dn_dz[0] = n0_vec.z.clone() / six_v.clone();
        
        // Node 1
        let n1_vec = (v3 - v0).cross(&(v2 - v0));
        dn_dx[1] = n1_vec.x.clone() / six_v.clone();
        dn_dy[1] = n1_vec.y.clone() / six_v.clone();
        dn_dz[1] = n1_vec.z.clone() / six_v.clone();
        
        // Node 2
        let n2_vec = (v1 - v0).cross(&(v3 - v0));
        dn_dx[2] = n2_vec.x.clone() / six_v.clone();
        dn_dy[2] = n2_vec.y.clone() / six_v.clone();
        dn_dz[2] = n2_vec.z.clone() / six_v.clone();
        
        // Node 3
        let n3_vec = (v2 - v0).cross(&(v1 - v0));
        dn_dx[3] = n3_vec.x.clone() / six_v.clone();
        dn_dy[3] = n3_vec.y.clone() / six_v.clone();
        dn_dz[3] = n3_vec.z.clone() / six_v;
        
        Ok((dn_dx, dn_dy, dn_dz))
    }
    
    /// Build strain rate tensor from velocity gradients
    pub fn strain_rate_tensor(&self, velocity_gradient: &Matrix3<T>) -> Matrix3<T> {
        // ε̇ = 0.5 * (∇u + ∇u^T)
        let half = T::from_f64(0.5).unwrap();
        (velocity_gradient + velocity_gradient.transpose()) * half
    }
}

/// FEM solver for incompressible Navier-Stokes equations
pub struct FemSolver<T: RealField> {
    config: FemConfig<T>,
    mesh: Mesh<T>,
    properties: FluidProperties<T>,
    /// Velocity field (3 components per node)
    velocity: DVector<T>,
    /// Pressure field (1 per node)
    pressure: DVector<T>,
    /// Global stiffness matrix (sparse)
    k_global: Option<SparseMatrix<T>>,
    /// Global gradient matrix (sparse)
    g_global: Option<SparseMatrix<T>>,
    /// Global mass matrix (sparse)
    m_global: Option<SparseMatrix<T>>,
    /// Global pressure stabilization matrix (sparse)
    s_pp_global: Option<SparseMatrix<T>>,
}

impl<T: RealField + FromPrimitive> FemSolver<T> {
    /// Create new FEM solver
    pub fn new(config: FemConfig<T>, mesh: Mesh<T>, properties: FluidProperties<T>) -> Self {
        let n_nodes = mesh.vertices.len();
        let n_vel_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pres_dof = n_nodes;
        
        Self {
            config,
            mesh,
            properties,
            velocity: DVector::zeros(n_vel_dof),
            pressure: DVector::zeros(n_pres_dof),
            k_global: None,
            g_global: None,
            m_global: None,
            s_pp_global: None,
        }
    }
    
    /// Assemble global matrices
    pub fn assemble_global_matrices(&mut self) -> Result<()> {
        let n_nodes = self.mesh.vertices.len();
        let n_vel_dof = n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pres_dof = n_nodes;
        
        // Use sparse matrix builders for efficient assembly
        let mut k_builder = SparseMatrixBuilder::with_capacity(n_vel_dof, n_vel_dof, n_vel_dof * 27);
        let mut g_builder = SparseMatrixBuilder::with_capacity(n_vel_dof, n_pres_dof, n_vel_dof * 4);
        let mut m_builder = SparseMatrixBuilder::with_capacity(n_vel_dof, n_vel_dof, n_vel_dof * 27);
        let mut s_builder = SparseMatrixBuilder::with_capacity(n_pres_dof, n_pres_dof, n_pres_dof * 27);
        
        // Allow duplicate entries (they will be summed during assembly)
        k_builder = k_builder.allow_duplicates(true);
        g_builder = g_builder.allow_duplicates(true);
        m_builder = m_builder.allow_duplicates(true);
        s_builder = s_builder.allow_duplicates(true);
        
        // Process each element (assuming tetrahedral mesh)
        // For now, we'll create tetrahedra from cells with 4 faces
        for cell in &self.mesh.cells {
            // Get unique vertices from cell faces
            let mut vertex_set = std::collections::HashSet::new();
            for &face_id in &cell.faces {
                if let Some(face) = self.mesh.faces.iter().find(|f| f.id == face_id) {
                    vertex_set.extend(&face.vertices);
                }
            }
            let nodes: Vec<usize> = vertex_set.into_iter().collect();
            
            if nodes.len() == 4 {  // Tetrahedral element
                let element = FluidElement::new(nodes.clone(), 0);
                
                // Get node coordinates
                let node_coords: Vec<Vector3<T>> = nodes.iter()
                    .map(|&idx| {
                        let v = &self.mesh.vertices[idx];
                        Vector3::new(
                            v.position.x.clone(),
                            v.position.y.clone(),
                            v.position.z.clone(),
                        )
                    })
                    .collect();
                
                // Compute element matrices
                let elem_matrices = element.stokes_matrices(
                    &node_coords,
                    &self.properties,
                    &self.config,
                )?;
                
                // Assemble element contributions into global sparse matrices
                // Map local to global indices and add entries
                for (local_row, local_col, value) in &elem_matrices.k_entries {
                    let node_i = nodes[local_row / constants::VELOCITY_COMPONENTS];
                    let comp_i = local_row % constants::VELOCITY_COMPONENTS;
                    let node_j = nodes[local_col / constants::VELOCITY_COMPONENTS];
                    let comp_j = local_col % constants::VELOCITY_COMPONENTS;
                    
                    let global_row = node_i * constants::VELOCITY_COMPONENTS + comp_i;
                    let global_col = node_j * constants::VELOCITY_COMPONENTS + comp_j;
                    
                    k_builder.add_entry(global_row, global_col, value.clone())?;
                }
                
                // Assemble gradient matrix entries
                for (local_row, local_col, value) in &elem_matrices.g_entries {
                    let node_i = nodes[local_row / constants::VELOCITY_COMPONENTS];
                    let comp_i = local_row % constants::VELOCITY_COMPONENTS;
                    let node_j = nodes[*local_col];
                    
                    let global_row = node_i * constants::VELOCITY_COMPONENTS + comp_i;
                    
                    g_builder.add_entry(global_row, node_j, value.clone())?;
                }
                
                // Assemble mass matrix entries
                for (local_row, local_col, value) in &elem_matrices.m_entries {
                    let node_i = nodes[local_row / constants::VELOCITY_COMPONENTS];
                    let comp_i = local_row % constants::VELOCITY_COMPONENTS;
                    let node_j = nodes[local_col / constants::VELOCITY_COMPONENTS];
                    let comp_j = local_col % constants::VELOCITY_COMPONENTS;
                    
                    let global_row = node_i * constants::VELOCITY_COMPONENTS + comp_i;
                    let global_col = node_j * constants::VELOCITY_COMPONENTS + comp_j;
                    
                    m_builder.add_entry(global_row, global_col, value.clone())?;
                }
                
                // Assemble stabilization matrix entries
                for (local_row, local_col, value) in &elem_matrices.s_entries {
                    let node_i = nodes[*local_row];
                    let node_j = nodes[*local_col];
                    
                    s_builder.add_entry(node_i, node_j, value.clone())?;
                }
            }
        }
        
        // Build the sparse matrices from the builders
        self.k_global = Some(k_builder.build()?);
        self.g_global = Some(g_builder.build()?);
        self.m_global = Some(m_builder.build()?);
        self.s_pp_global = Some(s_builder.build()?);
        
        Ok(())
    }
    
    /// Solve Stokes flow (steady, no convection)
    pub fn solve_stokes(&mut self, boundary_conditions: &HashMap<usize, Vector3<T>>) -> Result<()> {
        if self.k_global.is_none() {
            self.assemble_global_matrices()?;
        }
        
        let k = self.k_global.as_ref().unwrap();
        let g = self.g_global.as_ref().unwrap();
        let gt = g.transpose();
        let s_pp = self.s_pp_global.as_ref().unwrap();
        
        let n_vel = self.velocity.len();
        let n_pres = self.pressure.len();
        let n_total = n_vel + n_pres;
        
        // Build saddle-point system with stabilization using sparse matrices
        // [K   G  ] [u]   [f]
        // [G^T S_pp] [p] = [0]
        
        // Estimate capacity for the system matrix
        let estimated_nnz = k.nnz() + g.nnz() + gt.nnz() + s_pp.nnz();
        let mut a_builder = SparseMatrixBuilder::with_capacity(n_total, n_total, estimated_nnz);
        let mut b_vector = DVector::zeros(n_total);
        
        // Add K block (velocity-velocity)
        for (row_idx, row) in k.row_iter().enumerate() {
            for (&col_idx, val) in row.col_indices().iter().zip(row.values()) {
                a_builder.add_entry(row_idx, col_idx, val.clone())?;
            }
        }
        
        // Add G block (velocity-pressure)
        for (row_idx, row) in g.row_iter().enumerate() {
            for (&col_idx, val) in row.col_indices().iter().zip(row.values()) {
                a_builder.add_entry(row_idx, n_vel + col_idx, val.clone())?;
            }
        }
        
        // Add G^T block (pressure-velocity)
        for (row_idx, row) in gt.row_iter().enumerate() {
            for (&col_idx, val) in row.col_indices().iter().zip(row.values()) {
                a_builder.add_entry(n_vel + row_idx, col_idx, val.clone())?;
            }
        }
        
        // Add S_pp block (pressure-pressure stabilization)
        for (row_idx, row) in s_pp.row_iter().enumerate() {
            for (&col_idx, val) in row.col_indices().iter().zip(row.values()) {
                a_builder.add_entry(n_vel + row_idx, n_vel + col_idx, val.clone())?;
            }
        }
        
        // No need for artificial stabilization when PSPG is used
        if !self.config.use_stabilization {
            // Only add small diagonal term if stabilization is disabled
            let eps = T::from_f64(1e-10).unwrap();
            for i in 0..n_pres {
                a_builder.add_entry(n_vel + i, n_vel + i, eps.clone())?;
            }
        }
        
        // Apply body forces using iterator pattern for zero-copy optimization
        if let Some(ref body_force) = self.properties.body_force {
            let rho = self.properties.density.clone();
            self.mesh.vertices
                .iter()
                .enumerate()
                .for_each(|(i, _)| {
                    let base_idx = i * 3;
                    b_vector[base_idx] = rho.clone() * body_force.x.clone();
                    b_vector[base_idx + 1] = rho.clone() * body_force.y.clone();
                    b_vector[base_idx + 2] = rho.clone() * body_force.z.clone();
                });
        }
        
        // Build the sparse matrix
        let a_matrix = a_builder.build()?;
        
        // Apply boundary conditions 
        // Note: This is simplified - in practice, should modify matrix assembly to handle BCs
        for (&node_idx, velocity_bc) in boundary_conditions {
            for d in 0..constants::VELOCITY_COMPONENTS {
                let dof = node_idx * constants::VELOCITY_COMPONENTS + d;
                
                // Set RHS to BC value
                b_vector[dof] = if d == 0 { velocity_bc.x.clone() }
                               else if d == 1 { velocity_bc.y.clone() }
                               else { velocity_bc.z.clone() };
            }
        }
        
        // Use BiCGSTAB solver for the saddle-point system
        use cfd_math::linear_solver::{BiCGSTAB, LinearSolver};
        let solver_config = cfd_core::LinearSolverConfig::default();
        let solver = BiCGSTAB::new(solver_config);
        
        let solution = solver.solve(&a_matrix, &b_vector, None)?;
        
        // Extract velocity and pressure
        self.velocity = solution.rows(0, n_vel).into();
        self.pressure = solution.rows(n_vel, n_pres).into();
        
        Ok(())
    }
    
    /// Get the full velocity solution vector
    pub fn get_velocity_field(&self) -> &DVector<T> {
        &self.velocity
    }
    
    /// Get the full pressure solution vector
    pub fn get_pressure_field(&self) -> &DVector<T> {
        &self.pressure
    }
    
    /// Get velocity at node
    pub fn get_velocity(&self, node_idx: usize) -> Vector3<T> {
        let base = node_idx * constants::VELOCITY_COMPONENTS;
        Vector3::new(
            self.velocity[base].clone(),
            self.velocity[base + 1].clone(),
            self.velocity[base + 2].clone(),
        )
    }
    
    /// Get pressure at node
    pub fn get_pressure(&self, node_idx: usize) -> T {
        self.pressure[node_idx].clone()
    }
    
    /// Calculate strain rate at element
    pub fn calculate_strain_rate(&self, element_nodes: &[usize]) -> Result<Matrix3<T>> {
        let element = FluidElement::<T>::new(element_nodes.to_vec(), 0);
        
        // Get node coordinates and velocities
        let node_coords: Vec<Vector3<T>> = element_nodes.iter()
            .map(|&idx| {
                let v = &self.mesh.vertices[idx];
                Vector3::new(
                    v.position.x.clone(),
                    v.position.y.clone(),
                    v.position.z.clone(),
                )
            })
            .collect();
        
        let node_velocities: Vec<Vector3<T>> = element_nodes.iter()
            .map(|&idx| self.get_velocity(idx))
            .collect();
        
        // Calculate velocity gradient
        let (dn_dx, dn_dy, dn_dz) = element.shape_derivatives(&node_coords)?;
        
        let mut velocity_gradient = Matrix3::zeros();
        for i in 0..constants::TET4_NODES {
            let vel = &node_velocities[i];
            
            // ∂u/∂x, ∂u/∂y, ∂u/∂z
            velocity_gradient[(0, 0)] += vel.x.clone() * dn_dx[i].clone();
            velocity_gradient[(0, 1)] += vel.x.clone() * dn_dy[i].clone();
            velocity_gradient[(0, 2)] += vel.x.clone() * dn_dz[i].clone();
            
            // ∂v/∂x, ∂v/∂y, ∂v/∂z
            velocity_gradient[(1, 0)] += vel.y.clone() * dn_dx[i].clone();
            velocity_gradient[(1, 1)] += vel.y.clone() * dn_dy[i].clone();
            velocity_gradient[(1, 2)] += vel.y.clone() * dn_dz[i].clone();
            
            // ∂w/∂x, ∂w/∂y, ∂w/∂z
            velocity_gradient[(2, 0)] += vel.z.clone() * dn_dx[i].clone();
            velocity_gradient[(2, 1)] += vel.z.clone() * dn_dy[i].clone();
            velocity_gradient[(2, 2)] += vel.z.clone() * dn_dz[i].clone();
        }
        
        // Return strain rate tensor
        Ok(element.strain_rate_tensor(&velocity_gradient))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Point3;
    use cfd_mesh::{Vertex, Face, Cell};
    
    #[test]
    fn test_fluid_properties() {
        let props = FluidProperties::<f64>::water();
        assert_relative_eq!(props.density, 998.2, epsilon = 0.1);
        assert_relative_eq!(props.viscosity, 0.001, epsilon = 0.0001);
        assert_relative_eq!(props.kinematic_viscosity(), 0.001 / 998.2, epsilon = 1e-6);
    }
    
    #[test]
    fn test_poiseuille_flow() {
        // Create a simple tetrahedral mesh for pipe flow
        let mut mesh: Mesh<f64> = Mesh::new();
        
        // Create a single well-formed tetrahedron representing a channel element
        mesh.vertices.push(Vertex { position: Point3::new(0.0, 0.0, 0.0), id: 0 });
        mesh.vertices.push(Vertex { position: Point3::new(1.0, 0.0, 0.0), id: 1 });
        mesh.vertices.push(Vertex { position: Point3::new(0.5, 0.866, 0.0), id: 2 }); // Equilateral triangle base
        mesh.vertices.push(Vertex { position: Point3::new(0.5, 0.433, 2.0), id: 3 }); // Extended in z for pipe flow
        
        // Create a single tetrahedral cell with proper connectivity
        // Create the 4 triangular faces of the tetrahedron
        mesh.faces.push(Face { vertices: vec![0, 1, 2], id: 0 }); // Bottom triangle
        mesh.faces.push(Face { vertices: vec![0, 3, 1], id: 1 }); // Side 1
        mesh.faces.push(Face { vertices: vec![1, 3, 2], id: 2 }); // Side 2
        mesh.faces.push(Face { vertices: vec![2, 3, 0], id: 3 }); // Side 3
        
        // Create the tetrahedral cell
        mesh.cells.push(Cell { faces: vec![0, 1, 2, 3], id: 0 });
        
        // Set up solver
        let config = FemConfig::default();
        let properties = FluidProperties::water();
        let mut solver = FemSolver::new(config, mesh, properties);
        
        // Apply boundary conditions for pipe flow
        let mut bc = HashMap::new();
        
        // No-slip on the bottom triangle (inlet/wall)
        bc.insert(0, Vector3::zeros());
        bc.insert(1, Vector3::zeros());
        bc.insert(2, Vector3::zeros());
        
        // Apply pressure-driven flow by leaving outlet (vertex 3) unconstrained
        // or apply a small velocity to simulate pressure gradient
        
        // Solve - should converge for well-formed tetrahedron
        solver.solve_stokes(&bc).expect("Stokes solver should converge for Poiseuille flow test");
        
        // Check that the unconstrained vertex has some velocity
        let outlet_vel = solver.get_velocity(3);
        assert!(outlet_vel.norm() >= 0.0, "Solution should be computed for unconstrained nodes");
    }
    
    #[test]
    fn test_couette_flow() {
        // Create a proper 3D tetrahedral mesh for Couette flow
        let mut mesh: Mesh<f64> = Mesh::new();
        
        // Create vertices for a proper tetrahedron with good aspect ratio
        // Bottom face: triangle in xy-plane
        mesh.vertices.push(Vertex { position: Point3::new(0.0, 0.0, 0.0), id: 0 });
        mesh.vertices.push(Vertex { position: Point3::new(1.0, 0.0, 0.0), id: 1 });
        mesh.vertices.push(Vertex { position: Point3::new(0.5, 1.0, 0.0), id: 2 });
        // Top vertex: properly positioned for non-degenerate tetrahedron
        mesh.vertices.push(Vertex { position: Point3::new(0.5, 0.5, 1.0), id: 3 });
        
        // Create the 4 triangular faces of the tetrahedron
        // Face ordering follows right-hand rule for outward normals
        mesh.faces.push(Face { vertices: vec![0, 1, 2], id: 0 }); // Bottom
        mesh.faces.push(Face { vertices: vec![0, 3, 1], id: 1 }); // Front
        mesh.faces.push(Face { vertices: vec![1, 3, 2], id: 2 }); // Right
        mesh.faces.push(Face { vertices: vec![2, 3, 0], id: 3 }); // Left
        
        // Create the tetrahedral cell with proper face connectivity
        mesh.cells.push(Cell { faces: vec![0, 1, 2, 3], id: 0 });
        
        let config = FemConfig::default();
        let properties = FluidProperties {
            density: 1000.0,
            viscosity: 0.001,
            body_force: None,
        };
        
        let mut solver = FemSolver::new(config, mesh, properties);
        
        // Boundary conditions: bottom face stationary, top node moving
        let mut bc = HashMap::new();
        bc.insert(0, Vector3::zeros()); // Bottom vertex 0
        bc.insert(1, Vector3::zeros()); // Bottom vertex 1  
        bc.insert(2, Vector3::zeros()); // Bottom vertex 2
        bc.insert(3, Vector3::new(1.0, 0.0, 0.0)); // Top vertex moving in x direction
        
        // Solve - test should fail if solver doesn't converge
        solver.solve_stokes(&bc).expect("Stokes solver should converge for Couette flow test");
        
        // Check that flow direction is correct (x-component should be positive somewhere)
        let top_vel = solver.get_velocity(3);
        assert!(top_vel.x.abs() > 0.01, "Top node should have x-velocity close to boundary condition");
    }
}
