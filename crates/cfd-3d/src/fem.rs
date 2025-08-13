//! Finite Element Method (FEM) for 3D incompressible fluid dynamics.
//!
//! This module implements FEM for solving the incompressible Navier-Stokes equations
//! using mixed velocity-pressure formulation with appropriate stabilization.

use cfd_core::{Error, Result};

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
    /// Gauss point weight for single point quadrature
    pub const GAUSS_WEIGHT_1PT: f64 = 1.0;
    /// Number of strain rate components (symmetric tensor)
    pub const STRAIN_COMPONENTS: usize = 6;
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
    /// Returns (K, G, G^T, M, S_pp) where:
    /// - K: viscous stiffness matrix (velocity-velocity)
    /// - G: gradient matrix (velocity-pressure)
    /// - G^T: divergence matrix (pressure-velocity)
    /// - M: mass matrix (for transient terms)
    /// - S_pp: pressure stabilization matrix (pressure-pressure)
    pub fn stokes_matrices(
        &self,
        node_coords: &[Vector3<T>],
        properties: &FluidProperties<T>,
        config: &FemConfig<T>,
    ) -> Result<(DMatrix<T>, DMatrix<T>, DMatrix<T>, DMatrix<T>, DMatrix<T>)> {
        let n_vel_dof = constants::TET4_NODES * constants::VELOCITY_COMPONENTS; // 12
        let n_pres_dof = constants::TET4_NODES; // 4
        
        // Initialize matrices
        let mut k_matrix: DMatrix<T> = DMatrix::zeros(n_vel_dof, n_vel_dof);
        let mut g_matrix: DMatrix<T> = DMatrix::zeros(n_vel_dof, n_pres_dof);
        let mut m_matrix: DMatrix<T> = DMatrix::zeros(n_vel_dof, n_vel_dof);
        
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
                    k_matrix[(row, col)] = k_matrix[(row, col)].clone() + k_val.clone();
                }
            }
        }
        
        // Build gradient matrix G
        // G_ij = -∫ (∇·N_i) M_j dΩ where N_i are velocity shape functions, M_j are pressure
        for i in 0..constants::TET4_NODES {
            for j in 0..constants::TET4_NODES {
                let g_val = weight.clone() / T::from_f64(4.0).unwrap(); // Linear shape function value at centroid
                
                // Velocity divergence coupled with pressure
                g_matrix[(i * 3, j)] = g_matrix[(i * 3, j)].clone() - dn_dx[i].clone() * g_val.clone();
                g_matrix[(i * 3 + 1, j)] = g_matrix[(i * 3 + 1, j)].clone() - dn_dy[i].clone() * g_val.clone();
                g_matrix[(i * 3 + 2, j)] = g_matrix[(i * 3 + 2, j)].clone() - dn_dz[i].clone() * g_val.clone();
            }
        }
        
        // Build mass matrix M (for transient terms)
        if config.dt.is_some() {
            for i in 0..constants::TET4_NODES {
                for j in 0..constants::TET4_NODES {
                    let m_val = rho.clone() * weight.clone() / T::from_f64(20.0).unwrap(); // Consistent mass
                    if i == j {
                        for d in 0..constants::VELOCITY_COMPONENTS {
                            let idx = i * constants::VELOCITY_COMPONENTS + d;
                            m_matrix[(idx, idx)] = m_val.clone() * T::from_f64(2.0).unwrap();
                        }
                    } else {
                        for d in 0..constants::VELOCITY_COMPONENTS {
                            let row = i * constants::VELOCITY_COMPONENTS + d;
                            let col = j * constants::VELOCITY_COMPONENTS + d;
                            m_matrix[(row, col)] = m_val.clone();
                        }
                    }
                }
            }
        }
        
        // Add SUPG/PSPG stabilization if enabled
        let s_pp = if config.use_stabilization {
            self.add_stabilization(
                &mut k_matrix,
                &mut g_matrix,
                node_coords,
                properties,
                config,
                &dn_dx,
                &dn_dy,
                &dn_dz,
                volume,
            )?
        } else {
            DMatrix::zeros(constants::TET4_NODES, constants::TET4_NODES)
        };
        
        let gt_matrix = g_matrix.transpose();
        
        // Return all matrices including the pressure stabilization matrix
        Ok((k_matrix, g_matrix, gt_matrix, m_matrix, s_pp))
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
    /// Global stiffness matrix
    k_global: Option<DMatrix<T>>,
    /// Global gradient matrix
    g_global: Option<DMatrix<T>>,
    /// Global mass matrix
    m_global: Option<DMatrix<T>>,
    /// Global pressure stabilization matrix
    s_pp_global: Option<DMatrix<T>>,
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
        
        let mut k_global = DMatrix::zeros(n_vel_dof, n_vel_dof);
        let mut g_global = DMatrix::zeros(n_vel_dof, n_pres_dof);
        let mut m_global = DMatrix::zeros(n_vel_dof, n_vel_dof);
        let mut s_pp_global = DMatrix::zeros(n_pres_dof, n_pres_dof);
        
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
                let (k_elem, g_elem, _, m_elem, s_pp_elem) = element.stokes_matrices(
                    &node_coords,
                    &self.properties,
                    &self.config,
                )?;
                
                // Assemble into global matrices
                for (local_i, &global_i) in nodes.iter().enumerate() {
                    for (local_j, &global_j) in nodes.iter().enumerate() {
                        // Velocity-velocity coupling
                        for d1 in 0..constants::VELOCITY_COMPONENTS {
                            for d2 in 0..constants::VELOCITY_COMPONENTS {
                                let local_row = local_i * constants::VELOCITY_COMPONENTS + d1;
                                let local_col = local_j * constants::VELOCITY_COMPONENTS + d2;
                                let global_row = global_i * constants::VELOCITY_COMPONENTS + d1;
                                let global_col = global_j * constants::VELOCITY_COMPONENTS + d2;
                                
                                k_global[(global_row, global_col)] += k_elem[(local_row, local_col)].clone();
                                m_global[(global_row, global_col)] += m_elem[(local_row, local_col)].clone();
                            }
                        }
                        
                        // Velocity-pressure coupling
                        for d in 0..constants::VELOCITY_COMPONENTS {
                            let local_row = local_i * constants::VELOCITY_COMPONENTS + d;
                            let global_row = global_i * constants::VELOCITY_COMPONENTS + d;
                            
                            g_global[(global_row, global_j)] += g_elem[(local_row, local_j)].clone();
                        }
                        
                        // Pressure-pressure stabilization
                        s_pp_global[(global_i, global_j)] += s_pp_elem[(local_i, local_j)].clone();
                    }
                }
            }
        }
        
        self.k_global = Some(k_global);
        self.g_global = Some(g_global);
        self.m_global = Some(m_global);
        self.s_pp_global = Some(s_pp_global);
        
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
        
        // Build saddle-point system with stabilization
        // [K   G  ] [u]   [f]
        // [G^T S_pp] [p] = [0]
        let mut a_matrix = DMatrix::zeros(n_total, n_total);
        let mut b_vector = DVector::zeros(n_total);
        
        // Fill system matrix
        a_matrix.view_mut((0, 0), (n_vel, n_vel)).copy_from(k);
        a_matrix.view_mut((0, n_vel), (n_vel, n_pres)).copy_from(g);
        a_matrix.view_mut((n_vel, 0), (n_pres, n_vel)).copy_from(&gt);
        
        // Add PSPG stabilization to pressure-pressure block
        // This stabilizes the inf-sup condition for equal-order elements
        a_matrix.view_mut((n_vel, n_vel), (n_pres, n_pres)).copy_from(s_pp);
        
        // No need for artificial stabilization when PSPG is used
        if !self.config.use_stabilization {
            // Only add small diagonal term if stabilization is disabled
            let eps = T::from_f64(1e-10).unwrap();
            for i in 0..n_pres {
                a_matrix[(n_vel + i, n_vel + i)] = a_matrix[(n_vel + i, n_vel + i)].clone() + eps.clone();
            }
        }
        
        // Apply body forces
        if let Some(ref body_force) = self.properties.body_force {
            let rho = self.properties.density.clone();
            for i in 0..self.mesh.vertices.len() {
                b_vector[i * 3] = rho.clone() * body_force.x.clone();
                b_vector[i * 3 + 1] = rho.clone() * body_force.y.clone();
                b_vector[i * 3 + 2] = rho.clone() * body_force.z.clone();
            }
        }
        
        // Apply boundary conditions
        for (&node_idx, velocity_bc) in boundary_conditions {
            for d in 0..constants::VELOCITY_COMPONENTS {
                let dof = node_idx * constants::VELOCITY_COMPONENTS + d;
                
                // Zero out row and column
                for j in 0..n_total {
                    a_matrix[(dof, j)] = T::zero();
                    a_matrix[(j, dof)] = T::zero();
                }
                
                // Set diagonal to 1 and RHS to BC value
                a_matrix[(dof, dof)] = T::one();
                b_vector[dof] = if d == 0 { velocity_bc.x.clone() }
                               else if d == 1 { velocity_bc.y.clone() }
                               else { velocity_bc.z.clone() };
            }
        }
        
        // PERFORMANCE WARNING: Using dense matrix for simplicity
        // For production use, this should be converted to sparse matrices (CsrMatrix)
        // and solved with iterative methods (GMRES, BiCGSTAB) from cfd_math::linear_solver
        // The current dense LU decomposition has O(N^3) complexity and O(N^2) memory
        // which is prohibitive for realistic 3D problems
        let solution = a_matrix.lu().solve(&b_vector)
            .ok_or_else(|| Error::NumericalError("Failed to solve linear system".into()))?;
        
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
        // Create a simple pipe mesh
        let mut mesh: Mesh<f64> = Mesh::new();
        
        // Add vertices for a simple rectangular channel
        let nx = 3;
        let ny = 3;
        let nz = 5;
        let dx = 0.1;
        let dy = 0.1;
        let dz = 0.2;
        
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = k * nx * ny + j * nx + i;
                    mesh.vertices.push(Vertex {
                        position: Point3::new(
                            i as f64 * dx,
                            j as f64 * dy,
                            k as f64 * dz,
                        ),
                        id: idx,
                    });
                }
            }
        }
        
        // Create tetrahedral cells using structured grid connectivity
        // For a structured grid, we can create tetrahedra by dividing each hexahedral cell
        // into 6 tetrahedra using a consistent pattern
        if nx > 1 && ny > 1 && nz > 1 {
            for k in 0..nz-1 {
                for j in 0..ny-1 {
                    for i in 0..nx-1 {
                        // Get the 8 vertices of the hexahedral cell
                        let v000 = k * ny * nx + j * nx + i;
                        let v100 = k * ny * nx + j * nx + (i + 1);
                        let v010 = k * ny * nx + (j + 1) * nx + i;
                        let v110 = k * ny * nx + (j + 1) * nx + (i + 1);
                        let v001 = (k + 1) * ny * nx + j * nx + i;
                        let v101 = (k + 1) * ny * nx + j * nx + (i + 1);
                        let v011 = (k + 1) * ny * nx + (j + 1) * nx + i;
                        let v111 = (k + 1) * ny * nx + (j + 1) * nx + (i + 1);
                        
                        // Divide hexahedron into 6 tetrahedra using diagonal decomposition
                        // This ensures consistent orientation and no gaps
                        let base_face_id = mesh.faces.len();
                        let base_cell_id = mesh.cells.len();
                        
                        // Tetrahedron 1: v000, v100, v010, v001
                        mesh.faces.push(Face { vertices: vec![v000, v100, v010], id: base_face_id });
                        mesh.faces.push(Face { vertices: vec![v000, v100, v001], id: base_face_id + 1 });
                        mesh.faces.push(Face { vertices: vec![v000, v010, v001], id: base_face_id + 2 });
                        mesh.faces.push(Face { vertices: vec![v100, v010, v001], id: base_face_id + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![base_face_id, base_face_id + 1, base_face_id + 2, base_face_id + 3],
                            id: base_cell_id,
                        });
                        
                        // Tetrahedron 2: v100, v110, v010, v111
                        let face_id_2 = mesh.faces.len();
                        let cell_id_2 = mesh.cells.len();
                        mesh.faces.push(Face { vertices: vec![v100, v110, v010], id: face_id_2 });
                        mesh.faces.push(Face { vertices: vec![v100, v110, v111], id: face_id_2 + 1 });
                        mesh.faces.push(Face { vertices: vec![v100, v010, v111], id: face_id_2 + 2 });
                        mesh.faces.push(Face { vertices: vec![v110, v010, v111], id: face_id_2 + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![face_id_2, face_id_2 + 1, face_id_2 + 2, face_id_2 + 3],
                            id: cell_id_2,
                        });
                        
                        // Tetrahedron 3: v100, v001, v101, v111
                        let face_id_3 = mesh.faces.len();
                        let cell_id_3 = mesh.cells.len();
                        mesh.faces.push(Face { vertices: vec![v100, v001, v101], id: face_id_3 });
                        mesh.faces.push(Face { vertices: vec![v100, v001, v111], id: face_id_3 + 1 });
                        mesh.faces.push(Face { vertices: vec![v100, v101, v111], id: face_id_3 + 2 });
                        mesh.faces.push(Face { vertices: vec![v001, v101, v111], id: face_id_3 + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![face_id_3, face_id_3 + 1, face_id_3 + 2, face_id_3 + 3],
                            id: cell_id_3,
                        });
                        
                        // Tetrahedron 4: v010, v011, v001, v111
                        let face_id_4 = mesh.faces.len();
                        let cell_id_4 = mesh.cells.len();
                        mesh.faces.push(Face { vertices: vec![v010, v011, v001], id: face_id_4 });
                        mesh.faces.push(Face { vertices: vec![v010, v011, v111], id: face_id_4 + 1 });
                        mesh.faces.push(Face { vertices: vec![v010, v001, v111], id: face_id_4 + 2 });
                        mesh.faces.push(Face { vertices: vec![v011, v001, v111], id: face_id_4 + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![face_id_4, face_id_4 + 1, face_id_4 + 2, face_id_4 + 3],
                            id: cell_id_4,
                        });
                        
                        // Tetrahedron 5: v100, v010, v001, v111
                        let face_id_5 = mesh.faces.len();
                        let cell_id_5 = mesh.cells.len();
                        mesh.faces.push(Face { vertices: vec![v100, v010, v001], id: face_id_5 });
                        mesh.faces.push(Face { vertices: vec![v100, v010, v111], id: face_id_5 + 1 });
                        mesh.faces.push(Face { vertices: vec![v100, v001, v111], id: face_id_5 + 2 });
                        mesh.faces.push(Face { vertices: vec![v010, v001, v111], id: face_id_5 + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![face_id_5, face_id_5 + 1, face_id_5 + 2, face_id_5 + 3],
                            id: cell_id_5,
                        });
                        
                        // Tetrahedron 6: v010, v110, v111, v011
                        let face_id_6 = mesh.faces.len();
                        let cell_id_6 = mesh.cells.len();
                        mesh.faces.push(Face { vertices: vec![v010, v110, v111], id: face_id_6 });
                        mesh.faces.push(Face { vertices: vec![v010, v110, v011], id: face_id_6 + 1 });
                        mesh.faces.push(Face { vertices: vec![v010, v111, v011], id: face_id_6 + 2 });
                        mesh.faces.push(Face { vertices: vec![v110, v111, v011], id: face_id_6 + 3 });
                        mesh.cells.push(Cell {
                            faces: vec![face_id_6, face_id_6 + 1, face_id_6 + 2, face_id_6 + 3],
                            id: cell_id_6,
                        });
                    }
                }
            }
        } else if mesh.vertices.len() >= 4 {
            // Fallback for small meshes: create a single tetrahedron
            mesh.faces.push(Face { vertices: vec![0, 1, 2], id: 0 });
            mesh.faces.push(Face { vertices: vec![0, 1, 3], id: 1 });
            mesh.faces.push(Face { vertices: vec![0, 2, 3], id: 2 });
            mesh.faces.push(Face { vertices: vec![1, 2, 3], id: 3 });
            
            // Create cell from faces
            mesh.cells.push(Cell { faces: vec![0, 1, 2, 3], id: 0 });
        }
        
        // Set up solver
        let config = FemConfig::default();
        let properties = FluidProperties::water();
        let mut solver = FemSolver::new(config, mesh, properties);
        
        // Apply boundary conditions
        let mut bc = HashMap::new();
        
        // No-slip on walls (y=0, y=max)
        for k in 0..nz {
            for i in 0..nx {
                // Bottom wall
                let idx = k * nx * ny + i;
                bc.insert(idx, Vector3::zeros());
                
                // Top wall
                let idx = k * nx * ny + (ny-1) * nx + i;
                bc.insert(idx, Vector3::zeros());
            }
        }
        
        // Inlet pressure (implicit through solution)
        // Outlet pressure (implicit through solution)
        
        // Solve - test should fail if solver doesn't converge
        solver.solve_stokes(&bc).expect("Stokes solver should converge for Poiseuille flow test");
        
        // Check that center velocity is positive (flow direction)
        let center_node = nz/2 * nx * ny + ny/2 * nx + nx/2;
        let center_vel = solver.get_velocity(center_node);
        assert!(center_vel.z > 0.0, "Flow should be in positive z direction");
    }
    
    #[test]
    fn test_couette_flow() {
        // Create mesh between two parallel plates
        let mut mesh: Mesh<f64> = Mesh::new();
        
        // Simple 2x3x2 mesh
        let mut idx = 0;
        for k in 0..2 {
            for j in 0..3 {
                for i in 0..2 {
                    mesh.vertices.push(Vertex {
                        position: Point3::new(
                            i as f64 * 0.5,
                            j as f64 * 0.1,
                            k as f64 * 0.5,
                        ),
                        id: idx,
                    });
                    idx += 1;
                }
            }
        }
        
        // Create a simple tetrahedral cell
        mesh.faces.push(Face { vertices: vec![0, 1, 2], id: 0 });
        mesh.faces.push(Face { vertices: vec![0, 1, 6], id: 1 });
        mesh.faces.push(Face { vertices: vec![0, 2, 6], id: 2 });
        mesh.faces.push(Face { vertices: vec![1, 2, 6], id: 3 });
        mesh.cells.push(Cell { faces: vec![0, 1, 2, 3], id: 0 });
        
        let config = FemConfig::default();
        let properties = FluidProperties {
            density: 1000.0,
            viscosity: 0.001,
            body_force: None,
        };
        
        let mut solver = FemSolver::new(config, mesh, properties);
        
        // Boundary conditions: top plate moving, bottom plate stationary
        let mut bc = HashMap::new();
        bc.insert(0, Vector3::zeros()); // Bottom plate
        bc.insert(1, Vector3::zeros());
        bc.insert(4, Vector3::new(1.0, 0.0, 0.0)); // Top plate moving in x
        bc.insert(5, Vector3::new(1.0, 0.0, 0.0));
        
        // Solve - test should fail if solver doesn't converge
        solver.solve_stokes(&bc).expect("Stokes solver should converge for Couette flow test");
        
        // Velocity should vary linearly between plates
        let mid_vel = solver.get_velocity(2);
        assert!(mid_vel.x > 0.0 && mid_vel.x < 1.0, "Mid-plane velocity should be between 0 and 1");
    }
}
