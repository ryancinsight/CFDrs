//! Immersed Boundary Method (IBM) for complex geometries in 3D
//!
//! The IBM allows simulation of flow around complex objects without
//! body-fitted meshes by using forcing terms in the momentum equations.

use cfd_core::{Result, Error};
use nalgebra::{Vector3, DMatrix, DVector, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Named constants for IBM
const DEFAULT_SMOOTHING_WIDTH: f64 = 1.5;  // Delta function width in grid cells
const DEFAULT_FORCE_SCALE: f64 = 1.0;
const DEFAULT_MAX_ITERATIONS: usize = 100;
const DEFAULT_TOLERANCE: f64 = 1e-6;
const DELTA_FUNCTION_CUTOFF: f64 = 4.0;  // Support radius for delta function
const INTERPOLATION_STENCIL_SIZE: usize = 4;  // 4x4x4 stencil for interpolation

/// IBM configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IbmConfig {
    /// Width of smoothing kernel (in grid cells)
    pub smoothing_width: f64,
    /// Force scaling factor
    pub force_scale: f64,
    /// Maximum iterations for force calculation
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Use direct forcing method
    pub use_direct_forcing: bool,
}

impl Default for IbmConfig {
    fn default() -> Self {
        Self {
            smoothing_width: DEFAULT_SMOOTHING_WIDTH,
            force_scale: DEFAULT_FORCE_SCALE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            tolerance: DEFAULT_TOLERANCE,
            use_direct_forcing: true,
        }
    }
}

/// Lagrangian marker point for immersed boundary
#[derive(Debug, Clone)]
pub struct LagrangianPoint<T: RealField> {
    /// Position in 3D space
    pub position: Vector3<T>,
    /// Velocity at this point
    pub velocity: Vector3<T>,
    /// Force applied at this point
    pub force: Vector3<T>,
    /// Surface normal (if on boundary)
    pub normal: Option<Vector3<T>>,
    /// Reference position (for elastic boundaries)
    pub reference_position: Option<Vector3<T>>,
}

/// Immersed Boundary Method solver
pub struct IbmSolver<T: RealField + FromPrimitive> {
    config: IbmConfig,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    dz: T,
    /// Lagrangian points representing the immersed boundary
    lagrangian_points: Vec<LagrangianPoint<T>>,
    /// Eulerian force field
    force_field: Vec<Vector3<T>>,
    /// Velocity field (reference to external field)
    velocity_field: Vec<Vector3<T>>,
}

impl<T: RealField + FromPrimitive> IbmSolver<T> {
    /// Create a new IBM solver
    pub fn new(
        config: IbmConfig,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
    ) -> Self {
        let grid_size = nx * ny * nz;
        Self {
            config,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            lagrangian_points: Vec::new(),
            force_field: vec![Vector3::zeros(); grid_size],
            velocity_field: vec![Vector3::zeros(); grid_size],
        }
    }
    
    /// Add Lagrangian points from a surface mesh
    pub fn add_surface_mesh(&mut self, vertices: &[Vector3<T>], normals: Option<&[Vector3<T>]>) {
        self.lagrangian_points.clear();
        
        for (i, vertex) in vertices.iter().enumerate() {
            let normal = normals.map(|n| n[i].clone());
            self.lagrangian_points.push(LagrangianPoint {
                position: vertex.clone(),
                velocity: Vector3::zeros(),
                force: Vector3::zeros(),
                normal,
                reference_position: Some(vertex.clone()),
            });
        }
    }
    
    /// Discrete delta function (Roma et al., 1999)
    fn delta_function(&self, r: T) -> T {
        let h = self.dx.clone();  // Assuming uniform grid
        let r_norm = r.abs() / h.clone();
        
        if r_norm >= T::from_f64(2.0).unwrap() {
            T::zero()
        } else if r_norm < T::one() {
            let one = T::one();
            let three = T::from_f64(3.0).unwrap();
            (one.clone() + ComplexField::sqrt(one.clone() - three.clone() * r_norm.clone() * r_norm.clone())) / (three * h)
        } else {
            let one = T::one();
            let two = T::from_f64(2.0).unwrap();
            let three = T::from_f64(3.0).unwrap();
            let five = T::from_f64(5.0).unwrap();
            // For 1 <= r_norm < 2, use the correct formula
            let term = (one.clone() - r_norm.clone()) * (one - r_norm.clone());
            (five - three.clone() * r_norm - ComplexField::sqrt(three * term)) / (T::from_f64(6.0).unwrap() * h)
        }
    }
    
    /// Interpolate velocity from Eulerian grid to Lagrangian points
    pub fn interpolate_velocity(&mut self) {
        let dx = self.dx.clone();
        let dy = self.dy.clone();
        let dz = self.dz.clone();
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let velocity_field = self.velocity_field.clone();
        
        for point in &mut self.lagrangian_points {
            let mut vel = Vector3::zeros();
            
            // Find grid cell containing the point
            let i_center = if let Some(val) = (point.position[0].clone() / dx.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(nx - 1)
            } else {
                0
            };
            let j_center = if let Some(val) = (point.position[1].clone() / dy.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(ny - 1)
            } else {
                0
            };
            let k_center = if let Some(val) = (point.position[2].clone() / dz.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(nz - 1)
            } else {
                0
            };
            
            // Use stencil for interpolation
            let stencil_half = INTERPOLATION_STENCIL_SIZE / 2;
            
            for di in 0..INTERPOLATION_STENCIL_SIZE {
                for dj in 0..INTERPOLATION_STENCIL_SIZE {
                    for dk in 0..INTERPOLATION_STENCIL_SIZE {
                        let i = (i_center + di).saturating_sub(stencil_half).min(nx - 1);
                        let j = (j_center + dj).saturating_sub(stencil_half).min(ny - 1);
                        let k = (k_center + dk).saturating_sub(stencil_half).min(nz - 1);
                        
                        let grid_pos = Vector3::new(
                            T::from_usize(i).unwrap() * dx.clone(),
                            T::from_usize(j).unwrap() * dy.clone(),
                            T::from_usize(k).unwrap() * dz.clone(),
                        );
                        
                        let dx_val = point.position[0].clone() - grid_pos[0].clone();
                        let dy_val = point.position[1].clone() - grid_pos[1].clone();
                        let dz_val = point.position[2].clone() - grid_pos[2].clone();
                        
                        // Inline delta function calculation
                        let weight_x = Self::delta_function_static(dx_val, dx.clone());
                        let weight_y = Self::delta_function_static(dy_val, dy.clone());
                        let weight_z = Self::delta_function_static(dz_val, dz.clone());
                        
                        let weight = weight_x * weight_y * weight_z
                            * dx.clone() * dy.clone() * dz.clone();
                        
                        let idx = k * nx * ny + j * nx + i;
                        vel += velocity_field[idx].clone() * weight;
                    }
                }
            }
            
            point.velocity = vel;
        }
    }
    
    /// Static version of delta function for use in mutable contexts
    fn delta_function_static(r: T, h: T) -> T {
        let r_norm = r.abs() / h.clone();
        
        if r_norm >= T::from_f64(2.0).unwrap() {
            T::zero()
        } else if r_norm < T::one() {
            let one = T::one();
            let three = T::from_f64(3.0).unwrap();
            (one.clone() + ComplexField::sqrt(one.clone() - three.clone() * r_norm.clone() * r_norm.clone())) / (three * h)
        } else {
            let one = T::one();
            let three = T::from_f64(3.0).unwrap();
            let five = T::from_f64(5.0).unwrap();
            // For 1 <= r_norm < 2, use the correct formula  
            let term = (one.clone() - r_norm.clone()) * (one - r_norm.clone());
            (five - three.clone() * r_norm - ComplexField::sqrt(three * term)) / (T::from_f64(6.0).unwrap() * h)
        }
    }
    
    /// Calculate forces at Lagrangian points
    pub fn calculate_forces(&mut self, target_velocity: Option<&[Vector3<T>]>) {
        if self.config.use_direct_forcing {
            // Direct forcing: F = (u_target - u_interpolated) / dt
            for (i, point) in self.lagrangian_points.iter_mut().enumerate() {
                let u_target = if let Some(targets) = target_velocity {
                    targets[i].clone()
                } else {
                    Vector3::zeros()  // No-slip condition
                };
                
                // Force to enforce boundary condition
                point.force = (u_target - point.velocity.clone()) * T::from_f64(self.config.force_scale).unwrap();
            }
        } else {
            // Feedback forcing or elastic boundary
            for point in &mut self.lagrangian_points {
                if let Some(ref_pos) = &point.reference_position {
                    // Elastic force: F = -k * (X - X0)
                    let displacement = point.position.clone() - ref_pos.clone();
                    let stiffness = T::from_f64(100.0).unwrap();  // Spring stiffness
                    point.force = displacement * (-stiffness);
                }
            }
        }
    }
    
    /// Spread forces from Lagrangian points to Eulerian grid
    pub fn spread_forces(&mut self) {
        // Clear force field
        self.force_field.iter_mut().for_each(|f| *f = Vector3::zeros());
        
        for point in &self.lagrangian_points {
            // Find grid cell containing the point
            let i_center = if let Some(val) = (point.position[0].clone() / self.dx.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(self.nx - 1)
            } else {
                0
            };
            let j_center = if let Some(val) = (point.position[1].clone() / self.dy.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(self.ny - 1)
            } else {
                0
            };
            let k_center = if let Some(val) = (point.position[2].clone() / self.dz.clone()).floor().to_subset() {
                let val_f64: f64 = val;
                (val_f64 as usize).min(self.nz - 1)
            } else {
                0
            };
            
            // Spread force using delta function
            let stencil_half = INTERPOLATION_STENCIL_SIZE / 2;
            
            for di in 0..INTERPOLATION_STENCIL_SIZE {
                for dj in 0..INTERPOLATION_STENCIL_SIZE {
                    for dk in 0..INTERPOLATION_STENCIL_SIZE {
                        let i = (i_center + di).saturating_sub(stencil_half).min(self.nx - 1);
                        let j = (j_center + dj).saturating_sub(stencil_half).min(self.ny - 1);
                        let k = (k_center + dk).saturating_sub(stencil_half).min(self.nz - 1);
                        
                        let grid_pos = Vector3::new(
                            T::from_usize(i).unwrap() * self.dx.clone(),
                            T::from_usize(j).unwrap() * self.dy.clone(),
                            T::from_usize(k).unwrap() * self.dz.clone(),
                        );
                        
                        let dx = point.position[0].clone() - grid_pos[0].clone();
                        let dy = point.position[1].clone() - grid_pos[1].clone();
                        let dz = point.position[2].clone() - grid_pos[2].clone();
                        
                        let weight = self.delta_function(dx) 
                            * self.delta_function(dy)
                            * self.delta_function(dz);
                        
                        let idx = k * self.nx * self.ny + j * self.nx + i;
                        self.force_field[idx] = self.force_field[idx].clone() + point.force.clone() * weight;
                    }
                }
            }
        }
    }
    
    /// Update Lagrangian point positions (for moving boundaries)
    pub fn update_positions(&mut self, dt: T) {
        for point in &mut self.lagrangian_points {
            point.position += point.velocity.clone() * dt.clone();
        }
    }
    
    /// Main IBM step
    pub fn step(&mut self, velocity_field: &[Vector3<T>], dt: T) -> Result<()> {
        // Update velocity field reference
        self.velocity_field.clear();
        self.velocity_field.extend_from_slice(velocity_field);
        
        // Interpolate velocity to Lagrangian points
        self.interpolate_velocity();
        
        // Calculate forces at Lagrangian points
        self.calculate_forces(None);
        
        // Spread forces to Eulerian grid
        self.spread_forces();
        
        // Update positions if boundaries are moving
        if !self.config.use_direct_forcing {
            self.update_positions(dt);
        }
        
        Ok(())
    }
    
    /// Get the force field
    pub fn force_field(&self) -> &[Vector3<T>] {
        &self.force_field
    }
    
    /// Get Lagrangian points
    pub fn lagrangian_points(&self) -> &[LagrangianPoint<T>] {
        &self.lagrangian_points
    }
    
    /// Compute drag and lift forces on the immersed body
    pub fn compute_forces(&self) -> (T, T, T) {
        let mut drag = T::zero();
        let mut lift = T::zero();
        let mut side = T::zero();
        
        for point in &self.lagrangian_points {
            drag += point.force[0].clone();
            lift += point.force[1].clone();
            side += point.force[2].clone();
        }
        
        (drag, lift, side)
    }
}

use nalgebra::ComplexField;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_ibm_creation() {
        let config = IbmConfig::default();
        let solver: IbmSolver<f64> = IbmSolver::new(
            config,
            10, 10, 10,
            0.1, 0.1, 0.1,
        );
        
        assert_eq!(solver.force_field.len(), 1000);
    }
    
    #[test]
    fn test_delta_function() {
        let config = IbmConfig::default();
        let solver: IbmSolver<f64> = IbmSolver::new(
            config,
            10, 10, 10,
            1.0, 1.0, 1.0,
        );
        
        // Delta function should be maximum at r=0
        let d0 = solver.delta_function(0.0);
        let d1 = solver.delta_function(1.0);
        // The delta function should decrease as we move away from 0
        // For a grid spacing of 1.0, the value at r=1.0 should be less than at r=0
        assert!(d0 > d1, "d0={} should be > d1={}", d0, d1);
        
        // Should be zero beyond cutoff
        let d_far = solver.delta_function(3.0);
        assert_eq!(d_far, 0.0);
    }
    
    #[test]
    fn test_force_spreading() {
        let config = IbmConfig::default();
        let mut solver: IbmSolver<f64> = IbmSolver::new(
            config,
            10, 10, 10,
            1.0, 1.0, 1.0,
        );
        
        // Add a single Lagrangian point
        solver.lagrangian_points.push(LagrangianPoint {
            position: Vector3::new(5.0, 5.0, 5.0),
            velocity: Vector3::zeros(),
            force: Vector3::new(1.0, 0.0, 0.0),
            normal: None,
            reference_position: None,
        });
        
        solver.spread_forces();
        
        // Force should be spread to nearby grid points
        let center_idx = 5 * 100 + 5 * 10 + 5;
        assert!(solver.force_field[center_idx][0] > 0.0);
    }
}