//! Flow field operations
//!
//! Provides operations on flow fields including vorticity, divergence,
//! and other fluid mechanical quantities.

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use rayon::prelude::*;

use super::fields::VelocityField;

/// Operations on flow fields
pub struct FlowOperations;

impl FlowOperations {
    /// Calculate vorticity field (curl of velocity)
    pub fn vorticity<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity: &VelocityField<T>
    ) -> Vec<Vector3<T>> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut vorticity = vec![Vector3::zeros(); nx * ny * nz];
        
        // Use parallel iteration for better performance
        vorticity.par_iter_mut().enumerate().for_each(|(idx, vort)| {
            let k = idx / (nx * ny);
            let j = (idx % (nx * ny)) / nx;
            let i = idx % nx;
            
            *vort = vorticity_at_point(velocity, i, j, k, nx, ny, nz);
        });
        
        vorticity
    }
    
    /// Calculate divergence field
    pub fn divergence<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut divergence = vec![T::zero(); nx * ny * nz];
        
        // Use parallel iteration for better performance
        divergence.par_iter_mut().enumerate().for_each(|(idx, div)| {
            let k = idx / (nx * ny);
            let j = (idx % (nx * ny)) / nx;
            let i = idx % nx;
            
            let dx = T::one();
            let dy = T::one();
            let dz = T::one();
            
            // Central differences with boundary handling
            let dudx = if i > 0 && i < nx - 1 {
                let idx_plus = k * nx * ny + j * nx + (i + 1);
                let idx_minus = k * nx * ny + j * nx + (i - 1);
                (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / (T::from_f64(2.0).unwrap() * dx)
            } else {
                T::zero()
            };
            
            let dvdy = if j > 0 && j < ny - 1 {
                let idx_plus = k * nx * ny + (j + 1) * nx + i;
                let idx_minus = k * nx * ny + (j - 1) * nx + i;
                (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / (T::from_f64(2.0).unwrap() * dy)
            } else {
                T::zero()
            };
            
            let dwdz = if k > 0 && k < nz - 1 {
                let idx_plus = (k + 1) * nx * ny + j * nx + i;
                let idx_minus = (k - 1) * nx * ny + j * nx + i;
                (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / (T::from_f64(2.0).unwrap() * dz)
            } else {
                T::zero()
            };
            
            *div = dudx + dvdy + dwdz;
        });
        
        divergence
    }
    
    /// Calculate kinetic energy field
    pub fn kinetic_energy<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        // Use parallel iteration for better performance
        velocity.components
            .par_iter()
            .map(|v| T::from_f64(0.5).unwrap() * v.norm_squared())
            .collect()
    }
    
    /// Calculate enstrophy field
    pub fn enstrophy<T: RealField + Copy + FromPrimitive + Send + Sync>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        let vorticity = Self::vorticity(velocity);
        
        // Use parallel iteration for better performance
        vorticity
            .par_iter()
            .map(|w| T::from_f64(0.5).unwrap() * w.norm_squared())
            .collect()
    }
}

// Helper function for vorticity calculation at a point
fn vorticity_at_point<T: RealField + Copy + FromPrimitive>(
    velocity: &VelocityField<T>,
    i: usize, j: usize, k: usize,
    nx: usize, ny: usize, nz: usize
) -> Vector3<T> {
    let idx = k * nx * ny + j * nx + i;
    let dx = T::one();
    let dy = T::one();
    let dz = T::one();
    
    // Calculate velocity gradients using central differences
    let dvdz = if k > 0 && k < nz - 1 {
        let idx_plus = (k + 1) * nx * ny + j * nx + i;
        let idx_minus = (k - 1) * nx * ny + j * nx + i;
        (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / (T::from_f64(2.0).unwrap() * dz)
    } else {
        T::zero()
    };
    
    let dwdy = if j > 0 && j < ny - 1 {
        let idx_plus = k * nx * ny + (j + 1) * nx + i;
        let idx_minus = k * nx * ny + (j - 1) * nx + i;
        (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / (T::from_f64(2.0).unwrap() * dy)
    } else {
        T::zero()
    };
    
    let dudz = if k > 0 && k < nz - 1 {
        let idx_plus = (k + 1) * nx * ny + j * nx + i;
        let idx_minus = (k - 1) * nx * ny + j * nx + i;
        (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / (T::from_f64(2.0).unwrap() * dz)
    } else {
        T::zero()
    };
    
    let dwdx = if i > 0 && i < nx - 1 {
        let idx_plus = k * nx * ny + j * nx + (i + 1);
        let idx_minus = k * nx * ny + j * nx + (i - 1);
        (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / (T::from_f64(2.0).unwrap() * dx)
    } else {
        T::zero()
    };
    
    let dvdx = if i > 0 && i < nx - 1 {
        let idx_plus = k * nx * ny + j * nx + (i + 1);
        let idx_minus = k * nx * ny + j * nx + (i - 1);
        (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / (T::from_f64(2.0).unwrap() * dx)
    } else {
        T::zero()
    };
    
    let dudy = if j > 0 && j < ny - 1 {
        let idx_plus = k * nx * ny + (j + 1) * nx + i;
        let idx_minus = k * nx * ny + (j - 1) * nx + i;
        (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / (T::from_f64(2.0).unwrap() * dy)
    } else {
        T::zero()
    };
    
    // Vorticity = curl(velocity)
    Vector3::new(
        dwdy - dvdz,
        dudz - dwdx,
        dvdx - dudy
    )
}