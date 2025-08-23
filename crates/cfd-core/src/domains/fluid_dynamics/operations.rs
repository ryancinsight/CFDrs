//! Flow field operations
//!
//! Provides operations on flow fields including vorticity, divergence,
//! and other fluid mechanical quantities.

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use super::fields::VelocityField;

/// Flow field operations
pub struct FlowOperations;

impl FlowOperations {
    /// Calculate vorticity field: ω = ∇ × u
    pub fn vorticity<T: RealField + Copy + FromPrimitive>(
        velocity: &VelocityField<T>
    ) -> Vec<Vector3<T>> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut vorticity = Vec::with_capacity(velocity.components.len());
        
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let omega = Self::vorticity_at_point(velocity, i, j, k);
                    vorticity.push(omega);
                }
            }
        }
        
        vorticity
    }
    
    /// Calculate vorticity at a single point
    fn vorticity_at_point<T: RealField + Copy + FromPrimitive>(
        velocity: &VelocityField<T>,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let (nx, ny, nz) = velocity.dimensions;
        let dx = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let dy = T::from_f64(1.0 / ny as f64).unwrap_or_else(T::one);
        let dz = T::from_f64(1.0 / nz as f64).unwrap_or_else(T::one);
        
        let mut omega_x = T::zero();
        let mut omega_y = T::zero();
        let mut omega_z = T::zero();
        
        // ω_x = ∂w/∂y - ∂v/∂z
        if j > 0 && j < ny - 1 && k > 0 && k < nz - 1 {
            if let (Some(v_kp), Some(v_km), Some(w_jp), Some(w_jm)) = (
                velocity.get(i, j, k + 1),
                velocity.get(i, j, k - 1),
                velocity.get(i, j + 1, k),
                velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let dwdy = (w_jp.z - w_jm.z) / (two * dy);
                let dvdz = (v_kp.y - v_km.y) / (two * dz);
                omega_x = dwdy - dvdz;
            }
        }
        
        // ω_y = ∂u/∂z - ∂w/∂x
        if i > 0 && i < nx - 1 && k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km), Some(w_ip), Some(w_im)) = (
                velocity.get(i, j, k + 1),
                velocity.get(i, j, k - 1),
                velocity.get(i + 1, j, k),
                velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let dudz = (u_kp.x - u_km.x) / (two * dz);
                let dwdx = (w_ip.z - w_im.z) / (two * dx);
                omega_y = dudz - dwdx;
            }
        }
        
        // ω_z = ∂v/∂x - ∂u/∂y
        if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
            if let (Some(v_ip), Some(v_im), Some(u_jp), Some(u_jm)) = (
                velocity.get(i + 1, j, k),
                velocity.get(i - 1, j, k),
                velocity.get(i, j + 1, k),
                velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let dvdx = (v_ip.y - v_im.y) / (two * dx);
                let dudy = (u_jp.x - u_jm.x) / (two * dy);
                omega_z = dvdx - dudy;
            }
        }
        
        Vector3::new(omega_x, omega_y, omega_z)
    }
    
    /// Calculate divergence field: ∇·u
    pub fn divergence<T: RealField + Copy + FromPrimitive>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut divergence = Vec::with_capacity(velocity.components.len());
        
        let dx = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let dy = T::from_f64(1.0 / ny as f64).unwrap_or_else(T::one);
        let dz = T::from_f64(1.0 / nz as f64).unwrap_or_else(T::one);
        
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let mut div = T::zero();
                    
                    // ∂u/∂x
                    if i > 0 && i < nx - 1 {
                        if let (Some(u_ip), Some(u_im)) = (
                            velocity.get(i + 1, j, k),
                            velocity.get(i - 1, j, k)
                        ) {
                            let two = T::from_f64(2.0).unwrap_or_else(T::one);
                            div = div + (u_ip.x - u_im.x) / (two * dx);
                        }
                    }
                    
                    // ∂v/∂y
                    if j > 0 && j < ny - 1 {
                        if let (Some(v_jp), Some(v_jm)) = (
                            velocity.get(i, j + 1, k),
                            velocity.get(i, j - 1, k)
                        ) {
                            let two = T::from_f64(2.0).unwrap_or_else(T::one);
                            div = div + (v_jp.y - v_jm.y) / (two * dy);
                        }
                    }
                    
                    // ∂w/∂z
                    if k > 0 && k < nz - 1 {
                        if let (Some(w_kp), Some(w_km)) = (
                            velocity.get(i, j, k + 1),
                            velocity.get(i, j, k - 1)
                        ) {
                            let two = T::from_f64(2.0).unwrap_or_else(T::one);
                            div = div + (w_kp.z - w_km.z) / (two * dz);
                        }
                    }
                    
                    divergence.push(div);
                }
            }
        }
        
        divergence
    }
    
    /// Calculate kinetic energy field: KE = 0.5 * |u|²
    pub fn kinetic_energy<T: RealField + Copy + FromPrimitive>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        velocity.components
            .iter()
            .map(|v| half * v.dot(v))
            .collect()
    }
    
    /// Calculate enstrophy: E = 0.5 * |ω|²
    pub fn enstrophy<T: RealField + Copy + FromPrimitive>(
        velocity: &VelocityField<T>
    ) -> Vec<T> {
        let vorticity = Self::vorticity(velocity);
        let half = T::from_f64(0.5).unwrap_or_else(T::one);
        vorticity
            .iter()
            .map(|omega| half * omega.dot(omega))
            .collect()
    }
}