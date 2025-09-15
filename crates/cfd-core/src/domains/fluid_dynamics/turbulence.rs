//! Turbulence modeling abstractions
//!
//! Provides trait-based turbulence model implementations including
//! LES (Smagorinsky, Dynamic) and mixing length models.

use super::fields::{FlowField, VelocityField};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Turbulence model abstraction following Strategy pattern
pub trait TurbulenceModel<T: RealField + Copy>: Send + Sync {
    /// Calculate turbulent viscosity
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T>;

    /// Calculate turbulent kinetic energy
    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T>;

    /// Get model name
    fn name(&self) -> &str;
}

/// Smagorinsky subgrid-scale model for LES
#[derive(Debug, Clone)]
pub struct SmagorinskyModel<T: RealField + Copy> {
    /// Smagorinsky constant (typically 0.1-0.2)
    pub cs: T,
    /// Base Smagorinsky constant for dynamic model
    pub cs_base: T,
}

impl<T: RealField + Copy + FromPrimitive> SmagorinskyModel<T> {
    /// Create a new Smagorinsky model with standard constant
    pub fn new(cs: T) -> Self {
        Self { cs, cs_base: cs }
    }

    /// Calculate strain rate magnitude at a grid point
    fn calculate_strain_rate(&self, flow_field: &FlowField<T>, idx: usize) -> T {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let i = idx % nx;
        let j = (idx / nx) % ny;
        let k = idx / (nx * ny);

        // Grid spacing (assuming uniform)
        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);

        // Calculate strain rate tensor components
        let mut s11 = T::zero();
        let mut s22 = T::zero();
        let mut s33 = T::zero();

        // Central differences for velocity gradients
        if i > 0 && i < nx - 1 {
            if let (Some(u_plus), Some(u_minus)) = (
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s11 = (u_plus.x - u_minus.x) / (two * delta);
            }
        }

        if j > 0 && j < ny - 1 {
            if let (Some(v_plus), Some(v_minus)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s22 = (v_plus.y - v_minus.y) / (two * delta);
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(w_plus), Some(w_minus)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                s33 = (w_plus.z - w_minus.z) / (two * delta);
            }
        }

        // Off-diagonal components: Sij = 0.5 * (∂ui/∂xj + ∂uj/∂xi)
        let mut s12 = T::zero();
        let mut s13 = T::zero();
        let mut s23 = T::zero();

        // S12 = 0.5 * (∂u/∂y + ∂v/∂x)
        if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm), Some(v_ip), Some(v_im)) = (
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let du_dy = (u_jp.x - u_jm.x) / (two * delta);
                let dv_dx = (v_ip.y - v_im.y) / (two * delta);
                s12 = (du_dy + dv_dx) / two;
            }
        }

        // S13 = 0.5 * (∂u/∂z + ∂w/∂x)
        if i > 0 && i < nx - 1 && k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km), Some(w_ip), Some(w_im)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i + 1, j, k),
                flow_field.velocity.get(i - 1, j, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let du_dz = (u_kp.x - u_km.x) / (two * delta);
                let dw_dx = (w_ip.z - w_im.z) / (two * delta);
                s13 = (du_dz + dw_dx) / two;
            }
        }

        // S23 = 0.5 * (∂v/∂z + ∂w/∂y)
        if j > 0 && j < ny - 1 && k > 0 && k < nz - 1 {
            if let (Some(v_kp), Some(v_km), Some(w_jp), Some(w_jm)) = (
                flow_field.velocity.get(i, j, k + 1),
                flow_field.velocity.get(i, j, k - 1),
                flow_field.velocity.get(i, j + 1, k),
                flow_field.velocity.get(i, j - 1, k),
            ) {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                let dv_dz = (v_kp.y - v_km.y) / (two * delta);
                let dw_dy = (w_jp.z - w_jm.z) / (two * delta);
                s23 = (dv_dz + dw_dy) / two;
            }
        }

        // Strain rate magnitude: |S| = sqrt(2 * Sij * Sij)
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let s_mag_sq =
            two * (s11 * s11 + s22 * s22 + s33 * s33 + two * (s12 * s12 + s13 * s13 + s23 * s23));
        s_mag_sq.sqrt()
    }
}

impl<T: RealField + Copy + FromPrimitive> TurbulenceModel<T> for SmagorinskyModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, _, _) = flow_field.velocity.dimensions;
        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);

        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let strain_rate = self.calculate_strain_rate(flow_field, idx);
                // νₜ = (Cs * Δ)² * |S|
                self.cs * self.cs * delta * delta * strain_rate
            })
            .collect()
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // For Smagorinsky model, estimate TKE from strain rate
        // k ≈ (Cs * Δ * |S|)²
        let (nx, _, _) = flow_field.velocity.dimensions;
        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);

        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let strain_rate = self.calculate_strain_rate(flow_field, idx);
                let cs_delta_s = self.cs * delta * strain_rate;
                cs_delta_s * cs_delta_s
            })
            .collect()
    }

    fn name(&self) -> &str {
        "Smagorinsky"
    }
}

/// Mixing length turbulence model
#[derive(Debug, Clone)]
pub struct MixingLengthModel<T: RealField + Copy> {
    /// Mixing length scale
    pub length_scale: T,
    /// von Kármán constant
    pub kappa: T,
}

impl<T: RealField + Copy + FromPrimitive> MixingLengthModel<T> {
    /// Create a new mixing length model
    pub fn new(length_scale: T) -> Self {
        let kappa =
            T::from_f64(crate::constants::physics::fluid::VON_KARMAN).unwrap_or_else(T::one);
        Self {
            length_scale,
            kappa,
        }
    }

    /// Calculate velocity gradient magnitude
    fn calculate_velocity_gradient(&self, velocity: &VelocityField<T>, idx: usize) -> T {
        let (nx, ny, nz) = velocity.dimensions;
        let i = idx % nx;
        let j = (idx / nx) % ny;
        let k = idx / (nx * ny);

        let delta = T::from_f64(1.0 / nx as f64).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);

        let mut grad_u_sq = T::zero();

        // Calculate ∂u/∂y for wall-bounded flows (primary gradient)
        if j > 0 && j < ny - 1 {
            if let (Some(u_jp), Some(u_jm)) = (velocity.get(i, j + 1, k), velocity.get(i, j - 1, k))
            {
                let dudy = (u_jp.x - u_jm.x) / (two * delta);
                grad_u_sq = dudy * dudy;
            }
        }

        // Add other gradient components for 3D flows
        if i > 0 && i < nx - 1 {
            if let (Some(u_ip), Some(u_im)) = (velocity.get(i + 1, j, k), velocity.get(i - 1, j, k))
            {
                let u_gradient_x = (u_ip.x - u_im.x) / (two * delta);
                let v_gradient_x = (u_ip.y - u_im.y) / (two * delta);
                let w_gradient_x = (u_ip.z - u_im.z) / (two * delta);
                grad_u_sq = grad_u_sq + u_gradient_x * u_gradient_x + v_gradient_x * v_gradient_x + w_gradient_x * w_gradient_x;
            }
        }

        if k > 0 && k < nz - 1 {
            if let (Some(u_kp), Some(u_km)) = (velocity.get(i, j, k + 1), velocity.get(i, j, k - 1))
            {
                let u_gradient_z = (u_kp.x - u_km.x) / (two * delta);
                let v_gradient_z = (u_kp.y - u_km.y) / (two * delta);
                let w_gradient_z = (u_kp.z - u_km.z) / (two * delta);
                grad_u_sq = grad_u_sq + u_gradient_z * u_gradient_z + v_gradient_z * v_gradient_z + w_gradient_z * w_gradient_z;
            }
        }

        grad_u_sq.sqrt()
    }
}

impl<T: RealField + Copy + FromPrimitive> TurbulenceModel<T> for MixingLengthModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // νₜ = l² * |∂u/∂y| (Prandtl's mixing length hypothesis)
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                self.length_scale * self.length_scale * grad_u
            })
            .collect()
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // Estimate TKE from mixing length and velocity gradient
        // k ≈ (l * |∂u/∂y|)²
        flow_field
            .velocity
            .components
            .iter()
            .enumerate()
            .map(|(idx, _)| {
                let grad_u = self.calculate_velocity_gradient(&flow_field.velocity, idx);
                let l_grad_u = self.length_scale * grad_u;
                l_grad_u * l_grad_u
            })
            .collect()
    }

    fn name(&self) -> &str {
        "MixingLength"
    }
}
