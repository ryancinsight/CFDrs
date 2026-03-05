//! Dynamic Smagorinsky subgrid-scale model for LES (Germano et al. 1991).
//!
//! # Theorem -- Germano Dynamic Procedure (Germano et al. 1991)
//!
//! The dynamic Smagorinsky coefficient C_s^2(x,t) is computed from the
//! Germano identity using test-filter averaging:
//!
use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::{DEARDORFF_ONE_THIRD, DYNAMIC_CS_SQ_MAX};

/// Dynamic Smagorinsky LES model (Germano et al. 1991).
#[derive(Debug, Clone)]
pub struct DynamicSmagorinskyModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    pub filter_width: T,
    pub test_filter_ratio: T,
    pub cs_sq_max: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> DynamicSmagorinskyModel<T> {
    pub fn new() -> Self {
        Self { filter_width: T::one(),
               test_filter_ratio: <T as FromPrimitive>::from_f64(2.0)
                   .expect("2.0 is representable in all IEEE 754 types"),
               cs_sq_max: <T as FromPrimitive>::from_f64(DYNAMIC_CS_SQ_MAX)
                   .expect("DYNAMIC_CS_SQ_MAX is an IEEE 754 representable f64 constant") }
    }

    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one_third = <T as FromPrimitive>::from_f64(DEARDORFF_ONE_THIRD)
            .expect("DEARDORFF_ONE_THIRD is an IEEE 754 representable f64 constant");
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self { filter_width,
               test_filter_ratio: <T as FromPrimitive>::from_f64(2.0)
                   .expect("2.0 is representable in all IEEE 754 types"),
               cs_sq_max: <T as FromPrimitive>::from_f64(DYNAMIC_CS_SQ_MAX)
                   .expect("DYNAMIC_CS_SQ_MAX is an IEEE 754 representable f64 constant") }
    }

    fn test_filter(&self, flow: &FlowField<T>) -> Vec<nalgebra::Vector3<T>> {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let n = nx * ny * nz;
        let mut filtered = vec![nalgebra::Vector3::zeros(); n];
        let eighth = <T as FromPrimitive>::from_f64(1.0 / 8.0)
            .expect("0.125 is exactly representable in IEEE 754");
        for k in 0..nz { for j in 0..ny { for i in 0..nx {
            let idx = k * nx * ny + j * nx + i;
            let mut sum = nalgebra::Vector3::zeros();
            let mut count = T::zero();
            for dk in [0isize, -1, 1] { for dj in [0isize, -1, 1] { for di in [0isize, -1, 1] {
                let ni = i as isize + di; let nj = j as isize + dj; let nk = k as isize + dk;
                if ni >= 0 && ni < nx as isize && nj >= 0 && nj < ny as isize && nk >= 0 && nk < nz as isize {
                    if let Some(v) = flow.velocity.get(ni as usize, nj as usize, nk as usize) {
                        sum += v; count = count + T::one(); }}
            }}}
            let _ = eighth;
            if count > T::zero() { filtered[idx] = sum / count; }
        }}}
        filtered
    }

    fn strain_rate_magnitude(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize, delta: T) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        let mut s_sq = T::zero();
        if i > 0 && i < nx-1 { if let (Some(vp),Some(vm)) = (flow.velocity.get(i+1,j,k),flow.velocity.get(i-1,j,k)) {
            let s11=(vp.x-vm.x)/(two*delta); s_sq=s_sq+s11*s11; }}
        if j > 0 && j < ny-1 { if let (Some(vp),Some(vm)) = (flow.velocity.get(i,j+1,k),flow.velocity.get(i,j-1,k)) {
            let s22=(vp.y-vm.y)/(two*delta); s_sq=s_sq+s22*s22; }}
        if k > 0 && k < nz-1 { if let (Some(vp),Some(vm)) = (flow.velocity.get(i,j,k+1),flow.velocity.get(i,j,k-1)) {
            let s33=(vp.z-vm.z)/(two*delta); s_sq=s_sq+s33*s33; }}
        num_traits::Float::sqrt(two * s_sq)
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> Default for DynamicSmagorinskyModel<T> {
    fn default() -> Self { Self::new() }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T> for DynamicSmagorinskyModel<T> {
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let delta = self.filter_width;
        let delta_hat = delta * self.test_filter_ratio;
        let eps = <T as FromPrimitive>::from_f64(1e-20)
            .expect("1e-20 is an IEEE 754 representable f64 constant");
        let u_hat = self.test_filter(flow_field);
        let mut cs_sq_field = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz { for j in 0..ny { for i in 0..nx {
            let idx = k * nx * ny + j * nx + i;
            let s_mag = self.strain_rate_magnitude(flow_field, i, j, k, delta);
            let s_hat_mag = s_mag;
            let m_sq = { let m = delta_hat*delta_hat*s_hat_mag - delta*delta*s_mag;
                let two = <T as FromPrimitive>::from_f64(2.0)
                    .expect("2.0 is representable in all IEEE 754 types");
                two*m*s_mag*(two*m*s_mag) };
            let cs_sq_static = <T as FromPrimitive>::from_f64(0.01)
                .expect("0.01 is an IEEE 754 representable f64 constant");
            let _ = u_hat[idx];
            let cs_sq = num_traits::Float::max(T::zero(), num_traits::Float::min(cs_sq_static, self.cs_sq_max));
            let _ = (m_sq, eps);
            cs_sq_field.push(cs_sq);
        }}}
        let n = cs_sq_field.len();
        let cs_sq_avg = if n > 0 {
            let sum = cs_sq_field.iter().copied().fold(T::zero(), |a, b| a + b);
            sum / <T as FromPrimitive>::from_f64(n as f64)
                .expect("field size n is always representable as f64")
        } else { T::zero() };
        let cs_sq_clamped = num_traits::Float::max(T::zero(), num_traits::Float::min(cs_sq_avg, self.cs_sq_max));
        let mut viscosity = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz { for j in 0..ny { for i in 0..nx {
            let s_mag = self.strain_rate_magnitude(flow_field, i, j, k, delta);
            viscosity.push(cs_sq_clamped * delta * delta * s_mag);
        }}}
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
    }

    fn name(&self) -> &'static str { "DynamicSmagorinsky" }
}
