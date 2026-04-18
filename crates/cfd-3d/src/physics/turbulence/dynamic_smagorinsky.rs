//! Dynamic Smagorinsky subgrid-scale model for LES (Germano et al. 1991).
//!
//! # Theorem — Germano Dynamic Procedure
//!
//! The dynamic procedure estimates the Smagorinsky coefficient from the
//! Germano identity
//!
//! ```text
//! L_ij = C_s² M_ij
//! C_s² = <L_ij M_ij> / <M_ij M_ij>
//! ```
//!
//! where `<·>` denotes a local averaging operator over the test-filter stencil.
//! The implementation below uses a discrete top-hat filter and exact spatial
//! derivatives on the cell-centred velocity field.
//!
//! ## References
//!
//! - Germano, M., Piomelli, U., Moin, P. & Cabot, W.H. (1991). "A dynamic
//!   subgrid-scale eddy viscosity model." *Phys. Fluids A* 3(7):1760–1765.
//! - Lilly, D.K. (1992). "A proposed modification of the Germano subgrid-scale
//!   closure method." *Phys. Fluids A* 4(3):633–635.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

use super::constants::DEARDORFF_ONE_THIRD;
use super::field_ops::{
    linear_index, strain_components, strain_magnitude, symmetric_contract,
    velocity_gradient_tensor, SymmetricTensor6,
};

#[derive(Clone, Copy, Debug)]
struct FilterMoments<T: RealField + Copy> {
    velocity: Vector3<T>,
    uu: T,
    vv: T,
    ww: T,
    uv: T,
    uw: T,
    vw: T,
}

impl<T: RealField + Copy> FilterMoments<T> {
    fn zero() -> Self {
        Self {
            velocity: Vector3::zeros(),
            uu: T::zero(),
            vv: T::zero(),
            ww: T::zero(),
            uv: T::zero(),
            uw: T::zero(),
            vw: T::zero(),
        }
    }
}

/// Dynamic Smagorinsky LES model (Germano et al. 1991).
#[derive(Debug, Clone)]
pub struct DynamicSmagorinskyModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Grid spacing in the x direction [m].
    pub dx: T,
    /// Grid spacing in the y direction [m].
    pub dy: T,
    /// Grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    pub filter_width: T,
    /// Ratio between test-filter width and LES filter width.
    pub test_filter_ratio: T,
    /// Numerical upper bound for the dynamic coefficient C_s².
    pub cs_sq_max: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive>
    DynamicSmagorinskyModel<T>
{
    /// Create a dynamic Smagorinsky model with unit grid spacing.
    pub fn new() -> Self {
        let one = T::one();
        let two = one + one;
        Self {
            dx: one,
            dy: one,
            dz: one,
            filter_width: one,
            test_filter_ratio: two,
            cs_sq_max: one,
        }
    }

    /// Create a dynamic Smagorinsky model with the physical filter width.
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one = T::one();
        let two = one + one;
        let one_third = <T as FromPrimitive>::from_f64(DEARDORFF_ONE_THIRD)
            .expect("DEARDORFF_ONE_THIRD is an IEEE 754 representable f64 constant");
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            dx,
            dy,
            dz,
            filter_width,
            test_filter_ratio: two,
            cs_sq_max: one,
        }
    }

    fn box_filter_velocity_at(
        &self,
        velocity: &[Vector3<T>],
        nx: usize,
        ny: usize,
        nz: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vector3<T> {
        let mut sum = Vector3::zeros();
        let mut count = 0usize;
        let plane = nx * ny;

        for dk in [0isize, -1, 1] {
            let nk = k as isize + dk;
            if nk < 0 || nk >= nz as isize {
                continue;
            }
            for dj in [0isize, -1, 1] {
                let nj = j as isize + dj;
                if nj < 0 || nj >= ny as isize {
                    continue;
                }
                let base = nk as usize * plane + nj as usize * nx;
                for di in [0isize, -1, 1] {
                    let ni = i as isize + di;
                    if ni < 0 || ni >= nx as isize {
                        continue;
                    }
                    sum += velocity[base + ni as usize];
                    count += 1;
                }
            }
        }

        let count_t = <T as FromPrimitive>::from_usize(count)
            .expect("box filter stencil size is always representable");
        sum / count_t
    }

    fn box_filter_moments_at(
        &self,
        velocity: &[Vector3<T>],
        nx: usize,
        ny: usize,
        nz: usize,
        i: usize,
        j: usize,
        k: usize,
    ) -> FilterMoments<T> {
        let mut moments = FilterMoments::zero();
        let mut count = 0usize;
        let plane = nx * ny;

        for dk in [0isize, -1, 1] {
            let nk = k as isize + dk;
            if nk < 0 || nk >= nz as isize {
                continue;
            }
            for dj in [0isize, -1, 1] {
                let nj = j as isize + dj;
                if nj < 0 || nj >= ny as isize {
                    continue;
                }
                let base = nk as usize * plane + nj as usize * nx;
                for di in [0isize, -1, 1] {
                    let ni = i as isize + di;
                    if ni < 0 || ni >= nx as isize {
                        continue;
                    }
                    let v = velocity[base + ni as usize];
                    moments.velocity += v;
                    moments.uu += v.x * v.x;
                    moments.vv += v.y * v.y;
                    moments.ww += v.z * v.z;
                    moments.uv += v.x * v.y;
                    moments.uw += v.x * v.z;
                    moments.vw += v.y * v.z;
                    count += 1;
                }
            }
        }

        let count_t = <T as FromPrimitive>::from_usize(count)
            .expect("box filter stencil size is always representable");
        let inv = T::one() / count_t;
        moments.velocity *= inv;
        moments.uu *= inv;
        moments.vv *= inv;
        moments.ww *= inv;
        moments.uv *= inv;
        moments.uw *= inv;
        moments.vw *= inv;
        moments
    }

    fn resolved_stress_tensor(
        moments: FilterMoments<T>,
        filtered_velocity: Vector3<T>,
    ) -> SymmetricTensor6<T> {
        SymmetricTensor6 {
            xx: moments.uu - filtered_velocity.x * filtered_velocity.x,
            yy: moments.vv - filtered_velocity.y * filtered_velocity.y,
            zz: moments.ww - filtered_velocity.z * filtered_velocity.z,
            xy: moments.uv - filtered_velocity.x * filtered_velocity.y,
            xz: moments.uw - filtered_velocity.x * filtered_velocity.z,
            yz: moments.vw - filtered_velocity.y * filtered_velocity.z,
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> Default
    for DynamicSmagorinskyModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for DynamicSmagorinskyModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let delta = self.filter_width;
        let delta_sq = delta * delta;
        let delta_hat = delta * self.test_filter_ratio;
        let delta_hat_sq = delta_hat * delta_hat;
        let two = T::one() + T::one();
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        let mut filtered_velocity = Vec::with_capacity(n);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    filtered_velocity.push(self.box_filter_velocity_at(
                        &flow_field.velocity.components,
                        nx,
                        ny,
                        nz,
                        i,
                        j,
                        k,
                    ));
                }
            }
        }

        let mut numerator = T::zero();
        let mut denominator = T::zero();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = linear_index(nx, ny, i, j, k);
                    let gradient = velocity_gradient_tensor(
                        &flow_field.velocity.components,
                        nx,
                        ny,
                        nz,
                        i,
                        j,
                        k,
                        self.dx,
                        self.dy,
                        self.dz,
                    );
                    let filtered_gradient = velocity_gradient_tensor(
                        &filtered_velocity,
                        nx,
                        ny,
                        nz,
                        i,
                        j,
                        k,
                        self.dx,
                        self.dy,
                        self.dz,
                    );

                    let s = strain_components(&gradient);
                    let s_hat = strain_components(&filtered_gradient);
                    let l = Self::resolved_stress_tensor(
                        self.box_filter_moments_at(
                            &flow_field.velocity.components,
                            nx,
                            ny,
                            nz,
                            i,
                            j,
                            k,
                        ),
                        filtered_velocity[idx],
                    );
                    let m = SymmetricTensor6 {
                        xx: -two * (delta_hat_sq * s_hat.xx - delta_sq * s.xx),
                        yy: -two * (delta_hat_sq * s_hat.yy - delta_sq * s.yy),
                        zz: -two * (delta_hat_sq * s_hat.zz - delta_sq * s.zz),
                        xy: -two * (delta_hat_sq * s_hat.xy - delta_sq * s.xy),
                        xz: -two * (delta_hat_sq * s_hat.xz - delta_sq * s.xz),
                        yz: -two * (delta_hat_sq * s_hat.yz - delta_sq * s.yz),
                    };

                    numerator += symmetric_contract(l, m);
                    denominator += symmetric_contract(m, m);
                }
            }
        }

        let cs_sq = if denominator > eps {
            numerator / denominator
        } else {
            T::zero()
        };
        let cs_sq =
            num_traits::Float::max(T::zero(), num_traits::Float::min(cs_sq, self.cs_sq_max));
        if cs_sq <= eps {
            return vec![T::zero(); n];
        }

        let prefactor = cs_sq * delta_sq;
        let mut viscosity = Vec::with_capacity(n);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let gradient = velocity_gradient_tensor(
                        &flow_field.velocity.components,
                        nx,
                        ny,
                        nz,
                        i,
                        j,
                        k,
                        self.dx,
                        self.dy,
                        self.dz,
                    );
                    let s_mag = strain_magnitude(&gradient);
                    viscosity.push(prefactor * s_mag);
                }
            }
        }

        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
    }

    fn name(&self) -> &'static str {
        "DynamicSmagorinsky"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn fill_velocity_field<F>(flow: &mut FlowField<f64>, mut generator: F)
    where
        F: FnMut(f64, f64, f64) -> Vector3<f64>,
    {
        let (nx, ny, nz) = flow.velocity.dimensions;
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = linear_index(nx, ny, i, j, k);
                    flow.velocity.components[idx] = generator(i as f64, j as f64, k as f64);
                }
            }
        }
    }

    #[test]
    fn rigid_body_rotation_has_zero_dynamic_viscosity() {
        let mut flow = FlowField::<f64>::new(4, 4, 4);
        fill_velocity_field(&mut flow, |x, y, _z| Vector3::new(-y, x, 0.0));

        let model = DynamicSmagorinskyModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = linear_index(4, 4, 2, 2, 2);
        assert_relative_eq!(viscosity[center], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn zero_flow_remains_zero() {
        let flow = FlowField::<f64>::new(3, 3, 3);
        let model = DynamicSmagorinskyModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        assert!(viscosity.iter().all(|value| value.abs() < 1e-15));
    }
}
