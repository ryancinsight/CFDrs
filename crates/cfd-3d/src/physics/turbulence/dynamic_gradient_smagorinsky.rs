//! Dynamic gradient Smagorinsky LES model (Rozema, Bae & Verstappen 2022).
//!
//! # Theorem - Local Dynamic Gradient Smagorinsky Closure
//!
//! Rozema et al. show that the local dynamic Smagorinsky coefficient becomes
//! singular when it is built from the strain tensor alone. Replacing the strain
//! tensor in the Germano identity with the full velocity-gradient tensor yields
//! the dynamic gradient Smagorinsky model:
//!
//! ```text
//! C_GS = max( L_ij M_ij / (M_kl M_kl), 0 )
//! M_ij = -2 (Δ̂² |∇û| ∂_j û_i - Δ² |∇u| ∂_j u_i)
//! ν_t = C_GS Δ² |∇u|
//! ```
//!
//! where `Δ = (dx·dy·dz)^(1/3)` and `Δ̂ = r Δ` is the test-filter width.
//! The local dynamic coefficient is bounded at leading order and does not
//! require homogeneous-direction averaging for stability.
//!
//! ## References
//!
//! - Rozema, W., Bae, H. J. & Verstappen, R. W. C. P. (2022). "Local dynamic
//!   gradient Smagorinsky model for large-eddy simulation." *Phys. Rev. Fluids*
//!   7, 074604.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::field_ops::{
    gradient_magnitude, linear_index, matrix_contract, symmetric_matrix_contract,
    velocity_gradient_tensor,
};
use super::filter_ops::{box_filter_moments_at, box_filter_velocity_at, resolved_stress_tensor};
use super::sgs_energy::kinetic_energy_from_eddy_viscosity;

/// Local dynamic gradient Smagorinsky model.
#[derive(Debug, Clone)]
pub struct DynamicGradientSmagorinskyModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Grid spacing in the x direction [m].
    pub dx: T,
    /// Grid spacing in the y direction [m].
    pub dy: T,
    /// Grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width `Δ = (dx·dy·dz)^(1/3)` [m].
    pub filter_width: T,
    /// Ratio between the test filter and LES filter widths.
    pub test_filter_ratio: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    DynamicGradientSmagorinskyModel<T>
{
    /// Create a local DGSM with unit grid spacing.
    pub fn new() -> Self {
        let one = T::one();
        let two = one + one;
        Self {
            dx: one,
            dy: one,
            dz: one,
            filter_width: one,
            test_filter_ratio: two,
        }
    }

    /// Create a local DGSM with the physical filter width.
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one = T::one();
        let two = one + one;
        let one_third = <T as FromPrimitive>::from_f64(1.0 / 3.0)
            .expect("1/3 is an IEEE 754 representable f64 constant");
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            dx,
            dy,
            dz,
            filter_width,
            test_filter_ratio: two,
        }
    }

    #[inline]
    fn test_filter_width(&self) -> T {
        self.filter_width * self.test_filter_ratio
    }

    #[inline]
    fn gradient_model_tensor(
        &self,
        gradient: &[[T; 3]; 3],
        filtered_gradient: &[[T; 3]; 3],
        resolved_scale: T,
        filtered_scale: T,
    ) -> [[T; 3]; 3] {
        let two = T::one() + T::one();
        let delta_sq = self.filter_width * self.filter_width;
        let delta_hat = self.test_filter_width();
        let delta_hat_sq = delta_hat * delta_hat;
        let mut tensor = [[T::zero(); 3]; 3];

        for i in 0..3 {
            for j in 0..3 {
                tensor[i][j] = -two
                    * (delta_hat_sq * filtered_scale * filtered_gradient[i][j]
                        - delta_sq * resolved_scale * gradient[i][j]);
            }
        }

        tensor
    }

    #[inline]
    fn local_viscosity_at(
        &self,
        flow: &FlowField<T>,
        filtered_velocity: &[nalgebra::Vector3<T>],
        i: usize,
        j: usize,
        k: usize,
    ) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
        let idx = linear_index(nx, ny, i, j, k);
        let gradient = velocity_gradient_tensor(
            &flow.velocity.components,
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
            filtered_velocity,
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
        let resolved_scale = gradient_magnitude(&gradient);
        let filtered_scale = gradient_magnitude(&filtered_gradient);
        let moments = box_filter_moments_at(&flow.velocity.components, nx, ny, nz, i, j, k);
        let leonard = resolved_stress_tensor(moments, filtered_velocity[idx]);
        let m = self.gradient_model_tensor(
            &gradient,
            &filtered_gradient,
            resolved_scale,
            filtered_scale,
        );
        let numerator = symmetric_matrix_contract(leonard, &m);
        let denominator = matrix_contract(&m, &m);
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        if denominator <= eps {
            T::zero()
        } else {
            let coefficient = numerator / denominator;
            if coefficient > T::zero() {
                coefficient * self.filter_width * self.filter_width * resolved_scale
            } else {
                T::zero()
            }
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    Default for DynamicGradientSmagorinskyModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    TurbulenceModel<T> for DynamicGradientSmagorinskyModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let n = nx * ny * nz;
        let mut filtered_velocity = Vec::with_capacity(n);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    filtered_velocity.push(box_filter_velocity_at(
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

        let mut viscosity = Vec::with_capacity(n);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    viscosity.push(self.local_viscosity_at(
                        flow_field,
                        &filtered_velocity,
                        i,
                        j,
                        k,
                    ));
                }
            }
        }

        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        self.turbulent_viscosity(flow_field)
            .into_iter()
            .map(|nu_t| kinetic_energy_from_eddy_viscosity(nu_t, self.filter_width))
            .collect()
    }

    fn name(&self) -> &'static str {
        "DynamicGradientSmagorinsky"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

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
    fn zero_flow_has_zero_viscosity() {
        let flow = FlowField::<f64>::new(3, 3, 3);
        let model = DynamicGradientSmagorinskyModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);

        assert!(viscosity.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn rigid_body_rotation_remains_finite_without_singularities() {
        let mut flow = FlowField::<f64>::new(4, 4, 4);
        fill_velocity_field(&mut flow, |x, y, _| Vector3::new(-y, x, 0.0));

        let model = DynamicGradientSmagorinskyModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);

        assert!(viscosity
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0));
    }

    #[test]
    fn linear_field_matches_closed_form_value() {
        let dx = 0.5;
        let dy = 1.5;
        let dz = 2.0;
        let mut flow = FlowField::<f64>::new(5, 5, 5);
        let gradient = [[1.3, -0.4, 0.8], [0.9, -1.1, 0.2], [-0.6, 0.5, 0.7]];

        fill_velocity_field(&mut flow, |i, j, k| {
            let x = i * dx;
            let y = j * dy;
            let z = k * dz;
            Vector3::new(
                gradient[0][0] * x + gradient[0][1] * y + gradient[0][2] * z,
                gradient[1][0] * x + gradient[1][1] * y + gradient[1][2] * z,
                gradient[2][0] * x + gradient[2][1] * y + gradient[2][2] * z,
            )
        });

        let model = DynamicGradientSmagorinskyModel::<f64>::with_filter_width(dx, dy, dz);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = linear_index(5, 5, 2, 2, 2);

        let gradient_norm_sq: f64 = gradient.iter().flatten().map(|value| value * value).sum();
        let gradient_magnitude = (2.0 * gradient_norm_sq).sqrt();
        let delta = (dx * dy * dz).powf(1.0 / 3.0);
        let delta_hat = delta * model.test_filter_ratio;
        let scale = -2.0 * (delta_hat * delta_hat - delta * delta) * gradient_magnitude;
        let covariance = [
            2.0 * dx * dx / 3.0,
            2.0 * dy * dy / 3.0,
            2.0 * dz * dz / 3.0,
        ];

        let mut leonard = [[0.0; 3]; 3];
        for row in 0..3 {
            for col in 0..3 {
                leonard[row][col] = gradient[row][0] * gradient[col][0] * covariance[0]
                    + gradient[row][1] * gradient[col][1] * covariance[1]
                    + gradient[row][2] * gradient[col][2] * covariance[2];
            }
        }

        let mut m = [[0.0; 3]; 3];
        for row in 0..3 {
            for col in 0..3 {
                m[row][col] = scale * gradient[row][col];
            }
        }

        let numerator = (0..3)
            .flat_map(|row| (0..3).map(move |col| leonard[row][col] * m[row][col]))
            .sum::<f64>();
        let denominator = (0..3)
            .flat_map(|row| (0..3).map(move |col| m[row][col] * m[row][col]))
            .sum::<f64>();
        let coefficient = if denominator > 0.0 {
            (numerator / denominator).max(0.0)
        } else {
            0.0
        };
        let expected = coefficient * delta * delta * gradient_magnitude;

        assert_relative_eq!(viscosity[center], expected, epsilon = 1e-12);
    }
}
