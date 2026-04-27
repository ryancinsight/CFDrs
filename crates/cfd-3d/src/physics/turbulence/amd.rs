//! Anisotropic Minimum Dissipation (AMD) model for LES.
//!
//! # Theorem - AMD Closure (Rozema et al. 2015)
//!
//! The AMD model defines the eddy viscosity by the minimum dissipation needed
//! to dissipate the energy of the unresolved scales on anisotropic grids:
//!
//! ```text
//! ν_t = C_A * max(-B : S, 0) / |∇u|_{3×3}²
//! ```
//!
//! where:
//! - `S` is the resolved strain-rate tensor
//! - `B = (∇_Δ u)^T (∇_Δ u)` is the anisotropically scaled gradient metric
//! - `|∇u|_{3×3}²` is the squared Frobenius norm of the resolved gradient
//! - `C_A = 1/12` for the continuous derivation and `C_A = 1/3` for the
//!   second-order central-difference discretisation used in this codebase
//!
//! **Proof sketch.** The gradient-model expansion of the exact sub-grid stress
//! on a rectangular filter box gives `τ - tr(τ)I/3 = (1/12) B + O(Δ^4)`.
//! Matching the leading-order eddy dissipation with the minimum-dissipation
//! closure yields the AMD formula above.  On anisotropic grids this avoids the
//! filter-width approximation needed by isotropic QR-style models.
//!
//! ## References
//!
//! - Rozema, W., Bae, H. J., Moin, P. & Verstappen, R. (2015).
//!   "Minimum-dissipation models for large-eddy simulation."
//! - Rozema, W., Bae, H. J. & Verstappen, R. (2022).
//!   "Local dynamic gradient Smagorinsky model for large-eddy simulation."

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::AMD_C_A_SECOND_ORDER;
use super::field_ops::{
    strain_components, symmetric_contract, velocity_gradient_tensor, SymmetricTensor6,
};
use super::sgs_energy::kinetic_energy_from_eddy_viscosity;

/// Anisotropic Minimum Dissipation closure for large-eddy simulation.
#[derive(Debug, Clone)]
pub struct AnisotropicMinimumDissipationModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy>
{
    /// AMD coefficient `C_A`.
    pub c_a: T,
    /// Physical grid spacing in the x direction [m].
    pub dx: T,
    /// Physical grid spacing in the y direction [m].
    pub dy: T,
    /// Physical grid spacing in the z direction [m].
    pub dz: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive>
    AnisotropicMinimumDissipationModel<T>
{
    /// Create an AMD model on a unit Cartesian grid.
    ///
    /// Uses `C_A = 1/3`, the second-order central-difference correction from
    /// Rozema et al. (2015).
    pub fn new() -> Self {
        Self {
            c_a: <T as FromPrimitive>::from_f64(AMD_C_A_SECOND_ORDER)
                .expect("AMD_C_A_SECOND_ORDER is an IEEE 754 representable f64 constant"),
            dx: T::one(),
            dy: T::one(),
            dz: T::one(),
        }
    }

    /// Create an AMD model with a custom coefficient and physical spacings.
    pub fn with_constant(c_a: T, dx: T, dy: T, dz: T) -> Self {
        Self { c_a, dx, dy, dz }
    }

    /// Create an AMD model using the default coefficient and physical spacings.
    pub fn with_grid_spacing(dx: T, dy: T, dz: T) -> Self {
        Self {
            c_a: <T as FromPrimitive>::from_f64(AMD_C_A_SECOND_ORDER)
                .expect("AMD_C_A_SECOND_ORDER is an IEEE 754 representable f64 constant"),
            dx,
            dy,
            dz,
        }
    }

    #[inline]
    fn gradient_metric(&self, gradient: &[[T; 3]; 3]) -> SymmetricTensor6<T> {
        let spacing_sq = [self.dx * self.dx, self.dy * self.dy, self.dz * self.dz];
        let mut beta = [[T::zero(); 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let mut value = T::zero();
                for k in 0..3 {
                    value += spacing_sq[k] * gradient[i][k] * gradient[j][k];
                }
                beta[i][j] = value;
                beta[j][i] = value;
            }
        }

        SymmetricTensor6 {
            xx: beta[0][0],
            yy: beta[1][1],
            zz: beta[2][2],
            xy: beta[0][1],
            xz: beta[0][2],
            yz: beta[1][2],
        }
    }

    #[inline]
    fn gradient_frobenius_norm_sq(gradient: &[[T; 3]; 3]) -> T {
        let mut norm_sq = T::zero();
        for row in gradient {
            for &value in row {
                norm_sq += value * value;
            }
        }
        norm_sq
    }

    fn amd_viscosity_at(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (nx, ny, nz) = flow.velocity.dimensions;
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
        let strain = strain_components(&gradient);
        let metric = self.gradient_metric(&gradient);
        let numerator = -symmetric_contract(metric, strain);
        let numerator = if numerator > T::zero() {
            numerator
        } else {
            T::zero()
        };
        let denominator = Self::gradient_frobenius_norm_sq(&gradient);
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        if denominator <= eps || numerator <= T::zero() {
            T::zero()
        } else {
            self.c_a * numerator / denominator
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> Default
    for AnisotropicMinimumDissipationModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive> TurbulenceModel<T>
    for AnisotropicMinimumDissipationModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let mut viscosity = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    viscosity.push(self.amd_viscosity_at(flow_field, i, j, k));
                }
            }
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let length_scale = if self.dx >= self.dy && self.dx >= self.dz {
            self.dx
        } else if self.dy >= self.dz {
            self.dy
        } else {
            self.dz
        };
        self.turbulent_viscosity(flow_field)
            .into_iter()
            .map(|nu_t| kinetic_energy_from_eddy_viscosity(nu_t, length_scale))
            .collect()
    }

    fn name(&self) -> &'static str {
        "AMD"
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
                    let idx = k * nx * ny + j * nx + i;
                    flow.velocity.components[idx] = generator(i as f64, j as f64, k as f64);
                }
            }
        }
    }

    #[test]
    fn zero_flow_has_zero_viscosity() {
        let flow = FlowField::<f64>::new(3, 3, 3);
        let model = AnisotropicMinimumDissipationModel::<f64>::with_grid_spacing(0.5, 1.5, 2.0);
        let viscosity = model.turbulent_viscosity(&flow);

        assert!(viscosity.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn pure_shear_has_zero_viscosity() {
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |_, y, _| Vector3::new(y, 0.0, 0.0));

        let model = AnisotropicMinimumDissipationModel::<f64>::with_grid_spacing(0.5, 1.5, 2.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = 13;

        assert_relative_eq!(viscosity[center], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn anisotropic_linear_field_matches_analytic_amd_value() {
        let dx = 0.5;
        let dy = 1.5;
        let dz = 2.0;
        let mut flow = FlowField::<f64>::new(3, 3, 3);

        let a = [
            [1.613734946325263, 0.7232431980345768, -1.4146463799644549],
            [0.9674030241104554, -1.826607564819564, 0.9954169690130685],
            [0.9131006443117435, 1.7808987218488506, 0.21287261849430106],
        ];

        fill_velocity_field(&mut flow, |x, y, z| {
            let x = x * dx;
            let y = y * dy;
            let z = z * dz;
            Vector3::new(
                a[0][0] * x + a[0][1] * y + a[0][2] * z,
                a[1][0] * x + a[1][1] * y + a[1][2] * z,
                a[2][0] * x + a[2][1] * y + a[2][2] * z,
            )
        });

        let model = AnisotropicMinimumDissipationModel::<f64>::with_grid_spacing(dx, dy, dz);
        let viscosity = model.turbulent_viscosity(&flow);

        let strain = [
            [
                a[0][0],
                0.5 * (a[0][1] + a[1][0]),
                0.5 * (a[0][2] + a[2][0]),
            ],
            [
                0.5 * (a[1][0] + a[0][1]),
                a[1][1],
                0.5 * (a[1][2] + a[2][1]),
            ],
            [
                0.5 * (a[2][0] + a[0][2]),
                0.5 * (a[2][1] + a[1][2]),
                a[2][2],
            ],
        ];
        let spacing_sq = [dx * dx, dy * dy, dz * dz];
        let mut metric = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let mut value = 0.0;
                for k in 0..3 {
                    value += spacing_sq[k] * a[i][k] * a[j][k];
                }
                metric[i][j] = value;
                metric[j][i] = value;
            }
        }
        let b_dot_s = metric
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, value)| *value * strain[i][j])
                    .sum::<f64>()
            })
            .sum::<f64>();
        let grad_norm_sq = a.iter().flatten().map(|value| *value * *value).sum::<f64>();
        let expected = (1.0 / 3.0) * (-b_dot_s).max(0.0) / grad_norm_sq;

        assert_relative_eq!(viscosity[13], expected, epsilon = 1e-12);
        assert!(viscosity[13] > 0.0);
    }

    #[test]
    fn turbulent_kinetic_energy_uses_yoshizawa_relation_not_viscosity_alias() {
        let dx = 0.5;
        let dy = 1.5;
        let dz = 2.0;
        let mut flow = FlowField::<f64>::new(3, 3, 3);

        let a = [
            [1.613734946325263, 0.7232431980345768, -1.4146463799644549],
            [0.9674030241104554, -1.826607564819564, 0.9954169690130685],
            [0.9131006443117435, 1.7808987218488506, 0.21287261849430106],
        ];

        fill_velocity_field(&mut flow, |x, y, z| {
            let x = x * dx;
            let y = y * dy;
            let z = z * dz;
            Vector3::new(
                a[0][0] * x + a[0][1] * y + a[0][2] * z,
                a[1][0] * x + a[1][1] * y + a[1][2] * z,
                a[2][0] * x + a[2][1] * y + a[2][2] * z,
            )
        });

        let model = AnisotropicMinimumDissipationModel::<f64>::with_grid_spacing(dx, dy, dz);
        let viscosity = model.turbulent_viscosity(&flow);
        let kinetic_energy = model.turbulent_kinetic_energy(&flow);
        let center = 13;
        let expected = (viscosity[center] / (0.094 * dz)).powi(2);

        assert!(viscosity[center] > 0.0);
        assert_relative_eq!(kinetic_energy[center], expected, epsilon = 1e-12);
        assert!(
            (kinetic_energy[center] - viscosity[center]).abs() > 1e-9,
            "SGS kinetic energy must not alias eddy viscosity"
        );
    }
}
