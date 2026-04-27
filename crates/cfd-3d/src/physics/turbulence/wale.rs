//! Wall-Adapting Local Eddy-viscosity (WALE) model for LES.
//!
//! # Theorem - WALE SGS Model (Nicoud & Ducros 1999)
//!
//! The WALE closure defines the eddy viscosity from invariants of the
//! velocity-gradient tensor `G = ∇u`:
//!
//! ```text
//! nu_t = (C_w Delta)^2 * (S^d : S^d)^(3/2)
//!        / [ (S : S)^(5/2) + (S^d : S^d)^(5/4) ]
//! ```
//!
//! where `S = 1/2 (G + G^T)` is the resolved strain tensor and `S^d` is the
//! symmetric trace-free part of `G^2`.  The closure recovers the correct
//! `y^3` near-wall scaling without dynamic damping functions.
//!
//! **Proof sketch.** The square of the gradient tensor contains both strain
//! and rotation information.  Removing its trace yields a wall-sensitive
//! invariant that vanishes for simple shear, while the ratio above produces an
//! eddy viscosity with the proper asymptotic behaviour close to solid walls.
//!
//! ## References
//!
//! - Nicoud, F. & Ducros, F. (1999). "Subgrid-Scale Stress Modelling on the
//!   Square of the Velocity Gradient Tensor." Flow, Turbulence and Combustion.
//! - Recent comparative assessments continue to use WALE as a baseline LES
//!   model for wall-bounded flows.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::WALE_CW;
use super::field_ops::{
    matrix_square, strain_components, symmetric_contract, symmetric_trace_free_part,
    velocity_gradient_tensor,
};
use super::sgs_energy::kinetic_energy_from_eddy_viscosity;

/// Wall-Adapting Local Eddy-viscosity model.
#[derive(Debug, Clone)]
pub struct WaleModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// WALE constant `C_w`.
    pub c_w: T,
    /// Physical grid spacing in the x direction [m].
    pub dx: T,
    /// Physical grid spacing in the y direction [m].
    pub dy: T,
    /// Physical grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width `Delta = (dx dy dz)^(1/3)` [m].
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    WaleModel<T>
{
    /// Create a WALE model on a unit grid.
    pub fn new() -> Self {
        Self {
            c_w: <T as FromPrimitive>::from_f64(WALE_CW)
                .expect("WALE_CW is an IEEE 754 representable f64 constant"),
            dx: T::one(),
            dy: T::one(),
            dz: T::one(),
            filter_width: T::one(),
        }
    }

    /// Create a WALE model with a physically correct filter width.
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            c_w: <T as FromPrimitive>::from_f64(WALE_CW)
                .expect("WALE_CW is an IEEE 754 representable f64 constant"),
            dx,
            dy,
            dz,
            filter_width,
        }
    }

    /// Create a WALE model with a custom coefficient and grid spacing.
    pub fn with_constant(c_w: T, dx: T, dy: T, dz: T) -> Self {
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            c_w,
            dx,
            dy,
            dz,
            filter_width,
        }
    }

    #[inline]
    fn wale_invariants(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> (T, T) {
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
        let strain_energy = symmetric_contract(strain, strain);
        let gradient_square = matrix_square(&gradient);
        let sd = symmetric_trace_free_part(&gradient_square);
        let sd_energy = symmetric_contract(sd, sd);
        (strain_energy, sd_energy)
    }

    #[inline]
    fn wale_viscosity_at(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        let (strain_energy, sd_energy) = self.wale_invariants(flow, i, j, k);
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        if sd_energy <= eps {
            return T::zero();
        }

        let three_over_two =
            <T as FromPrimitive>::from_f64(1.5).expect("1.5 is representable in IEEE 754");
        let five_over_two =
            <T as FromPrimitive>::from_f64(2.5).expect("2.5 is representable in IEEE 754");
        let five_over_four =
            <T as FromPrimitive>::from_f64(1.25).expect("1.25 is representable in IEEE 754");
        let c_delta = self.c_w * self.filter_width;
        let numerator = num_traits::Float::powf(sd_energy, three_over_two);
        let denominator = num_traits::Float::powf(strain_energy, five_over_two)
            + num_traits::Float::powf(sd_energy, five_over_four);

        if denominator <= eps {
            T::zero()
        } else {
            c_delta * c_delta * numerator / denominator
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    Default for WaleModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    TurbulenceModel<T> for WaleModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let mut viscosity = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    viscosity.push(self.wale_viscosity_at(flow_field, i, j, k));
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
        "WALE"
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
        let model = WaleModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);

        assert!(viscosity.iter().all(|&value| value == 0.0));
    }

    #[test]
    fn simple_shear_vanishes() {
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |_, y, _| Vector3::new(y, 0.0, 0.0));

        let model = WaleModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = 13;

        assert_relative_eq!(viscosity[center], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn linear_field_matches_closed_form_value() {
        let dx = 0.5;
        let dy = 1.5;
        let dz = 2.0;
        let mut flow = FlowField::<f64>::new(3, 3, 3);
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

        let model = WaleModel::<f64>::with_filter_width(dx, dy, dz);
        let viscosity = model.turbulent_viscosity(&flow);

        let strain = [
            [
                gradient[0][0],
                0.5 * (gradient[0][1] + gradient[1][0]),
                0.5 * (gradient[0][2] + gradient[2][0]),
            ],
            [
                0.5 * (gradient[1][0] + gradient[0][1]),
                gradient[1][1],
                0.5 * (gradient[1][2] + gradient[2][1]),
            ],
            [
                0.5 * (gradient[2][0] + gradient[0][2]),
                0.5 * (gradient[2][1] + gradient[1][2]),
                gradient[2][2],
            ],
        ];
        let strain_energy = {
            let mut sum = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let value = strain[i][j];
                    sum += value * value;
                }
            }
            sum
        };

        let gradient_square = {
            let mut square = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    let mut value = 0.0;
                    for k in 0..3 {
                        value += gradient[i][k] * gradient[k][j];
                    }
                    square[i][j] = value;
                }
            }
            square
        };

        let trace = gradient_square[0][0] + gradient_square[1][1] + gradient_square[2][2];
        let sd = [
            [
                gradient_square[0][0] - trace / 3.0,
                0.5 * (gradient_square[0][1] + gradient_square[1][0]),
                0.5 * (gradient_square[0][2] + gradient_square[2][0]),
            ],
            [
                0.5 * (gradient_square[1][0] + gradient_square[0][1]),
                gradient_square[1][1] - trace / 3.0,
                0.5 * (gradient_square[1][2] + gradient_square[2][1]),
            ],
            [
                0.5 * (gradient_square[2][0] + gradient_square[0][2]),
                0.5 * (gradient_square[2][1] + gradient_square[1][2]),
                gradient_square[2][2] - trace / 3.0,
            ],
        ];
        let sd_energy = {
            let mut sum = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    let value = sd[i][j];
                    sum += value * value;
                }
            }
            sum
        };

        let c_delta = model.c_w * model.filter_width;
        let expected = if strain_energy > 0.0 && sd_energy > 0.0 {
            let numerator = sd_energy.powf(1.5);
            let denominator = strain_energy.powf(2.5) + sd_energy.powf(1.25);
            c_delta * c_delta * numerator / denominator
        } else {
            0.0
        };

        assert_relative_eq!(viscosity[13], expected, epsilon = 1e-12);
    }
}
