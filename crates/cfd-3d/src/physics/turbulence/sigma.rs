//! Sigma subgrid-scale model for LES.
//!
//! # Theorem — Sigma SGS Model (Nicoud et al. 2011)
//!
//! The eddy viscosity is computed from the singular values σ₁ ≥ σ₂ ≥ σ₃ ≥ 0
//! of the velocity gradient tensor **g** = ∇**u**:
//!
//! ```text
//! νₜ = (C_σ Δ)² · σ₃(σ₁ − σ₂)(σ₂ − σ₃) / σ₁²
//! ```
//!
//! **Properties** (proved in Nicoud et al. §II):
//! - νₜ = 0 in solid-body rotation and irrotational flow (no SGS dissipation in these regimes)
//! - νₜ = 0 at solid walls (correct near-wall scaling without damping functions)
//! - Strictly non-negative for any velocity gradient
//!
//! ## Algorithm
//!
//! ```text
//! 1. For each grid point (i,j,k):
//!    a. Compute velocity gradient tensor G_ij = ∂u_i/∂x_j  (central differences)
//!    b. Compute G^T G and its eigenvalues λ₁ ≥ λ₂ ≥ λ₃ ≥ 0
//!    c. σ_k = sqrt(λ_k)
//!    d. νₜ = (C_σ·Δ)² · σ₃(σ₁−σ₂)(σ₂−σ₃) / max(σ₁², ε)
//! 2. Return νₜ field over all points.
//! ```
//!
//! ## References
//!
//! - Nicoud, F., Baya Toda, H., Cabrit, O., Bose, S. & Lee, J. (2011).
//!   "Using singular values to build a subgrid-scale model for LES."
//!   *Phys. Fluids* 23(8):085106.

use cfd_core::physics::fluid_dynamics::fields::FlowField;
use cfd_core::physics::fluid_dynamics::turbulence::TurbulenceModel;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::constants::SIGMA_C;
use super::field_ops::{strain_magnitude, velocity_gradient_tensor};

/// Sigma subgrid-scale model for LES (Nicoud et al. 2011).
///
/// Computes eddy viscosity from singular values of the velocity gradient tensor.
/// Automatically satisfies νₜ = 0 at walls without ad hoc damping.
#[derive(Debug, Clone)]
pub struct SigmaModel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Sigma model constant C_σ = 1.35 (Nicoud et al. 2011).
    pub c_sigma: T,
    /// Physical grid spacing in the x direction [m].
    pub dx: T,
    /// Physical grid spacing in the y direction [m].
    pub dy: T,
    /// Physical grid spacing in the z direction [m].
    pub dz: T,
    /// Physical LES filter width Δ = (dx·dy·dz)^(1/3) [m].
    pub filter_width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    SigmaModel<T>
{
    /// Create a Sigma model with default constant and unit filter width.
    pub fn new() -> Self {
        Self {
            c_sigma: <T as FromPrimitive>::from_f64(SIGMA_C)
                .expect("SIGMA_C is an IEEE 754 representable f64 constant"),
            dx: T::one(),
            dy: T::one(),
            dz: T::one(),
            filter_width: T::one(),
        }
    }

    /// Create a Sigma model with physically correct filter width Δ = (dx·dy·dz)^(1/3).
    ///
    /// # Arguments
    /// * `dx`, `dy`, `dz` — physical cell dimensions [m]
    pub fn with_filter_width(dx: T, dy: T, dz: T) -> Self {
        let one_third = T::one() / (T::one() + T::one() + T::one());
        let filter_width = num_traits::Float::powf(dx * dy * dz, one_third);
        Self {
            c_sigma: <T as FromPrimitive>::from_f64(SIGMA_C)
                .expect("SIGMA_C is an IEEE 754 representable f64 constant"),
            dx,
            dy,
            dz,
            filter_width,
        }
    }

    /// Compute the exact eigenvalues of a real symmetric 3×3 matrix.
    ///
    /// # Proof sketch
    /// A real symmetric matrix has three real eigenvalues and an orthogonal
    /// eigenbasis. For `G^T G`, those eigenvalues are exactly `σ_i²`, where
    /// `σ_i` are the singular values of `G`. The closed-form 3×3 solver below
    /// evaluates the characteristic polynomial without heap allocation, so the
    /// singular values derived from it are exact up to floating-point roundoff.
    fn symmetric_eigenvalues_3x3(matrix: &[[T; 3]; 3]) -> [T; 3] {
        let two = T::one() + T::one();
        let three = two + T::one();
        let six = three + three;
        let trace = matrix[0][0] + matrix[1][1] + matrix[2][2];
        let mean = trace / three;
        let b00 = matrix[0][0] - mean;
        let b11 = matrix[1][1] - mean;
        let b22 = matrix[2][2] - mean;
        let b01 = matrix[0][1];
        let b02 = matrix[0][2];
        let b12 = matrix[1][2];
        let p2 =
            (b00 * b00 + b11 * b11 + b22 * b22 + two * (b01 * b01 + b02 * b02 + b12 * b12)) / six;
        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");

        if p2 <= eps {
            return [mean, mean, mean];
        }

        let p = <T as num_traits::Float>::sqrt(p2);
        let p3 = p * p2;
        if p3 <= eps {
            return [mean, mean, mean];
        }

        let det_b = b00 * (b11 * b22 - b12 * b12) - b01 * (b01 * b22 - b12 * b02)
            + b02 * (b01 * b12 - b11 * b02);
        let denom = two * p3;
        let mut r = det_b / denom;
        let one = T::one();
        if r > one {
            r = one;
        } else if r < -one {
            r = -one;
        }

        let phi = <T as num_traits::Float>::acos(r) / three;
        let pi = <T as FromPrimitive>::from_f64(std::f64::consts::PI)
            .expect("PI is an IEEE 754 representable f64 constant");
        let two_pi_over_three = two * pi / three;
        let four_pi_over_three = two_pi_over_three * two;
        let mut eigenvalues = [
            mean + two * p * num_traits::Float::cos(phi),
            mean + two * p * num_traits::Float::cos(phi + two_pi_over_three),
            mean + two * p * num_traits::Float::cos(phi + four_pi_over_three),
        ];

        eigenvalues.sort_by(|a, b| b.partial_cmp(a).expect("Sigma eigenvalues must be finite"));
        eigenvalues
    }

    /// Compute the local Sigma-model viscosity and strain magnitude.
    fn point_metrics(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> (T, T) {
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
        let strain_mag = strain_magnitude(&gradient);

        let mut gtg = [[T::zero(); 3]; 3];
        for ii in 0..3 {
            for jj in 0..3 {
                for kk in 0..3 {
                    gtg[ii][jj] += gradient[kk][ii] * gradient[kk][jj];
                }
            }
        }

        let eigenvalues = Self::symmetric_eigenvalues_3x3(&gtg);
        let sigma1 = if eigenvalues[0] > T::zero() {
            <T as num_traits::Float>::sqrt(eigenvalues[0])
        } else {
            T::zero()
        };
        let sigma2 = if eigenvalues[1] > T::zero() {
            <T as num_traits::Float>::sqrt(eigenvalues[1])
        } else {
            T::zero()
        };
        let sigma3 = if eigenvalues[2] > T::zero() {
            <T as num_traits::Float>::sqrt(eigenvalues[2])
        } else {
            T::zero()
        };

        let eps = <T as FromPrimitive>::from_f64(1e-30)
            .expect("1e-30 is an IEEE 754 representable f64 constant");
        let denominator = sigma1 * sigma1;
        let sigma_product = if denominator > eps {
            let product = sigma3 * (sigma1 - sigma2) * (sigma2 - sigma3) / denominator;
            if product > T::zero() {
                product
            } else {
                T::zero()
            }
        } else {
            T::zero()
        };

        let c_delta = self.c_sigma * self.filter_width;
        let viscosity = c_delta * c_delta * sigma_product;
        (viscosity, strain_mag)
    }

    /// Compute the Sigma eddy viscosity at a single grid point.
    fn sigma_viscosity_at(&self, flow: &FlowField<T>, i: usize, j: usize, k: usize) -> T {
        self.point_metrics(flow, i, j, k).0
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    Default for SigmaModel<T>
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + num_traits::Float>
    TurbulenceModel<T> for SigmaModel<T>
{
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T> {
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let mut viscosity = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    viscosity.push(self.sigma_viscosity_at(flow_field, i, j, k));
                }
            }
        }
        viscosity
    }

    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T> {
        // Estimate SGS TKE from the resolved eddy viscosity and strain magnitude.
        let (nx, ny, nz) = flow_field.velocity.dimensions;
        let mut kinetic_energy = Vec::with_capacity(nx * ny * nz);
        let delta = if self.filter_width > T::zero() {
            self.filter_width
        } else {
            <T as FromPrimitive>::from_f64(1e-30)
                .expect("1e-30 is an IEEE 754 representable f64 constant")
        };

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let (viscosity, strain_mag) = self.point_metrics(flow_field, i, j, k);
                    kinetic_energy.push(viscosity * strain_mag / delta);
                }
            }
        }

        kinetic_energy
    }

    fn name(&self) -> &'static str {
        "Sigma"
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
    fn diagonal_strain_produces_exact_sigma_viscosity() {
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |x, y, z| Vector3::new(x, 2.0 * y, 3.0 * z));

        let model = SigmaModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = 13;
        let expected = model.c_sigma * model.c_sigma / 9.0;

        assert_relative_eq!(viscosity[center], expected, epsilon = 1e-12);
    }

    #[test]
    fn pure_rotation_vanishes_for_viscosity_and_tke() {
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |x, y, _z| Vector3::new(-y, x, 0.0));

        let model = SigmaModel::<f64>::with_filter_width(1.0, 1.0, 1.0);
        let viscosity = model.turbulent_viscosity(&flow);
        let kinetic_energy = model.turbulent_kinetic_energy(&flow);
        let center = 13;

        assert_relative_eq!(viscosity[center], 0.0, epsilon = 1e-12);
        assert_relative_eq!(kinetic_energy[center], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn anisotropic_spacing_matches_linear_strain_reference() {
        let dx = 1.0;
        let dy = 2.0;
        let dz = 4.0;
        let mut flow = FlowField::<f64>::new(3, 3, 3);
        fill_velocity_field(&mut flow, |i, j, k| {
            let x = i * dx;
            let y = j * dy;
            let z = k * dz;
            Vector3::new(x, 2.0 * y, 3.0 * z)
        });

        let model = SigmaModel::<f64>::with_filter_width(dx, dy, dz);
        let viscosity = model.turbulent_viscosity(&flow);
        let center = 13;
        let expected =
            model.c_sigma * model.c_sigma * model.filter_width * model.filter_width / 9.0;

        assert_relative_eq!(viscosity[center], expected, epsilon = 1e-12);
    }
}
