//! Momentum equation coefficients

use super::solver::MomentumComponent;
use crate::fields::{Field2D, SimulationFields};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Coefficients for momentum discretization
#[derive(Debug, Clone)]
pub struct MomentumCoefficients<T: RealField + Copy> {
    /// Central coefficient (aP)
    pub ap: Field2D<T>,
    /// East coefficient (aE)
    pub ae: Field2D<T>,
    /// West coefficient (aW)
    pub aw: Field2D<T>,
    /// North coefficient (aN)
    pub an: Field2D<T>,
    /// South coefficient (aS)
    pub as_: Field2D<T>,
    /// Source term
    pub source: Field2D<T>,
}

impl<T: RealField + Copy + FromPrimitive> MomentumCoefficients<T> {
    /// Compute momentum equation coefficients
    pub fn compute(
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        dt: T,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
    ) -> cfd_core::Result<Self> {
        let mut coeffs = Self {
            ap: Field2D::new(nx, ny, T::zero()),
            ae: Field2D::new(nx, ny, T::zero()),
            aw: Field2D::new(nx, ny, T::zero()),
            an: Field2D::new(nx, ny, T::zero()),
            as_: Field2D::new(nx, ny, T::zero()),
            source: Field2D::new(nx, ny, T::zero()),
        };

        // Compute diffusion coefficients
        let dx2 = dx * dx;
        let dy2 = dy * dy;

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let mu = fields.viscosity.at(i, j);

                // Diffusion coefficients
                coeffs.ae.set(i, j, mu / dx2);
                coeffs.aw.set(i, j, mu / dx2);
                coeffs.an.set(i, j, mu / dy2);
                coeffs.as_.set(i, j, mu / dy2);

                // Convection coefficients (using upwind)
                let (u, v) = match component {
                    MomentumComponent::U => {
                        let u = fields.u.at(i, j);
                        let v = (fields.v.at(i, j) + fields.v.at(i + 1, j)) / (T::one() + T::one());
                        (u, v)
                    }
                    MomentumComponent::V => {
                        let u = (fields.u.at(i, j) + fields.u.at(i, j + 1)) / (T::one() + T::one());
                        let v = fields.v.at(i, j);
                        (u, v)
                    }
                };

                // Add convection to coefficients (upwind scheme)
                if u > T::zero() {
                    coeffs.ae.set(i, j, coeffs.ae.at(i, j) + u / dx);
                    coeffs.ap.set(i, j, coeffs.ap.at(i, j) + u / dx);
                } else {
                    coeffs.aw.set(i, j, coeffs.aw.at(i, j) - u / dx);
                    coeffs.ap.set(i, j, coeffs.ap.at(i, j) - u / dx);
                }

                if v > T::zero() {
                    coeffs.an.set(i, j, coeffs.an.at(i, j) + v / dy);
                    coeffs.ap.set(i, j, coeffs.ap.at(i, j) + v / dy);
                } else {
                    coeffs.as_.set(i, j, coeffs.as_.at(i, j) - v / dy);
                    coeffs.ap.set(i, j, coeffs.ap.at(i, j) - v / dy);
                }

                // Central coefficient (including time term)
                let ap_sum = coeffs.ae.at(i, j)
                    + coeffs.aw.at(i, j)
                    + coeffs.an.at(i, j)
                    + coeffs.as_.at(i, j);
                coeffs.ap.set(i, j, ap_sum + fields.density.at(i, j) / dt);

                // Source term (including previous time step)
                let old_vel = match component {
                    MomentumComponent::U => fields.u.at(i, j),
                    MomentumComponent::V => fields.v.at(i, j),
                };
                coeffs
                    .source
                    .set(i, j, fields.density.at(i, j) * old_vel / dt);
            }
        }

        Ok(coeffs)
    }
}
