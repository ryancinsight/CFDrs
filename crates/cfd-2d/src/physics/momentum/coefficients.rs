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
    ) -> cfd_core::error::Result<Self> {
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
                if let Some(ae) = coeffs.ae.at_mut(i, j) {
                    *ae = mu / dx2;
                }
                if let Some(aw) = coeffs.aw.at_mut(i, j) {
                    *aw = mu / dx2;
                }
                if let Some(an) = coeffs.an.at_mut(i, j) {
                    *an = mu / dy2;
                }
                if let Some(as_) = coeffs.as_.at_mut(i, j) {
                    *as_ = mu / dy2;
                }

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
                    let ae_val = coeffs.ae.at(i, j);
                    let ap_val = coeffs.ap.at(i, j);
                    if let Some(ae) = coeffs.ae.at_mut(i, j) {
                        *ae = ae_val + u / dx;
                    }
                    if let Some(ap) = coeffs.ap.at_mut(i, j) {
                        *ap = ap_val + u / dx;
                    }
                } else {
                    let aw_val = coeffs.aw.at(i, j);
                    let ap_val = coeffs.ap.at(i, j);
                    if let Some(aw) = coeffs.aw.at_mut(i, j) {
                        *aw = aw_val - u / dx;
                    }
                    if let Some(ap) = coeffs.ap.at_mut(i, j) {
                        *ap = ap_val - u / dx;
                    }
                }

                if v > T::zero() {
                    let an_val = coeffs.an.at(i, j);
                    let ap_val = coeffs.ap.at(i, j);
                    if let Some(an) = coeffs.an.at_mut(i, j) {
                        *an = an_val + v / dy;
                    }
                    if let Some(ap) = coeffs.ap.at_mut(i, j) {
                        *ap = ap_val + v / dy;
                    }
                } else {
                    let as_val = coeffs.as_.at(i, j);
                    let ap_val = coeffs.ap.at(i, j);
                    if let Some(as_) = coeffs.as_.at_mut(i, j) {
                        *as_ = as_val - v / dy;
                    }
                    if let Some(ap) = coeffs.ap.at_mut(i, j) {
                        *ap = ap_val - v / dy;
                    }
                }

                // Central coefficient (including time term)
                let ap_sum = coeffs.ae.at(i, j)
                    + coeffs.aw.at(i, j)
                    + coeffs.an.at(i, j)
                    + coeffs.as_.at(i, j);
                if let Some(ap) = coeffs.ap.at_mut(i, j) {
                    *ap = ap_sum + fields.density.at(i, j) / dt;
                }

                // Source term (including previous time step and pressure gradient)
                let previous_velocity = match component {
                    MomentumComponent::U => fields.u.at(i, j),
                    MomentumComponent::V => fields.v.at(i, j),
                };

                // Calculate pressure gradient term for momentum equation
                // Momentum: ρ ∂u/∂t = μ ∇²u - ∂p/∂x
                // Rearranging: (ρ/dt + diffusion) * u = (ρ/dt) * u_old - ∂p/∂x
                // So RHS source term is: ρ * u_old / dt - ∂p/∂x
                let pressure_gradient = match component {
                    MomentumComponent::U => {
                        // -∂p/∂x using central difference
                        if i > 0 && i < nx - 1 {
                            -(fields.p.at(i + 1, j) - fields.p.at(i - 1, j))
                                / (T::one() + T::one())
                                / dx
                        } else {
                            T::zero()
                        }
                    }
                    MomentumComponent::V => {
                        // -∂p/∂y using central difference
                        if j > 0 && j < ny - 1 {
                            -(fields.p.at(i, j + 1) - fields.p.at(i, j - 1))
                                / (T::one() + T::one())
                                / dy
                        } else {
                            T::zero()
                        }
                    }
                };

                if let Some(source) = coeffs.source.at_mut(i, j) {
                    // RHS = ρ * u_old / dt + pressure_gradient_term
                    // where pressure_gradient_term = -∂p/∂x (already computed above)
                    *source = fields.density.at(i, j) * previous_velocity / dt + pressure_gradient;
                }
            }
        }

        Ok(coeffs)
    }
}
