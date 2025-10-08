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

        // Compute diffusion coefficients with proper finite volume scaling
        // Diffusion flux through a face = -μ * A * ∂φ/∂n / Δn
        // For a control volume, coefficient in discretized equation:
        // coeff_e = μ * A_e / Δx_e = μ * dy / dx (for east/west faces)
        // coeff_n = μ * A_n / Δy_n = μ * dx / dy (for north/south faces)
        let diff_coeff_ew = dy / dx;  // East-West diffusion scaling
        let diff_coeff_ns = dx / dy;  // North-South diffusion scaling

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let mu = fields.viscosity.at(i, j);

                // Diffusion coefficients with finite volume scaling
                // East/West: μ * (dy / dx)
                // North/South: μ * (dx / dy)
                if let Some(ae) = coeffs.ae.at_mut(i, j) {
                    *ae = mu * diff_coeff_ew;
                }
                if let Some(aw) = coeffs.aw.at_mut(i, j) {
                    *aw = mu * diff_coeff_ew;
                }
                if let Some(an) = coeffs.an.at_mut(i, j) {
                    *an = mu * diff_coeff_ns;
                }
                if let Some(as_) = coeffs.as_.at_mut(i, j) {
                    *as_ = mu * diff_coeff_ns;
                }

                // Convection coefficients for momentum equation
                // The nonlinear convective term is: ρ * (u * ∂φ/∂x + v * ∂φ/∂y)
                // where φ is the velocity component being solved (u or v)
                // 
                // Face velocity (for determining flow direction) vs transported quantity:
                // - Face velocity determines upwind direction
                // - Transported quantity (φ) is taken from upwind side
                //
                // For u-momentum: convective flux through east face = ρ * u_e * φ_e * A_e
                //   where u_e is velocity at east face, φ is u-velocity being transported
                // For v-momentum: similar but with appropriate velocity components

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
                // Using Patankar's formulation: a_W = D_w + max(F_w, 0)
                // where F_w = ρ * u_w * A_w is the mass flow rate through west face
                // and D_w = μ * A_w / δx is the diffusion conductance
                let rho = fields.density.at(i, j);
                
                // X-direction convection
                // Mass flux through face = ρ * u * A = ρ * u * dy
                let mass_flux_x = rho * u * dy;
                
                if u > T::zero() {
                    // Flow to right (W→P): add mass_flux to a_W
                    let aw_val = coeffs.aw.at(i, j);
                    if let Some(aw) = coeffs.aw.at_mut(i, j) {
                        *aw = aw_val + mass_flux_x;
                    }
                } else {
                    // Flow to left (E→P): add |mass_flux| to a_E
                    let ae_val = coeffs.ae.at(i, j);
                    if let Some(ae) = coeffs.ae.at_mut(i, j) {
                        *ae = ae_val - mass_flux_x;  // mass_flux_x is negative
                    }
                }

                // Y-direction convection  
                // Mass flux through face = ρ * v * A = ρ * v * dx
                let mass_flux_y = rho * v * dx;
                
                if v > T::zero() {
                    // Flow up (S→P): add mass_flux to a_S
                    let as_val = coeffs.as_.at(i, j);
                    if let Some(as_) = coeffs.as_.at_mut(i, j) {
                        *as_ = as_val + mass_flux_y;
                    }
                } else {
                    // Flow down (N→P): add |mass_flux| to a_N
                    let an_val = coeffs.an.at(i, j);
                    if let Some(an) = coeffs.an.at_mut(i, j) {
                        *an = an_val - mass_flux_y;  // mass_flux_y is negative
                    }
                }

                // Central coefficient (including time term)
                // a_P = sum of neighbor coefficients + ρ * V / dt
                // For incompressible flow, ΔF = 0, so no additional convection term needed
                let volume = dx * dy;
                let rho = fields.density.at(i, j);
                
                let ap_sum = coeffs.ae.at(i, j)
                    + coeffs.aw.at(i, j)
                    + coeffs.an.at(i, j)
                    + coeffs.as_.at(i, j);
                if let Some(ap) = coeffs.ap.at_mut(i, j) {
                    *ap = ap_sum + rho * volume / dt;
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
                    // RHS = ρ * V * u_old / dt + pressure_force
                    // where pressure_force = -∂p/∂x * V (total force on control volume)
                    // and V = dx * dy (cell volume in 2D, per unit depth)
                    let volume = dx * dy;
                    *source = fields.density.at(i, j) * volume * previous_velocity / dt 
                            + pressure_gradient * volume;
                }
            }
        }

        Ok(coeffs)
    }
}
