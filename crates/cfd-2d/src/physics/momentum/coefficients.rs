//! Momentum equation coefficients with advanced convection schemes
//!
//! This module implements the finite volume discretization of the momentum equations
//! with support for multiple convection schemes:
//!
//! * **Upwind** (First-order): Unconditionally stable but dissipative (Pe > 2)
//! * **Deferred Correction**: Combines upwind stability with QUICK accuracy
//!
//! # References
//! * Patankar, S.V. (1980). "Numerical Heat Transfer and Fluid Flow", §5.4
//! * Leonard, B.P. (1979). "A stable and accurate convective modelling procedure"
//!
//! # High-Peclet Number Flows
//!
//! For flows with Peclet number Pe = ρuL/μ >> 2, convection schemes introduce
//! numerical diffusion. This is a fundamental CFD challenge, not a solver bug.
//!
//! **Mitigation Strategies:**
//! 1. Use deferred correction with relaxation factor 0.7-0.9
//! 2. Apply velocity under-relaxation (0.5-0.8)
//! 3. Consider TVD limiters for shock-like flows
//! 4. For fully-developed flows (∂u/∂x ≈ 0), use pure diffusion solver
//!
//! # Example
//! ```ignore
//! use cfd_2d::physics::momentum::{ConvectionScheme, MomentumSolver};
//!
//! // Default: Deferred correction with α=0.7
//! let mut solver = MomentumSolver::new(&grid);
//!
//! // For high-Peclet flows, increase relaxation
//! solver.set_convection_scheme(ConvectionScheme::DeferredCorrectionQuick {
//!     relaxation_factor: 0.9
//! });
//! solver.set_velocity_relaxation(0.8);
//!
//! // For debugging or comparison, use pure upwind
//! solver.set_convection_scheme(ConvectionScheme::Upwind);
//! ```

use super::coefficient_corrections::{
    compute_quick_correction_x, compute_quick_correction_y, compute_tvd_correction_x,
    compute_tvd_correction_y,
};
use super::solver::MomentumComponent;
use super::tvd_limiters::{Minmod, Superbee, VanLeer};
use crate::fields::{Field2D, SimulationFields};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Coefficients for momentum discretization
#[derive(Debug, Clone)]
pub struct MomentumCoefficients<T: RealField + Copy> {
    /// Central coefficient (aP)
    pub ap: Field2D<T>,
    /// Consistent diagonal coefficient (aP - sum(aNb)) for SIMPLEC
    pub ap_consistent: Field2D<T>,
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

/// Convection discretization strategy for momentum equations
#[derive(Debug, Clone, Copy)]
pub enum ConvectionScheme {
    /// First-order upwind (stable, dissipative)
    Upwind,
    /// Deferred correction with QUICK explicit part (Patankar 1980 §5.4.3)
    /// Combines stability of upwind with accuracy of QUICK
    DeferredCorrectionQuick {
        /// Under-relaxation factor (0.5-0.8 recommended, default 0.7)
        relaxation_factor: f64,
    },
    /// TVD scheme with Superbee limiter (Roe 1986)
    ///
    /// Excellent for high-Peclet flows with shock-like features.
    /// Most compressive limiter, best shock resolution.
    ///
    /// # References
    /// * Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
    /// * Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"
    TvdSuperbee {
        /// Under-relaxation factor (0.7-0.9 recommended for high-Pe flows)
        relaxation_factor: f64,
    },
    /// TVD scheme with Van Leer limiter (Van Leer 1974)
    ///
    /// Smooth limiter, good balance of accuracy and stability.
    /// Recommended for general high-Peclet flows.
    ///
    /// # Reference
    /// Van Leer, B. (1974). "Towards the Ultimate Conservative Difference Scheme"
    TvdVanLeer {
        /// Under-relaxation factor (0.7-0.9 recommended for high-Pe flows)
        relaxation_factor: f64,
    },
    /// TVD scheme with Minmod limiter (Roe 1986)
    ///
    /// Most diffusive TVD limiter, extremely stable.
    /// Use for very difficult high-Peclet flows.
    ///
    /// # Reference
    /// Roe, P.L. (1986). "Characteristic-Based Schemes for the Euler Equations"
    TvdMinmod {
        /// Under-relaxation factor (0.7-0.9 recommended for high-Pe flows)
        relaxation_factor: f64,
    },
}

impl Default for ConvectionScheme {
    fn default() -> Self {
        Self::DeferredCorrectionQuick {
            relaxation_factor: 0.7,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> MomentumCoefficients<T> {
    /// Compute momentum equation coefficients with advanced convection scheme
    ///
    /// # Arguments
    /// * `scheme` - Convection discretization scheme (upwind or deferred correction)
    ///
    /// # References
    /// * Patankar (1980) "Numerical Heat Transfer and Fluid Flow" §5.4.3
    /// * Leonard (1979) "A stable and accurate convective modelling procedure"
    pub fn compute(
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        dt: T,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        scheme: ConvectionScheme,
    ) -> cfd_core::error::Result<Self> {
        let mut coeffs = Self {
            ap: Field2D::new(nx, ny, T::zero()),
            ap_consistent: Field2D::new(nx, ny, T::zero()),
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
        let diff_coeff_ew = dy / dx; // East-West diffusion scaling
        let diff_coeff_ns = dx / dy; // North-South diffusion scaling

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
                // Deferred correction approach (Patankar 1980 §5.4.3):
                // - Implicit part: Use first-order upwind for unconditional stability
                // - Explicit part: Add QUICK correction to source term for accuracy
                // - Combined: Stable iteration with high-order accuracy
                //
                // Face velocity (for determining flow direction) vs transported quantity:
                // - Face velocity determines upwind direction
                // - Transported quantity (φ) is taken from upwind side or QUICK interpolation

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

                let rho = fields.density.at(i, j);

                // Implicit part: Always use upwind for stability (added to coefficients)
                let mass_flux_x = rho * u * dy;
                let mass_flux_y = rho * v * dx;

                if u > T::zero() {
                    let aw_val = coeffs.aw.at(i, j);
                    if let Some(aw) = coeffs.aw.at_mut(i, j) {
                        *aw = aw_val + mass_flux_x;
                    }
                } else {
                    let ae_val = coeffs.ae.at(i, j);
                    if let Some(ae) = coeffs.ae.at_mut(i, j) {
                        *ae = ae_val - mass_flux_x;
                    }
                }

                if v > T::zero() {
                    let as_val = coeffs.as_.at(i, j);
                    if let Some(as_) = coeffs.as_.at_mut(i, j) {
                        *as_ = as_val + mass_flux_y;
                    }
                } else {
                    let an_val = coeffs.an.at(i, j);
                    if let Some(an) = coeffs.an.at_mut(i, j) {
                        *an = an_val - mass_flux_y;
                    }
                }

                // Explicit part: Deferred correction with QUICK scheme
                // This is added to the source term after upwind is applied implicitly
                match scheme {
                    ConvectionScheme::Upwind => {
                        // Pure upwind - no correction needed
                    }
                    ConvectionScheme::DeferredCorrectionQuick { relaxation_factor } => {
                        let alpha = T::from_f64(relaxation_factor)
                            .unwrap_or_else(|| T::from_f64(0.7).unwrap_or_else(T::one));

                        // Compute QUICK-upwind correction for X-direction
                        let quick_correction_x = if i >= 2 && i < nx - 2 {
                            compute_quick_correction_x(i, j, u, rho, dy, fields, component)
                        } else {
                            T::zero()
                        };

                        // Compute QUICK-upwind correction for Y-direction
                        let quick_correction_y = if j >= 2 && j < ny - 2 {
                            compute_quick_correction_y(i, j, v, rho, dx, fields, component)
                        } else {
                            T::zero()
                        };

                        // Add deferred correction to source with under-relaxation
                        // Source correction = α * (QUICK_flux - upwind_flux)
                        if let Some(source) = coeffs.source.at_mut(i, j) {
                            let correction = alpha * (quick_correction_x + quick_correction_y);
                            *source += correction;
                        }
                    }
                    ConvectionScheme::TvdSuperbee { relaxation_factor } => {
                        let alpha = T::from_f64(relaxation_factor)
                            .unwrap_or_else(|| T::from_f64(0.8).unwrap_or_else(T::one));
                        let limiter = Superbee;

                        let tvd_correction_x = if i >= 1 && i < nx - 1 {
                            compute_tvd_correction_x(i, j, u, rho, dy, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        let tvd_correction_y = if j >= 1 && j < ny - 1 {
                            compute_tvd_correction_y(i, j, v, rho, dx, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        if let Some(source) = coeffs.source.at_mut(i, j) {
                            let correction = alpha * (tvd_correction_x + tvd_correction_y);
                            *source += correction;
                        }
                    }
                    ConvectionScheme::TvdVanLeer { relaxation_factor } => {
                        let alpha = T::from_f64(relaxation_factor)
                            .unwrap_or_else(|| T::from_f64(0.8).unwrap_or_else(T::one));
                        let limiter = VanLeer;

                        let tvd_correction_x = if i >= 1 && i < nx - 1 {
                            compute_tvd_correction_x(i, j, u, rho, dy, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        let tvd_correction_y = if j >= 1 && j < ny - 1 {
                            compute_tvd_correction_y(i, j, v, rho, dx, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        if let Some(source) = coeffs.source.at_mut(i, j) {
                            let correction = alpha * (tvd_correction_x + tvd_correction_y);
                            *source += correction;
                        }
                    }
                    ConvectionScheme::TvdMinmod { relaxation_factor } => {
                        let alpha = T::from_f64(relaxation_factor)
                            .unwrap_or_else(|| T::from_f64(0.8).unwrap_or_else(T::one));
                        let limiter = Minmod;

                        let tvd_correction_x = if i >= 1 && i < nx - 1 {
                            compute_tvd_correction_x(i, j, u, rho, dy, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        let tvd_correction_y = if j >= 1 && j < ny - 1 {
                            compute_tvd_correction_y(i, j, v, rho, dx, fields, component, &limiter)
                        } else {
                            T::zero()
                        };

                        if let Some(source) = coeffs.source.at_mut(i, j) {
                            let correction = alpha * (tvd_correction_x + tvd_correction_y);
                            *source += correction;
                        }
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

                // SIMPLEC consistent diagonal: aP - sum(aNb)
                // This represents the part of aP that does not depend on neighbor velocities.
                // For SIMPLEC consistency, we use this in the velocity correction equation.
                if let Some(ap_c) = coeffs.ap_consistent.at_mut(i, j) {
                    // Start with the time term
                    let time_term = rho * volume / dt;

                    // Add a small epsilon to ensure we never divide by zero in steady state
                    // This is standard practice in SIMPLEC implementations
                    let eps = T::from_f64(1e-10).unwrap_or_else(T::zero);

                    *ap_c = time_term + eps;
                }

                // Source term (including previous time step and pressure gradient)
                let previous_velocity = match component {
                    MomentumComponent::U => fields.u.at(i, j),
                    MomentumComponent::V => fields.v.at(i, j),
                };

                // Body force term
                let body_force = match component {
                    MomentumComponent::U => fields.force_u.at(i, j),
                    MomentumComponent::V => fields.force_v.at(i, j),
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
                    // RHS = ρ * V * u_old / dt + pressure_force + body_force
                    // where pressure_force = -∂p/∂x * V (total force on control volume)
                    // and V = dx * dy (cell volume in 2D, per unit depth)
                    //
                    // NOTE: Use += to ACCUMULATE with any deferred correction (QUICK/TVD)
                    // that was already added to the source term above. For Upwind scheme,
                    // source is still zero here, so += is equivalent to =.
                    let volume = dx * dy;
                    *source += fields.density.at(i, j) * volume * previous_velocity / dt
                        + pressure_gradient * volume
                        + body_force * volume;
                }
            }
        }

        Ok(coeffs)
    }
}
