//! Momentum equation coefficients

use super::solver::MomentumComponent;
use crate::discretization::extended_stencil::{ExtendedStencilScheme, QuickScheme};
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
                        let alpha = T::from_f64(relaxation_factor).unwrap_or_else(|| {
                            T::from_f64(0.7).unwrap_or_else(T::one)
                        });

                        // Compute QUICK-upwind correction for X-direction
                        let quick_correction_x = if i >= 2 && i < nx - 2 {
                            Self::compute_quick_correction_x(
                                i, j, u, rho, dy, &fields, component,
                            )
                        } else {
                            T::zero()
                        };

                        // Compute QUICK-upwind correction for Y-direction
                        let quick_correction_y = if j >= 2 && j < ny - 2 {
                            Self::compute_quick_correction_y(
                                i, j, v, rho, dx, &fields, component,
                            )
                        } else {
                            T::zero()
                        };

                        // Add deferred correction to source with under-relaxation
                        // Source correction = α * (QUICK_flux - upwind_flux)
                        if let Some(source) = coeffs.source.at_mut(i, j) {
                            let correction = alpha * (quick_correction_x + quick_correction_y);
                            
                            // Debug logging for first interior cell
                            #[cfg(debug_assertions)]
                            if i == 2 && j == ny / 2 {
                                tracing::debug!(
                                    "Deferred correction at center: x_corr={:.3e}, y_corr={:.3e}, total={:.3e}, alpha={:.3}",
                                    quick_correction_x, quick_correction_y, correction, alpha
                                );
                            }
                            
                            *source = *source + correction;
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

    /// Compute QUICK correction for X-direction convection
    ///
    /// Returns the difference between QUICK flux and upwind flux
    /// Reference: Leonard (1979), Patankar (1980) §5.4.3
    fn compute_quick_correction_x(
        i: usize,
        j: usize,
        u: T,
        rho: T,
        dy: T,
        fields: &SimulationFields<T>,
        component: MomentumComponent,
    ) -> T {
        let quick = QuickScheme;

        // Get velocity values for 5-point stencil [i-2, i-1, i, i+1, i+2]
        let phi_values: [T; 5] = match component {
            MomentumComponent::U => [
                fields.u.at(i - 2, j),
                fields.u.at(i - 1, j),
                fields.u.at(i, j),
                fields.u.at(i + 1, j),
                fields.u.at(i + 2, j),
            ],
            MomentumComponent::V => [
                fields.v.at(i - 2, j),
                fields.v.at(i - 1, j),
                fields.v.at(i, j),
                fields.v.at(i + 1, j),
                fields.v.at(i + 2, j),
            ],
        };

        // QUICK face value at east face (between i and i+1)
        let phi_e_quick = quick.face_value(&phi_values, u, None);

        // Upwind face value at east face
        let phi_e_upwind = if u > T::zero() {
            phi_values[2] // From cell i (current)
        } else {
            phi_values[3] // From cell i+1 (downstream)
        };

        // Flux correction = mass_flux * (φ_QUICK - φ_upwind)
        let mass_flux = rho * u * dy;
        mass_flux * (phi_e_quick - phi_e_upwind)
    }

    /// Compute QUICK correction for Y-direction convection
    ///
    /// Returns the difference between QUICK flux and upwind flux
    /// Reference: Leonard (1979), Patankar (1980) §5.4.3
    fn compute_quick_correction_y(
        i: usize,
        j: usize,
        v: T,
        rho: T,
        dx: T,
        fields: &SimulationFields<T>,
        component: MomentumComponent,
    ) -> T {
        let quick = QuickScheme;

        // Get velocity values for 5-point stencil [j-2, j-1, j, j+1, j+2]
        let phi_values: [T; 5] = match component {
            MomentumComponent::U => [
                fields.u.at(i, j - 2),
                fields.u.at(i, j - 1),
                fields.u.at(i, j),
                fields.u.at(i, j + 1),
                fields.u.at(i, j + 2),
            ],
            MomentumComponent::V => [
                fields.v.at(i, j - 2),
                fields.v.at(i, j - 1),
                fields.v.at(i, j),
                fields.v.at(i, j + 1),
                fields.v.at(i, j + 2),
            ],
        };

        // QUICK face value at north face (between j and j+1)
        let phi_n_quick = quick.face_value(&phi_values, v, None);

        // Upwind face value at north face
        let phi_n_upwind = if v > T::zero() {
            phi_values[2] // From cell j (current)
        } else {
            phi_values[3] // From cell j+1 (downstream)
        };

        // Flux correction = mass_flux * (φ_QUICK - φ_upwind)
        let mass_flux = rho * v * dx;
        mass_flux * (phi_n_quick - phi_n_upwind)
    }
}
