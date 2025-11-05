//! Central difference schemes with stability analysis
//!
//! ## Von Neumann Stability Analysis
//!
//! ### Second-Order Central Difference
//!
//! For the 1D advection equation: ∂u/∂t + c ∂u/∂x = 0
//!
//! The central difference scheme: u_j^{n+1} = u_j^n - (cΔt/(2Δx)) (u_{j+1}^n - u_{j-1}^n)
//!
//! **Amplification factor**: G(k) = 1 - i (cΔt/Δx) sin(kΔx)
//!
//! **Stability condition**: |G(k)| ≤ 1 for all k
//!
//! The scheme is unconditionally unstable for pure advection. For diffusion-dominated
//! problems, stability requires: cΔt/Δx ≤ 1 (CFL condition for central schemes)
//!
//! ### Fourth-Order Central Difference
//!
//! Higher-order central schemes have better dispersion properties but remain
//! unstable for advection-dominated problems without artificial dissipation.
//!
//! ## Local Truncation Error (LTE) Bounds
//!
//! ### Second-Order Central Difference
//!
//! For the general PDE ∂u/∂t + c ∂u/∂x = 0, the LTE is:
//!
//! τ = - (c Δt / 6) ∂³u/∂x³ (Δx)² + O(Δt Δx² + Δx⁴)
//!
//! **LTE bound**: |τ| ≤ C (c Δt / 6) ||∂³u/∂x³||_∞ Δx²
//!
//! ### Fourth-Order Central Difference
//!
//! LTE = - (c Δt / 30) ∂⁵u/∂x⁵ (Δx)⁴ + O(Δt Δx⁴ + Δx⁶)
//!
//! **LTE bound**: |τ| ≤ C (c Δt / 30) ||∂⁵u/∂x⁵||_∞ Δx⁴
//!
//! ## Stability Regions
//!
//! ### Second-Order Central Difference
//!
//! **Stability region**: |G(ξ)| ≤ 1 where ξ = c Δt / Δx, G(ξ) = 1 - i ξ sin(θ)
//!
//! The scheme is unstable for pure advection but stable for:
//! - Diffusion-dominated problems (Re < 1)
//! - With artificial dissipation
//! - Implicit time integration
//!
//! ### Fourth-Order Central Difference
//!
//! More accurate but narrower stability region than second-order schemes.
//! Requires smaller CFL numbers for stability.
//!
//! ## CFL Conditions
//!
//! ### For Advection: cΔt/Δx ≤ 1/2 (second-order), cΔt/Δx ≤ 1 (fourth-order)
//! ### For Diffusion: Unconditionally stable (for explicit schemes)
//! ### For Navier-Stokes: CFL ≤ min(CFL_advection, CFL_diffusion, CFL_pressure)
//!
//! ## References
//!
//! - Hirsch, C. (1990). *Numerical Computation of Internal and External Flows*.
//!   Wiley. Chapter 6: Stability Theory.
//! - LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*.
//!   Cambridge University Press. Chapter 2: Linear Advection.
//! - Morton, K. W., & Mayers, D. F. (2005). *Numerical solution of partial differential equations*.
//!   Cambridge University Press.

use super::{Grid2D, SpatialDiscretization};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Second-order central difference scheme
pub struct CentralDifference<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for CentralDifference<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> CentralDifference<T> {
    /// Create new central difference scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for CentralDifference<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let divisor = T::from_f64(super::constants::CENTRAL_DIFF_DIVISOR).unwrap_or_else(T::zero);
        (grid.data[(i + 1, j)] - grid.data[(i - 1, j)]) / (divisor * grid.dx)
    }

    fn order(&self) -> usize {
        2
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for second-order central difference: c*dt/dx ≤ 1
    /// The scheme is unstable for pure advection, but stable for diffusion-dominated problems
    fn cfl_limit(&self) -> f64 {
        1.0 // Conservative limit for advection-diffusion problems
    }
}

/// Fourth-order central difference scheme
pub struct FourthOrderCentral<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for FourthOrderCentral<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> FourthOrderCentral<T> {
    /// Create new fourth-order central scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T>
    for FourthOrderCentral<T>
{
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let eight = T::from_f64(8.0).unwrap_or_else(T::zero);
        let twelve = T::from_f64(12.0).unwrap_or_else(T::zero);

        (-grid.data[(i + 2, j)] + eight * grid.data[(i + 1, j)] - eight * grid.data[(i - 1, j)]
            + grid.data[(i - 2, j)])
            / (twelve * grid.dx)
    }

    fn order(&self) -> usize {
        4
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for fourth-order central difference: c*dt/dx ≤ 2
    /// Higher-order schemes allow larger time steps than second-order schemes
    fn cfl_limit(&self) -> f64 {
        2.0 // Fourth-order schemes are more stable than second-order
    }
}
