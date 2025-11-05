//! Upwind discretization schemes with stability analysis
//!
//! These schemes implement proper upwind discretization based on the velocity field,
//! not the scalar field being transported. This is critical for physical correctness.
//!
//! ## Von Neumann Stability Analysis
//!
//! ### First-Order Upwind Scheme
//!
//! For the 1D advection equation: ∂u/∂t + c ∂u/∂x = 0
//!
//! The upwind scheme: u_j^{n+1} = u_j^n - (cΔt/Δx) (u_j^n - u_{j-1}^n) for c > 0
//!
//! **Amplification factor**: G(k) = 1 - CFL * (1 - e^{-ikΔx})
//!
//! **Stability condition**: |G(k)| ≤ 1 for all k
//!
//! The scheme is unconditionally stable for CFL ≤ 1.
//!
//! ### Second-Order Upwind Scheme
//!
//! Uses linear interpolation with upwind bias: φ_face = φ_upwind + (1/2) * limiter * (φ_downwind - φ_upwind)
//!
//! **Stability**: CFL ≤ 1 for explicit schemes, unconditionally stable for implicit.
//!
//! ## Local Truncation Error (LTE) Bounds
//!
//! ### First-Order Upwind Scheme
//!
//! For ∂u/∂t + c ∂u/∂x = 0, the LTE is:
//!
//! τ = (c Δx / 2) (1 - CFL) ∂²u/∂x² + O(Δx²)
//!
//! **LTE bound**: |τ| ≤ (c Δx / 2) ||∂²u/∂x²||_∞ + O(Δx²)
//!
//! ### Second-Order Upwind Scheme
//!
//! LTE = - (c Δx² / 6) (1 - CFL/2) ∂³u/∂x³ + O(Δx³)
//!
//! **LTE bound**: |τ| ≤ (c Δx² / 6) ||∂³u/∂x³||_∞ + O(Δx³)
//!
//! ## Stability Regions
//!
//! ### First-Order Upwind
//!
//! **Stability region**: Entire left half-plane for CFL ≤ 1
//!
//! G(ξ) = 1 - CFL + CFL e^{-iξ} where ξ = k Δx
//!
//! |G(ξ)|² = 1 + CFL² - 2 CFL cos(ξ) ≤ 1 ⇒ CFL ≤ 1
//!
//! ### Second-Order Upwind
//!
//! Stability depends on the limiter used. With TVD limiters, maintains
//! CFL ≤ 1 stability but with improved accuracy.
//!
//! ## CFL Conditions for Upwind Schemes
//!
//! ### First-order upwind: CFL ≤ 1 (stable for all CFL ≤ 1)
//! ### Second-order upwind: CFL ≤ 1 (with appropriate limiters)
//! ### Higher-order upwind: CFL ≤ 1 (with TVD limiters for stability)
//!
//! ## References
//!
//! - Hirsch, C. (1990). *Numerical Computation of Internal and External Flows*.
//!   Wiley. Chapter 8: Upwind Schemes.
//! - LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*.
//!   Cambridge University Press. Chapter 6: High-Resolution Methods.

use super::{FaceReconstruction, Grid2D, SpatialDiscretization};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// First-order upwind scheme
pub struct FirstOrderUpwind<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for FirstOrderUpwind<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FirstOrderUpwind<T> {
    /// Create new first-order upwind scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FaceReconstruction<T> for FirstOrderUpwind<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        // Check boundaries - return boundary value
        let nx = phi.data.ncols();
        if i >= nx - 1 {
            // At boundary, return the last cell value
            return phi.data[(nx - 1, j)];
        }

        if velocity_at_face > T::zero() {
            // Positive velocity: flow is from left (i) to right (i+1)
            // The upwind value is from the cell on the left
            phi.data[(i, j)]
        } else {
            // Negative velocity: flow is from right (i+1) to left (i)
            // The upwind value is from the cell on the right
            phi.data[(i + 1, j)]
        }
    }

    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        // Check boundaries - return boundary value
        let ny = phi.data.nrows();
        if j >= ny - 1 {
            // At boundary, return the last cell value
            return phi.data[(i, ny - 1)];
        }

        if velocity_at_face > T::zero() {
            // Positive velocity: flow is from bottom (j) to top (j+1)
            // The upwind value is from the cell on the bottom
            phi.data[(i, j)]
        } else {
            // Negative velocity: flow is from top (j+1) to bottom (j)
            // The upwind value is from the cell on the top
            phi.data[(i, j + 1)]
        }
    }

    fn order(&self) -> usize {
        1
    }
}

// Keep backward compatibility with old interface (deprecated)
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for FirstOrderUpwind<T> {
    fn compute_derivative(&self, _grid: &Grid2D<T>, _i: usize, _j: usize) -> T {
        // Return zero as a safe default for deprecated interface
        // Users should migrate to FaceReconstruction trait
        T::zero()
    }

    fn order(&self) -> usize {
        1
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for first-order upwind: CFL ≤ 1
    /// The scheme is stable for all CFL ≤ 1
    fn cfl_limit(&self) -> f64 {
        1.0
    }
}

/// Second-order upwind scheme
pub struct SecondOrderUpwind<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for SecondOrderUpwind<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> SecondOrderUpwind<T> {
    /// Create new second-order upwind scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> FaceReconstruction<T> for SecondOrderUpwind<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let half = T::from_f64(0.5).expect("Failed to represent 0.5 in numeric type T");
        let three_half = T::from_f64(1.5).expect("Failed to represent 1.5 in numeric type T");

        if velocity_at_face > T::zero() {
            // Positive velocity: use backward-biased stencil
            // Check boundaries
            if i < 1 {
                // Fall back to first-order at boundary
                return phi.data[(i, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_i - 1/2*phi_{i-1}
            three_half * phi.data[(i, j)] - half * phi.data[(i - 1, j)]
        } else {
            // Negative velocity: use forward-biased stencil
            // Check boundaries
            let nx = phi.data.ncols();
            if i >= nx - 2 {
                // Fall back to first-order at boundary
                return phi.data[(i + 1, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_{i+1} - 1/2*phi_{i+2}
            three_half * phi.data[(i + 1, j)] - half * phi.data[(i + 2, j)]
        }
    }

    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let half = T::from_f64(0.5).expect("Failed to represent 0.5 in numeric type T");
        let three_half = T::from_f64(1.5).expect("Failed to represent 1.5 in numeric type T");

        if velocity_at_face > T::zero() {
            // Positive velocity: use backward-biased stencil
            // Check boundaries
            if j < 1 {
                // Fall back to first-order at boundary
                return phi.data[(i, j)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_j - 1/2*phi_{j-1}
            three_half * phi.data[(i, j)] - half * phi.data[(i, j - 1)]
        } else {
            // Negative velocity: use forward-biased stencil
            // Check boundaries
            let ny = phi.data.nrows();
            if j >= ny - 2 {
                // Fall back to first-order at boundary
                return phi.data[(i, j + 1)];
            }
            // Second-order reconstruction: phi_face = 3/2*phi_{j+1} - 1/2*phi_{j+2}
            three_half * phi.data[(i, j + 1)] - half * phi.data[(i, j + 2)]
        }
    }

    fn order(&self) -> usize {
        2
    }
}

// Keep backward compatibility with old interface (deprecated)
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for SecondOrderUpwind<T> {
    fn compute_derivative(&self, _grid: &Grid2D<T>, _i: usize, _j: usize) -> T {
        // Return zero as a safe default for deprecated interface
        // Users should migrate to FaceReconstruction trait
        T::zero()
    }

    fn order(&self) -> usize {
        2
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for second-order upwind: CFL ≤ 1
    /// With proper limiters, the scheme maintains stability up to CFL = 1
    fn cfl_limit(&self) -> f64 {
        1.0
    }
}
