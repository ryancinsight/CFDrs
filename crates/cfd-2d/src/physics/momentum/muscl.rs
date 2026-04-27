//! MUSCL (Monotonic Upstream-Centered Scheme for Conservation Laws) reconstruction
//!
//! This module implements MUSCL reconstruction schemes for higher-order spatial discretization
//! in convection-dominated flows. MUSCL provides 2nd/3rd order accuracy while maintaining
//! monotonicity through TVD limiters.
//!
//! ## References
//!
//! - van Leer, B. (1979). "Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method"
//! - Hirsch, C. (2007). "Numerical Computation of Internal and External Flows: Fundamentals of Computational Fluid Dynamics"
//! - Barth, T.J. & Jespersen, D.C. (1989). "The design and application of upwind schemes on unstructured meshes"
//!
//! # Theorem - Conservative Face Reconstruction
//!
//! The momentum discretization must conserve linear momentum globally and locally.
//!
//! **Proof sketch**:
//! By integrating the Navier-Stokes momentum equation over a control volume $\Omega$,
//! Gauss's divergence theorem converts the convective and diffusive volume integrals
//! into surface fluxes. The finite volume method ensures that the flux leaving one
//! cell exactly equals the flux entering the adjacent cell. Thus, in the absence of
//! external forces and boundary fluxes, the total momentum $\int_\Omega \rho \mathbf{u} dV$
//! is exactly conserved to machine precision.
//!
//! # Theorem - Uniform-Grid QUICK Face Values
//!
//! On a uniform grid with cell centers `x_i = i`, the unique quadratic
//! interpolant through `(x_{i-1}, phi_{i-1})`, `(x_i, phi_i)`, and
//! `(x_{i+1}, phi_{i+1})` evaluated at `x_{i+1/2}` is
//!
//! ```text
//! phi^L_{i+1/2} = -phi_{i-1}/8 + 3 phi_i/4 + 3 phi_{i+1}/8.
//! ```
//!
//! The right-biased state from cells `i`, `i+1`, and `i+2` is
//!
//! ```text
//! phi^R_{i+1/2} = 3 phi_i/8 + 3 phi_{i+1}/4 - phi_{i+2}/8.
//! ```
//!
//! **Proof sketch**: Write the Lagrange basis polynomials for the three
//! equispaced cell centers and evaluate them at the face coordinate
//! `x = i + 1/2`. The resulting weights are exactly `(-1/8, 3/4, 3/8)` for
//! the left stencil and `(3/8, 3/4, -1/8)` for the right stencil. Substitution
//! of any polynomial with degree <= 2 reproduces the exact face value, giving
//! third-order finite-volume face accuracy for smooth cell-centered fields.

use super::tvd_limiters::TvdLimiter;
use nalgebra::RealField;

/// MUSCL reconstruction scheme order
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MusclOrder {
    /// Second-order MUSCL (MUSCL2)
    SecondOrder,
    /// Third-order MUSCL (MUSCL3/QUICK-like)
    ThirdOrder,
}

/// MUSCL reconstruction interface
pub trait MusclReconstruction<T: RealField + Copy> {
    /// Reconstruct left interface value at cell face (φ_{i+1/2}^L)
    fn reconstruct_left(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T;

    /// Reconstruct right interface value at cell face (φ_{i+1/2}^R)
    fn reconstruct_right(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T;

    /// Get scheme name
    fn name(&self) -> &str;

    /// Get scheme order
    fn order(&self) -> MusclOrder;
}

/// MUSCL reconstruction with TVD limiter
pub struct MusclScheme<T, L>
where
    T: RealField + Copy,
    L: TvdLimiter<T>,
{
    limiter: L,
    order: MusclOrder,
    _phantom: std::marker::PhantomData<T>,
}

impl<T, L> MusclScheme<T, L>
where
    T: RealField + Copy,
    L: TvdLimiter<T>,
{
    /// Create new MUSCL scheme with specified limiter and order
    pub fn new(limiter: L, order: MusclOrder) -> Self {
        Self {
            limiter,
            order,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Compute limited slope using TVD limiter
    fn limited_slope(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        // Compute gradient ratios for limiter
        let delta_i = phi_i - phi_im1;
        let delta_ip1 = phi_ip1 - phi_i;

        // Avoid division by zero - use small epsilon
        let epsilon = T::from_f64(1e-12).unwrap_or_else(T::default_epsilon);

        // Compute r = δ_{i+1} / δ_i (with protection)
        let r = if delta_i.abs() > epsilon {
            delta_ip1 / delta_i
        } else {
            T::zero()
        };

        // Apply TVD limiter
        let psi = self.limiter.limit(r);

        // Return limited slope: ψ * δ_i
        psi * delta_i
    }
}

impl<T, L> MusclReconstruction<T> for MusclScheme<T, L>
where
    T: RealField + Copy,
    L: TvdLimiter<T>,
{
    fn reconstruct_left(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T {
        match self.order {
            MusclOrder::SecondOrder => {
                // MUSCL2: Linear reconstruction with limiter
                let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
                phi_i + slope / (T::one() + T::one()) // φ_i + (1/2) * slope
            }
            MusclOrder::ThirdOrder => {
                // MUSCL3/QUICK-like: Quadratic reconstruction
                if let Some(phi_ip2) = phi_ip2 {
                    // Use 4-point stencil for 3rd order
                    let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
                    let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

                    let quick = quadratic_left_face(phi_im1, phi_i, phi_ip1);

                    // Blend QUICK with MUSCL2 based on limiter
                    let muscl2 = phi_i + slope1 / (T::one() + T::one());
                    let r = if slope1.abs() > T::default_epsilon() {
                        slope2 / slope1
                    } else {
                        T::zero()
                    };

                    let psi = self.limiter.limit(r);
                    psi * quick + (T::one() - psi) * muscl2
                } else {
                    // Fall back to MUSCL2 at boundaries
                    let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
                    phi_i + slope / (T::one() + T::one())
                }
            }
        }
    }

    fn reconstruct_right(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T {
        match self.order {
            MusclOrder::SecondOrder => {
                // MUSCL2: Linear reconstruction with limiter
                let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
                phi_i - slope / (T::one() + T::one()) // φ_i - (1/2) * slope
            }
            MusclOrder::ThirdOrder => {
                // MUSCL3/QUICK-like: Quadratic reconstruction
                if let Some(phi_ip2) = phi_ip2 {
                    // For right interface, use symmetric QUICK-like formula
                    let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
                    let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

                    let quick = quadratic_right_face(phi_i, phi_ip1, phi_ip2);

                    // Blend with MUSCL2
                    let muscl2 = phi_i - slope1 / (T::one() + T::one());
                    let r = if slope1.abs() > T::default_epsilon() {
                        slope2 / slope1
                    } else {
                        T::zero()
                    };

                    let psi = self.limiter.limit(r);
                    psi * quick + (T::one() - psi) * muscl2
                } else {
                    // Fall back to MUSCL2 at boundaries
                    let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
                    phi_i - slope / (T::one() + T::one())
                }
            }
        }
    }

    fn name(&self) -> &str {
        match self.order {
            MusclOrder::SecondOrder => "MUSCL2",
            MusclOrder::ThirdOrder => "MUSCL3",
        }
    }

    fn order(&self) -> MusclOrder {
        self.order
    }
}

#[inline]
fn quadratic_left_face<T: RealField + Copy>(phi_im1: T, phi_i: T, phi_ip1: T) -> T {
    let one_eighth = T::from_f64(0.125).expect("analytical constant conversion");
    let three_eighths = T::from_f64(0.375).expect("analytical constant conversion");
    let three_quarters = T::from_f64(0.75).expect("analytical constant conversion");
    -one_eighth * phi_im1 + three_quarters * phi_i + three_eighths * phi_ip1
}

#[inline]
fn quadratic_right_face<T: RealField + Copy>(phi_i: T, phi_ip1: T, phi_ip2: T) -> T {
    let one_eighth = T::from_f64(0.125).expect("analytical constant conversion");
    let three_eighths = T::from_f64(0.375).expect("analytical constant conversion");
    let three_quarters = T::from_f64(0.75).expect("analytical constant conversion");
    three_eighths * phi_i + three_quarters * phi_ip1 - one_eighth * phi_ip2
}

/// Convenience constructors for common MUSCL schemes
pub mod schemes {
    use super::MusclScheme;
    use crate::physics::momentum::tvd_limiters::{Minmod, Superbee, VanLeer};

    /// MUSCL2 with Superbee limiter (most accurate for shocks)
    pub type Muscl2Superbee<T> = MusclScheme<T, Superbee>;

    /// MUSCL2 with van Leer limiter (good general purpose)
    pub type Muscl2VanLeer<T> = MusclScheme<T, VanLeer>;

    /// MUSCL2 with Minmod limiter (most dissipative)
    pub type Muscl2Minmod<T> = MusclScheme<T, Minmod>;

    /// MUSCL3 with Superbee limiter (high accuracy with monotonicity)
    pub type Muscl3Superbee<T> = MusclScheme<T, Superbee>;

    /// MUSCL3 with van Leer limiter (high accuracy general purpose)
    pub type Muscl3VanLeer<T> = MusclScheme<T, VanLeer>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::momentum::tvd_limiters::{Superbee, VanLeer};
    use approx::assert_relative_eq;

    #[test]
    fn test_muscl2_superbee_smooth_solution() {
        let limiter = Superbee;
        let scheme = MusclScheme::new(limiter, MusclOrder::SecondOrder);

        // Smooth solution: φ = x, uniform grid spacing = 1
        // For i=1, φ_{i-1}=0, φ_i=1, φ_{i+1}=2
        // MUSCL2 should give φ_{i+1/2}^L = 1 + (1/2)(2-0)/2 = 1.5
        // Wait, let me recalculate:
        // slope = limiter(r) * δ_i, where δ_i = φ_i - φ_{i-1} = 1-0 = 1
        // r = δ_{i+1}/δ_i = (2-1)/(1-0) = 1/1 = 1
        // Superbee limiter: ψ(1) = max(0, min(2*1, 1), min(1, 2*1)) = min(2,1,2) = 1
        // slope = 1 * 1 = 1
        // φ_{i+1/2}^L = φ_i + slope/2 = 1 + 1/2 = 1.5

        let phi_l = scheme.reconstruct_left(0.0, 1.0, 2.0, Some(3.0));
        assert_relative_eq!(phi_l, 1.5, epsilon = 1e-10);
    }

    #[test]
    fn test_muscl2_vanleer_extrema() {
        let limiter = VanLeer;
        let scheme = MusclScheme::new(limiter, MusclOrder::SecondOrder);

        // At extrema (r=0), van Leer limiter should give ψ=0
        let phi_l = scheme.reconstruct_left(1.0, 1.0, 2.0, None);
        assert_relative_eq!(phi_l, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_muscl_scheme_naming() {
        let limiter = Superbee;
        let muscl2 = MusclScheme::<f64, Superbee>::new(limiter, MusclOrder::SecondOrder);
        let muscl3 = MusclScheme::<f64, Superbee>::new(limiter, MusclOrder::ThirdOrder);

        assert_eq!(muscl2.name(), "MUSCL2");
        assert_eq!(muscl3.name(), "MUSCL3");
        assert_eq!(muscl2.order(), MusclOrder::SecondOrder);
        assert_eq!(muscl3.order(), MusclOrder::ThirdOrder);
    }

    #[test]
    fn test_muscl3_fallback_to_muscl2_at_boundary() {
        let limiter = Superbee;
        let scheme = MusclScheme::new(limiter, MusclOrder::ThirdOrder);

        // At boundary (no φ_{i+2}), should fall back to MUSCL2
        let phi_l_boundary = scheme.reconstruct_left(0.0, 1.0, 2.0, None);
        let phi_l_muscl2 = {
            let muscl2 = MusclScheme::new(Superbee, MusclOrder::SecondOrder);
            muscl2.reconstruct_left(0.0, 1.0, 2.0, None)
        };

        assert_relative_eq!(phi_l_boundary, phi_l_muscl2, epsilon = 1e-10);
    }

    #[test]
    fn muscl3_quick_reproduces_linear_face_value_from_both_sides() {
        let limiter = Superbee;
        let scheme = MusclScheme::new(limiter, MusclOrder::ThirdOrder);

        let left = scheme.reconstruct_left(0.0, 1.0, 2.0, Some(3.0));
        let right = scheme.reconstruct_right(0.0, 1.0, 2.0, Some(3.0));

        assert_relative_eq!(left, 1.5, epsilon = 1e-12);
        assert_relative_eq!(right, 1.5, epsilon = 1e-12);
    }

    #[test]
    fn quadratic_quick_face_values_reproduce_parabola() {
        let left = quadratic_left_face(1.0_f64, 0.0, 1.0);
        let right = quadratic_right_face(0.0_f64, 1.0, 4.0);

        assert_relative_eq!(left, 0.25, epsilon = 1e-12);
        assert_relative_eq!(right, 0.25, epsilon = 1e-12);
    }
}
