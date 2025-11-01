//! Total Variation Diminishing (TVD) schemes and MUSCL reconstruction
//!
//! This module provides both the traditional TVD limiter scheme interface
//! for backward compatibility and complete MUSCL (Monotonic Upstream-centered
//! Scheme for Conservation Laws) face reconstruction implementations per
//! Barth & Jespersen (1989).
//!
//! ## References
//!
//! - Barth, T.J. & Jespersen, D.C. (1989). "The design and application of upwind schemes on unstructured meshes"
//! - van Leer, B. (1979). "Towards the ultimate conservative difference scheme. V. A second-order sequel to Godunov's method"
//! - Sweby, P.K. (1984). "High Resolution Schemes Using Flux Limiters"

// Direct access to TVD components
use crate::schemes::grid::Grid2D;
use crate::schemes::FaceReconstruction;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// Flux limiter for TVD schemes
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FluxLimiter {
    /// No limiter (unlimited)
    None,
    /// Van Leer limiter
    VanLeer,
    /// Van Albada limiter
    VanAlbada,
    /// Superbee limiter
    Superbee,
    /// MC (Monotonized Central) limiter
    MC,
    /// Minmod limiter
    Minmod,
}

impl FluxLimiter {
    /// Apply flux limiter function
    pub fn apply<T: RealField + Copy + FromPrimitive + Copy>(&self, r: T) -> T {
        match self {
            FluxLimiter::None => T::one(),
            FluxLimiter::VanLeer => {
                if r > T::zero() {
                    let _two =
                        super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                    (r + r.abs()) / (T::one() + r.abs())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::VanAlbada => {
                if r > T::zero() {
                    (r * r + r) / (r * r + T::one())
                } else {
                    T::zero()
                }
            }
            FluxLimiter::Superbee => {
                let two = super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                // Superbee limiter: φ(r) = max(0, min(1, 2r), min(2, r))
                T::zero().max(T::one().min(two * r).max(two.min(r)))
            }
            FluxLimiter::MC => {
                let two = super::constants::to_realfield::<T>(super::constants::FLUX_LIMITER_TWO);
                T::zero().max(((T::one() + r) / two).min(two.min(two * r)))
            }
            FluxLimiter::Minmod => {
                if r > T::zero() {
                    T::one().min(r)
                } else {
                    T::zero()
                }
            }
        }
    }
}

/// MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
pub struct MUSCLScheme<T: RealField + Copy> {
    /// MUSCL reconstruction order
    order: MUSCLOrder,
    /// TVD flux limiter
    limiter: FluxLimiter,
    _phantom: std::marker::PhantomData<T>,
}

/// MUSCL reconstruction order
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MUSCLOrder {
    /// Second-order MUSCL (MUSCL2)
    SecondOrder,
    /// Third-order MUSCL (MUSCL3/QUICK-like)
    ThirdOrder,
}

impl<T: RealField + Copy + ToPrimitive> MUSCLScheme<T> {
    /// Create new MUSCL scheme with specified order and limiter
    #[must_use]
    pub fn new(order: MUSCLOrder, limiter: FluxLimiter) -> Self {
        Self {
            order,
            limiter,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create MUSCL2 scheme with Superbee limiter
    #[must_use]
    pub fn muscl2_superbee() -> Self {
        Self::new(MUSCLOrder::SecondOrder, FluxLimiter::Superbee)
    }

    /// Create MUSCL2 scheme with van Leer limiter
    #[must_use]
    pub fn muscl2_van_leer() -> Self {
        Self::new(MUSCLOrder::SecondOrder, FluxLimiter::VanLeer)
    }

    /// Create MUSCL2 scheme with Minmod limiter
    #[must_use]
    pub fn muscl2_minmod() -> Self {
        Self::new(MUSCLOrder::SecondOrder, FluxLimiter::Minmod)
    }

    /// Create MUSCL3 scheme with Superbee limiter
    #[must_use]
    pub fn muscl3_superbee() -> Self {
        Self::new(MUSCLOrder::ThirdOrder, FluxLimiter::Superbee)
    }

    /// Compute limited slope using TVD limiter
    fn limited_slope(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T
    where
        T: FromPrimitive + Copy,
    {
        // Compute gradient ratios for limiter
        let delta_i = phi_i - phi_im1;
        let delta_ip1 = phi_ip1 - phi_i;

        // Avoid division by zero - use small epsilon
        let epsilon = T::default_epsilon();

        // Compute r = δ_{i+1} / δ_i (with protection)
        let r = if delta_i.abs() > epsilon {
            delta_ip1 / delta_i
        } else {
            T::zero()
        };

        // Convert r to f64 for limiter function, then back to T
        let r_f64 = r.to_f64().unwrap_or(0.0);
        let psi_f64 = self.limiter.apply(r_f64);
        let psi = T::from_f64(psi_f64).unwrap_or_else(|| T::from_i32(1).unwrap_or(T::one()));

        // Return limited slope: ψ * δ_i
        psi * delta_i
    }

    /// Reconstruct left interface value at cell face for MUSCL2
    fn reconstruct_left_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i + slope * T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) // φ_i + (1/2) * slope
    }

    /// Reconstruct right interface value at cell face for MUSCL2
    fn reconstruct_right_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i - slope * T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) // φ_i - (1/2) * slope
    }

    /// Reconstruct left interface value at cell face for MUSCL3
    fn reconstruct_left_muscl3(&self, _phi_im2: T, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T {
        if let Some(phi_ip2) = phi_ip2 {
            // Use 4-point stencil for 3rd order
            let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
            let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

            // QUICK-like scheme with limiter blending
            let quick = (T::from_f64(6.0).unwrap() * phi_i
                       - T::from_f64(2.0).unwrap() * phi_im1
                       + T::from_f64(8.0).unwrap() * phi_ip1
                       - phi_ip2)
                      / T::from_f64(12.0).unwrap();

            // Blend QUICK with MUSCL2 based on limiter
            let muscl2 = self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if slope1.abs() > T::default_epsilon() {
                slope2 / slope1
            } else {
                T::zero()
            };

            // Convert r to f64 for limiter function, then back to T
            let r_f64 = r.to_f64().unwrap_or(0.0);
            let psi_f64 = self.limiter.apply(r_f64);
            let psi_field = T::from_f64(psi_f64).unwrap_or(T::from_i32(1).unwrap_or(T::one()));
            psi_field * quick + (T::one() - psi_field) * muscl2
        } else {
            // Fall back to MUSCL2 at boundaries where φ_{i+2} is not available
            self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1)
        }
    }

    /// Reconstruct right interface value at cell face for MUSCL3
    fn reconstruct_right_muscl3(&self, _phi_im2: T, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: Option<T>) -> T {
        if let Some(phi_ip2) = phi_ip2 {
            // For right interface, use symmetric QUICK-like formula
            let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
            let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

            // Symmetric QUICK-like for right interface
            let quick = (T::from_f64(6.0).unwrap() * phi_ip1
                       - T::from_f64(2.0).unwrap() * phi_i
                       + T::from_f64(8.0).unwrap() * phi_ip2
                       - phi_ip2) // Approximation - full QUICK would need more points
                      / T::from_f64(12.0).unwrap();

            // Blend with MUSCL2
            let muscl2 = self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if slope1.abs() > T::default_epsilon() {
                slope2 / slope1
            } else {
                T::zero()
            };

            // Convert r to f64 for limiter function, then back to T
            let r_f64 = r.to_f64().unwrap_or(0.0);
            let psi_f64 = self.limiter.apply(r_f64);
            let psi_field = T::from_f64(psi_f64).unwrap_or(T::from_i32(1).unwrap_or(T::one()));
            psi_field * quick + (T::one() - psi_field) * muscl2
        } else {
            // Fall back to MUSCL2 at boundaries
            self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1)
        }
    }
}

/// QUICK (Quadratic Upstream Interpolation for Convective Kinematics)
///
/// Third-order accurate scheme from Leonard (1979) "A stable and accurate convective modelling procedure"
/// based on quadratic upstream interpolation. Combines second-order central differencing with
/// third-order upwind-biased interpolation to minimize numerical diffusion.
pub struct QUICKScheme<T: RealField + Copy> {
    /// Courant number for stability analysis
    courant_max: T,
}

impl<T: RealField + Copy> Default for QUICKScheme<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> QUICKScheme<T> {
    /// Create new QUICK scheme with default Courant limit of 0.8
    #[must_use]
    pub fn new() -> Self {
        Self {
            courant_max: super::constants::to_realfield::<T>(super::constants::CFL_QUICK_SCHEME),
        }
    }

    /// Compute face value using QUICK interpolation
    ///
    /// For uniform grid with flow from left to right:
    /// `φ_face` = 6/8 * `φ_U` + 3/8 * `φ_C` - 1/8 * `φ_UU`
    /// where U = upstream, C = central, UU = far upstream
    pub fn interpolate_face(
        &self,
        phi_uu: T, // Far upstream value
        phi_u: T,  // Upstream value
        phi_c: T,  // Central value
        phi_d: T,  // Downstream value
        velocity: T,
    ) -> T {
        if velocity > T::zero() {
            // Flow from left to right
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_u + three_eighths * phi_c - one_eighth * phi_uu
        } else {
            // Flow from right to left (use downstream values)
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_c + three_eighths * phi_u - one_eighth * phi_d
        }
    }

    /// Check Courant number for stability
    pub fn is_stable(&self, courant: T) -> bool {
        courant <= self.courant_max
    }
}

// FaceReconstruction implementation for MUSCL schemes
impl<T: RealField + Copy + FromPrimitive + ToPrimitive + Copy> FaceReconstruction<T> for MUSCLScheme<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let nx = phi.data.ncols();
        let _ny = phi.data.nrows();

        // Handle boundary cases with one-sided reconstruction
        if i == 0 {
            // Left boundary - use one-sided reconstruction
            match self.order {
                MUSCLOrder::SecondOrder => {
                    if nx > 1usize {
                        let phi_0 = phi.data[(0, j)];
                        let phi_1 = phi.data[(1, j)];

                        // For outflow, use first-order upwind; for inflow, one-sided MUSCL2
                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            // One-sided backward reconstruction: φ_{1/2}^L = φ_0 - (1/2)δφ for inflow
                            phi_0 - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(0, j)]
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    // Fall back to MUSCL2 at boundary (not enough points for 3rd order)
                    if nx > 1usize {
                        let phi_0 = phi.data[(0, j)];
                        let phi_1 = phi.data[(1, j)];

                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            phi_0 - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(0, j)]
                    }
                }
            }
        } else if i >= nx - 1 {
            // Right boundary - use one-sided reconstruction
            match self.order {
                MUSCLOrder::SecondOrder => {
                    if nx > 1usize {
                        let phi_nm1 = phi.data[(nx - 1, j)];
                        let phi_nm2 = phi.data[(nx - 2, j)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            // One-sided forward reconstruction: φ_{N-1/2}^L = φ_{N-1} + (1/2)δφ for outflow
                            phi_nm1 + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(nx - 1, j)]
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    // Fall back to MUSCL2 at boundary
                    if nx > 1usize {
                        let phi_nm1 = phi.data[(nx - 1, j)];
                        let phi_nm2 = phi.data[(nx - 2, j)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            phi_nm1 + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(nx - 1, j)]
                    }
                }
            }
        } else {
            // Interior points - use centered MUSCL reconstruction
            let phi_im1 = phi.data[(i - 1, j)];
            let phi_i = phi.data[(i, j)];
            let phi_ip1 = phi.data[(i + 1, j)];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= T::zero() {
                        // Flow from left to right - reconstruct φ_{i+1/2}^L
                        self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1)
                    } else {
                        // Flow from right to left - reconstruct φ_{i+1/2}^R
                        self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1)
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    // MUSCL3 needs additional points
                    let phi_im2 = if i > 1 { phi.data[(i - 2, j)] } else { phi_im1 };
                    let phi_ip2 = if i + 2 < nx { Some(phi.data[(i + 2, j)]) } else { None };

                    if velocity_at_face >= T::zero() {
                        // Flow from left to right - reconstruct φ_{i+1/2}^L
                        self.reconstruct_left_muscl3(phi_im2, phi_im1, phi_i, phi_ip1, phi_ip2)
                    } else {
                        // Flow from right to left - reconstruct φ_{i+1/2}^R
                        self.reconstruct_right_muscl3(phi_im2, phi_im1, phi_i, phi_ip1, phi_ip2)
                    }
                }
            }
        }
    }

    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let _nx = phi.data.ncols();
        let ny = phi.data.nrows();

        // Handle boundary cases with one-sided reconstruction (analogous to x-direction)
        if j == 0 {
            // Bottom boundary
            match self.order {
                MUSCLOrder::SecondOrder => {
                    if ny > 1usize {
                        let phi_0 = phi.data[(i, 0)];
                        let phi_1 = phi.data[(i, 1)];

                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            phi_0 - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(i, 0)]
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_0 = phi.data[(i, 0)];
                        let phi_1 = phi.data[(i, 1)];

                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            phi_0 - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(i, 0)]
                    }
                }
            }
        } else if j >= ny - 1 {
            // Top boundary
            match self.order {
                MUSCLOrder::SecondOrder => {
                    if ny > 1usize {
                        let phi_nm1 = phi.data[(i, ny - 1)];
                        let phi_nm2 = phi.data[(i, ny - 2)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            phi_nm1 + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(i, ny - 1)]
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_nm1 = phi.data[(i, ny - 1)];
                        let phi_nm2 = phi.data[(i, ny - 2)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            phi_nm1 + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one())) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(i, ny - 1)]
                    }
                }
            }
        } else {
            // Interior points - use centered MUSCL reconstruction
            let phi_jm1 = phi.data[(i, j - 1)];
            let phi_j = phi.data[(i, j)];
            let phi_jp1 = phi.data[(i, j + 1)];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= T::zero() {
                        self.reconstruct_left_muscl2(phi_jm1, phi_j, phi_jp1)
                    } else {
                        self.reconstruct_right_muscl2(phi_jm1, phi_j, phi_jp1)
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    let phi_jm2 = if j > 1 { phi.data[(i, j - 2)] } else { phi_jm1 };
                    let phi_jp2 = if j + 2 < ny { Some(phi.data[(i, j + 2)]) } else { None };

                    if velocity_at_face >= T::zero() {
                        self.reconstruct_left_muscl3(phi_jm2, phi_jm1, phi_j, phi_jp1, phi_jp2)
                    } else {
                        self.reconstruct_right_muscl3(phi_jm2, phi_jm1, phi_j, phi_jp1, phi_jp2)
                    }
                }
            }
        }
    }

    fn order(&self) -> usize {
        match self.order {
            MUSCLOrder::SecondOrder => 2,
            MUSCLOrder::ThirdOrder => 3,
        }
    }
}
