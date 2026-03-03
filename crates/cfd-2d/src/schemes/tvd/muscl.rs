//! MUSCL reconstruction scheme with TVD flux limiters.

use super::FluxLimiter;
use crate::schemes::grid::Grid2D;
use crate::schemes::FaceReconstruction;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

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
        let delta_i = phi_i - phi_im1;
        let delta_ip1 = phi_ip1 - phi_i;

        let epsilon = T::default_epsilon();

        let r = if delta_i.abs() > epsilon {
            delta_ip1 / delta_i
        } else {
            T::zero()
        };

        let r_f64 = r.to_f64().unwrap_or(0.0);
        let psi_f64 = self.limiter.apply(r_f64);
        let psi = T::from_f64(psi_f64).unwrap_or_else(|| T::from_i32(1).unwrap_or(T::one()));

        psi * delta_i
    }

    /// Reconstruct left interface value at cell face for MUSCL2
    fn reconstruct_left_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i + slope * T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
    }

    /// Reconstruct right interface value at cell face for MUSCL2
    fn reconstruct_right_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i - slope * T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
    }

    /// Reconstruct left interface value at cell face for MUSCL3
    fn reconstruct_left_muscl3(
        &self,
        _phi_im2: T,
        phi_im1: T,
        phi_i: T,
        phi_ip1: T,
        phi_ip2: Option<T>,
    ) -> T {
        if let Some(phi_ip2) = phi_ip2 {
            let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
            let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

            let quick = (T::from_f64(6.0).unwrap_or_else(num_traits::Zero::zero) * phi_i - T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero) * phi_im1
                + T::from_f64(8.0).unwrap_or_else(num_traits::Zero::zero) * phi_ip1
                - phi_ip2)
                / T::from_f64(12.0).unwrap_or_else(num_traits::Zero::zero);

            let muscl2 = self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if slope1.abs() > T::default_epsilon() {
                slope2 / slope1
            } else {
                T::zero()
            };

            let r_f64 = r.to_f64().unwrap_or(0.0);
            let psi_f64 = self.limiter.apply(r_f64);
            let psi_field = T::from_f64(psi_f64).unwrap_or(T::from_i32(1).unwrap_or(T::one()));
            psi_field * quick + (T::one() - psi_field) * muscl2
        } else {
            self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1)
        }
    }

    /// Reconstruct right interface value at cell face for MUSCL3
    fn reconstruct_right_muscl3(
        &self,
        _phi_im2: T,
        phi_im1: T,
        phi_i: T,
        phi_ip1: T,
        phi_ip2: Option<T>,
    ) -> T {
        if let Some(phi_ip2) = phi_ip2 {
            let slope1 = self.limited_slope(phi_im1, phi_i, phi_ip1);
            let slope2 = self.limited_slope(phi_i, phi_ip1, phi_ip2);

            let quick = (-phi_i
                + T::from_f64(5.0).unwrap_or_else(num_traits::Zero::zero) * phi_ip1
                + T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero) * phi_ip2)
                / T::from_f64(6.0).unwrap_or_else(num_traits::Zero::zero);

            let muscl2 = self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if slope1.abs() > T::default_epsilon() {
                slope2 / slope1
            } else {
                T::zero()
            };

            let r_f64 = r.to_f64().unwrap_or(0.0);
            let psi_f64 = self.limiter.apply(r_f64);
            let psi_field = T::from_f64(psi_f64).unwrap_or(T::from_i32(1).unwrap_or(T::one()));
            psi_field * quick + (T::one() - psi_field) * muscl2
        } else {
            self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1)
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + Copy> FaceReconstruction<T>
    for MUSCLScheme<T>
{
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let nx = phi.data.nrows();
        let _ny = phi.data.ncols();

        if i == 0 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if nx > 1usize {
                        let phi_0 = phi.data[(0, j)];
                        let phi_1 = phi.data[(1, j)];

                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            phi_0
                                - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
                                    * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(0, j)]
                    }
                }
            }
        } else if i >= nx - 1 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if nx > 1usize {
                        let phi_nm1 = phi.data[(nx - 1, j)];
                        let phi_nm2 = phi.data[(nx - 2, j)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            phi_nm1
                                + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
                                    * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(nx - 1, j)]
                    }
                }
            }
        } else {
            let phi_im1 = phi.data[(i - 1, j)];
            let phi_i = phi.data[(i, j)];
            let phi_ip1 = phi.data[(i + 1, j)];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= T::zero() {
                        self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1)
                    } else if i + 2 < nx {
                        let phi_ip2 = phi.data[(i + 2, j)];
                        self.reconstruct_right_muscl2(phi_i, phi_ip1, phi_ip2)
                    } else {
                        phi_ip1
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    let phi_im2 = if i > 1 { phi.data[(i - 2, j)] } else { phi_im1 };
                    let phi_ip2 = if i + 2 < nx {
                        Some(phi.data[(i + 2, j)])
                    } else {
                        None
                    };

                    if velocity_at_face >= T::zero() {
                        self.reconstruct_left_muscl3(phi_im2, phi_im1, phi_i, phi_ip1, phi_ip2)
                    } else {
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
        let _nx = phi.data.nrows();
        let ny = phi.data.ncols();

        if j == 0 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_0 = phi.data[(i, 0)];
                        let phi_1 = phi.data[(i, 1)];

                        if velocity_at_face >= T::zero() {
                            phi_0
                        } else {
                            phi_0
                                - T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
                                    * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[(i, 0)]
                    }
                }
            }
        } else if j >= ny - 1 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_nm1 = phi.data[(i, ny - 1)];
                        let phi_nm2 = phi.data[(i, ny - 2)];

                        if velocity_at_face <= T::zero() {
                            phi_nm1
                        } else {
                            phi_nm1
                                + T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()))
                                    * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[(i, ny - 1)]
                    }
                }
            }
        } else {
            let phi_jm1 = phi.data[(i, j - 1)];
            let phi_j = phi.data[(i, j)];
            let phi_jp1 = phi.data[(i, j + 1)];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= T::zero() {
                        self.reconstruct_left_muscl2(phi_jm1, phi_j, phi_jp1)
                    } else if j + 2 < ny {
                        let phi_jp2 = phi.data[(i, j + 2)];
                        self.reconstruct_right_muscl2(phi_j, phi_jp1, phi_jp2)
                    } else {
                        phi_jp1
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    let phi_jm2 = if j > 1 { phi.data[(i, j - 2)] } else { phi_jm1 };
                    let phi_jp2 = if j + 2 < ny {
                        Some(phi.data[(i, j + 2)])
                    } else {
                        None
                    };

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
