//! MUSCL reconstruction scheme with TVD flux limiters.

use super::FluxLimiter;
use crate::scalar::Cfd2dScalar;
use crate::scalar::{from_f64, one, zero};
use crate::schemes::grid::Grid2D;
use crate::schemes::FaceReconstruction;
use eunomia::{FloatElement, NumericElement};

/// MUSCL (Monotonic Upstream-centered Scheme for Conservation Laws)
pub struct MUSCLScheme<T: Cfd2dScalar + Copy> {
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

impl<T: Cfd2dScalar + Copy + FloatElement> MUSCLScheme<T> {
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
    fn limited_slope(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let delta_i = phi_i - phi_im1;
        let delta_ip1 = phi_ip1 - phi_i;

        let epsilon = T::default_epsilon();

        let r = if <T as NumericElement>::abs(delta_i) > epsilon {
            delta_ip1 / delta_i
        } else {
            zero()
        };

        let psi = self.limiter.apply(r);

        psi * delta_i
    }

    /// Reconstruct left interface value at cell face for MUSCL2
    fn reconstruct_left_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i + slope * from_f64::<T>(0.5)
    }

    /// Reconstruct right interface value at cell face for MUSCL2
    fn reconstruct_right_muscl2(&self, phi_im1: T, phi_i: T, phi_ip1: T) -> T {
        let slope = self.limited_slope(phi_im1, phi_i, phi_ip1);
        phi_i - slope * from_f64::<T>(0.5)
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

            let quick = (from_f64::<T>(6.0) * phi_i - from_f64::<T>(2.0) * phi_im1
                + from_f64::<T>(8.0) * phi_ip1
                - phi_ip2)
                / from_f64::<T>(12.0);

            let muscl2 = self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if <T as NumericElement>::abs(slope1) > T::default_epsilon() {
                slope2 / slope1
            } else {
                zero()
            };

            let psi = self.limiter.apply(r);
            psi * quick + (one::<T>() - psi) * muscl2
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

            let quick = (-phi_i + from_f64::<T>(5.0) * phi_ip1 + from_f64::<T>(2.0) * phi_ip2)
                / from_f64::<T>(6.0);

            let muscl2 = self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1);
            let r = if <T as NumericElement>::abs(slope1) > T::default_epsilon() {
                slope2 / slope1
            } else {
                zero()
            };

            let psi = self.limiter.apply(r);
            psi * quick + (one::<T>() - psi) * muscl2
        } else {
            self.reconstruct_right_muscl2(phi_im1, phi_i, phi_ip1)
        }
    }
}

impl<T: Cfd2dScalar + Copy + FloatElement> FaceReconstruction<T> for MUSCLScheme<T> {
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,
        velocity_at_face: T,
        i: usize,
        j: usize,
    ) -> T {
        let [nx, _ny] = phi.data.shape();

        if i == 0 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if nx > 1usize {
                        let phi_0 = phi.data[[0, j]];
                        let phi_1 = phi.data[[1, j]];

                        if velocity_at_face >= zero() {
                            phi_0
                        } else {
                            phi_0 - from_f64::<T>(0.5) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[[0, j]]
                    }
                }
            }
        } else if i >= nx - 1 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if nx > 1usize {
                        let phi_nm1 = phi.data[[nx - 1, j]];
                        let phi_nm2 = phi.data[[nx - 2, j]];

                        if velocity_at_face <= zero() {
                            phi_nm1
                        } else {
                            phi_nm1 + from_f64::<T>(0.5) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[[nx - 1, j]]
                    }
                }
            }
        } else {
            let phi_im1 = phi.data[[i - 1, j]];
            let phi_i = phi.data[[i, j]];
            let phi_ip1 = phi.data[[i + 1, j]];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= zero() {
                        self.reconstruct_left_muscl2(phi_im1, phi_i, phi_ip1)
                    } else if i + 2 < nx {
                        let phi_ip2 = phi.data[[i + 2, j]];
                        self.reconstruct_right_muscl2(phi_i, phi_ip1, phi_ip2)
                    } else {
                        phi_ip1
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    let phi_im2 = if i > 1 { phi.data[[i - 2, j]] } else { phi_im1 };
                    let phi_ip2 = if i + 2 < nx {
                        Some(phi.data[[i + 2, j]])
                    } else {
                        None
                    };

                    if velocity_at_face >= zero() {
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
        let [_nx, ny] = phi.data.shape();

        if j == 0 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_0 = phi.data[[i, 0]];
                        let phi_1 = phi.data[[i, 1]];

                        if velocity_at_face >= zero() {
                            phi_0
                        } else {
                            phi_0 - from_f64::<T>(0.5) * (phi_1 - phi_0)
                        }
                    } else {
                        phi.data[[i, 0]]
                    }
                }
            }
        } else if j >= ny - 1 {
            match self.order {
                MUSCLOrder::SecondOrder | MUSCLOrder::ThirdOrder => {
                    if ny > 1usize {
                        let phi_nm1 = phi.data[[i, ny - 1]];
                        let phi_nm2 = phi.data[[i, ny - 2]];

                        if velocity_at_face <= zero() {
                            phi_nm1
                        } else {
                            phi_nm1 + from_f64::<T>(0.5) * (phi_nm1 - phi_nm2)
                        }
                    } else {
                        phi.data[[i, ny - 1]]
                    }
                }
            }
        } else {
            let phi_jm1 = phi.data[[i, j - 1]];
            let phi_j = phi.data[[i, j]];
            let phi_jp1 = phi.data[[i, j + 1]];

            match self.order {
                MUSCLOrder::SecondOrder => {
                    if velocity_at_face >= zero() {
                        self.reconstruct_left_muscl2(phi_jm1, phi_j, phi_jp1)
                    } else if j + 2 < ny {
                        let phi_jp2 = phi.data[[i, j + 2]];
                        self.reconstruct_right_muscl2(phi_j, phi_jp1, phi_jp2)
                    } else {
                        phi_jp1
                    }
                }
                MUSCLOrder::ThirdOrder => {
                    let phi_jm2 = if j > 1 { phi.data[[i, j - 2]] } else { phi_jm1 };
                    let phi_jp2 = if j + 2 < ny {
                        Some(phi.data[[i, j + 2]])
                    } else {
                        None
                    };

                    if velocity_at_face >= zero() {
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
