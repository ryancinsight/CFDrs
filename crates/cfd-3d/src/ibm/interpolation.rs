//! Interpolation kernels and delta functions for IBM

use nalgebra::RealField;
use num_traits::FromPrimitive;

// Numerical constants
const TWO: f64 = 2.0;
const THREE: f64 = 3.0;
const FIVE: f64 = 5.0;
const INTERPOLATION_STENCIL_SIZE: usize = 4;

/// Delta function types for IBM
#[derive(Debug, Clone, Copy)]
pub enum DeltaFunction {
    /// Roma & Peskin (2000) 3-point delta function
    RomaPeskin3,
    /// Roma & Peskin (2000) 4-point delta function
    RomaPeskin4,
    /// Peskin (2002) 4-point delta function
    Peskin4,
}

/// Interpolation kernel for transferring between Eulerian and Lagrangian grids
pub struct InterpolationKernel<T: RealField + Copy> {
    delta_type: DeltaFunction,
    #[allow(dead_code)]
    width: T,
}

impl<T: RealField + FromPrimitive + Copy> InterpolationKernel<T> {
    /// Create a new interpolation kernel
    pub fn new(delta_type: DeltaFunction, width: T) -> Self {
        Self { delta_type, width }
    }

    /// Evaluate the delta function
    pub fn delta(&self, r: T) -> T {
        let r_abs = r.abs();

        match self.delta_type {
            DeltaFunction::RomaPeskin3 => self.roma_peskin_3(r_abs),
            DeltaFunction::RomaPeskin4 => self.roma_peskin_4(r_abs),
            DeltaFunction::Peskin4 => self.peskin_4(r_abs),
        }
    }

    /// Roma & Peskin 3-point delta function
    fn roma_peskin_3(&self, r: T) -> T {
        if r >= T::from_f64(1.5).unwrap_or_else(T::zero) {
            return T::zero();
        }

        let one = T::one();
        let half = T::from_f64(0.5).unwrap_or_else(T::zero);

        if r <= half {
            (one + (T::from_f64(TWO).unwrap_or_else(T::one) * r
                - T::from_f64(THREE).unwrap_or_else(T::one))
                * r
                * r)
                / T::from_f64(THREE).unwrap_or_else(T::one)
        } else {
            (T::from_f64(FIVE).unwrap_or_else(T::one)
                - T::from_f64(THREE).unwrap_or_else(T::one) * r
                - ((T::from_f64(4.0).unwrap_or_else(T::one)
                    - T::from_f64(TWO).unwrap_or_else(T::one) * r)
                    * r)
                    .sqrt())
                / T::from_f64(6.0).unwrap_or_else(T::one)
        }
    }

    /// Roma & Peskin 4-point delta function
    fn roma_peskin_4(&self, r: T) -> T {
        if r >= T::from_f64(TWO).unwrap_or_else(T::zero) {
            return T::zero();
        }

        let one = T::one();

        if r <= one {
            (T::from_f64(THREE).unwrap_or_else(T::one)
                - T::from_f64(TWO).unwrap_or_else(T::one) * r
                + ((one + T::from_f64(4.0).unwrap_or_else(T::one) * r
                    - T::from_f64(4.0).unwrap_or_else(T::one) * r * r)
                    .sqrt()))
                / T::from_f64(8.0).unwrap_or_else(T::one)
        } else {
            (T::from_f64(FIVE).unwrap_or_else(T::one)
                - T::from_f64(TWO).unwrap_or_else(T::one) * r
                - ((-T::from_f64(7.0).unwrap_or_else(T::one)
                    + T::from_f64(12.0).unwrap_or_else(T::one) * r
                    - T::from_f64(4.0).unwrap_or_else(T::one) * r * r)
                    .sqrt()))
                / T::from_f64(8.0).unwrap_or_else(T::one)
        }
    }

    /// Peskin 4-point delta function
    fn peskin_4(&self, r: T) -> T {
        if r >= T::from_f64(TWO).unwrap_or_else(T::zero) {
            return T::zero();
        }

        let one = T::one();
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one);

        if r <= one {
            (one + (pi * r).cos()) / T::from_f64(4.0).unwrap_or_else(T::one)
        } else {
            (one + (pi * (T::from_f64(TWO).unwrap_or_else(T::one) - r)).cos())
                / T::from_f64(4.0).unwrap_or_else(T::one)
        }
    }

    /// Get stencil size for interpolation
    pub fn stencil_size(&self) -> usize {
        match self.delta_type {
            DeltaFunction::RomaPeskin3 => 3,
            DeltaFunction::RomaPeskin4 | DeltaFunction::Peskin4 => INTERPOLATION_STENCIL_SIZE,
        }
    }
}
