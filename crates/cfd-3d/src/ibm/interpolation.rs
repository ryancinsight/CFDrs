//! Interpolation kernels and delta functions for IBM

use nalgebra::RealField;
use num_traits::FromPrimitive;

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
pub struct InterpolationKernel<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    delta_type: DeltaFunction,
    width: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> InterpolationKernel<T> {
    /// Create a new interpolation kernel
    pub fn new(delta_type: DeltaFunction, width: T) -> Self {
        Self { delta_type, width }
    }

    /// Kernel support half-width in grid-spacing units.
    pub fn width(&self) -> T {
        self.width
    }

    /// Delta function type used by this kernel.
    pub fn delta_type(&self) -> &DeltaFunction {
        &self.delta_type
    }

    /// Evaluate the delta function
    pub fn delta(&self, r: T) -> T {
        let r_abs = num_traits::Float::abs(r);

        match self.delta_type {
            DeltaFunction::RomaPeskin3 => self.roma_peskin_3(r_abs),
            DeltaFunction::RomaPeskin4 => self.roma_peskin_4(r_abs),
            DeltaFunction::Peskin4 => self.peskin_4(r_abs),
        }
    }

    /// Roma, Peskin & Berger (1999) 3-point delta function φ₃(r):
    ///
    /// # Theorem
    ///
    /// The 3-point regularized delta function satisfies:
    /// 1. Σ φ₃(r - j) = 1 for all r (partition of unity)
    /// 2. φ₃(r) ≥ 0 for all r (non-negativity)
    /// 3. φ₃(r) = 0 for |r| ≥ 3/2 (compact support)
    ///
    /// **Proof sketch**: Direct substitution into the piecewise formula
    /// and summation over the stencil j ∈ {-1, 0, 1}.
    fn roma_peskin_3(&self, r: T) -> T {
        let three_half = <T as FromPrimitive>::from_f64(1.5)
            .expect("1.5 is an IEEE 754 representable f64 constant");
        if r >= three_half {
            return T::zero();
        }

        let one = T::one();
        let three = <T as FromPrimitive>::from_f64(3.0)
            .expect("3.0 is representable in all IEEE 754 types");
        let five = <T as FromPrimitive>::from_f64(5.0)
            .expect("5.0 is representable in all IEEE 754 types");
        let six = <T as FromPrimitive>::from_f64(6.0)
            .expect("6.0 is representable in all IEEE 754 types");
        let half = <T as FromPrimitive>::from_f64(0.5)
            .expect("0.5 is exactly representable in IEEE 754");

        if r <= half {
            // φ₃(r) = (1 + √(1 - 3r²)) / 3
            let disc = one - three * r * r;
            (one + num_traits::Float::sqrt(num_traits::Float::max(disc, T::zero()))) / three
        } else {
            // φ₃(r) = (5 - 3r - √(-3(1-r)² + 1)) / 6
            let one_minus_r = one - r;
            let disc = one - three * one_minus_r * one_minus_r;
            (five - three * r - num_traits::Float::sqrt(num_traits::Float::max(disc, T::zero())))
                / six
        }
    }

    /// Roma, Peskin & Berger (1999) 4-point delta function φ₄(r):
    ///
    /// # Theorem
    ///
    /// The 4-point regularized delta function satisfies:
    /// 1. Σ φ₄(r - j) = 1 for all r (partition of unity)
    /// 2. φ₄(r) ≥ 0 for all r (non-negativity)  
    /// 3. φ₄(r) = 0 for |r| ≥ 2 (compact support)
    /// 4. Σ (r - j) φ₄(r - j) = 0 (first moment condition)
    ///
    /// **Proof sketch**: The piecewise formula is derived by requiring conditions
    /// 1–4 plus continuity at r = 0, 1, 2.
    fn roma_peskin_4(&self, r: T) -> T {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        if r >= two {
            return T::zero();
        }

        let one = T::one();
        let three = <T as FromPrimitive>::from_f64(3.0)
            .expect("3.0 is representable in all IEEE 754 types");
        let four = <T as FromPrimitive>::from_f64(4.0)
            .expect("4.0 is representable in all IEEE 754 types");
        let five = <T as FromPrimitive>::from_f64(5.0)
            .expect("5.0 is representable in all IEEE 754 types");
        let seven = <T as FromPrimitive>::from_f64(7.0)
            .expect("7.0 is representable in all IEEE 754 types");
        let eight = <T as FromPrimitive>::from_f64(8.0)
            .expect("8.0 is representable in all IEEE 754 types");
        let twelve = <T as FromPrimitive>::from_f64(12.0)
            .expect("12.0 is representable in all IEEE 754 types");

        if r <= one {
            // φ₄(r) = (3 - 2r + √(1 + 4r - 4r²)) / 8
            let disc = one + four * r - four * r * r;
            (three - two * r + num_traits::Float::sqrt(num_traits::Float::max(disc, T::zero())))
                / eight
        } else {
            // φ₄(r) = (5 - 2r - √(-7 + 12r - 4r²)) / 8
            let disc = -seven + twelve * r - four * r * r;
            (five - two * r - num_traits::Float::sqrt(num_traits::Float::max(disc, T::zero())))
                / eight
        }
    }

    /// Peskin (2002) cosine 4-point delta function:
    ///
    /// φ(r) = (1 + cos(π r/2)) / 4 for |r| < 2, else 0
    ///
    /// # Theorem
    ///
    /// This kernel satisfies the partition-of-unity property.
    fn peskin_4(&self, r: T) -> T {
        let two = <T as FromPrimitive>::from_f64(2.0)
            .expect("2.0 is representable in all IEEE 754 types");
        if r >= two {
            return T::zero();
        }

        let one = T::one();
        let four = <T as FromPrimitive>::from_f64(4.0)
            .expect("4.0 is representable in all IEEE 754 types");
        let pi = <T as FromPrimitive>::from_f64(std::f64::consts::PI)
            .expect("π is an IEEE 754 representable f64 constant");

        if r <= one {
            (one + num_traits::Float::cos(pi * r)) / four
        } else {
            (one + num_traits::Float::cos(pi * (two - r))) / four
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
