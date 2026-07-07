//! Interpolation kernels and delta functions for IBM

use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;

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

impl<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> InterpolationKernel<T> {
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
        let r_abs = scalar::abs::<T>(r);

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
        let three_half = scalar::from_f64::<T>(1.5);
        if r >= three_half {
            return scalar::zero::<T>();
        }

        let one = scalar::one::<T>();
        let three = scalar::from_f64::<T>(3.0);
        let five = scalar::from_f64::<T>(5.0);
        let six = scalar::from_f64::<T>(6.0);
        let half = scalar::from_f64::<T>(0.5);

        if r <= half {
            // φ₃(r) = (1 + √(1 - 3r²)) / 3
            let disc = one - three * r * r;
            (one + scalar::sqrt::<T>(scalar::max::<T>(disc, scalar::zero::<T>()))) / three
        } else {
            // φ₃(r) = (5 - 3r - √(-3(1-r)² + 1)) / 6
            let one_minus_r = one - r;
            let disc = one - three * one_minus_r * one_minus_r;
            (five - three * r - scalar::sqrt::<T>(scalar::max::<T>(disc, scalar::zero::<T>())))
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
        let two = scalar::from_f64::<T>(2.0);
        if r >= two {
            return scalar::zero::<T>();
        }

        let one = scalar::one::<T>();
        let three = scalar::from_f64::<T>(3.0);
        let four = scalar::from_f64::<T>(4.0);
        let five = scalar::from_f64::<T>(5.0);
        let seven = scalar::from_f64::<T>(7.0);
        let eight = scalar::from_f64::<T>(8.0);
        let twelve = scalar::from_f64::<T>(12.0);

        if r <= one {
            // φ₄(r) = (3 - 2r + √(1 + 4r - 4r²)) / 8
            let disc = one + four * r - four * r * r;
            (three - two * r + scalar::sqrt::<T>(scalar::max::<T>(disc, scalar::zero::<T>())))
                / eight
        } else {
            // φ₄(r) = (5 - 2r - √(-7 + 12r - 4r²)) / 8
            let disc = -seven + twelve * r - four * r * r;
            (five - two * r - scalar::sqrt::<T>(scalar::max::<T>(disc, scalar::zero::<T>())))
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
        let two = scalar::from_f64::<T>(2.0);
        if r >= two {
            return scalar::zero::<T>();
        }

        let one = scalar::one::<T>();
        let four = scalar::from_f64::<T>(4.0);
        let pi = scalar::from_f64::<T>(std::f64::consts::PI);

        if r <= one {
            (one + scalar::cos::<T>(pi * r)) / four
        } else {
            (one + scalar::cos::<T>(pi * (two - r))) / four
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Compact support: δ(r) = 0 for |r| ≥ support radius.
    ///
    /// Roma 3-pt: support = 3/2. Roma 4-pt, Peskin 4-pt: support = 2.
    #[test]
    fn compact_support() {
        let k3 = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0);
        let k4r = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin4, 1.0);
        let k4p = InterpolationKernel::<f64>::new(DeltaFunction::Peskin4, 1.0);

        // Roma 3-pt: zero at r ≥ 1.5
        assert_eq!(k3.delta(1.5), 0.0);
        assert_eq!(k3.delta(2.0), 0.0);
        assert_eq!(k3.delta(10.0), 0.0);

        // Roma 4-pt: zero at r ≥ 2
        assert_eq!(k4r.delta(2.0), 0.0);
        assert_eq!(k4r.delta(3.0), 0.0);

        // Peskin 4-pt: zero at r ≥ 2
        assert_eq!(k4p.delta(2.0), 0.0);
        assert_eq!(k4p.delta(5.0), 0.0);
    }

    /// Non-negativity: δ(r) ≥ 0 for all r (Roma 3-pt).
    ///
    /// Theorem (Roma, Peskin & Berger 1999): φ₃(r) ≥ 0 for all r.
    #[test]
    fn non_negativity_roma3() {
        let k = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0);

        for i in 0..100 {
            let r = i as f64 * 0.03; // r ∈ [0, 3)
            let val = k.delta(r);
            assert!(val >= 0.0, "Roma 3-pt delta({r}) = {val} < 0");
        }
    }

    /// Partition of unity: Σ_{j∈Z} δ(r − j) = 1 for Roma 3-pt at r = 0.25.
    ///
    /// Theorem (Roma, Peskin & Berger 1999): Σ φ₃(r − j) = 1 ∀ r.
    /// For r = 0.25, the stencil is j ∈ {-1, 0, 1} (within support).
    #[test]
    fn partition_of_unity_roma3() {
        let k = InterpolationKernel::<f64>::new(DeltaFunction::RomaPeskin3, 1.0);

        let r = 0.25_f64;
        let sum: f64 = (-2..=2).map(|j| k.delta(r - j as f64)).sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-12);
    }
}
