//! The scalar seam for the whole CFD suite.
//!
//! Every solver, discretization, and validation kernel in this workspace is
//! generic over one real scalar type. [`CfdScalar`] is that single contract:
//! there is no per-crate, per-solver, or per-kernel variant of it.

use core::fmt::Display;
use core::ops::DivAssign;

/// Real scalar over which every CFD kernel in this workspace is generic.
///
/// The canonical stack seam is [`eunomia::RealField`], which already supplies
/// the numeric surface (`FloatElement`, `NumericElement`, `Copy`, `Default`,
/// `Debug`, `Send`, `Sync`, `'static`, `PartialOrd`, `Neg`, and the arithmetic
/// operators). This trait adds exactly the bounds the CFD domain needs on top
/// of it, and nothing else:
///
/// - [`leto_ops::RealScalar`] — every field in this workspace is stored in a
///   Leto array, and the dense norm/reduction kernels the linear solvers call
///   are methods on this trait. It supertraits [`leto_ops::Scalar`], which
///   supplies `from_usize` and the slice-level kernels used by the sparse and
///   stencil code.
/// - [`leto::ScalarOperand`] — required to place a scalar on the right-hand
///   side of Leto array arithmetic (`&field * dt`), which the time integrators
///   and relaxation loops rely on.
/// - [`DivAssign`] — `NumericElement` provides `Div` together with
///   `AddAssign`/`SubAssign`/`MulAssign` but not `DivAssign`; in-place
///   normalization of residuals and conserved fields needs it.
/// - [`Display`] — solver diagnostics, convergence traces, and validation
///   reports format scalar values directly.
///
/// The blanket implementation below is the only implementation: any type
/// satisfying the bounds is a `CfdScalar`, so admitting a new precision is a
/// change to `eunomia`/`leto`, never to this workspace.
///
/// Mesh-facing signatures additionally bind `cfd_mesh::domain::core::Scalar`
/// at the use site, because that contract belongs to the geometry substrate
/// rather than to the numerics; it is deliberately not a supertrait here so
/// that `cfd-core` stays free of a mesh dependency.
pub trait CfdScalar:
    eunomia::RealField + leto_ops::RealScalar + leto::ScalarOperand + DivAssign + Display
{
}

impl<T> CfdScalar for T where
    T: eunomia::RealField + leto_ops::RealScalar + leto::ScalarOperand + DivAssign + Display
{
}

#[cfg(test)]
mod tests {
    use super::CfdScalar;
    use eunomia::{FloatElement, NumericElement};

    fn round_trip<T: CfdScalar>(value: f64) -> f64 {
        let mut scalar = <T as FloatElement>::from_f64(value);
        scalar /= <T as FloatElement>::from_f64(2.0);
        <T as NumericElement>::to_f64(scalar)
    }

    #[test]
    fn seam_is_implemented_for_both_shipped_precisions() {
        assert!((round_trip::<f64>(3.5) - 1.75).abs() < f64::EPSILON);
        assert!((round_trip::<f32>(3.5) - 1.75).abs() < f64::from(f32::EPSILON));
    }
}
