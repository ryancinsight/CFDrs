//! Scalar contract for the 1D network solver.
//!
//! `cfd-1d` uses Eunomia and Leto for all linear-system operations.
//! This trait is the single scalar seam that must satisfy all Atlas provider contracts.

use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};
use leto_ops::RealScalar;

/// Real scalar supported by the 1D network solver.
pub trait Cfd1dScalar:
    EunomiaRealField
    + FloatElement
    + NumericElement
    + RealScalar
    + Copy
    + SafeFromF64
    + SafeFromUsize
    + std::fmt::Debug
    + std::fmt::Display
    + std::ops::DivAssign
    + Send
    + Sync
    + 'static
{
    /// Additive identity from the Eunomia numeric contract.
    #[inline]
    fn zero() -> Self {
        <Self as NumericElement>::ZERO
    }

    /// Multiplicative identity from the Eunomia numeric contract.
    #[inline]
    fn one() -> Self {
        <Self as NumericElement>::ONE
    }
}

impl<T> Cfd1dScalar for T where
    T: EunomiaRealField
        + FloatElement
        + NumericElement
        + RealScalar
        + Copy
        + SafeFromF64
        + SafeFromUsize
        + std::fmt::Debug
        + std::fmt::Display
        + std::ops::DivAssign
        + Send
        + Sync
        + 'static
{
}
