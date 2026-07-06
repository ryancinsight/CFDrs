//! Scalar contract for the 1D network solver.
//!
//! `cfd-1d` still uses nalgebra for its sparse/dense linear-system boundary
//! while CFDrs core physics and geometry contracts are migrating to Eunomia.
//! This trait is the single scalar seam for that transition: every solver scalar
//! must satisfy both the remaining nalgebra matrix requirements and the Eunomia
//! scalar provider contract.

use cfd_core::conversion::{SafeFromF64, SafeFromUsize};
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};

/// Real scalar supported by the 1D network solver.
pub trait Cfd1dScalar:
    nalgebra::RealField
    + EunomiaRealField
    + FloatElement
    + NumericElement
    + Copy
    + SafeFromF64
    + SafeFromUsize
    + std::fmt::Debug
    + Send
    + Sync
    + 'static
{
}

impl<T> Cfd1dScalar for T where
    T: nalgebra::RealField
        + EunomiaRealField
        + FloatElement
        + NumericElement
        + Copy
        + SafeFromF64
        + SafeFromUsize
        + std::fmt::Debug
        + Send
        + Sync
        + 'static
{
}
