//! Finite-difference scheme enumeration — re-export from the Atlas SSOT.
//!
//! The canonical definition lives in `leto-ops::application::diff`; this
//! module re-exports it under the `cfd-math` vocabulary so that higher
//! layers depend on one import path.

pub use leto_ops::FiniteDifferenceScheme;