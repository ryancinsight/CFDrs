//! Slope limiters for Discontinuous Galerkin methods.
//!
//! This module provides various slope limiters that can be used to control
//! oscillations in DG solutions, especially in the presence of shocks or discontinuities.

//! This module root is a passthrough facade: each limiter lives in its own leaf
//! behind the shared [`Limiter`] trait, so adding a limiter means adding a file
//! rather than editing one thousand-line root.

mod factory;
mod minmod;
mod moment;
mod none;
mod params;
mod traits;
mod tvb;
mod weno;

pub use factory::LimiterFactory;
pub use minmod::MinmodLimiter;
pub use moment::MomentLimiter;
pub use none::NoLimiter;
pub use params::{LimiterParams, LimiterType};
pub use traits::Limiter;
pub use tvb::TVBLimiter;
pub use weno::WENOLimiter;

#[cfg(test)]
mod tests;
