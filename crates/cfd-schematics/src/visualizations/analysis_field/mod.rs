//! Typed CFD analysis overlays.
//!
//! CFD owns field interpretation and sparse channel/node association. Iris owns
//! normalized color laws. Overlay construction validates scalar values and
//! computes ranges once, so rendering performs constant-time color lookup
//! without temporary range allocations.

mod color;
mod field;
mod overlay;
mod range;

pub use color::colorize;
pub use field::AnalysisField;
pub use overlay::AnalysisOverlay;

#[cfg(test)]
mod tests;
