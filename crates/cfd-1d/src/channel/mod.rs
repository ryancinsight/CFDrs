//! Channel geometry and flow characteristics for 1D CFD.
//!
//! This module provides advanced channel modeling capabilities including
//! complex geometries, surface effects, and flow regime transitions.

mod cross_section;
mod flow;
mod geometry;
mod solver;
mod surface;

pub use cross_section::CrossSection;
pub use flow::{Channel, FlowRegime, FlowState, NumericalParameters};
pub use geometry::{ChannelGeometry, ChannelType, GeometricVariation};
pub use surface::{SurfaceProperties, Wettability};
