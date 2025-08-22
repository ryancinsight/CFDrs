//! Channel geometry and flow characteristics for 1D CFD.
//!
//! This module provides advanced channel modeling capabilities including
//! complex geometries, surface effects, and flow regime transitions.

mod geometry;
mod cross_section;
mod surface;
mod flow;
mod solver;

pub use geometry::{ChannelGeometry, ChannelType, GeometricVariation};
pub use cross_section::CrossSection;
pub use surface::{SurfaceProperties, Wettability};
pub use flow::{Channel, FlowState, FlowRegime, NumericalParameters};