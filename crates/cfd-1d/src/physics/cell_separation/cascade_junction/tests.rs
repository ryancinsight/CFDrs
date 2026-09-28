//! Unit tests for the Zweifach-Fung cascade junction model.
//!
//! Split by behaviour family: each leaf covers one routing concern, so a
//! failure names the physics it broke instead of pointing at a single
//! thousand-line file. The shared prelude below is re-exported so every leaf
//! reaches the model through one `use super::*;`.

pub(crate) use super::cascade_routing::*;
pub(crate) use super::incremental_filtration::*;
pub(crate) use super::routing_probability::*;
pub(crate) use super::*;
pub(crate) use aequitas::systems::si::quantities::{Length, Velocity, VolumetricFlowRate};

pub(crate) fn length(value: f64) -> Length {
    Length::from_base(value)
}

pub(crate) fn flow(value: f64) -> VolumetricFlowRate {
    VolumetricFlowRate::from_base(value)
}

mod cascade_split;
mod core_probability;
mod fahraeus;
mod kappa_aware;
mod peripheral_recovery;
