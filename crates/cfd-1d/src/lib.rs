//! 1D CFD simulations and microfluidic network solvers.
//! 
//! This crate provides 1D computational fluid dynamics solvers specifically
//! designed for microfluidic and millifluidic networks. It includes support
//! for:
//! - Channel-based flow networks
//! - Electrical circuit analogy solvers
//! - Pressure-driven flow simulations
//! - Component-based microfluidic devices
//! 
//! ## 2D Schematic Support
//! 
//! This crate is designed to integrate with the `scheme` library for
//! 2D schematic representation of 1D microfluidic networks, similar to
//! how `cfd-3d` integrates with `csgrs` for 3D mesh handling.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod network;
pub mod channel;
pub mod components;
pub mod solver;
pub mod resistance;

#[cfg(feature = "scheme-integration")]
pub mod scheme_integration;

// When scheme integration is not enabled, provide a stub module
#[cfg(not(feature = "scheme-integration"))]
pub mod scheme_integration;

// TODO: Implement these exports
// pub use network::{Network, Node, NetworkBuilder};
// pub use channel::{Channel, ChannelGeometry};
// pub use components::{Pump, Valve, Sensor, Component};
// pub use solver::{NetworkSolver, ElectricalAnalogySolver, HagenPoiseuilleSolver};
// pub use resistance::{ResistanceModel, RectangularChannel, CircularChannel};

/// Common 1D CFD types and traits
pub mod prelude {
    pub use crate::scheme_integration::{SchemeConversion, ComponentType, SchematicLayout};
    
    #[cfg(feature = "scheme-integration")]
    pub use crate::scheme_integration::{JunctionType, ChannelPathType, helpers};
    
    // TODO: Add exports when implemented
    // pub use crate::{
    //     network::{Network, NetworkBuilder},
    //     channel::Channel,
    //     components::{Pump, Component},
    //     solver::NetworkSolver,
    //     resistance::ResistanceModel,
    // };
}