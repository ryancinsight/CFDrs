//! 1D CFD simulations and microfluidic network solvers.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]

pub mod network;
pub mod channel;
pub mod components;
pub mod solver;
pub mod resistance;

// TODO: Implement these exports
// pub use network::{Network, Node, NetworkBuilder};
// pub use channel::{Channel, ChannelGeometry};
// pub use components::{Pump, Valve, Sensor, Component};
// pub use solver::{NetworkSolver, ElectricalAnalogySolver, HagenPoiseuilleSolver};
// pub use resistance::{ResistanceModel, RectangularChannel, CircularChannel};

/// Common 1D CFD types and traits
pub mod prelude {
    // TODO: Add exports when implemented
    // pub use crate::{
    //     network::{Network, NetworkBuilder},
    //     channel::Channel,
    //     components::{Pump, Component},
    //     solver::NetworkSolver,
    //     resistance::ResistanceModel,
    // };
}