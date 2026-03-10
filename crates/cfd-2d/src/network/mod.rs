//! Full-network 2D solving driven by authoritative `cfd-schematics` blueprints.
//!
//! The network builder validates a [`NetworkBlueprint`], computes a normalized
//! `cfd-1d` reference solve, and then configures one 2D Navier-Stokes solve per
//! blueprint channel using that trace as the inlet/outlet contract.
//!
//! Every channel result reports both field-resolved wall shear and a
//! field-integrated outlet-flow error against the corresponding 1D reference
//! channel flow so cross-fidelity audits stay tied to one canonical topology
//! solve path.

mod build;
mod postprocess;
mod reference;
mod solve;
#[cfg(test)]
mod tests;
mod types;
mod validation;

pub use build::Network2dBuilderSink;
pub use reference::{
    solve_reference_trace, ChannelReferenceTrace, NetworkReferenceTrace, NodeReferenceTrace,
};
pub use types::{Channel2dResult, Network2DSolver, Network2dResult};
pub use validation::validate_blueprint_for_2d_projection;
