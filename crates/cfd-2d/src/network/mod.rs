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
//!
//! The optional projected solve path also returns per-channel projection
//! summaries, including solver-domain extents and fluid-cell occupancy, so the
//! schematics-to-grid mapping is observable in tests and benchmarks.

mod build;
mod coupled;
mod channel;
mod postprocess;
mod projection;
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
pub use projection::NetworkProjectionSummary;
pub use types::{
    Channel2dResult, ChannelProjectionSummary, CoupledNetwork2dResult, Network2DSolver,
    Network2dResult, ProjectedNetwork2dResult,
};
pub use validation::validate_blueprint_for_2d_projection;
