//! Bridge between `scheme` 2D schematics and `cfd-1d` simulation networks.
//!
//! This module converts a [`scheme::geometry::ChannelSystem`] into a
//! [`crate::network::Network`] that the 1D solver can consume directly.
//!
//! # Pipeline
//!
//! ```text
//! scheme::ChannelSystem ──► SchemeNetworkConverter ──► Network<f64>
//! ```
//!
//! The converter:
//! 1. Infers node types (inlet / outlet / junction) from channel connectivity.
//! 2. Derives channel length from the 2D path polyline.
//! 3. Maps scheme `Channel` width/height to a rectangular cross-section
//!    with the appropriate hydraulic resistance.
//! 4. Handles frustum (tapered) channels by using average cross-section
//!    dimensions for hydraulic calculations.
//! 5. Computes initial resistance coefficients using the
//!    `ResistanceCalculator` when a fluid model is provided.
//!
//! # Example
//!
//! ```rust,no_run
//! use scheme::geometry::generator::create_geometry;
//! use scheme::geometry::SplitType;
//! use scheme::config::{GeometryConfig, ChannelTypeConfig};
//! use cfd_1d::scheme_bridge::SchemeNetworkConverter;
//!
//! let system = create_geometry(
//!     (200.0, 100.0),
//!     &[SplitType::Bifurcation],
//!     &GeometryConfig::default(),
//!     &ChannelTypeConfig::AllStraight,
//! );
//!
//! let converter = SchemeNetworkConverter::new(&system);
//! let network = converter.build_network_with_water().unwrap();
//!
//! println!("Nodes: {}, Edges: {}", network.node_count(), network.edge_count());
//! ```

mod converter;
mod error;

pub use converter::SchemeNetworkConverter;
pub use error::BridgeError;
