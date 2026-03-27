//! Network builder for constructing fluid networks.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`blueprint_conversion`] | `network_from_blueprint()` — canonical entry-point |
//! | [`venturi_coefficients`] | Venturi-specific (R, K) computation |
//! | [`network_builder`] | `NetworkBuilder<T>` — low-level graph construction |

pub mod blueprint_conversion;
pub mod network_builder;
pub(crate) mod venturi_coefficients;

pub use blueprint_conversion::network_from_blueprint;
pub use network_builder::NetworkBuilder;
