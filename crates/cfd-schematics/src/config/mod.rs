//! config/mod.rs - 2D Schematic Configuration (Hierarchical)
//!
//! This module provides comprehensive configuration management for 2D microfluidic
//! schematic generation. It centralizes all configuration logic and eliminates
//! Single Source of Truth (SSOT) violations by providing validated configuration
//! types with clear constraints and relationships.
//!
//! # Design Principles
//!
//! - **Single Source of Truth**: All configuration values are defined once
//! - **Validation**: All configurations are validated at creation time
//! - **Immutability**: Configurations are immutable after creation
//! - **Composability**: Complex configurations are built from simpler ones
//! - **Discoverability**: Presets and builders make common configurations easy

pub mod adaptive;
pub mod channel;
pub mod constants;
pub mod geometry;
pub mod presets;

// Re-export constants module content
pub use constants::*;

// Re-export adaptive config
pub use adaptive::*;

// Re-export channel configs (flat)
pub use channel::*;

// Re-export geometry configs (flat)
pub use geometry::*;

// Re-export presets as a module?
// Original config.rs had `pub mod presets`.
// So we keep `pub mod presets`, which is already declared above.
