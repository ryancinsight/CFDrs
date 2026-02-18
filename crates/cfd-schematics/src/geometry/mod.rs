//! geometry/mod.rs - 2D Microfluidic Schematic Geometry
//!
//! This module provides 2D geometry types and generation functions for
//! microfluidic schematic design, including bifurcation and trifurcation patterns.
//!
//! # Architecture
//!
//! The geometry module is organized into several submodules:
//! - `types`: Core geometric types and data structures
//! - `strategies`: Channel type generation strategies (Strategy pattern)
//! - `generator`: Main geometry generation logic with optional metadata support
//! - `metadata`: Extensible metadata system for tracking additional information
//! - `builders`: Builder pattern implementations for nodes and channels
//! - `optimization`: Optimization algorithms for serpentine channels
//!
//! # Design Patterns
//!
//! This module implements several design patterns:
//! - **Strategy Pattern**: For channel type generation (`strategies` module)
//! - **Factory Pattern**: For creating appropriate strategies
//! - **Builder Pattern**: For constructing complex geometries and metadata

pub mod adaptive_collision;
pub mod builders;
pub mod collision_detection;
pub mod generator;
pub mod metadata;
pub mod optimization;
pub mod state_integration;
pub mod strategies;
pub mod types;

pub use self::{
    generator::{create_geometry, create_geometry_with_metadata, MetadataConfig},
    types::{
        Channel, ChannelSystem, ChannelType, ChannelTypeCategory, InterchangeChannel,
        InterchangeChannelProfile, InterchangeChannelSystem, InterchangeNode, Node, Point2D,
        SplitType,
    },
};

/// Alias for `Point2D` for backward compatibility and convenience
pub type Point = Point2D;
