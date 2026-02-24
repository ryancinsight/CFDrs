//! Ruppert's refinement algorithm and supporting quality metrics.
//!
//! # Overview
//!
//! After constructing a CDT, Ruppert's algorithm inserts Steiner points to
//! guarantee a minimum angle bound (typically > 20.7° for a radius-edge ratio
//! bound of √2).

pub mod circumcenter;
pub mod encroachment;
pub mod quality;
pub mod ruppert;

pub use quality::TriangleQuality;
pub use ruppert::RuppertRefiner;
