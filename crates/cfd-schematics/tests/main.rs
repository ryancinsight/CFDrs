//! Hierarchical integration harness for schematic contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (6 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/blueprint_render_parity.rs"]
mod blueprint_render_parity;
#[path = "main/canonical_geometry_authoring.rs"]
mod canonical_geometry_authoring;
#[path = "main/overlap_reconstruction.rs"]
mod overlap_reconstruction;
#[path = "main/plotters_drawer.rs"]
mod plotters_drawer;
#[path = "main/preset_autolayout.rs"]
mod preset_autolayout;
#[path = "main/schematic_annotations.rs"]
mod schematic_annotations;
