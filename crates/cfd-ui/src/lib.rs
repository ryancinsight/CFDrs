//! cfd-ui — CAD-like GUI for millifluidic CFD simulation.
//!
//! This crate provides a desktop application for creating and visualizing
//! meshes (via `cfd-mesh`), designing and running CFD simulations (via `cfd-3d`),
//! and exporting engineering drawings on standard templates.
//!
//! # Architecture
//!
//! Four Clean Architecture layers (inner to outer):
//! - **Domain**: Pure data models — scene graph, camera, drawing specs, document
//! - **Application**: Use-case orchestrators — mesh creation, simulation, export
//! - **Infrastructure**: GPU rendering (wgpu), file I/O, SVG/PDF export
//! - **Presentation**: UI views, panels, toolbars, dialogs, menus

#![warn(clippy::all)]

pub mod domain;
pub mod application;
pub mod infrastructure;
pub mod presentation;

#[cfg(target_arch = "wasm32")]
pub mod wasm;
