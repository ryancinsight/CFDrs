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
// ── Pedantic suppressions ────────────────────────────────────────────────────
// GPU/WGPU rendering pipelines and UI coordinate geometry inherently use
// index-to-float casts for clip-space transforms and pixel addressing.
// Render-pass functions require many arguments to configure pipeline state
// without indirection. Complex GPU pipeline structs carry many orthogonal
// boolean mode flags; bitset refactors add complexity without safety gains.
// These allow-attrs suppress warnings raised by external `-- -W pedantic`
// invocations that bypass this crate's intentional lint profile.
#![allow(
    clippy::cast_possible_truncation, // pixel coords bounded by window dimensions
    clippy::cast_precision_loss,      // usize→f32 for clip-space; values < 2^24
    clippy::cast_sign_loss,           // non-negative component index arithmetic
    clippy::cast_possible_wrap,       // modular screen-tile arithmetic
    clippy::similar_names,            // pos_x/pos_y, uv_x/uv_y — GPU shader convention
    clippy::too_many_lines,           // wgpu pipeline builders are inherently verbose
    clippy::too_many_arguments,       // render-pass fns must accept all pipeline params
    clippy::missing_panics_doc,       // GPU/wgpu helpers panic only on device loss (fatal)
    clippy::missing_errors_doc,       // error paths documented inline in struct variants
    clippy::missing_fields_in_debug,  // large GPU buffer fields omitted from Debug output
    clippy::struct_excessive_bools,   // pipeline mode flags are orthogonal bit fields
    clippy::needless_pass_by_value,   // mesh/scene ownership model passes by value
    clippy::trivially_copy_pass_by_ref, // Copy coords passed by ref for pipeline symmetry
    clippy::many_single_char_names,  // ray intersection / geometry uses t,a,b,c,d — math convention
)]

pub mod application;
pub mod domain;
pub mod infrastructure;
pub mod presentation;

#[cfg(target_arch = "wasm32")]
pub mod wasm;
