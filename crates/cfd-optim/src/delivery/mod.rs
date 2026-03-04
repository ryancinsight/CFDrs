//! Artifact delivery: pipeline, JSON/SVG export, and CSV export.
//!
//! All output-concern code lives here (SoC).  Nothing in this module
//! computes metrics or scores candidates.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`csv_export`] | Per-channel hemolysis CSV |
//! [`export`]     | JSON and SVG report artefacts |
//! [`pipeline`]   | End-to-end mesh + schematic pipeline (mesh-export feature) |

pub mod csv_export;
pub mod export;
#[cfg(feature = "mesh-export")]
pub mod pipeline;

pub use csv_export::save_per_channel_csv;
pub use export::{
    save_all_modes_json, save_annotated_cif_svg, save_comparison_svg, save_schematic_svg,
    save_top5_json,
};
#[cfg(feature = "mesh-export")]
pub use pipeline::{DesignArtifacts, DesignPipeline};
