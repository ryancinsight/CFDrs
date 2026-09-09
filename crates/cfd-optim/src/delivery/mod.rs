//! Artifact delivery: blueprint-native JSON/SVG export.
//!
//! All output-concern code lives here (SoC). Nothing in this module
//! computes metrics or scores candidates.
//!
//! | Sub-module | Responsibility |
//! |------------|----------------|
//! [`export`]   | JSON and SVG report artefacts |

mod export;

pub use export::{
    load_pareto_points, load_top5_report_json, save_blueprint_schematic_svg, save_json_pretty,
    save_pareto_points, save_top5_report_json,
};
