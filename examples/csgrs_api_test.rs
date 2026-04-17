//! Canonical schematic export example.
//!
//! The historical example name is retained for compatibility, but the
//! authoritative 2D design API in this workspace is `cfd-schematics`, not an
//! external `csgrs` dependency. This example exercises the current blueprint
//! generation and annotated SVG export path end-to-end.

use std::fs;
use std::path::PathBuf;

use cfd_1d::validate_blueprint_for_1d_solve;
use cfd_schematics::visualizations::RenderConfig;
use cfd_schematics::{plot_blueprint_auto_annotated, venturi_rect};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_dir = PathBuf::from("outputs/examples/csgrs_api_test");
    fs::create_dir_all(&output_dir)?;

    let blueprint = venturi_rect("ui_venturi_rect", 2.0e-3, 0.8e-3, 1.0e-3, 4.0e-3);
    blueprint.validate()?;
    validate_blueprint_for_1d_solve(&blueprint)?;

    let json_path = output_dir.join("ui_venturi_rect.json");
    fs::write(&json_path, blueprint.to_json_pretty()?)?;

    let svg_path = output_dir.join("ui_venturi_rect.svg");
    let mut render = RenderConfig::well_plate_96_report_annotated();
    render.title = "Canonical venturi schematic export".to_string();
    plot_blueprint_auto_annotated(&blueprint, &svg_path.to_string_lossy(), &render)?;

    println!("Blueprint: {}", blueprint.name);
    println!("Geometry-authored: {}", blueprint.is_geometry_authored());
    println!("Channels: {}", blueprint.channels.len());
    println!("JSON: {}", json_path.display());
    println!("SVG: {}", svg_path.display());

    Ok(())
}
