use cfd_schematics::domain::model::NetworkBlueprint;
use cfd_schematics::plot_geometry;
use std::fs;
use std::path::PathBuf;

/// Standardized output generator for all schematic examples.
/// This ensures a unified directory structure: `output/examples/<example_name>/`
/// containing both the `blueprint.json` and `geometry.svg` for each run.
pub fn save_example_output(blueprint: &NetworkBlueprint, example_name: &str) {
    let current_dir = std::env::current_dir().expect("Failed to get current directory");
    let workspace_root = current_dir
        .ancestors()
        .find(|p| p.join("Cargo.toml").exists() && p.join("crates").exists())
        .expect("Must be run from within the CFDrs workspace");

    // Output directory: CFDrs/crates/cfd-schematics/output/examples/<name>/
    let output_dir = workspace_root
        .join("crates")
        .join("cfd-schematics")
        .join("output")
        .join("examples")
        .join(example_name);

    fs::create_dir_all(&output_dir)
        .unwrap_or_else(|e| panic!("Failed to create output directory {:?}: {}", output_dir, e));

    println!("Saving outputs to {:?}", output_dir);

    // 1. Save JSON Blueprint
    let json_path = output_dir.join(format!("{}.json", example_name));
    let json_data =
        serde_json::to_string_pretty(blueprint).expect("Failed to serialize blueprint to JSON");
    fs::write(&json_path, json_data)
        .unwrap_or_else(|e| panic!("Failed to write JSON path {:?}: {}", json_path, e));
    println!("  - JSON: {:?}", json_path.file_name().unwrap());

    // 2. Save SVG Visualization
    let svg_path = output_dir.join(format!("{}.svg", example_name));
    plot_geometry(blueprint, svg_path.to_str().unwrap())
        .map_err(|e| e.to_string())
        .unwrap_or_else(|e| panic!("Failed to plot geometry {:?}: {}", svg_path, e));
    println!("  - SVG : {:?}", svg_path.file_name().unwrap());

    // Print statistics
    println!("\nBlueprint Statistics:");
    println!("  - Nodes: {}", blueprint.nodes.len());
    println!("  - Channels: {}", blueprint.channels.len());
}

/// Standardized output generator that allows specifying a custom filename
/// within the example's directory. Useful for examples that generate multiple permutations.
pub fn save_example_output_with_name(
    blueprint: &NetworkBlueprint,
    example_name: &str,
    file_name: &str,
) {
    let current_dir = std::env::current_dir().expect("Failed to get current directory");
    let workspace_root = current_dir
        .ancestors()
        .find(|p| p.join("Cargo.toml").exists() && p.join("crates").exists())
        .expect("Must be run from within the CFDrs workspace");

    // Output directory: CFDrs/crates/cfd-schematics/output/examples/<name>/
    let output_dir = workspace_root
        .join("crates")
        .join("cfd-schematics")
        .join("output")
        .join("examples")
        .join(example_name);

    fs::create_dir_all(&output_dir)
        .unwrap_or_else(|e| panic!("Failed to create output directory {:?}: {}", output_dir, e));

    // 1. Save JSON Blueprint
    let json_path = output_dir.join(format!("{}.json", file_name));
    let json_data =
        serde_json::to_string_pretty(blueprint).expect("Failed to serialize blueprint to JSON");
    fs::write(&json_path, json_data)
        .unwrap_or_else(|e| panic!("Failed to write JSON path {:?}: {}", json_path, e));

    // 2. Save SVG Visualization
    let svg_path = output_dir.join(format!("{}.svg", file_name));
    plot_geometry(blueprint, svg_path.to_str().unwrap())
        .map_err(|e| e.to_string())
        .unwrap_or_else(|e| panic!("Failed to plot geometry {:?}: {}", svg_path, e));

    println!("Saved -> {:?}", svg_path.file_name().unwrap());
}
