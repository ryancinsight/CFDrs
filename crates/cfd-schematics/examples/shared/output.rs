use cfd_schematics::domain::model::NetworkBlueprint;
use cfd_schematics::plot_geometry;
use std::fs;
use std::io::{self, Write};
use std::path::PathBuf;

fn workspace_root() -> PathBuf {
    std::env::current_dir()
        .expect("Failed to get current directory")
        .ancestors()
        .find(|path| path.join("Cargo.toml").exists() && path.join("crates").exists())
        .expect("Must be run from within the CFDrs workspace")
        .to_path_buf()
}

fn example_output_dir(example_name: &str) -> PathBuf {
    workspace_root()
        .join("crates")
        .join("cfd-schematics")
        .join("output")
        .join("examples")
        .join(example_name)
}

/// Standardized output generator for all schematic examples.
/// This ensures a unified directory structure: `output/examples/<example_name>/`
/// containing both the `blueprint.json` and `geometry.svg` for each run.
#[allow(dead_code)]
pub fn save_example_output(blueprint: &NetworkBlueprint, example_name: &str) {
    let output_dir = example_output_dir(example_name);

    fs::create_dir_all(&output_dir).unwrap_or_else(|error| {
        panic!(
            "failed to create output directory {}: {error}",
            output_dir.display()
        )
    });

    writeln!(io::stdout(), "Saving outputs to {}", output_dir.display())
        .expect("failed to report example output directory");

    // 1. Save JSON Blueprint
    let json_path = output_dir.join(format!("{example_name}.json"));
    let json_data =
        serde_json::to_string_pretty(blueprint).expect("Failed to serialize blueprint to JSON");
    fs::write(&json_path, json_data).unwrap_or_else(|error| {
        panic!("failed to write JSON path {}: {error}", json_path.display())
    });
    writeln!(io::stdout(), "  - JSON: {}", json_path.display())
        .expect("failed to report JSON output");

    // 2. Save SVG Visualization
    let svg_path = output_dir.join(format!("{example_name}.svg"));
    plot_geometry(blueprint, &svg_path)
        .unwrap_or_else(|error| panic!("failed to plot geometry {}: {error}", svg_path.display()));
    writeln!(io::stdout(), "  - SVG : {}", svg_path.display())
        .expect("failed to report SVG output");

    // Print statistics
    writeln!(io::stdout(), "\nBlueprint Statistics:")
        .expect("failed to report blueprint statistics");
    writeln!(io::stdout(), "  - Nodes: {}", blueprint.nodes.len())
        .expect("failed to report node statistics");
    writeln!(io::stdout(), "  - Channels: {}", blueprint.channels.len())
        .expect("failed to report channel statistics");
}

/// Standardized output generator that allows specifying a custom filename
/// within the example's directory. Useful for examples that generate multiple permutations.
#[allow(dead_code)]
pub fn save_example_output_with_name(
    blueprint: &NetworkBlueprint,
    example_name: &str,
    file_name: &str,
) {
    let output_dir = example_output_dir(example_name);

    fs::create_dir_all(&output_dir).unwrap_or_else(|error| {
        panic!(
            "failed to create output directory {}: {error}",
            output_dir.display()
        )
    });

    // 1. Save JSON Blueprint
    let json_path = output_dir.join(format!("{file_name}.json"));
    let json_data =
        serde_json::to_string_pretty(blueprint).expect("Failed to serialize blueprint to JSON");
    fs::write(&json_path, json_data).unwrap_or_else(|error| {
        panic!("failed to write JSON path {}: {error}", json_path.display())
    });

    // 2. Save SVG Visualization
    let svg_path = output_dir.join(format!("{file_name}.svg"));
    plot_geometry(blueprint, &svg_path)
        .unwrap_or_else(|error| panic!("failed to plot geometry {}: {error}", svg_path.display()));

    writeln!(io::stdout(), "Saved -> {}", svg_path.display()).expect("failed to report SVG output");
}
