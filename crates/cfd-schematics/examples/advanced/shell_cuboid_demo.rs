//! `shell_cuboid_demo` — TPMS-filled Shell Millifluidic Cuboid
//!
//! Demonstrates building a channelless millifluidic cuboid (shell only)
//! on a standard 96-well plate footprint (ANSI/SLAS 1-2004: 127.76 × 85.47 mm)
//! with a 2 mm wall thickness and a Gyroid TPMS fill.  The schematic is rendered
//! to an SVG showing the outer/inner rectangles, port stubs, and a cross-hatch
//! pattern indicating the TPMS lattice fill.
//!
//! # Run
//!
//! ```bash
//! cargo run -p cfd-schematics --example shell_cuboid_demo
//! ```
//!
//! # Outputs
//!
//! - `cfd-schematics/output/examples/shell_cuboid_demo/shell_cuboid_demo.svg`  — schematic SVG
//! - Stdout — interchange JSON for downstream mesh/CFD toolchains

use cfd_schematics::geometry::{ShellCuboid, TpmsFillSpec, TpmsSurfaceKind};
use cfd_schematics::visualizations::{PlottersRenderer, RenderConfig};
use std::fs;

fn main() {
    // ── 1. Geometry ──────────────────────────────────────────────────────────
    // 96-well plate footprint: 127.76 × 85.47 mm  |  Shell thickness: 2 mm
    let fill = TpmsFillSpec {
        surface: TpmsSurfaceKind::Gyroid,
        period_mm: 5.0,
        iso_value: 0.0,
        resolution: 64,
        gradient: None,
    };

    let shell = ShellCuboid::well_plate_96(2.0)
        .expect("failed to create shell cuboid")
        .with_tpms_fill(fill);

    tracing::info!("Shell cuboid created (96-well plate footprint):");
    tracing::info!(
        "  outer_dims : {:.2} mm × {:.2} mm",
        shell.outer_dims.0,
        shell.outer_dims.1
    );
    tracing::info!(
        "  inner_dims : {:.2} mm × {:.2} mm",
        shell.inner_dims.0,
        shell.inner_dims.1
    );
    tracing::info!("  shell_thick: {} mm", shell.shell_thickness_mm);
    if let Some(ref f) = shell.tpms_fill {
        tracing::info!(
            "  tpms_fill  : {} λ={:.1}mm iso={} res={}",
            f.surface.label(),
            f.period_mm,
            f.iso_value,
            f.resolution
        );
    }
    tracing::info!(
        "  inlet port : outer ({:.2}, {:.3})  →  inner ({:.2}, {:.3})",
        shell.inlet_port.0,
        shell.inlet_port.1,
        shell.inlet_port.0 + shell.shell_thickness_mm,
        shell.inlet_port.1,
    );
    tracing::info!(
        "  outlet port: inner ({:.2}, {:.3})  →  outer ({:.2}, {:.3})",
        shell.outlet_port.0 - shell.shell_thickness_mm,
        shell.outlet_port.1,
        shell.outlet_port.0,
        shell.outlet_port.1,
    );
    tracing::info!("  outline segs: {}", shell.box_outline.len());

    // ── 2. Validate ──────────────────────────────────────────────────────────
    shell.validate().expect("shell cuboid invariants must hold");

    // ── 3. Render to SVG ─────────────────────────────────────────────────────
    let current_dir = std::env::current_dir().expect("Failed to get current directory");
    let workspace_root = current_dir
        .ancestors()
        .find(|p| p.join("Cargo.toml").exists() && p.join("crates").exists())
        .expect("Must be run from within the CFDrs workspace");

    // Output directory: CFDrs/crates/cfd-schematics/output/examples/shell_cuboid_demo/
    let output_dir = workspace_root
        .join("crates")
        .join("cfd-schematics")
        .join("output")
        .join("examples")
        .join("shell_cuboid_demo");

    fs::create_dir_all(&output_dir).expect("could not create outputs dir");
    let output_path = output_dir.join("shell_cuboid_demo.svg");
    let output_path_str = output_path.to_str().unwrap().to_string();

    let config = RenderConfig {
        title: "TPMS-Filled Shell Cuboid — 96-well (Gyroid λ=5mm)".to_string(),
        width: 1000,
        height: 700,
        ..RenderConfig::default()
    };

    let renderer = PlottersRenderer;
    renderer
        .render_shell_cuboid(&shell, &output_path_str, &config)
        .expect("rendering failed");

    // ── 4. Interchange JSON ──────────────────────────────────────────────────
    let json = shell
        .to_interchange_json()
        .expect("interchange serialisation failed");

    tracing::info!("\n── Interchange JSON ──────────────────────────────────────────");
    tracing::info!("{json}");

    let json_path = output_dir.join("shell_cuboid_demo.json");
    fs::write(&json_path, &json).expect("failed to write json");
    tracing::info!("Saved files to {}", output_dir.display());
}
