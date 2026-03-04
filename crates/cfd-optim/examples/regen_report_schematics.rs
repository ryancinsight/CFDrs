//! Regenerate Milestone 12 report schematics from current generated rankings.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example regen_report_schematics --no-default-features
//! ```

use cfd_optim::{save_schematic_svg, DesignCandidate};
use std::path::{Path, PathBuf};

fn load_top_candidate(
    path: &Path,
    label: &str,
) -> Result<DesignCandidate, Box<dyn std::error::Error>> {
    let data = std::fs::read_to_string(path)
        .map_err(|e| format!("failed reading {label} source {}: {e}", path.display()))?;
    let raw: serde_json::Value = serde_json::from_str(&data)
        .map_err(|e| format!("failed parsing {label} source {}: {e}", path.display()))?;
    let first = raw
        .as_array()
        .and_then(|a| a.first())
        .ok_or_else(|| format!("{label} source {} contains no candidates", path.display()))?;
    let candidate = first
        .get("candidate")
        .ok_or_else(|| format!("{label} source {} missing candidate field", path.display()))?;
    let decoded: DesignCandidate = serde_json::from_value(candidate.clone()).map_err(|e| {
        format!(
            "failed decoding candidate from {label} source {}: {e}",
            path.display()
        )
    })?;
    Ok(decoded)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate parent")
        .parent()
        .expect("workspace root")
        .to_path_buf();
    let figures_dir = workspace_root.join("report").join("figures");
    let milestone12_dir = workspace_root.join("report").join("milestone12");
    std::fs::create_dir_all(&figures_dir)?;

    let option1 = load_top_candidate(
        &milestone12_dir.join("two_concept_option1_ultrasound_top5.json"),
        "option1_ultrasound",
    )?;
    let option2 = load_top_candidate(
        &milestone12_dir.join("two_concept_option2_venturi_top5.json"),
        "option2_oncology_cif_cct",
    )?;
    let option2_rbc = load_top_candidate(
        &milestone12_dir.join("top5_rbc_protected.json"),
        "option2_rbc_protected",
    )?;

    save_schematic_svg(&option1, &figures_dir.join("selected_ga_schematic.svg"))?;
    save_schematic_svg(
        &option2,
        &figures_dir.join("selected_cifx_combined_schematic.svg"),
    )?;
    save_schematic_svg(
        &option2_rbc,
        &figures_dir.join("selected_cif_schematic.svg"),
    )?;

    println!("Generated report schematics:");
    println!("  Figure 4: {}", option1.id);
    println!("  Figure 5: {}", option2.id);
    println!("  Figure 6: {}", option2_rbc.id);

    Ok(())
}
