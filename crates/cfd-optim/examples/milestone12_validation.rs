//! Milestone 12 multi-fidelity venturi validation.
//!
//! This companion example consumes the selected Option 2 and GA designs written by
//! `milestone12_report` and performs the expensive 1D/2D/3D pressure-drop
//! confirmation separately.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_validation --no-default-features
//! ```

use cfd_optim::{load_top5_report_json, run_milestone12_validation};
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate has a parent")
        .parent()
        .expect("crates/ has a workspace root")
        .to_path_buf();
    let out_dir = workspace_root.join("report").join("milestone12");
    std::fs::create_dir_all(&out_dir)?;

    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let ga_ranked = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    let option2 = option2_ranked
        .first()
        .ok_or("two_concept_option2_venturi_top5.json is empty")?;
    let ga = ga_ranked.first().ok_or("ga_hydrosdt_top5.json is empty")?;

    run_milestone12_validation(&out_dir, option2, ga)?;
    Ok(())
}
