//! Primitive selective mesh generation pipeline for top-ranked primitive split-sequence designs.
//!
//! Processes GA-evolved primitive selective designs through the [`DesignPipeline`]
//! to produce watertight surface meshes (STL) and volume meshes (OpenFOAM
//! polyMesh) for downstream 3D CFD validation and fabrication.
//!
//! Output: `crates/cfd-optim/outputs/primitive_selective_meshes/`
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example sdt_primitive_selective_mesh --features mesh-export --no-default-features
//! ```

use cfd_optim::{DesignPipeline, GeneticOptimizer, OptimMode, RankedDesign, SdtWeights};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = Path::new("crates/cfd-optim/outputs/primitive_selective_meshes");
    let pipeline = DesignPipeline::new(out);
    let w = SdtWeights::default();

    println!("=== Primitive Selective Mesh Generation Pipeline ===\n");
    println!("Output: {}\n", out.display());

    let modes: [(&str, OptimMode); 3] = [
        ("cell_sep", OptimMode::CellSeparation),
        ("three_pop", OptimMode::ThreePopSeparation),
        (
            "sdt_leuka",
            OptimMode::CombinedSdtLeukapheresis {
                leuka_weight: 0.5,
                sdt_weight: 0.5,
                patient_weight_kg: 3.0,
            },
        ),
    ];

    let mut submissions: Vec<(&str, RankedDesign)> = Vec::new();

    for (slug, mode) in &modes {
        println!("[GA] {slug}: running 100 pop × 120 gen …");
        let res = GeneticOptimizer::new(*mode, w)
            .with_population(100)
            .with_max_generations(120)
            .with_top_k(2)
            .with_rng_seed(42)
            .run()?;

        for (i, d) in res.top_designs.into_iter().enumerate() {
            println!(
                "      rank {}: {}  ({})  score={:.4}",
                i + 1,
                d.candidate.id,
                d.candidate.topology.short(),
                d.score,
            );
            submissions.push((slug, d));
        }
    }

    println!(
        "\n=== Mesh Generation ({} designs) ===\n",
        submissions.len()
    );
    println!(
        "{:<12} {:<40} {:>8} {:>8} {:>10}",
        "Mode", "Design", "Verts", "Faces", "Watertight"
    );
    println!("{}", "-".repeat(82));

    let mut successes = 0_usize;
    let mut failures = 0_usize;

    for (slug, design) in &submissions {
        match pipeline.export_design(design) {
            Ok(a) => {
                successes += 1;
                println!(
                    "{:<12} {:<40} {:>8} {:>8} {:>10}",
                    slug,
                    a.design_id,
                    a.vertex_count,
                    a.face_count,
                    if a.watertight { "YES" } else { "NO" },
                );
            }
            Err(e) => {
                failures += 1;
                println!(
                    "{:<12} {:<40} {:>8} {:>8} {:>10}",
                    slug, design.candidate.id, "—", "—", "FAIL",
                );
                eprintln!("      Error: {e}");
            }
        }
    }

    println!(
        "\n{successes} meshes generated, {failures} failed  →  {}\n",
        out.display()
    );
    Ok(())
}
