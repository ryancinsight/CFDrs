//! CIF mesh generation pipeline for top-ranked CIF/CCT designs.
//!
//! Processes GA-evolved CIF and CCT designs through the [`DesignPipeline`]
//! to produce watertight surface meshes (STL) and volume meshes (OpenFOAM
//! polyMesh) for downstream 3D CFD validation and SLA fabrication
//! (Section 8.24 of the Milestone 12 report).
//!
//! Submits the top-2 designs from each of three GA modes (6 total);
//! expects 5 successful mesh constructions (one CSG failure is a known
//! limitation of the current Boolean engine for topologies with closely
//! spaced branch junctions).
//!
//! Output: `crates/cfd-optim/outputs/cif_meshes/`
//! (5 design directories with STL + OpenFOAM polyMesh + schematic SVG)
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example sdt_cif_mesh --features mesh-export --no-default-features
//! ```

use cfd_optim::{DesignPipeline, GeneticOptimizer, OptimMode, RankedDesign, SdtWeights};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = Path::new("crates/cfd-optim/outputs/cif_meshes");
    let pipeline = DesignPipeline::new(out);
    let w = SdtWeights::default();

    println!("=== CIF Mesh Generation Pipeline ===\n");
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

    // Collect top-2 designs per GA mode (6 designs total).
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

    // Process each design through the mesh pipeline.
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
