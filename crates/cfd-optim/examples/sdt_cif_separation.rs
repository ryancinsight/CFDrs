//! CIF-dedicated 6-part pipeline for optimal design determination.
//!
//! Executes the full CIF/CCT determination strategy reported in
//! Sections 8.18–8.23 of the Milestone 12 report:
//!
//! 1. **Parametric sweep** — exhaustive CIF/CCT exploration across three
//!    cell-separation objectives.
//! 2. **GA evolutionary search** — 100-individual × 120-generation runs for
//!    `CellSeparation`, `ThreePopSeparation`, and `CombinedSdtLeukapheresis`.
//! 3. **Pareto front isolation** — non-dominated front across cancer
//!    separation, WBC recovery, and lysis risk.
//! 4. **Venturi throat analysis** — per-design throat metrics summary.
//! 5. **SVG schematics** — 6 comparison charts + 5 annotated CIF diagrams.
//! 6. **CSV export** — per-channel hemolysis decomposition.
//!
//! Output: `crates/cfd-optim/outputs/cif_separation/`
//! (6 JSON, 6 SVG comparison charts, 5 CIF schematic SVGs, 1 CSV = 18 files)
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example sdt_cif_separation --no-default-features
//! ```

use cfd_optim::delivery::csv_export::save_per_channel_csv;
use cfd_optim::{
    save_annotated_cif_svg, save_comparison_svg, save_top5_json, DesignTopology, GeneticOptimizer,
    OptimMode, RankedDesign, SdtOptimizer, SdtWeights,
};
use std::path::Path;

/// Returns `true` for CIF and CCT topology families.
fn is_cif_or_cct(t: &DesignTopology) -> bool {
    matches!(
        t,
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
            | DesignTopology::CascadeCenterTrifurcationSeparator { .. }
    )
}

/// Extract (cancer_sep↑, wbc_recovery↑, lysis_risk↓) objectives.
fn objectives(d: &RankedDesign) -> (f64, f64, f64) {
    (
        d.metrics.cell_separation_efficiency,
        d.metrics.wbc_recovery,
        d.metrics.hemolysis_index_per_pass,
    )
}

/// Non-dominated Pareto front across three objectives.
fn pareto_front(pool: &[RankedDesign]) -> Vec<RankedDesign> {
    pool.iter()
        .filter(|d| {
            !pool.iter().any(|o| {
                let (ds, dw, dl) = objectives(d);
                let (os, ow, ol) = objectives(o);
                // `o` dominates `d` if better-or-equal in all, strictly better in ≥1
                os >= ds && ow >= dw && ol <= dl && (os > ds || ow > dw || ol < dl)
            })
        })
        .cloned()
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out = Path::new("crates/cfd-optim/outputs/cif_separation");
    std::fs::create_dir_all(out)?;
    let w = SdtWeights::default();

    println!("=== CIF Separation Pipeline ===\n");
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

    // ── Part 1: Parametric sweep ──────────────────────────────────────────
    println!("[1/6] Parametric sweep across 3 modes …");
    let mut param_cif: Vec<Vec<RankedDesign>> = Vec::new();

    for (slug, mode) in &modes {
        let all = SdtOptimizer::new(*mode, w).all_ranked()?;
        let cif: Vec<RankedDesign> = all
            .into_iter()
            .filter(|d| is_cif_or_cct(&d.candidate.topology))
            .collect();
        let top5: Vec<_> = cif.iter().take(5).cloned().collect();

        save_top5_json(&top5, &out.join(format!("parametric_{slug}_top5.json")))?;
        save_comparison_svg(
            &top5,
            &out.join(format!("parametric_{slug}_comparison.svg")),
            *mode,
        )?;
        println!(
            "      {slug}: {} CIF/CCT candidates, top-5 saved",
            cif.len()
        );
        param_cif.push(cif);
    }

    // ── Part 2: GA evolutionary search ────────────────────────────────────
    println!("\n[2/6] GA evolution (100 pop × 120 gen) …");
    let mut ga_top: Vec<Vec<RankedDesign>> = Vec::new();

    for (slug, mode) in &modes {
        let res = GeneticOptimizer::new(*mode, w)
            .with_population(100)
            .with_max_generations(120)
            .with_top_k(7)
            .with_rng_seed(42)
            .run()?;

        save_top5_json(&res.top_designs, &out.join(format!("ga_{slug}.json")))?;
        save_comparison_svg(
            &res.top_designs,
            &out.join(format!("ga_{slug}_comparison.svg")),
            *mode,
        )?;
        println!(
            "      ga_{slug}: best {:.4} → {:.4} ({} gen)",
            res.best_per_gen.first().copied().unwrap_or(0.0),
            res.best_per_gen.last().copied().unwrap_or(0.0),
            res.best_per_gen.len(),
        );
        ga_top.push(res.top_designs);
    }

    // ── Part 3: Pareto front isolation ────────────────────────────────────
    println!("\n[3/6] Pareto front isolation …");
    let pool: Vec<RankedDesign> = param_cif
        .iter()
        .flatten()
        .chain(ga_top.iter().flatten())
        .cloned()
        .collect();
    let front = pareto_front(&pool);
    println!("      {} non-dominated designs", front.len());
    println!(
        "      {:<40} {:>8} {:>8} {:>10}",
        "Design", "CancSep", "WBCrec", "LysisRisk"
    );
    for d in &front {
        let (s, wr, l) = objectives(d);
        println!(
            "      {:<40} {:>8.3} {:>8.3} {:>10.2e}",
            d.candidate.id, s, wr, l,
        );
    }

    // ── Part 4: Venturi throat analysis ───────────────────────────────────
    println!("\n[4/6] Venturi throat analysis …");
    let mut venturi_pool: Vec<&RankedDesign> = pool
        .iter()
        .filter(|d| d.candidate.topology.has_venturi())
        .collect();
    venturi_pool.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    venturi_pool.dedup_by(|a, b| a.candidate.id == b.candidate.id);

    println!(
        "      {:<40} {:>8} {:>10} {:>10} {:>10}",
        "Design", "σ", "τ_thr(Pa)", "CancCtr%", "HI/pass"
    );
    for d in venturi_pool.iter().take(6) {
        let m = &d.metrics;
        println!(
            "      {:<40} {:>8.3} {:>10.0} {:>10.1} {:>10.2e}",
            d.candidate.id,
            m.cavitation_number,
            m.throat_shear_pa,
            m.cancer_center_fraction * 100.0,
            m.hemolysis_index_per_pass,
        );
    }

    // ── Part 5: CIF schematic SVGs ────────────────────────────────────────
    println!("\n[5/6] Annotated CIF schematics …");
    let mut best: Vec<RankedDesign> = pool;
    best.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    best.dedup_by(|a, b| a.candidate.id == b.candidate.id);

    for (i, d) in best.iter().take(5).enumerate() {
        let path = out.join(format!("cif_schematic_{}.svg", i + 1));
        save_annotated_cif_svg(d, &path)?;
        println!("      #{}: {}", i + 1, d.candidate.id);
    }

    // ── Part 6: CSV export ────────────────────────────────────────────────
    println!("\n[6/6] Per-channel hemolysis CSV …");
    let csv_set: Vec<_> = best.into_iter().take(5).collect();
    save_per_channel_csv(&csv_set, &out.join("per_channel_hemolysis.csv"))?;
    println!("      CSV: {} designs", csv_set.len());

    println!("\n=== Complete: {} ===", out.display());
    Ok(())
}
