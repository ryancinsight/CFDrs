//! Formatted console output and file export for the leukapheresis example.
//!
//! Provides table-printing helpers for the three pipeline phases and a
//! convenience wrapper that saves JSON + SVG + per-design schematics.

use std::path::Path;

use cfd_optim::{save_comparison_svg, save_schematic_svg, save_top5_json, OptimMode, RankedDesign};

// ── Phase 1 — Benchmark comparison ───────────────────────────────────────────

/// Print the Phase 1 comparison table header.
pub fn print_benchmark_header() {
    println!(
        "  {:<32}  {:>9}  {:>9}  |  {:>9}  {:>9}",
        "Design", "WBC-rec", "RBC-rem", "Ref-WBC", "Ref-RBC"
    );
    println!(
        "  {:<32}  {:>9}  {:>9}  |  {:>9}  {:>9}",
        "", "Model", "Model", "Paper", "Paper"
    );
    println!("  {}", "-".repeat(76));
}

/// Print one row comparing model predictions to paper-reported values.
///
/// `model_rbc_rem` = 1 − rbc_pass_fraction (fraction removed to periphery).
pub fn print_benchmark_comparison(
    design_id: &str,
    model_wbc_rec: f64,
    model_rbc_rem: f64,
    paper_wbc_rec: f64,
    paper_rbc_rem: f64,
) {
    println!(
        "  {:<32}  {:>8.1}%  {:>8.1}%  |  {:>8.1}%  {:>8.1}%",
        design_id,
        model_wbc_rec * 100.0,
        model_rbc_rem * 100.0,
        paper_wbc_rec * 100.0,
        paper_rbc_rem * 100.0,
    );
}

// ── Phase 3 — GA results table ────────────────────────────────────────────────

/// Print the leukapheresis GA results table header.
pub fn print_leuka_table_header() {
    println!(
        "  {:>4}  {:<38}  {:>7}  {:>8}  {:>8}  {:>8}  {:>9}",
        "#", "Candidate ID", "Score", "WBC-rec%", "RBC-rem%", "Purity%", "ECV mL"
    );
    println!("  {}", "-".repeat(95));
}

/// Print one leukapheresis result row.
pub fn print_leuka_row(d: &RankedDesign) {
    let m = &d.metrics;
    println!(
        "  {:>4}  {:<38}  {:>7.4}  {:>7.1}%  {:>7.1}%  {:>7.1}%  {:>8.2}",
        d.rank,
        truncate(&d.candidate.id, 38),
        d.score,
        m.wbc_recovery * 100.0,
        (1.0 - m.rbc_pass_fraction) * 100.0,
        m.wbc_purity * 100.0,
        m.total_ecv_ml,
    );
}

// ── File export ───────────────────────────────────────────────────────────────

/// Save JSON result file, comparison SVG, and per-design schematic SVGs.
///
/// Files are named `{tag}_top5.json`, `{tag}_comparison.svg`, and
/// `{tag}_rank{N}_{id}.svg` respectively.
pub fn save_leuka_results(
    designs: &[RankedDesign],
    tag: &str,
    out_dir: &Path,
    patient_weight_kg: f64,
) {
    let mode = OptimMode::PediatricLeukapheresis { patient_weight_kg };

    // JSON
    let json_path = out_dir.join(format!("{tag}_top5.json"));
    match save_top5_json(designs, &json_path) {
        Ok(_) => println!("  Saved: {}", json_path.display()),
        Err(e) => eprintln!("  WARN (JSON): {e}"),
    }

    // Comparison SVG (all designs side-by-side bar chart)
    if !designs.is_empty() {
        let svg_path = out_dir.join(format!("{tag}_comparison.svg"));
        match save_comparison_svg(designs, &svg_path, mode) {
            Ok(_) => println!("  Saved: {}", svg_path.display()),
            Err(e) => eprintln!("  WARN (comparison SVG): {e}"),
        }
    }

    // Per-design channel schematics
    for d in designs {
        let sch = out_dir.join(format!(
            "{}_rank{:02}_{}.svg",
            tag,
            d.rank,
            safe_id(&d.candidate.id)
        ));
        match save_schematic_svg(&d.candidate, &sch) {
            Ok(_) => println!("  Saved: {}", sch.display()),
            Err(e) => eprintln!("  WARN (schematic {}): {e}", d.candidate.id),
        }
    }
}

// ── Utilities ─────────────────────────────────────────────────────────────────

fn truncate(s: &str, n: usize) -> &str {
    if s.len() <= n {
        s
    } else {
        &s[..n]
    }
}

fn safe_id(id: &str) -> String {
    id.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
