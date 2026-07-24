//! check-figures — SSOT drift verification for docs/book/{SUMMARY,README}.md.
//!
//! Compares figure links parsed from SUMMARY.md + README.md against the
//! canonical `FIGURE_SPECS` list declared in `prebook.rs`.  Drift in either
//! direction is reported; CI-friendly exit codes: 0 = in sync, 1 = drift.

use anyhow::{anyhow, Context, Result};
use serde::Serialize;
use std::collections::BTreeSet;
use std::fs;
use std::path::Path;

/// One figure link parsed from SUMMARY.md or README.md.
#[derive(Debug, Clone, Serialize)]
pub struct DocsFigureRef {
    /// "SUMMARY.md" or "README.md" (or any other docs source file passed in).
    pub source_file: String,
    /// 1-indexed line number in the source file.
    pub line_no: usize,
    /// Filename portion (e.g. "poiseuille_parabolic_profile.svg"), without
    /// the "figures/" prefix.
    pub figure_file: String,
}

/// Cross-check report — drift is any inconsistency between the two
/// sets (docs-referenced figures vs FIGURE_SPECS-declared figures).
#[derive(Debug, Clone)]
pub struct CheckFiguresReport {
    /// Unique figure filenames referenced in any docs source, sorted.
    pub docs_figures: BTreeSet<String>,
    /// Unique figure filenames declared in `FIGURE_SPECS`, sorted.
    pub specs_figures: BTreeSet<String>,
    /// Drifted docs references: figure files used in a docs source that
    /// are NOT in `FIGURE_SPECS`.
    pub orphan_docs_refs: Vec<DocsFigureRef>,
    /// Drifted specs: figure files in `FIGURE_SPECS` that are NOT
    /// referenced in any docs source.
    pub orphan_specs: Vec<String>,
}

/// Parse a markdown file for figure link references of the form
/// `[Title](figures/<name>.svg)`.  Pure byte-scan — no regex dep —
/// which lets us avoid pulling `regex` into the `xtask` workspace graph.
fn parse_figure_refs(path: &Path) -> Result<(String, Vec<DocsFigureRef>)> {
    let source_file = path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| anyhow!("path has no filename: {}", path.display()))?
        .to_owned();
    let content = fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let mut refs: Vec<DocsFigureRef> = Vec::new();

    for (lineno, line) in content.lines().enumerate() {
        // For each line, scan for occurrences of `figures/` followed by a
        // path-segment terminator (close-paren or whitespace).  This is
        // lexically precise enough for the small set of valid figure
        // references used by this book.
        let bytes = line.as_bytes();
        let mut idx: usize = 0;
        while let Some(pos) = line[idx..].find("figures/") {
            let abs_start = idx + pos;
            let mut end = abs_start + "figures/".len();
            while end < bytes.len() && bytes[end] != b')' && !bytes[end].is_ascii_whitespace() {
                end += 1;
            }
            let candidate = &line[abs_start..end];
            if let Some(stripped) = candidate.strip_prefix("figures/") {
                if stripped.ends_with(".svg") && !stripped.contains('/') {
                    refs.push(DocsFigureRef {
                        source_file: source_file.clone(),
                        line_no: lineno + 1,
                        figure_file: stripped.to_owned(),
                    });
                }
            }
            // Advance past this occurrence to look for more on the same line.
            idx = end.max(abs_start + 1);
        }
    }

    Ok((source_file, refs))
}

/// Run `check-figures` against `workspace_root`.  Verifies that every
/// figure link in `docs/book/SUMMARY.md` and `docs/book/README.md` is
/// listed in `super::prebook::FIGURE_SPECS`, and reports drift in
/// either direction (docs referencing unknown figures, or
/// `FIGURE_SPECS` declaring figures never referenced).
pub fn run_check_figures(workspace_root: &Path) -> Result<CheckFiguresReport> {
    let docs_dir = workspace_root.join("docs/book");
    let summary_path = docs_dir.join("SUMMARY.md");
    let readme_path = docs_dir.join("README.md");

    let (_, summary_refs) = parse_figure_refs(&summary_path)
        .with_context(|| format!("parsing {}", summary_path.display()))?;
    let (_, readme_refs) = parse_figure_refs(&readme_path)
        .with_context(|| format!("parsing {}", readme_path.display()))?;

    let all_docs_refs: Vec<DocsFigureRef> = summary_refs
        .into_iter()
        .chain(readme_refs)
        .collect();

    let mut docs_figures: BTreeSet<String> = BTreeSet::new();
    for r in &all_docs_refs {
        docs_figures.insert(r.figure_file.clone());
    }

    let specs_figures: BTreeSet<String> = super::prebook::FIGURE_SPECS
        .iter()
        .map(|s| s.name.to_owned())
        .collect();

    let orphan_docs_refs: Vec<DocsFigureRef> = all_docs_refs
        .into_iter()
        .filter(|r| !specs_figures.contains(&r.figure_file))
        .collect();

    let orphan_specs: Vec<String> = specs_figures.difference(&docs_figures).cloned().collect();

    Ok(CheckFiguresReport {
        docs_figures,
        specs_figures,
        orphan_docs_refs,
        orphan_specs,
    })
}

/// Print the report to stdout.  Returns the desired process exit code:
/// 0 = in sync, 1 = drift detected.  The caller (`main.rs`) is expected
/// to map a `1` to `std::process::exit(1)` so CI sees the failure.
pub fn print_check_figures_report(report: &CheckFiguresReport) -> i32 {
    println!(
        "check-figures: SSOT drift verification (SUMMARY.md + README.md vs FIGURE_SPECS)"
    );
    println!();
    println!("  FIGURE_SPECS entries  : {}", report.specs_figures.len());
    println!("  Docs references       : {}", report.docs_figures.len());
    println!(
        "  In sync (intersection): {}",
        report
            .specs_figures
            .intersection(&report.docs_figures)
            .count()
    );
    println!();

    if report.orphan_docs_refs.is_empty() && report.orphan_specs.is_empty() {
        println!("SSOT_IN_SYNC: every docs figure link is listed in FIGURE_SPECS.");
        return 0;
    }

    let exit_code: i32 = 1;

    if !report.orphan_docs_refs.is_empty() {
        println!(
            "DRIFT_DOCS_NOT_IN_SPECS: {} docs figure link(s) missing from FIGURE_SPECS:",
            report.orphan_docs_refs.len()
        );
        for r in &report.orphan_docs_refs {
            println!("  - {}:L{}  ({})", r.source_file, r.line_no, r.figure_file);
        }
        println!();
    }

    if !report.orphan_specs.is_empty() {
        println!(
            "DRIFT_SPECS_NOT_IN_DOCS: {} FIGURE_SPECS entry(ies) not referenced in any docs file:",
            report.orphan_specs.len()
        );
        for name in &report.orphan_specs {
            println!("  - {name}");
        }
        println!();
    }

    exit_code
}
