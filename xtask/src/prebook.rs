//! Prebook / figure-management for the CFDrs mdbook.
//!
//! Build the deterministic figure set committed at `docs/book/figures/`
//! and emit a `MANIFEST.json` manifest that downstream tooling
//! (link-checker pre-flight, CI evidence chain) can fingerprint.
//!
//! # Contract
//!
//! - `FIGURE_SPECS` is the **single source of truth** (SSOT) for the
//!   figure set.  `CFDrs/docs/book/SUMMARY.md` and
//!   `CFDrs/docs/book/README.md` stay manually consistent with this
//!   list.
//! - Each figure is `HandAuthored`: the `.svg` is committed by hand to
//!   mirror the deterministic output of a specific example.  Both the
//!   SVG and the example's printed values are deterministic, so two
//!   `prebook` invocations on unchanged inputs produce byte-identical
//!   output.
//! - The manifest hash is the first 16 hex chars of SHA-256 over the
//!   figure file content; 64 bits is plenty for a book-scale figure
//!   set (collision probability for 7 entries is roughly 2.7 × 10⁻¹⁵).
//! - `crate_name` mirrors the `**Crate**:` header on the linked
//!   example `.md` page; the canonical run form is
//!   `cargo run -p <crate_name> --example <example_name>` for crate
//!   members, or `cargo run --example <example_name>` (workspace
//!   root) when `crate_name == "cfd-suite"`.

use anyhow::{anyhow, Context, Result};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};

/// A figure referenced by the CFDrs mdbook.
///
/// `source_example` is `Some(_)` for figures that mirror a runnable
/// example's deterministic stdout / staged PNG; `None` for purely
/// authored schematics (workflow / architecture diagrams).
#[derive(Debug, Clone, Copy, Serialize)]
pub struct FigureSpec {
    /// File name (relative to `docs/book/figures/`).
    pub name: &'static str,
    /// Optional runnable example this figure mirrors.
    pub source_example: Option<ExampleRef>,
    /// One-line summary used in the README figure index.
    pub summary: &'static str,
}

/// A runnable example reference:
/// `cargo run [-p <crate_name>] --example <example_name>`.
///
/// `crate_name == "cfd-suite"` (workspace root) means the example is
/// declared at the root `Cargo.toml` [[example]] table and is run
/// without `-p`; all other values map 1:1 to a workspace member crate.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct ExampleRef {
    pub crate_name: &'static str,
    pub example_name: &'static str,
}

/// Authoritative figure list — kept synchronised with `SUMMARY.md` and
/// `README.md` figure indexes by hand.
///
/// `crate_name` values were verified against the `**Crate**:` header
/// of each linked example `.md` page (post BOOK-DETERMINISTIC-FIGURES-1
/// validation pass).  All entries resolve to either:
///
/// - `cfd-suite` — workspace-root example (run as
///   `cargo run --example <name>`), or
/// - `cfd-validation` — `cfd-validation` member crate (run as
///   `cargo run -p cfd-validation --example <name>`).
pub const FIGURE_SPECS: &[FigureSpec] = &[
    FigureSpec {
        name: "poiseuille_parabolic_profile.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-suite",
            example_name: "pipe_flow_validation",
        }),
        summary: "Plane Poiseuille parabolic profile u(y) = u_max (1 - (2y/H)^2), the 1-D analytical solution used as the laminar benchmark.",
    },
    FigureSpec {
        name: "cavity_streamfunction_contour.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-suite",
            example_name: "cavity_validation",
        }),
        summary: "Lid-driven cavity primary-vortex streamfunction contours (Re ~ 100, square domain).",
    },
    FigureSpec {
        name: "residual_convergence_semilog.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-suite",
            example_name: "turbulent_channel_flow",
        }),
        summary: "Semi-log L2 residual vs iteration, SIMPLE vs PISO (two corrector steps).",
    },
    FigureSpec {
        name: "channel_mesh_layout.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-suite",
            example_name: "mesh_3d_integration",
        }),
        summary: "16 x 4 structured rectangular mesh layout with dy / dx labelled; boundary cells highlighted.",
    },
    FigureSpec {
        name: "reynolds_regime_map.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-suite",
            example_name: "pipe_flow_validation",
        }),
        summary: "Pipe-flow Reynolds regime map: laminar (Re < 2300), transitional (2300 - 4000), turbulent (>= 4000).",
    },
    FigureSpec {
        name: "richardson_loglog.svg",
        source_example: Some(ExampleRef {
            crate_name: "cfd-validation",
            example_name: "richardson_convergence",
        }),
        summary: "Richardson extrapolation: log L2 error vs log grid spacing, with 2nd-order reference slope.",
    },
    FigureSpec {
        name: "architecture_stack.svg",
        source_example: None,
        summary: "CFDrs layered architecture on the Atlas stack: validation / optim / python top, dimensional schemes mid, core / math / io core, Atlas leto / hephaestus / coeus substrate below.",
    },
];

/// Single manifest entry — produced by [`run_prebook`] and serialised as
/// JSON to `docs/book/figures/MANIFEST.json`. Byte-deterministic across
/// runs of identical inputs (no timestamps, no ordering dependence).
#[derive(Debug, Clone, Serialize)]
pub struct ManifestEntry {
    pub name: String,
    pub source_example: Option<ExampleRef>,
    pub sha256_16: String,
    pub bytes: usize,
    pub summary: String,
}

/// Aggregated report returned to the CLI surface.
#[derive(Debug, Clone)]
pub struct PrebookReport {
    pub entries: Vec<ManifestEntry>,
    pub manifest_path: PathBuf,
}

/// Run prebook against `workspace_root`. For each spec, verifies the
/// figure file exists at `<workspace_root>/docs/book/figures/<name>`,
/// hashes it, and writes a deterministic `MANIFEST.json` next to the
/// figures. The returned report lists every verified entry.
pub fn run_prebook(workspace_root: &Path) -> Result<PrebookReport> {
    let figs_dir = workspace_root.join("docs/book/figures");
    if !figs_dir.is_dir() {
        return Err(anyhow!(
            "figures directory not found: {} (prebook expects <workspace>/docs/book/figures/)",
            figs_dir.display()
        ));
    }

    // Iterate specs in declaration order — stable, JSON-ordered manifest.
    let mut entries: Vec<ManifestEntry> = Vec::with_capacity(FIGURE_SPECS.len());
    for spec in FIGURE_SPECS {
        let path = figs_dir.join(spec.name);
        let bytes = fs::read(&path)
            .with_context(|| format!("reading figure file {}", path.display()))?;
        let sha = sha256_hex_first_16(&bytes);
        entries.push(ManifestEntry {
            name: spec.name.to_owned(),
            source_example: spec.source_example,
            sha256_16: sha,
            bytes: bytes.len(),
            summary: spec.summary.to_owned(),
        });
    }

    // Serialise with sorted keys (serde_json default) + no pretty-print
    // trailing whitespace; deterministic across runs and machines.
    let manifest_path = figs_dir.join("MANIFEST.json");
    let json = serde_json::to_string(&entries)
        .context("serialising MANIFEST.json entries")?;
    fs::write(&manifest_path, format!("{json}\n"))
        .with_context(|| format!("writing manifest {}", manifest_path.display()))?;

    Ok(PrebookReport {
        entries,
        manifest_path,
    })
}

/// SHA-256 hex digest; only the first 16 hex chars are surfaced in the
/// manifest (16 hex chars at ~6.5·10^7 distinct values is plenty for a
/// book-scale figure set and the file content is small + committed, so
/// anyone needing the full hash can recompute).
pub fn sha256_hex_first_16(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let hex = format!("{:x}", digest);
    hex[..16].to_owned()
}
