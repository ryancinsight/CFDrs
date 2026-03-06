//! Consolidated diagnostic audit for `cfd-schematics` + `cfd-optim`.
//!
//! | Part | Content |
//! |------|---------|
//! | 1 | Preset catalog — all `cfd-schematics` presets with channel counts |
//! | 2 | Leukapheresis infeasibility audit — CE/SP/PM topology diagnosis |
//! | 3 | Per-topology/mode feasibility summary table |
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example audit --no-default-features
//! ```

use cfd_optim::{
    build_candidate_space, compute_metrics, score_candidate, DesignCandidate, DesignTopology,
    OptimMode, SdtMetrics, SdtWeights,
};
use cfd_schematics::interface::presets::{
    asymmetric_bifurcation_serpentine_rect, bifurcation_rect, bifurcation_serpentine_rect,
    bifurcation_trifurcation_venturi_rect, bifurcation_venturi_rect,
    cascade_center_trifurcation_rect, cell_separation_rect, constriction_expansion_array_rect,
    double_bifurcation_serpentine_rect, double_bifurcation_venturi_rect,
    double_trifurcation_venturi_rect, incremental_filtration_tri_bi_rect,
    parallel_microchannel_array_rect, quad_trifurcation_venturi_rect, serial_double_venturi_rect,
    serpentine_chain, serpentine_rect, spiral_channel_rect, symmetric_bifurcation,
    symmetric_trifurcation, trifurcation_bifurcation_bifurcation_venturi_rect,
    trifurcation_bifurcation_venturi_rect, trifurcation_rect, trifurcation_serpentine_rect,
    trifurcation_venturi_rect, triple_bifurcation_venturi_rect, triple_trifurcation_venturi_rect,
    venturi_chain, venturi_rect, venturi_serpentine_rect,
};
use cfd_schematics::NetworkBlueprint;
use serde::Serialize;
use std::fmt::Write as FmtWrite;
use std::path::Path;

// ── Preset catalog types ──────────────────────────────────────────────────────

struct PresetDef {
    id: &'static str,
    family: &'static str,
    support: &'static str,
    mapped_to: &'static str,
    build: fn(&str) -> NetworkBlueprint,
}

#[derive(Debug, Clone, Serialize)]
struct PresetRow {
    preset_id: String,
    family: String,
    support: String,
    mapped_to: String,
    nodes: usize,
    channels: usize,
    inlets: usize,
    outlets: usize,
    venturi_channels: usize,
    serpentine_channels: usize,
    total_length_mm: f64,
}

// ── Leukapheresis infeasibility types ────────────────────────────────────────

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
enum PrimaryReason {
    MetricsError,
    PressureFail,
    FdaFail,
    PlateFail,
    EcvZero,
    WbcZero,
    ScoreOther,
    Feasible,
}

#[derive(Debug, Default, Clone)]
struct ModeCounters {
    positive_score: usize,
    metrics_error: usize,
    pressure_fail: usize,
    fda_fail: usize,
    plate_fail: usize,
    ecv_zero: usize,
    wbc_zero: usize,
    score_other: usize,
}

impl ModeCounters {
    fn record(&mut self, r: PrimaryReason) {
        match r {
            PrimaryReason::Feasible => self.positive_score += 1,
            PrimaryReason::MetricsError => self.metrics_error += 1,
            PrimaryReason::PressureFail => self.pressure_fail += 1,
            PrimaryReason::FdaFail => self.fda_fail += 1,
            PrimaryReason::PlateFail => self.plate_fail += 1,
            PrimaryReason::EcvZero => self.ecv_zero += 1,
            PrimaryReason::WbcZero => self.wbc_zero += 1,
            PrimaryReason::ScoreOther => self.score_other += 1,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct TopologyAudit {
    short: String,
    name: String,
    total: usize,
    pediatric: ModeCounters,
    combined: ModeCounters,
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("sdt_audit");
    std::fs::create_dir_all(&out_dir).expect("create output dir");
    println!("=== SDT Audit (3 parts) ===\n");
    println!("Outputs: {}\n", out_dir.display());

    // ── Part 1: Preset catalog ────────────────────────────────────────────────
    println!("{}", "=".repeat(100));
    println!("  PART 1 — Preset Catalog  (all cfd-schematics presets)");
    println!("{}", "=".repeat(100));

    let presets = preset_defs();
    let catalog: Vec<PresetRow> = presets.iter().map(inspect_preset).collect();

    println!(
        "  {:>30}  {:>16}  {:>10}  {:>5}  {:>5}  {:>5}  {:>5}  {:>6}",
        "Preset ID", "Support", "Family", "Nodes", "Chan", "Vent", "Serp", "Len mm"
    );
    println!("  {}", "-".repeat(100));
    for r in &catalog {
        println!(
            "  {:>30}  {:>16}  {:>10}  {:>5}  {:>5}  {:>5}  {:>5}  {:>6.0}",
            r.preset_id,
            r.support,
            &r.family[..r.family.len().min(10)],
            r.nodes,
            r.channels,
            r.venturi_channels,
            r.serpentine_channels,
            r.total_length_mm,
        );
    }
    println!("\n  Total presets: {}", catalog.len());

    let catalog_csv = out_dir.join("preset_catalog.csv");
    let mut csv = String::from("preset_id,family,support,mapped_to,nodes,channels,inlets,outlets,venturi,serpentine,length_mm\n");
    for r in &catalog {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{},{:.3}",
            r.preset_id,
            r.family,
            r.support,
            r.mapped_to,
            r.nodes,
            r.channels,
            r.inlets,
            r.outlets,
            r.venturi_channels,
            r.serpentine_channels,
            r.total_length_mm,
        );
    }
    write_file(&catalog_csv, &csv, "preset_catalog.csv");

    // ── Part 2: Leukapheresis infeasibility audit ─────────────────────────────
    println!("\n{}", "=".repeat(100));
    println!("  PART 2 — Leukapheresis Infeasibility Audit  (CE/SP/PM topologies)");
    println!("{}", "=".repeat(100));
    println!("  Classification order: metrics error → pressure → FDA → plate → ECV → WBC → other");
    println!("{}", "-".repeat(100));

    let pediatric_mode = OptimMode::PediatricLeukapheresis {
        patient_weight_kg: 3.0,
    };
    let combined_mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let weights = SdtWeights::default();

    let mut audits: Vec<TopologyAudit> = vec![
        TopologyAudit {
            short: "CE".into(),
            name: "Constriction-Expansion Array".into(),
            ..Default::default()
        },
        TopologyAudit {
            short: "SP".into(),
            name: "Spiral Serpentine".into(),
            ..Default::default()
        },
        TopologyAudit {
            short: "PM".into(),
            name: "Parallel Microchannel Array".into(),
            ..Default::default()
        },
    ];

    let mut all_candidates = build_candidate_space();
    println!("  Total parametric candidates: {}", all_candidates.len());
    all_candidates.retain(|c| {
        matches!(
            c.topology,
            DesignTopology::ConstrictionExpansionArray { .. }
                | DesignTopology::SpiralSerpentine { .. }
                | DesignTopology::ParallelMicrochannelArray { .. }
        )
    });
    println!("  CE/SP/PM subset: {} candidates", all_candidates.len());

    for c in &all_candidates {
        let slot = audits
            .iter_mut()
            .find(|a| a.short == c.topology.short())
            .expect("slot exists");
        slot.total += 1;
        slot.pediatric
            .record(classify_zero(c, pediatric_mode, &weights));
        slot.combined
            .record(classify_zero(c, combined_mode, &weights));
    }

    println!(
        "\n  {:20}  {:>6}  {:>12}  {:>8}  {:>8}  {:>8}",
        "Topology", "Total", "Mode", "Pos", "PFail", "FdaFail"
    );
    println!("  {}", "-".repeat(75));
    for a in &audits {
        for (mlabel, cnt) in [
            ("Pediatric(3kg)", &a.pediatric),
            ("CombinedLeuka", &a.combined),
        ] {
            println!(
                "  {:<20}  {:>6}  {:>12}  {:>8}  {:>8}  {:>8}",
                a.name, a.total, mlabel, cnt.positive_score, cnt.pressure_fail, cnt.fda_fail,
            );
        }
    }

    let leuka_csv = out_dir.join("leukapheresis_infeasibility.csv");
    let mut csv = String::from("topology,total,mode,pos_score,metrics_err,pressure_fail,fda_fail,plate_fail,ecv_zero,wbc_zero,score_other\n");
    for a in &audits {
        for (mkey, c) in [
            ("pediatric_3kg", &a.pediatric),
            ("combined_sdt_leuka", &a.combined),
        ] {
            let _ = writeln!(
                csv,
                "{},{},{},{},{},{},{},{},{},{},{}",
                a.short,
                a.total,
                mkey,
                c.positive_score,
                c.metrics_error,
                c.pressure_fail,
                c.fda_fail,
                c.plate_fail,
                c.ecv_zero,
                c.wbc_zero,
                c.score_other,
            );
        }
    }
    write_file(&leuka_csv, &csv, "leukapheresis_infeasibility.csv");

    let leuka_md = out_dir.join("leukapheresis_infeasibility.md");
    let mut md = String::new();
    let _ = writeln!(md, "# Leukapheresis Infeasibility Audit\n");
    let _ = writeln!(md, "| Topology | Total | Mode | Score>0 | PressureFail | FdaFail | PlateFail | ECV=0 | WBC=0 | Other |");
    let _ = writeln!(md, "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|");
    for a in &audits {
        for (ml, c) in [
            ("Pediatric (3 kg)", &a.pediatric),
            ("Combined SDT+Leuka", &a.combined),
        ] {
            let _ = writeln!(
                md,
                "| {} ({}) | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                a.name,
                a.short,
                a.total,
                ml,
                c.positive_score,
                c.pressure_fail,
                c.fda_fail,
                c.plate_fail,
                c.ecv_zero,
                c.wbc_zero,
                c.score_other,
            );
        }
    }
    write_file(&leuka_md, &md, "leukapheresis_infeasibility.md");

    // ── Part 3: Per-topology/mode feasibility summary ─────────────────────────
    println!("\n{}", "=".repeat(100));
    println!("  PART 3 — Per-Topology/Mode Feasibility Summary");
    println!("{}", "=".repeat(100));
    println!("  Modes: SdtCavitation, ThreePopSep, HydroSDT, RbcProtected, CombinedLeuka");
    println!("{}", "-".repeat(100));

    let modes: &[(&str, OptimMode)] = &[
        ("SdtCav", OptimMode::SdtCavitation),
        ("ThreePop", OptimMode::ThreePopSeparation),
        ("HydroSDT", OptimMode::HydrodynamicCavitationSDT),
        ("RbcProt", OptimMode::RbcProtectedSdt),
        ("Combined", combined_mode),
    ];

    let full_space = build_candidate_space();
    println!("  Full parametric space: {} candidates\n", full_space.len());

    // Pre-compute metrics for all candidates once.
    let evaluated: Vec<(DesignCandidate, Option<SdtMetrics>)> = full_space
        .into_iter()
        .map(|c| {
            let m = compute_metrics(&c).ok();
            (c, m)
        })
        .collect();

    println!(
        "  {:>8}  {:>28}  {:>8}  {:>8}  {:>8}  {:>8}",
        "Mode", "Topology", "Total", "Eval'd", "Feasible", "BestScr"
    );
    println!("  {}", "-".repeat(80));

    let mut summary_csv =
        String::from("mode,topology,total,evaluated,feasible,feasible_pct,best_score\n");

    for (mode_key, mode) in modes {
        // Group by topology short name.
        let mut by_topo: std::collections::BTreeMap<String, (usize, usize, usize, f64)> =
            std::collections::BTreeMap::new();

        for (c, m_opt) in &evaluated {
            let entry = by_topo
                .entry(c.topology.short().to_string())
                .or_insert((0, 0, 0, 0.0f64));
            entry.0 += 1; // total
            if let Some(m) = m_opt {
                entry.1 += 1; // evaluated
                let score = score_candidate(m, *mode, &weights);
                if score > 0.0 {
                    entry.2 += 1; // feasible
                    if score > entry.3 {
                        entry.3 = score;
                    }
                }
            }
        }

        for (topo_short, (total, evaluated_n, feasible, best_score)) in &by_topo {
            let feasible_pct = if *evaluated_n > 0 {
                *feasible as f64 / *evaluated_n as f64 * 100.0
            } else {
                0.0
            };
            println!(
                "  {:>8}  {:>28}  {:>8}  {:>8}  {:>6} ({:>4.0}%)  {:>8.4}",
                mode_key, topo_short, total, evaluated_n, feasible, feasible_pct, best_score,
            );
            let _ = writeln!(
                summary_csv,
                "{},{},{},{},{},{:.2},{:.6}",
                mode_key, topo_short, total, evaluated_n, feasible, feasible_pct, best_score,
            );
        }
        println!("  {}", "-".repeat(80));
    }

    let summary_csv_path = out_dir.join("topology_mode_feasibility.csv");
    write_file(
        &summary_csv_path,
        &summary_csv,
        "topology_mode_feasibility.csv",
    );

    println!("\n=== Done.  Outputs: {} ===", out_dir.display());
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn write_file(path: &Path, content: &str, label: &str) {
    match std::fs::write(path, content) {
        Ok(()) => println!("  Saved: {}", path.display()),
        Err(e) => eprintln!("  WARN: {label} write failed: {e}"),
    }
}

fn inspect_preset(def: &PresetDef) -> PresetRow {
    let bp = (def.build)(def.id);
    let venturi_channels = bp
        .channels
        .iter()
        .filter(|c| c.id.as_str().contains("throat"))
        .count();
    let serpentine_channels = bp
        .channels
        .iter()
        .filter(|c| {
            let id = c.id.as_str();
            id.contains("segment")
                || id.contains("serp")
                || id.contains("spiral")
                || id.starts_with("wide_")
                || id.starts_with("narrow_")
                || id.starts_with("ch_")
        })
        .count();
    PresetRow {
        preset_id: def.id.to_string(),
        family: def.family.to_string(),
        support: def.support.to_string(),
        mapped_to: def.mapped_to.to_string(),
        nodes: bp.nodes.len(),
        channels: bp.channels.len(),
        inlets: bp.inlet_count(),
        outlets: bp.outlet_count(),
        venturi_channels,
        serpentine_channels,
        total_length_mm: bp.total_length_m() * 1e3,
    }
}

fn classify_zero(c: &DesignCandidate, mode: OptimMode, weights: &SdtWeights) -> PrimaryReason {
    let Ok(m) = compute_metrics(c) else {
        return PrimaryReason::MetricsError;
    };
    let score = score_candidate(&m, mode, weights);
    if score > 0.0 {
        return PrimaryReason::Feasible;
    }
    if !m.pressure_feasible {
        return PrimaryReason::PressureFail;
    }
    if !m.fda_main_compliant {
        return PrimaryReason::FdaFail;
    }
    if !m.plate_fits {
        return PrimaryReason::PlateFail;
    }
    if m.total_ecv_ml <= 0.0 {
        return PrimaryReason::EcvZero;
    }
    if m.wbc_recovery <= 1e-9 {
        return PrimaryReason::WbcZero;
    }
    PrimaryReason::ScoreOther
}

fn preset_defs() -> Vec<PresetDef> {
    vec![
        PresetDef {
            id: "venturi_chain",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| venturi_chain(n, 45e-3, 4e-3, 150e-6),
        },
        PresetDef {
            id: "venturi_rect",
            family: "base_rectangular",
            support: "parametric+ga",
            mapped_to: "SV",
            build: |n| venturi_rect(n, 2e-3, 100e-6, 1e-3, 300e-6),
        },
        PresetDef {
            id: "symmetric_bifurcation",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| symmetric_bifurcation(n, 12e-3, 6e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "bifurcation_rect",
            family: "base_rectangular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| bifurcation_rect(n, 12e-3, 6e-3, 2e-3, 1e-3, 1e-3),
        },
        PresetDef {
            id: "symmetric_trifurcation",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| symmetric_trifurcation(n, 12e-3, 6e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "trifurcation_rect",
            family: "base_rectangular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| trifurcation_rect(n, 12e-3, 6e-3, 2e-3, 1e-3, 1e-3),
        },
        PresetDef {
            id: "serpentine_chain",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: |n| serpentine_chain(n, 6, 7.5e-3, 2e-3),
        },
        PresetDef {
            id: "serpentine_rect",
            family: "base_rectangular",
            support: "parametric+ga",
            mapped_to: "SG",
            build: |n| serpentine_rect(n, 6, 7.5e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "venturi_serpentine_rect",
            family: "series_composite",
            support: "parametric+ga",
            mapped_to: "VS",
            build: |n| venturi_serpentine_rect(n, 2e-3, 100e-6, 1e-3, 300e-6, 6, 7.5e-3),
        },
        PresetDef {
            id: "serial_double_venturi_rect",
            family: "series_composite",
            support: "parametric+ga",
            mapped_to: "S2",
            build: |n| serial_double_venturi_rect(n, 2e-3, 100e-6, 1e-3, 300e-6, 7.5e-3),
        },
        PresetDef {
            id: "bifurcation_venturi_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "BV",
            build: |n| bifurcation_venturi_rect(n, 15e-3, 2e-3, 100e-6, 1e-3, 300e-6),
        },
        PresetDef {
            id: "bifurcation_serpentine_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "BS",
            build: |n| bifurcation_serpentine_rect(n, 15e-3, 6, 7.5e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "double_bifurcation_serpentine_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "D4S",
            build: |n| double_bifurcation_serpentine_rect(n, 15e-3, 6, 7.5e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "trifurcation_venturi_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "TV",
            build: |n| trifurcation_venturi_rect(n, 15e-3, 2e-3, 100e-6, 1e-3, 300e-6),
        },
        PresetDef {
            id: "trifurcation_serpentine_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "TS",
            build: |n| trifurcation_serpentine_rect(n, 15e-3, 6, 7.5e-3, 2e-3, 1e-3),
        },
        PresetDef {
            id: "double_bifurcation_venturi_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "DBV",
            build: |n| double_bifurcation_venturi_rect(n, 12e-3, 9e-3, 2e-3, 100e-6, 1e-3, 300e-6),
        },
        PresetDef {
            id: "triple_bifurcation_venturi_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "3BV",
            build: |n| {
                triple_bifurcation_venturi_rect(n, 10e-3, 8e-3, 6e-3, 2e-3, 100e-6, 1e-3, 300e-6)
            },
        },
        PresetDef {
            id: "double_trifurcation_venturi_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "2TV",
            build: |n| double_trifurcation_venturi_rect(n, 12e-3, 8e-3, 2e-3, 100e-6, 1e-3, 300e-6),
        },
        PresetDef {
            id: "triple_trifurcation_venturi_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "3TV",
            build: |n| {
                triple_trifurcation_venturi_rect(
                    n, 9e-3, 7e-3, 5e-3, 2e-3, 0.45, 100e-6, 1e-3, 300e-6,
                )
            },
        },
        PresetDef {
            id: "quad_trifurcation_venturi_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "4TV",
            build: |n| {
                quad_trifurcation_venturi_rect(
                    n, 8e-3, 6e-3, 5e-3, 4e-3, 2e-3, 0.45, 100e-6, 1e-3, 300e-6,
                )
            },
        },
        PresetDef {
            id: "bifurcation_trifurcation_venturi_rect",
            family: "mixed_composite",
            support: "parametric+ga",
            mapped_to: "BTV",
            build: |n| {
                bifurcation_trifurcation_venturi_rect(n, 12e-3, 8e-3, 2e-3, 100e-6, 1e-3, 300e-6)
            },
        },
        PresetDef {
            id: "trifurcation_bifurcation_venturi_rect",
            family: "mixed_composite",
            support: "parametric+ga",
            mapped_to: "TBV",
            build: |n| {
                trifurcation_bifurcation_venturi_rect(n, 12e-3, 8e-3, 2e-3, 100e-6, 1e-3, 300e-6)
            },
        },
        PresetDef {
            id: "trifurcation_bifurcation_bifurcation_venturi_rect",
            family: "mixed_composite",
            support: "parametric+ga",
            mapped_to: "TBBV",
            build: |n| {
                trifurcation_bifurcation_bifurcation_venturi_rect(
                    n, 10e-3, 7e-3, 5e-3, 2e-3, 0.45, 100e-6, 1e-3, 300e-6,
                )
            },
        },
        PresetDef {
            id: "cascade_center_trifurcation_rect",
            family: "cascade_composite",
            support: "parametric+ga",
            mapped_to: "PST",
            build: |n| {
                cascade_center_trifurcation_rect(
                    n, 12e-3, 8e-3, 2, 2e-3, 0.45, 100e-6, 300e-6, 1e-3, true, None,
                )
            },
        },
        PresetDef {
            id: "incremental_filtration_tri_bi_rect",
            family: "cascade_composite",
            support: "parametric+ga",
            mapped_to: "CIF",
            build: |n| {
                incremental_filtration_tri_bi_rect(
                    n, 12e-3, 8e-3, 6e-3, 2, 2e-3, 0.45, 0.68, 100e-6, 300e-6, 1e-3,
                )
            },
        },
        PresetDef {
            id: "cell_separation_rect",
            family: "cell_separation",
            support: "parametric+ga",
            mapped_to: "CS",
            build: |n| cell_separation_rect(n, 22.5e-3, 2e-3, 100e-6, 1e-3, 300e-6, 22.5e-3),
        },
        PresetDef {
            id: "constriction_expansion_array_rect",
            family: "constriction_array",
            support: "parametric+ga",
            mapped_to: "CE",
            build: |n| constriction_expansion_array_rect(n, 10, 3e-3, 1.5e-3, 200e-6, 80e-6, 60e-6),
        },
        PresetDef {
            id: "spiral_channel_rect",
            family: "spiral",
            support: "parametric+ga",
            mapped_to: "SP",
            build: |n| spiral_channel_rect(n, 8, 5e-3, 200e-6, 60e-6),
        },
        PresetDef {
            id: "parallel_microchannel_array_rect",
            family: "parallel_array",
            support: "parametric+ga",
            mapped_to: "PM",
            build: |n| parallel_microchannel_array_rect(n, 64, 45e-3, 200e-6, 60e-6),
        },
        PresetDef {
            id: "asymmetric_bifurcation_serpentine_rect",
            family: "asymmetric",
            support: "parametric+ga",
            mapped_to: "ABS",
            build: |n| asymmetric_bifurcation_serpentine_rect(n, 15e-3, 6, 7.5e-3, 0.5, 2e-3, 1e-3),
        },
    ]
}
