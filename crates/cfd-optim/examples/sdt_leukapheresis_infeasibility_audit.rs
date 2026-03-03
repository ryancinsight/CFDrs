//! Targeted infeasibility audit for leukapheresis-focused topologies.
//!
//! This diagnostic isolates `CE`, `SP`, and `PM` candidates from the full
//! parametric space and reports why pediatric and combined scores remain zero.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example sdt_leukapheresis_infeasibility_audit --no-default-features
//! ```

use cfd_optim::{
    build_candidate_space, compute_metrics, score_candidate, DesignCandidate, DesignTopology,
    OptimMode, SdtWeights,
};
use std::fmt::Write as _;
use std::path::Path;

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
    fn record_reason(&mut self, reason: PrimaryReason) {
        match reason {
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
    topology_short: String,
    topology_name: String,
    total_candidates: usize,
    pediatric: ModeCounters,
    combined: ModeCounters,
}

fn is_target_topology(topology: DesignTopology) -> bool {
    matches!(
        topology,
        DesignTopology::ConstrictionExpansionArray { .. }
            | DesignTopology::SpiralSerpentine { .. }
            | DesignTopology::ParallelMicrochannelArray { .. }
    )
}

fn classify_score_zero(
    candidate: &DesignCandidate,
    mode: OptimMode,
    weights: &SdtWeights,
) -> PrimaryReason {
    let metrics = match compute_metrics(candidate) {
        Ok(m) => m,
        Err(_) => return PrimaryReason::MetricsError,
    };

    let score = score_candidate(&metrics, mode, weights);
    if score > 0.0 {
        return PrimaryReason::Feasible;
    }

    if !metrics.pressure_feasible {
        return PrimaryReason::PressureFail;
    }
    if !metrics.fda_main_compliant {
        return PrimaryReason::FdaFail;
    }
    if !metrics.plate_fits {
        return PrimaryReason::PlateFail;
    }
    if metrics.total_ecv_ml <= 0.0 {
        return PrimaryReason::EcvZero;
    }
    if metrics.wbc_recovery <= 1.0e-9 {
        return PrimaryReason::WbcZero;
    }

    PrimaryReason::ScoreOther
}

fn main() {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("sdt_scenario_audit");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    let report_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map(|p| p.join("report"))
        .expect("failed to resolve report directory");
    std::fs::create_dir_all(&report_dir).expect("failed to create report directory");

    let pediatric_mode = OptimMode::PediatricLeukapheresis {
        patient_weight_kg: 3.0,
    };
    let combined_mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let weights = SdtWeights::default();

    let mut audits: Vec<TopologyAudit> = Vec::new();
    for (short, name) in [
        ("CE", "Constriction-Expansion Array"),
        ("SP", "Spiral Serpentine"),
        ("PM", "Parallel Microchannel Array"),
    ] {
        audits.push(TopologyAudit {
            topology_short: short.to_string(),
            topology_name: name.to_string(),
            ..TopologyAudit::default()
        });
    }

    let mut candidates = build_candidate_space();
    candidates.retain(|c| is_target_topology(c.topology));

    for c in &candidates {
        let slot = audits
            .iter_mut()
            .find(|a| a.topology_short == c.topology.short())
            .expect("target topology must exist in audit table");
        slot.total_candidates += 1;

        let pedi_reason = classify_score_zero(c, pediatric_mode, &weights);
        slot.pediatric.record_reason(pedi_reason);

        let comb_reason = classify_score_zero(c, combined_mode, &weights);
        slot.combined.record_reason(comb_reason);
    }

    // CSV output
    let csv_path = out_dir.join("leukapheresis_infeasibility_audit.csv");
    let mut csv = String::from(
        "topology_short,topology_name,total,mode,positive_score,metrics_error,pressure_fail,fda_fail,plate_fail,ecv_zero,wbc_zero,score_other\n",
    );
    for a in &audits {
        for (mode_key, c) in [
            ("pediatric_3kg", &a.pediatric),
            ("combined_sdt_leuka", &a.combined),
        ] {
            let _ = writeln!(
                csv,
                "{},{},{},{},{},{},{},{},{},{},{},{}",
                a.topology_short,
                a.topology_name,
                a.total_candidates,
                mode_key,
                c.positive_score,
                c.metrics_error,
                c.pressure_fail,
                c.fda_fail,
                c.plate_fail,
                c.ecv_zero,
                c.wbc_zero,
                c.score_other
            );
        }
    }
    std::fs::write(&csv_path, csv).expect("failed to write audit csv");

    // Markdown summary output
    let md_path = report_dir.join("milestone12_leukapheresis_infeasibility.md");
    let mut md = String::new();
    md.push_str("# Milestone 12 Leukapheresis Infeasibility Audit\n\n");
    md.push_str("Target topologies: `CE`, `SP`, `PM` from `build_candidate_space()`.\n\n");
    md.push_str("| Topology | Total | Mode | Score>0 | MetricsErr | PressureFail | FdaFail | PlateFail | ECV=0 | WBC=0 | Other |\n");
    md.push_str("|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    for a in &audits {
        for (mode_label, c) in [("Pediatric (3kg)", &a.pediatric), ("Combined", &a.combined)] {
            let _ = writeln!(
                md,
                "| {} ({}) | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
                a.topology_name,
                a.topology_short,
                a.total_candidates,
                mode_label,
                c.positive_score,
                c.metrics_error,
                c.pressure_fail,
                c.fda_fail,
                c.plate_fail,
                c.ecv_zero,
                c.wbc_zero,
                c.score_other
            );
        }
    }

    md.push_str("\n## Primary Reason Priority\n\n");
    md.push_str("Classification order for score=0 was: metrics error -> pressure -> FDA main shear -> plate fit -> ECV -> WBC -> other.\n");
    std::fs::write(&md_path, md).expect("failed to write markdown summary");

    println!("Saved: {}", csv_path.display());
    println!("Saved: {}", md_path.display());
}
