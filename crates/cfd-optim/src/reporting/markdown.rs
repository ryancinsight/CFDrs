//! Canonical Milestone 12 markdown report writer.
//!
//! [`write_milestone12_results`] writes `report/milestone12_results.md`.
//! [`ValidationRow`] stores one row of multi-fidelity validation output.

use std::fmt::Write as FmtWrite;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::analysis::RobustnessReport;
use crate::constraints::{
    M12_GA_HYDRO_SEED, MILESTONE_TREATMENT_DURATION_MIN,
    PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG,
};
use crate::reporting::Milestone12ReportDesign;

use super::ranking::oncology_priority_score;

fn workspace_geometry_penalty_weight() -> f64 {
    std::env::var("M12_GA_GEOMETRY_PENALTY_WEIGHT")
        .ok()
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value >= 0.0)
        .unwrap_or(0.00043)
}

fn workspace_cavitation_gain_margin() -> f64 {
    std::env::var("M12_GA_CAVITATION_GAIN_MARGIN")
        .ok()
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value >= 0.0)
        .unwrap_or(0.01)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationRow {
    pub track: String,
    pub id: String,
    pub topology: String,
    /// 1D total pressure-loss coefficient `K = ΔP / (1/2 ρ u_in^2)`.
    #[serde(default)]
    pub k_loss_1d: f64,
    pub dp_1d_bernoulli_pa: f64,
    /// 2D total pressure-loss coefficient `K = ΔP / (1/2 ρ u_in^2)`.
    #[serde(default)]
    pub k_loss_2d: f64,
    pub dp_2d_fvm_pa: f64,
    /// 3D total pressure-loss coefficient `K = ΔP / (1/2 ρ u_in^2)`.
    #[serde(default)]
    pub k_loss_3d: f64,
    pub dp_3d_fem_pa: f64,
    pub agreement_1d_2d_pct: f64,
    pub agreement_2d_3d_pct: f64,
    pub mass_error_3d_pct: f64,
    pub sigma_1d: f64,
    pub sigma_2d: f64,
    pub score: f64,
    /// Whether the 2D FVM SIMPLE solver converged and the velocity field is physical.
    #[serde(default)]
    pub two_d_converged: bool,
    /// Whether the throat Reynolds number exceeds the laminar Navier-Stokes validity limit (Re > 2000).
    /// When true the 3D Laminar FEM ΔP will be << Bernoulli (inertia-dominated turbulent flow).
    #[serde(default)]
    pub high_re_laminar_mismatch: bool,
}

pub fn write_milestone12_results(
    total_candidates: usize,
    option1_pool_len: usize,
    option2_pool_len: usize,
    option1_ranked: &[Milestone12ReportDesign],
    option2_ranked: &[Milestone12ReportDesign],
    ga_top: &Milestone12ReportDesign,
    _validation_rows: &[ValidationRow],
    option2_robustness: &[RobustnessReport],
    authoritative_run: bool,
    canonical_source: &str,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut md = String::new();

    let option1 = option1_ranked.first();
    let option2 = option2_ranked
        .first()
        .ok_or("option2_ranked is empty; strict-core Option 2 selection failed")?;

    writeln!(md, "# Milestone 12 Results - Canonical\n")?;
    writeln!(
        md,
        "> Generated artifact. Do not edit manually.\n> Source: `cargo run -p cfd-optim --example milestone12_report --no-default-features`.\n> Artifact class: **{}**.\n> Published location: `{canonical_source}`.\n",
        if authoritative_run {
            "authoritative report"
        } else {
            "fast-mode report"
        }
    )?;

    writeln!(md, "## Dataset And Strict Eligibility")?;
    writeln!(md, "- Total candidates in the canonical parameter lattice: **{total_candidates}**")?;
    writeln!(
        md,
        "- Option 1 eligibility pool (selective acoustic + pressure/plate/FDA-main): **{option1_pool_len}**"
    )?;
    writeln!(
        md,
        "- Option 2 eligibility pool (selective venturi + pressure/plate/FDA-main + sigma<1): **{option2_pool_len}**\n"
    )?;

    writeln!(md, "## Selected Designs")?;
    writeln!(
        md,
        "ECV is reported explicitly as `ECV = Σ(L_i A_i) = Q_in t_res` and benchmarked against the 3 kg neonatal 10% circuit-volume limit of **25.5 mL**.\n"
    )?;
    writeln!(
        md,
        "| Track | Candidate | Topology | Mode | Active venturi throats | Score | sigma | K_loss | Cumulative cavitation dose | WBC treatment exposure | RBC treatment exposure | HI/pass | P95 shear (Pa) | ECV (mL) | ECV / 3kg limit (%) |"
    )?;
    writeln!(
        md,
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    )?;
    if let Some(option1) = option1 {
        writeln!(
            md,
            "| Option 1 (Selective acoustic center treatment) | `{}` | {} | {} | {} | {:.4} | n/a | n/a | {:.4} | {:.4} | {:.4} | {:.2} | {:.3} | {:.1} |",
            option1.candidate.id,
            option1.topology_display_name(),
            option1.metrics.treatment_zone_mode,
            option1.metrics.active_venturi_throat_count,
            option1.score,
            option1.metrics.wbc_recovery,
            option1.metrics.rbc_pass_fraction,
            option1.metrics.hemolysis_index_per_pass,
            option1.metrics.wall_shear_p95_pa,
            option1.metrics.total_ecv_ml,
            100.0 * option1.metrics.total_ecv_ml / 25.5,
        )?;
    } else {
        writeln!(
            md,
            "| Option 1 (Selective acoustic center treatment) | _No eligible design under current physics regime_ | n/a | Unavailable | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |"
        )?;
    }
    writeln!(
        md,
        "| Option 2 (Selective venturi cavitation selectivity) | `{}` | {} | {} | {} | {:.4} | {:.4} | {:.3} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2} | {:.3} | {:.1} |",
        option2.candidate.id,
        option2.topology_display_name(),
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.score,
        option2.metrics.cavitation_number,
        option2.metrics.venturi_total_loss_coefficient,
        option2.metrics.serial_cavitation_dose_fraction,
        option2.metrics.wbc_recovery,
        option2.metrics.rbc_venturi_exposure_fraction,
        option2.metrics.hemolysis_index_per_pass,
        option2.metrics.wall_shear_p95_pa,
        option2.metrics.total_ecv_ml,
        100.0 * option2.metrics.total_ecv_ml / 25.5,
    )?;
    writeln!(md)?;

    writeln!(md, "## Option 2 Ranked Top-5 (Deterministic Tie-Break)")?;
    writeln!(
        md,
        "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, clot risk asc, candidate id asc. `sigma` remains a per-throat severity metric; `K_loss` is the inlet-normalized venturi pressure-loss coefficient; `Cumulative cavitation dose` reports serial exposure across the treatment path.\n"
    )?;
    writeln!(
        md,
        "| Rank | Candidate | Mode | Active venturi throats | Score | Oncology priority | RBC venturi exposure | Clot risk | sigma | K_loss | Cumulative cavitation dose |"
    )?;
    writeln!(md, "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")?;
    for d in option2_ranked.iter().take(5) {
        writeln!(
            md,
            "| {} | `{}` | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.3} | {:.4} |",
            d.rank,
            d.candidate.id,
            d.metrics.treatment_zone_mode,
            d.metrics.active_venturi_throat_count,
            d.score,
            oncology_priority_score(&d.metrics),
            d.metrics.rbc_venturi_exposure_fraction,
            d.metrics.clotting_risk_index,
            d.metrics.cavitation_number,
            d.metrics.venturi_total_loss_coefficient,
            d.metrics.serial_cavitation_dose_fraction,
        )?;
    }
    writeln!(md)?;

    writeln!(
        md,
        "## Strict-Core Gate Evidence (Selected Option 2 Venturi Design)"
    )?;
    writeln!(
        md,
        "| Track | Candidate | Pressure feasible | Plate fits | FDA main | sigma finite | sigma<1 |"
    )?;
    writeln!(md, "|---|---|---|---|---|---|---|")?;
    {
        let (label, d) = ("Option 2 Combined Selective Venturi", option2);
        let m = &d.metrics;
        writeln!(
            md,
            "| {} | `{}` | {} | {} | {} | {} | {} |",
            label,
            d.candidate.id,
            if m.pressure_feasible { "PASS" } else { "FAIL" },
            if m.plate_fits { "PASS" } else { "FAIL" },
            if m.fda_main_compliant { "PASS" } else { "FAIL" },
            if m.cavitation_number.is_finite() {
                "PASS"
            } else {
                "FAIL"
            },
            if m.cavitation_number.is_finite() && m.cavitation_number < 1.0 {
                "PASS"
            } else {
                "FAIL"
            },
        )?;
    }
    writeln!(md)?;

    if !option2_robustness.is_empty() {
        writeln!(
            md,
            "## Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)"
        )?;
        let robust_count = option2_robustness.iter().filter(|r| r.is_robust).count();
        writeln!(
            md,
            "- Robust pass-rate: **{}/{}** candidates",
            robust_count,
            option2_robustness.len()
        )?;
        writeln!(
            md,
            "| Candidate | Nominal score | Min score | Max score | CV | Robust | Worst-case parameter |"
        )?;
        writeln!(md, "|---|---:|---:|---:|---:|---|---|")?;
        for r in option2_robustness {
            writeln!(
                md,
                "| `{}` | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} |",
                r.candidate_id,
                r.score_nominal,
                r.score_min,
                r.score_max,
                r.score_cv,
                if r.is_robust { "YES" } else { "NO" },
                r.worst_case_param
            )?;
        }
        writeln!(md)?;
    } else {
        writeln!(md, "## Robustness Screening")?;
        writeln!(
            md,
            "No standalone Option 2 perturbation sweep was emitted for this canonical run, so no robustness table is included in this artifact."
        )?;
        writeln!(md)?;
    }

    writeln!(md, "## GA Reproducibility")?;
    writeln!(md, "- HydroSDT GA RNG seed: `{M12_GA_HYDRO_SEED}`")?;
    writeln!(
        md,
        "- GA rank-1: `{}` (score {:.4}, sigma {:.4})\n",
        ga_top.candidate.id, ga_top.score, ga_top.metrics.cavitation_number
    )?;

    writeln!(md, "## Workspace-Level Configuration")?;
    writeln!(
        md,
        "- `M12_GA_GEOMETRY_PENALTY_WEIGHT={:.5}`: lineage-concentration penalty applied during HydroSDT GA shortlist ranking.",
        workspace_geometry_penalty_weight()
    )?;
    writeln!(
        md,
        "- `M12_GA_CAVITATION_GAIN_MARGIN={:.3}`: minimum cumulative cavitation-dose gain over the selected Option 2 baseline required for a GA venturi design to remain in the report-ranked HydroSDT subset.\n",
        workspace_cavitation_gain_margin()
    )?;

    writeln!(
        md,
        "## Limits Of Usage (From Selected Option 2 Oncology-Directed Design)"
    )?;
    writeln!(
        md,
        "- Operating point: **{:.1} mL/min**, **{:.0} kPa gauge**",
        option2.flow_rate_ml_min(),
        option2.inlet_gauge_kpa()
    )?;
    writeln!(
        md,
        "- Venturi throat: **{:.0} um**, throat length **{:.0} um**",
        option2.throat_width_um().unwrap_or(0.0),
        option2.throat_length_um().unwrap_or(0.0)
    )?;
    writeln!(
        md,
        "- Treatment mode: **{}**, active venturi throats: **{}** (serial stages per path: **{}**)",
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.metrics.serial_venturi_stages_per_path
    )?;
    writeln!(
        md,
        "- FDA main-channel compliance: **{}**",
        if option2.metrics.fda_main_compliant {
            "PASS"
        } else {
            "FAIL"
        }
    )?;
    writeln!(
        md,
        "- Clot-risk indicators: `clotting_risk_index={:.4}`, `clotting_risk_index_10ml_s={:.4}`",
        option2.metrics.clotting_risk_index, option2.metrics.clotting_risk_index_10ml_s
    )?;
    writeln!(
        md,
        "- Flow caution flags: `Q>=200={}`, `Q>=600={}`",
        if option2.metrics.clotting_flow_compliant {
            "PASS"
        } else {
            "FAIL"
        },
        if option2.metrics.clotting_flow_compliant_10ml_s {
            "PASS"
        } else {
            "FAIL"
        }
    )?;

    let q_ml_min = option2.flow_rate_ml_min().max(1.0e-9);
    let q_ml_s = q_ml_min / 60.0;
    let feed_hematocrit = option2
        .candidate
        .operating_point
        .feed_hematocrit
        .clamp(0.0, 0.95);
    let plasma_fraction = (1.0 - feed_hematocrit).max(0.05);
    let blood_volume_ml = option2
        .candidate
        .operating_point
        .patient_context
        .as_ref()
        .map_or(
            PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG,
            |ctx| ctx.blood_volume_ml.max(1.0),
        );
    let plasma_volume_ml = blood_volume_ml * plasma_fraction;
    let circuit_turnover_s = option2.metrics.total_ecv_ml / q_ml_s;
    let time_to_one_pv_min = plasma_volume_ml / (q_ml_min * plasma_fraction);
    let processed_plasma_15_ml = q_ml_min * MILESTONE_TREATMENT_DURATION_MIN * plasma_fraction;
    let pv_15 = processed_plasma_15_ml / plasma_volume_ml;

    writeln!(md, "\n## Treatment-Time Analysis")?;
    writeln!(
        md,
        "- Feed hematocrit: **{feed_hematocrit:.2}** → plasma fraction **{plasma_fraction:.2}**"
    )?;
    writeln!(
        md,
        "- Reference blood / plasma volume: **{blood_volume_ml:.1} mL / {plasma_volume_ml:.1} mL**"
    )?;
    writeln!(
        md,
        "- Millifluidic ECV: **{:.3} mL** ({:.2}% of plasma volume)",
        option2.metrics.total_ecv_ml,
        100.0 * option2.metrics.total_ecv_ml / plasma_volume_ml
    )?;
    writeln!(
        md,
        "- Circuit turnover time (`ECV / Q`): **{circuit_turnover_s:.2} s**"
    )?;
    writeln!(
        md,
        "- Treatment-zone residence time per pass: **{:.2} ms**",
        option2.metrics.mean_residence_time_s * 1.0e3
    )?;
    writeln!(
        md,
        "- Plasma processed in {MILESTONE_TREATMENT_DURATION_MIN:.0} min: **{processed_plasma_15_ml:.1} mL** ({pv_15:.2} plasma-volume equivalents)"
    )?;
    writeln!(
        md,
        "- Single-pass transit time for one plasma-volume equivalent: **{time_to_one_pv_min:.1} min** (not the clinical session duration)"
    )?;

    std::fs::write(path, md)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use crate::domain::fixtures::{
        canonical_option1_candidate, canonical_option2_candidate, operating_point,
    };
    use crate::reporting::{compute_blueprint_report_metrics, Milestone12ReportDesign};

    use super::write_milestone12_results;

    #[test]
    fn canonical_markdown_uses_treatment_exposure_language() {
        let operating_point = operating_point(2.0e-6, 30_000.0, 0.18);
        let option1_candidate = canonical_option1_candidate("option1", operating_point.clone());
        let option2_candidate = canonical_option2_candidate("option2", operating_point);
        let option1 = Milestone12ReportDesign::new(
            1,
            option1_candidate.clone(),
            compute_blueprint_report_metrics(&option1_candidate).expect("option1 metrics"),
            0.61,
        );
        let option2 = Milestone12ReportDesign::new(
            1,
            option2_candidate.clone(),
            compute_blueprint_report_metrics(&option2_candidate).expect("option2 metrics"),
            0.74,
        );

        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("m12-markdown-{unique}.md"));
        write_milestone12_results(
            2,
            1,
            1,
            std::slice::from_ref(&option1),
            std::slice::from_ref(&option2),
            &option2,
            &[],
            &[],
            true,
            "report/milestone12_results.md",
            &path,
        )
        .expect("canonical markdown should render");

        let rendered = std::fs::read_to_string(path).expect("rendered markdown should exist");
        assert!(rendered.contains("WBC treatment exposure"));
        assert!(rendered.contains("Treatment-Time Analysis"));
        assert!(rendered.contains("K_loss"));
        assert!(!rendered.contains("WBC recovery"));
        assert!(!rendered.to_ascii_lowercase().contains("leukapheresis"));
    }
}
