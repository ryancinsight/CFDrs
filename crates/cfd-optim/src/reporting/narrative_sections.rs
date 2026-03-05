//! Section builders for Milestone 12 narrative markdown.

use std::fmt::Write as _;

use crate::analysis::RobustnessReport;
use crate::reporting::figures::NarrativeFigureSpec;
use crate::reporting::ValidationRow;
use crate::RankedDesign;

pub(super) fn build_selected_table(option1: &RankedDesign, option2: &RankedDesign) -> String {
    let mut out = String::new();
    out.push_str(
        "| Track | Candidate | Topology | Mode | Active venturi throats | Score | sigma | WBC recovery | RBC venturi exposure | HI/pass | P95 shear (Pa) | ECV (mL) |\n",
    );
    out.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    let _ = writeln!(
        out,
        "| Option 1 (Selective acoustic center treatment) | `{}` | {} | {} | {} | {:.4} | n/a | {:.4} | {:.4} | {:.2e} | {:.2} | {:.3} |",
        option1.candidate.id,
        option1.candidate.topology.name(),
        option1.metrics.treatment_zone_mode,
        option1.metrics.active_venturi_throat_count,
        option1.score,
        option1.metrics.wbc_recovery,
        option1.metrics.rbc_venturi_exposure_fraction,
        option1.metrics.hemolysis_index_per_pass,
        option1.metrics.wall_shear_p95_pa,
        option1.metrics.total_ecv_ml
    );
    let _ = writeln!(
        out,
        "| Option 2 (Selective venturi, combined cancer + leukapheresis score) | `{}` | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2e} | {:.2} | {:.3} |",
        option2.candidate.id,
        option2.candidate.topology.name(),
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.score,
        option2.metrics.cavitation_number,
        option2.metrics.wbc_recovery,
        option2.metrics.rbc_venturi_exposure_fraction,
        option2.metrics.hemolysis_index_per_pass,
        option2.metrics.wall_shear_p95_pa,
        option2.metrics.total_ecv_ml
    );
    out
}

pub(super) fn build_option2_top5_table(option2_ranked: &[RankedDesign]) -> String {
    let mut out = String::new();
    out.push_str("| Rank | Candidate | Mode | Active venturi throats | Score | RBC venturi exposure | Clot risk | sigma |\n");
    out.push_str("|---:|---|---|---:|---:|---:|---:|---:|\n");
    for d in option2_ranked.iter().take(5) {
        let _ = writeln!(
            out,
            "| {} | `{}` | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} |",
            d.rank,
            d.candidate.id,
            d.metrics.treatment_zone_mode,
            d.metrics.active_venturi_throat_count,
            d.score,
            d.metrics.rbc_venturi_exposure_fraction,
            d.metrics.clotting_risk_index,
            d.metrics.cavitation_number
        );
    }
    out
}

pub(super) fn build_tri_cell_table(
    option2_ranked: &[RankedDesign],
    ga_ranked: &[RankedDesign],
) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | cancer_center_fraction | wbc_recovery | rbc_venturi_exposure | three_pop_sep_efficiency |\n");
    out.push_str("|---|---|---:|---:|---:|---:|\n");
    for (label, d_opt) in [
        ("Option 2 Combined", option2_ranked.first()),
        ("GA HydroSDT", ga_ranked.first()),
    ] {
        if let Some(d) = d_opt {
            let _ = writeln!(
                out,
                "| {} | `{}` | {:.4} | {:.4} | {:.4} | {:.4} |",
                label,
                d.candidate.id,
                d.metrics.cancer_center_fraction,
                d.metrics.wbc_recovery,
                d.metrics.rbc_venturi_exposure_fraction,
                d.metrics.three_pop_sep_efficiency
            );
        }
    }
    out
}

pub(super) fn build_strict_core_table(
    option2_ranked: &[RankedDesign],
    ga_ranked: &[RankedDesign],
) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | Pressure feasible | Plate fits | FDA main | sigma finite | sigma<1 |\n");
    out.push_str("|---|---|---|---|---|---|---|\n");
    for (label, d_opt) in [
        ("Option 2 Combined", option2_ranked.first()),
        ("GA HydroSDT", ga_ranked.first()),
    ] {
        if let Some(d) = d_opt {
            let m = &d.metrics;
            let _ = writeln!(
                out,
                "| {} | `{}` | {} | {} | {} | {} | {} |",
                label,
                d.candidate.id,
                pass_fail(m.pressure_feasible),
                pass_fail(m.plate_fits),
                pass_fail(m.fda_main_compliant),
                pass_fail(m.cavitation_number.is_finite()),
                pass_fail(m.cavitation_number.is_finite() && m.cavitation_number < 1.0)
            );
        }
    }
    out
}

pub(super) fn build_robustness_section(robustness: &[RobustnessReport], fast_mode: bool) -> String {
    if robustness.is_empty() {
        return String::new();
    }
    let _ = fast_mode;
    let mut out = String::new();
    out.push_str("### Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)\n\n");
    let robust_count = robustness.iter().filter(|r| r.is_robust).count();
    let _ = writeln!(
        out,
        "- Robust pass-rate: **{}/{}** candidates\n",
        robust_count,
        robustness.len()
    );
    out.push_str(
        "| Candidate | Nominal score | Min score | Max score | CV | Robust | Worst-case parameter |\n",
    );
    out.push_str("|---|---:|---:|---:|---:|---|---|\n");
    for r in robustness {
        let _ = writeln!(
            out,
            "| `{}` | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} |",
            r.candidate_id,
            r.score_nominal,
            r.score_min,
            r.score_max,
            r.score_cv,
            if r.is_robust { "YES" } else { "NO" },
            r.worst_case_param
        );
    }
    out
}

pub(super) fn build_validation_section(rows: &[ValidationRow], fast_mode: bool) -> String {
    if rows.is_empty() {
        return String::new();
    }
    let _ = fast_mode;
    let mut out = String::new();
    out.push_str("### Multi-Fidelity Pressure-Drop Validation (Selected Venturi Designs)\n\n");
    out.push_str("| Track | Candidate | dp1D Bernoulli (Pa) | dp2D FVM (Pa) | dp3D FEM (Pa) | 1D-2D diff (%) | 2D-3D diff (%) | Mass error (%) |\n");
    out.push_str("|---|---|---:|---:|---:|---:|---:|---:|\n");
    for row in rows {
        let _ = writeln!(
            out,
            "| {} | `{}` | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} |",
            row.track,
            row.id,
            row.dp_1d_bernoulli_pa,
            row.dp_2d_fvm_pa,
            row.dp_3d_fem_pa,
            row.agreement_1d_2d_pct,
            row.agreement_2d_3d_pct,
            row.mass_error_3d_pct
        );
    }
    out
}

pub(super) fn build_limits_of_usage(option2: &RankedDesign) -> String {
    format!(
        "- Operating point: **{:.1} mL/min**, **{:.0} kPa gauge**\n- Venturi throat: **{:.0} um**, throat length **{:.0} um**\n- Treatment mode: **{}**, active venturi throats: **{}** (serial stages per path: **{}**)\n- FDA main-channel compliance: **{}**\n- Clot-risk indicators: `clotting_risk_index={:.4}`, `clotting_risk_index_10ml_s={:.4}`\n- Flow caution flags: `Q>=200={}`, `Q>=600={}`",
        option2.candidate.flow_rate_m3_s * 6.0e7,
        option2.candidate.inlet_gauge_pa * 1.0e-3,
        option2.candidate.throat_diameter_m * 1.0e6,
        option2.candidate.throat_length_m * 1.0e6,
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.metrics.serial_venturi_stages_per_path,
        pass_fail(option2.metrics.fda_main_compliant),
        option2.metrics.clotting_risk_index,
        option2.metrics.clotting_risk_index_10ml_s,
        pass_fail(option2.metrics.clotting_flow_compliant),
        pass_fail(option2.metrics.clotting_flow_compliant_10ml_s)
    )
}

pub(super) fn build_figure_toc_rows(specs: &[NarrativeFigureSpec]) -> String {
    let mut out = String::new();
    for spec in specs {
        let _ = writeln!(
            out,
            "| {} | [{}](#fig-{}) |",
            spec.number, spec.title, spec.number
        );
    }
    out
}

pub(super) fn build_figure_sections(specs: &[NarrativeFigureSpec]) -> String {
    let mut out = String::new();
    for spec in specs {
        let _ = writeln!(out, "<a id=\"fig-{}\"></a>\n", spec.number);
        let _ = writeln!(
            out,
            "<p align=\"center\"><img src=\"{}\" alt=\"{}\" style=\"max-width:100%;width:6.5in;\" /></p>\n",
            spec.path, spec.alt
        );
        let _ = writeln!(out, "*Figure {}. {}*\n", spec.number, spec.caption);
    }
    out
}

pub(super) fn build_storage_policy_section() -> String {
    "\
To ensure integrity, security, and accessibility of Milestone 12 computational data, SonALAsense applies policy-aligned storage controls across its managed engineering systems.\n\n\
**Data Governance and Access Control:** versioned simulation/report assets are maintained in the SonALAsense GitHub repository with controlled branch/merge workflows and auditable change history.\n\n\
**Secure Collaboration and Operational Storage:** project data packages and controlled shared artifacts are managed in SonALAsense Egnyte with role-based access and organization-managed retention policies.\n\n\
**Backup, Recovery, and Traceability:** canonical milestone outputs are regenerated deterministically from source code and preserved alongside run artifacts and trace matrices to support audit reproducibility.\n"
        .to_string()
}

pub(super) fn build_storage_artifact_index() -> String {
    "\
- Canonical results: Appendix A (embedded below)\n\
- Contract trace matrix: Appendix B (embedded below)\n\
- Figure manifest: Appendix C (embedded below)\n\
- Figures: inline throughout §5 Results\n\
- Generation artifacts: top-5 JSON, validation summaries, robustness outputs, and GA artifacts are included in Appendix A canonical data"
        .to_string()
}

fn pass_fail(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}
