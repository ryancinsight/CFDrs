//! Cross-fidelity validation traceability helpers for Milestone 12 reports.
//!
//! # Theorem
//! The validation summary table is a deterministic projection of the supplied
//! [`ValidationRow`] slice: every row appears exactly once, and no computed
//! value can depend on hidden state outside the provided diagnostics.
//!
//! **Proof sketch**
//! The builder performs a single forward pass over the row slice. Each row is
//! rendered once, and all aggregate values are derived only from the row
//! contents through min/max/count reductions. Therefore the output is a pure
//! function of the supplied diagnostics.

use std::fmt::Write as _;

use super::ValidationRow;

#[derive(Debug, Clone, Copy, PartialEq)]
struct ValidationSummary {
    row_count: usize,
    converged_count: usize,
    high_re_mismatch_count: usize,
    min_agreement_1d_2d_pct: f64,
    max_agreement_1d_2d_pct: f64,
    min_agreement_2d_3d_pct: f64,
    max_agreement_2d_3d_pct: f64,
}

fn validation_summary(validation_rows: &[ValidationRow]) -> Option<ValidationSummary> {
    let first = validation_rows.first()?;
    let mut summary = ValidationSummary {
        row_count: validation_rows.len(),
        converged_count: usize::from(first.two_d_converged),
        high_re_mismatch_count: usize::from(first.high_re_laminar_mismatch),
        min_agreement_1d_2d_pct: first.agreement_1d_2d_pct,
        max_agreement_1d_2d_pct: first.agreement_1d_2d_pct,
        min_agreement_2d_3d_pct: first.agreement_2d_3d_pct,
        max_agreement_2d_3d_pct: first.agreement_2d_3d_pct,
    };

    for row in validation_rows.iter().skip(1) {
        if row.two_d_converged {
            summary.converged_count += 1;
        }
        if row.high_re_laminar_mismatch {
            summary.high_re_mismatch_count += 1;
        }
        summary.min_agreement_1d_2d_pct =
            summary.min_agreement_1d_2d_pct.min(row.agreement_1d_2d_pct);
        summary.max_agreement_1d_2d_pct =
            summary.max_agreement_1d_2d_pct.max(row.agreement_1d_2d_pct);
        summary.min_agreement_2d_3d_pct =
            summary.min_agreement_2d_3d_pct.min(row.agreement_2d_3d_pct);
        summary.max_agreement_2d_3d_pct =
            summary.max_agreement_2d_3d_pct.max(row.agreement_2d_3d_pct);
    }

    Some(summary)
}

fn pass_fail(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

/// Render the traceability note for the cross-fidelity validation section.
pub(crate) fn build_validation_traceability_intro(
    canonical_source: &str,
    figure_count: Option<usize>,
    figure_manifest_path: &str,
    asset_review_manifest_path: &str,
    validation_rows_path: &str,
    validation_rows: &[ValidationRow],
) -> String {
    let Some(summary) = validation_summary(validation_rows) else {
        return format!(
            "The canonical report source is `{canonical_source}`. The figure manifest at `{figure_manifest_path}` and the asset-review manifest at `{asset_review_manifest_path}` remain the authoritative artifact registers for the report package. No cross-fidelity validation rows were emitted for this run, so the validation subsection is intentionally empty rather than fabricated."
        );
    };

    let figure_manifest_clause = match figure_count {
        Some(count) => format!("The current figure manifest at `{figure_manifest_path}` enumerates **{count}** figures."),
        None => format!(
            "The current figure manifest at `{figure_manifest_path}` enumerates the report figures."
        ),
    };

    format!(
        "The canonical report source is `{canonical_source}`. {figure_manifest_clause} The asset-review manifest at `{asset_review_manifest_path}` records the review state for the generated narrative and figures, and `{validation_rows_path}` stores the cross-fidelity venturi rows used below. Across **{row_count}** validation rows, 1D→2D agreement spans **{min_1d_2d:.3}%–{max_1d_2d:.3}%**, 2D→3D agreement spans **{min_2d_3d:.3}%–{max_2d_3d:.3}%**, and the 2D solver converges on **{converged_count}/{row_count}** rows. The `mass_error_3d_pct` field is reproduced verbatim from the solver diagnostics and is not renormalized in the report.",
        row_count = summary.row_count,
        min_1d_2d = summary.min_agreement_1d_2d_pct,
        max_1d_2d = summary.max_agreement_1d_2d_pct,
        min_2d_3d = summary.min_agreement_2d_3d_pct,
        max_2d_3d = summary.max_agreement_2d_3d_pct,
        converged_count = summary.converged_count,
    )
}

/// Build the cross-fidelity validation table.
pub(crate) fn build_validation_summary_table(validation_rows: &[ValidationRow]) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | Topology | 1D→2D agreement (%) | 2D→3D agreement (%) | 3D mass error (%) | σ_1d | σ_2d | 2D converged | High-Re mismatch |\n");
    out.push_str("|---|---|---|---:|---:|---:|---:|---:|---|---|\n");

    if validation_rows.is_empty() {
        out.push_str("| _No validation rows emitted_ | _No validation rows emitted_ | _No validation rows emitted_ | n/a | n/a | n/a | n/a | n/a | n/a | n/a |\n");
        return out;
    }

    for row in validation_rows {
        let _ = writeln!(
            out,
            "| {} | `{}` | {} | {:.3} | {:.3} | {:.3} | {:.4} | {:.4} | {} | {} |",
            row.track,
            row.id,
            row.topology,
            row.agreement_1d_2d_pct,
            row.agreement_2d_3d_pct,
            row.mass_error_3d_pct,
            row.sigma_1d,
            row.sigma_2d,
            pass_fail(row.two_d_converged),
            pass_fail(row.high_re_laminar_mismatch)
        );
    }

    out
}

#[cfg(test)]
mod tests {
    use super::{build_validation_summary_table, build_validation_traceability_intro};
    use crate::reporting::ValidationRow;

    fn sample_rows() -> Vec<ValidationRow> {
        vec![
            ValidationRow {
                track: "Option 2 Combined Selective Venturi".to_string(),
                id: "9013012-PST-PentaTri-X-pcf330-tcf450-btf680-vt1-q300ml-g300kPa-dt45um-tl10-w1500um-h1000-n5".to_string(),
                topology: "PTV".to_string(),
                k_loss_1d: 2_010_849.8626242229,
                dp_1d_bernoulli_pa: 20_905_527.14365864,
                k_loss_2d: 17_905.791782043456,
                dp_2d_fvm_pa: 186_155.12927439384,
                k_loss_3d: 261.5663519097687,
                dp_3d_fem_pa: 2_719.3390075285315,
                agreement_1d_2d_pct: 99.1095410893245,
                agreement_2d_3d_pct: 98.53920812274788,
                mass_error_3d_pct: 99.9884485379281,
                sigma_1d: -0.8561223008286081,
                sigma_2d: 18.08333023458977,
                score: 0.7412766731882545,
                two_d_converged: true,
                high_re_laminar_mismatch: true,
            },
            ValidationRow {
                track: "GA HydroSDT Venturi".to_string(),
                id: "9017012-PST-PentaTri-X-pcf330-tcf450-btf680-vt2-q200ml-g50kPa-dt35um-tl10-w1500um-h1000-n5-ga-dean-seed-ga-ssd-stage_0_center-x-blend-refine-q85-p115-refine-q85-p115-refine-q85-p85-refine-q85-p95-refine-q85-p85-refine-q85-p85".to_string(),
                topology: "PTV".to_string(),
                k_loss_1d: 10_128_448.033079874,
                dp_1d_bernoulli_pa: 6_655_521.7735658465,
                k_loss_2d: 27_359.314775358067,
                dp_2d_fvm_pa: 17_978.126027059945,
                k_loss_3d: 1_044.542879907451,
                dp_3d_fem_pa: 686.3813545709871,
                agreement_1d_2d_pct: 99.72987653502292,
                agreement_2d_3d_pct: 96.18213069850621,
                mass_error_3d_pct: 99.9962273390378,
                sigma_1d: -1.0427830973550918,
                sigma_2d: 87.82608005894014,
                score: 0.8416012846182056,
                two_d_converged: true,
                high_re_laminar_mismatch: false,
            },
        ]
    }

    #[test]
    fn validation_intro_reports_traceability_and_spans() {
        let rows = sample_rows();
        let intro = build_validation_traceability_intro(
            "report/milestone12_results.md",
            Some(14),
            "report/milestone12/figure_manifest.json",
            "report/milestone12/asset_review_manifest.json",
            "report/milestone12/milestone12_validation_rows.json",
            &rows,
        );

        assert!(intro.contains("report/milestone12_results.md"));
        assert!(intro.contains("figure manifest"));
        assert!(intro.contains("14"));
        assert!(intro.contains("99.110%–99.730%"));
        assert!(intro.contains("96.182%–98.539%"));
        assert!(intro.contains("**2/2** rows"));
        assert!(intro.contains("mass_error_3d_pct"));
    }

    #[test]
    fn validation_table_renders_all_rows() {
        let rows = sample_rows();
        let table = build_validation_summary_table(&rows);

        assert!(table.contains("Option 2 Combined Selective Venturi"));
        assert!(table.contains("GA HydroSDT Venturi"));
        assert!(table.contains("PASS"));
        assert!(table.contains("FAIL"));
        assert!(table.contains("99.730"));
        assert!(table.contains("99.996"));
    }
}
