//! Section builders for Milestone 12 narrative markdown.

use std::fmt::Write as _;

use crate::analysis::RobustnessReport;
use crate::constraints::{EXPANSION_RATIO_LOW_RISK, VENTURI_EXPANSION_RATIO_HIGH_RISK};
use crate::reporting::ranking::oncology_priority_score;
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

/// Format a boolean gate as a PASS/FAIL string for markdown tables.
fn pass_fail(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

/// Pediatric 3 kg neonatal ECV as a percentage of the 25.5 mL limit.
fn pediatric_limit_pct(ecv_ml: f64) -> f64 {
    100.0 * ecv_ml / 25.5_f64
}

fn cavitation_regime_clause(sigma: f64) -> String {
    if sigma < 0.0 {
        format!(
            "σ = {sigma:.4} < 0 places the throat below the blood vapor-pressure reference, indicating an actively cavitating hydrodynamic regime"
        )
    } else if sigma < 1.0 {
        format!(
            "σ = {sigma:.4} lies in the inception-capable hydrodynamic cavitation window (0 < σ < 1), so the throat is cavitation-eligible but not yet below vapor pressure"
        )
    } else {
        format!(
            "σ = {sigma:.4} is above the hydrodynamic cavitation threshold, so the venturi is not operating in a cavitation-capable regime"
        )
    }
}

pub(super) fn build_selected_table(
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) -> String {
    let mut out = String::new();
    out.push_str(
        "| Track | Candidate | Topology | Mode | Active venturi throats | Score | sigma | Cumulative cavitation dose | WBC treatment exposure | RBC treatment exposure | HI/pass | P95 shear (Pa) | ECV (mL) | ECV / 3kg limit (%) |\n",
    );
    out.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    if let Some(option1) = option1 {
        let _ = writeln!(
            out,
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
                pediatric_limit_pct(option1.metrics.total_ecv_ml)
        );
    } else {
        out.push_str("| Option 1 (Selective acoustic center treatment) | _No eligible design under current physics regime_ | n/a | Unavailable | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |\n");
    }
    let _ = writeln!(
        out,
            "| Option 2 (Selective venturi cavitation selectivity) | `{}` | {} | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2} | {:.3} | {:.1} |",
        option2.candidate.id,
        option2.topology_display_name(),
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.score,
        option2.metrics.cavitation_number,
        option2.metrics.serial_cavitation_dose_fraction,
        option2.metrics.wbc_recovery,
        option2.metrics.rbc_venturi_exposure_fraction,
        option2.metrics.hemolysis_index_per_pass,
        option2.metrics.wall_shear_p95_pa,
            option2.metrics.total_ecv_ml,
            pediatric_limit_pct(option2.metrics.total_ecv_ml)
    );
    out
}

pub(super) fn build_option2_top5_table(option2_ranked: &[Milestone12ReportDesign]) -> String {
    let mut out = String::new();
    out.push_str("| Rank | Candidate | Mode | Active venturi throats | Score | Oncology priority | RBC venturi exposure | Clot risk | sigma | Cumulative cavitation dose |\n");
    out.push_str("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|\n");
    for d in option2_ranked.iter().take(5) {
        let _ = writeln!(
            out,
            "| {} | `{}` | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} |",
            d.rank,
            d.candidate.id,
            d.metrics.treatment_zone_mode,
            d.metrics.active_venturi_throat_count,
            d.score,
            oncology_priority_score(&d.metrics),
            d.metrics.rbc_venturi_exposure_fraction,
            d.metrics.clotting_risk_index,
            d.metrics.cavitation_number,
            d.metrics.serial_cavitation_dose_fraction
        );
    }
    out
}

pub(super) fn build_tri_cell_table(
    option2_ranked: &[Milestone12ReportDesign],
    ga_ranked: &[Milestone12ReportDesign],
) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | cancer_center_fraction | wbc_treatment_exposure | rbc_treatment_exposure | three_pop_sep_efficiency |\n");
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
    option2_ranked: &[Milestone12ReportDesign],
    ga_ranked: &[Milestone12ReportDesign],
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
        let _ = fast_mode;
        return "\
### Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)\n\n\
*Robustness screening was not computed in this run (fast mode). \
Re-run without `M12_FAST=1` to populate.*\n"
            .to_string();
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
        let _ = fast_mode;
        return "\
### Multi-Fidelity Pressure-Drop Validation (Selected Venturi Designs)\n\n\
*Multi-fidelity validation is produced by the companion example \
`cargo run -p cfd-optim --example milestone12_validation --no-default-features`. \
This report run did not embed 2D/3D validation rows.*\n"
            .to_string();
    }
    let _ = fast_mode;
    let mut out = String::new();
    out.push_str("### Multi-Fidelity Pressure-Drop Validation (Selected Venturi Designs)\n\n");
    let unconverged_count = rows.iter().filter(|row| !row.two_d_converged).count();
    let inertial_count = rows
        .iter()
        .filter(|row| row.high_re_stokes_mismatch)
        .count();
    if unconverged_count > 0 || inertial_count > 0 {
        let _ = writeln!(
            out,
            "These rows are diagnostic rather than confirmatory wherever the companion solver falls outside its validity range. Unconverged 2D rows: **{}**. Inertial rows where the 3D Stokes model is not valid: **{}**.\n",
            unconverged_count,
            inertial_count
        );
    } else {
        out.push_str(
            "All embedded rows lie within the nominal validity envelope of the companion solvers and may be interpreted as confirmation-grade cross-fidelity checks.\n\n",
        );
    }
    out.push_str("| Track | Candidate | dp1D Bernoulli+K+f (Pa) | dp2D FVM (Pa) | dp3D FEM (Pa) | 1D\u{2013}2D diff (%) | 2D\u{2013}3D diff (%) | Mass error (%) | 2D Conv? | Re regime |\n");
    out.push_str("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |\n");
    for row in rows {
        let conv_flag = if row.two_d_converged { "YES" } else { "NO" };
        let re_label = if row.high_re_stokes_mismatch {
            "inertial (3D Stokes invalid)"
        } else {
            "creeping (Stokes valid)"
        };
        let _ = writeln!(
            out,
            "| {} | `{}` | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} | {} | {} |",
            row.track,
            row.id,
            row.dp_1d_bernoulli_pa,
            row.dp_2d_fvm_pa,
            row.dp_3d_fem_pa,
            row.agreement_1d_2d_pct,
            row.agreement_2d_3d_pct,
            row.mass_error_3d_pct,
            conv_flag,
            re_label
        );
    }

    // Add explanatory notes when any solver is outside its validity range.
    let any_high_re = rows.iter().any(|r| r.high_re_stokes_mismatch);
    let any_2d_nonconv = rows.iter().any(|r| !r.two_d_converged);
    if any_high_re || any_2d_nonconv {
        out.push_str("\n**Validation Notes:**\n\n");
    }
    if any_high_re {
        out.push_str(
            "- **Re regime \u{2014} inertial**: The 3D Taylor\u{2013}Hood Stokes FEM is skipped for \
Re\u{2009}>\u{2009}50 (the solver omits \u{03c1}(u\u{00b7}\u{2207})u and cannot represent inertia-dominated \
flows). The 1D model includes Bernoulli dynamic pressure, Idelchik contraction / \
Borda\u{2013}Carnot expansion loss coefficients, and Darcy\u{2013}Weisbach pipe friction in the throat \
(f\u{2009}=\u{2009}64/Re laminar, Haaland turbulent). For high contraction-ratio venturis where the \
throat length-to-hydraulic-diameter ratio L/D\u{2095}\u{2009}\u{226b}\u{2009}1, viscous friction dominates.\n",
        );
    }
    if any_2d_nonconv {
        out.push_str(
            "- **2D FVM \u{2014} skipped/non-converged**: For contraction ratios CR\u{2009}>\u{2009}80, the 2D \
Cartesian immersed-boundary solver cannot resolve the throat (sub-cell geometry). \
In these cases dp2D reports the 1D estimate. For moderate CR designs that run the 2D \
SIMPLE solver, the velocity field is initialised from 1D mass-continuity (parabolic profile) \
which gives ~33% agreement with 1D for converged cases.\n",
        );
    }
    out
}

pub(super) fn build_limits_of_usage(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let thermal_line = format!(
        "- Throat viscous heating: **{} K** rise (FDA 5 K / 42 °C ceiling): **{}**",
        if m.throat_temperature_rise_k < 0.01 && m.throat_temperature_rise_k > 0.0 {
            "< 0.01".to_string()
        } else {
            format!("{:.2}", m.throat_temperature_rise_k)
        },
        pass_fail(m.fda_thermal_compliant)
    );
    format!(
        "- Operating point: **{:.1} mL/min**, **{:.0} kPa gauge**\n\
         - Venturi throat: **{:.0} µm**, throat length **{:.0} µm**\n\
         - Treatment mode: **{}**, active venturi throats: **{}** (serial stages per path: **{}**)\n\
         - FDA main-channel compliance: **{}**\n\
         {}\n\
         - Clot-risk indicators: `clotting_risk_index={:.4}`, `clotting_risk_index_10ml_s={:.4}`\n\
         - Flow caution flags: `Q>=200={}`, `Q>=600={}`",
        option2.flow_rate_ml_min(),
        option2.inlet_gauge_kpa(),
        option2.throat_width_um().unwrap_or(0.0),
        option2.throat_length_um().unwrap_or(0.0),
        m.treatment_zone_mode,
        m.active_venturi_throat_count,
        m.serial_venturi_stages_per_path,
        pass_fail(m.fda_main_compliant),
        thermal_line,
        m.clotting_risk_index,
        m.clotting_risk_index_10ml_s,
        pass_fail(m.clotting_flow_compliant),
        pass_fail(m.clotting_flow_compliant_10ml_s)
    )
}

// ── Narrative paragraph builders ─────────────────────────────────────────────

pub(super) fn build_results_intro(
    total_candidates: usize,
    opt1_pool: usize,
    opt2_pool: usize,
) -> String {
    format!(
        "From {total_candidates} total candidates generated from canonical Milestone 12 split-sequence scaffolds, \
strict eligibility gating produced {opt1_pool} Option 1 qualified designs \
(AsymmetricSplitResidenceSeparation track) and {opt2_pool} Option 2 qualified designs \
(AsymmetricSplitVenturiCavitationSelectivity track). The ranked pool therefore reflects only physically admissible \
designs that preserve selective split-width partitioning, treatment-lane residence time, and \
healthy-cell protection. The following sub-sections present the selected designs (§5.1), gate \
evidence (§5.2), robustness and multi-fidelity validation (§5.3), design visualizations (§5.4), \
derived metric formulas (§5.5), and operating limits (§5.6). Extracorporeal circuit volume is \
reported explicitly as ECV = Σ(L_i A_i) = Q_in t_res, and each selected design is benchmarked \
against the 3 kg neonatal reference limit of 25.5 mL (10% of 3 × 85 mL blood volume).",
    )
}

pub(super) fn build_selected_table_intro(
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) -> String {
    let m2 = &option2.metrics;
    let ccf_pct = m2.cancer_center_fraction * 100.0;
    let wbc2_exposure_pct = m2.wbc_recovery * 100.0;
    let option1_summary = if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        format!(
            "**Option 1 — Selective Acoustic Center Treatment** \
(`{}`, score {:.4}): branch-diameter biasing keeps {:.0}% of WBCs out of the treatment lane \
per pass while concentrating CTC-rich flow for ultrasound exposure; HI/pass = {:.4}% \
(FDA 0.1% non-therapeutic limit); ECV = {:.3} mL \
({:.1}% of the 3 kg neonatal circuit-volume limit). \
This design carries zero active venturi throats — treatment relies entirely on externally \
applied 412 kHz ultrasound acting on cells concentrated in the center lane.",
            option1.candidate.id,
            option1.score,
            (1.0 - m1.wbc_recovery) * 100.0,
            m1.hemolysis_index_per_pass,
            m1.total_ecv_ml,
            pediatric_limit_pct(m1.total_ecv_ml),
        )
    } else {
        "**Option 1 — Selective Acoustic Center Treatment**: no design satisfied the strict acoustic eligibility gates under the current physics regime, so the report records Option 1 explicitly as an empty shortlist.".to_string()
    };
    let cavitation_clause = cavitation_regime_clause(m2.cavitation_number);
    format!(
        "{}\n\n\
**Option 2 — Selective Venturi Hydrodynamic Cavitation** \
(`{}`, score {:.4}): {} at {} serial venturi throat(s) per path; {:.0}% of CTCs route through the venturi \
treatment lane (cancer_center_fraction); WBC treatment exposure = {:.0}%; therapeutic window \
score = {:.3}; HI/pass = {:.4}%; ECV = {:.3} mL ({:.1}% of the same neonatal limit). \
**Note:** Option 1 and Option 2 scores use different objective functions \
(AsymmetricSplitResidenceSeparation vs AsymmetricSplitVenturiCavitationSelectivity) and are not comparable across tracks. \
This report is emitted from structured canonical data and should be regenerated rather than edited manually.",
        option1_summary,
        option2.candidate.id,
        option2.score,
        cavitation_clause,
        m2.serial_venturi_stages_per_path,
        ccf_pct,
        wbc2_exposure_pct,
        m2.therapeutic_window_score,
        m2.hemolysis_index_per_pass,
        m2.total_ecv_ml,
        pediatric_limit_pct(m2.total_ecv_ml),
    )
}

pub(super) fn build_top5_intro(option2_ranked: &[Milestone12ReportDesign]) -> String {
    let n = option2_ranked.len().min(5);
    if n < 2 {
        return "The top-ranked Option 2 design was selected by deterministic tie-break \
(score desc, oncology-priority desc, RBC venturi exposure asc, clot risk asc, candidate id asc)."
            .to_string();
    }
    let score_spread =
        if let (Some(top), Some(rank5)) = (option2_ranked.first(), option2_ranked.get(n - 1)) {
            (top.score - rank5.score).abs()
        } else {
            0.0
        };

    // When all top-N scores are identical, the narrative must honestly report
    // that ranking is determined by tie-breaks, not score differentiation.
    if score_spread < 1.0e-12 {
        let score_display = option2_ranked.first().map_or(0.0, |d| d.score);
        return format!(
            "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, \
clot risk asc, candidate id asc. All {n} top-ranked designs share the same composite score \
({score_display:.2e}); ranking among them is determined entirely by the deterministic \
tie-break sequence. Differentiation arises from non-scored parameters (gauge pressure, \
serial throat count) rather than from the objective function."
        );
    }

    format!(
        "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, \
clot risk asc, candidate id asc. Ranks 1–{n} span a score range of {score_spread:.2e}, \
demonstrating that the shortlist is resolved by measurable physics rather than by a zero-score \
collapse. Adjacent high-ranked designs differ in asymmetric daughter-width partitioning, \
pressure budget, venturi stage count, and throat geometry; the top-ranked lineages therefore \
represent a family of nearby but not redundant hydrodynamic operating points."
    )
}

pub(super) fn build_tri_cell_intro(
    option2: &Milestone12ReportDesign,
    ga_best: Option<&Milestone12ReportDesign>,
) -> String {
    let m = &option2.metrics;
    let ccf_pct = m.cancer_center_fraction * 100.0;
    let wbc_pct = m.wbc_recovery * 100.0;
    let rbc_pct = m.rbc_venturi_exposure_fraction * 100.0;
    let mut s = format!(
        "Three-population routing provides the clinical selectivity basis for Option 2. \
Cell equilibrium positions are computed from inertial lift correlations; \
the Amini (2014) confinement-dependent lift correction refines the standard Di Carlo \
model for particle-to-channel size ratios a/D_h > 0.1 typical of CTCs in millifluidic \
geometries. Blood rheology incorporates the Fahraeus–Lindqvist effect (Pries 1992) for \
apparent viscosity reduction in channels below ~300 µm, plasma skimming (Pries 1989) for \
hematocrit partitioning at asymmetric bifurcations, and the Quemada (1978) rouleaux \
aggregation viscosity model for low-shear zones where RBC aggregation elevates the \
effective viscosity above the Casson baseline. \
For the selected design: {ccf_pct:.0}% of CTCs enter the venturi treatment lane \
(cancer_center_fraction); {wbc_pct:.0}% of WBCs enter the treatment lane, so {:.0}% are kept out \
of the treatment region before all flowpaths remerge upstream of the outlet; {rbc_pct:.0}% of RBCs \
pass through venturi throats (rbc_venturi_exposure) — reducing this fraction is the primary \
hemocompatibility lever. \
Composite three-population separation efficiency = {:.4}.",
        (1.0 - m.wbc_recovery) * 100.0,
        m.three_pop_sep_efficiency
    );
    if let Some(ga) = ga_best {
        let gm = &ga.metrics;
        let _ = write!(
            s,
            " GA rank-1 (`{}`) achieves cancer_center_fraction={:.4}, \
wbc_treatment_exposure={:.4}, rbc_venturi_exposure={:.4} — a narrower, more cavitation-intensive \
configuration than the parametric Option 2 selection.",
            ga.candidate.id,
            gm.cancer_center_fraction,
            gm.wbc_recovery,
            gm.rbc_venturi_exposure_fraction
        );
    }
    s
}

pub(super) fn build_strict_core_intro() -> String {
    "These five binary gates represent hard physical constraints derived from FDA guidance, \
ANSI/SLAS plate standards, and Hagen–Poiseuille network feasibility. A design failing any \
gate receives the infeasibility floor score (0.001) and is excluded from the eligible pool \
regardless of objective metrics.\n\n\
1. **Pressure feasible** — total pressure drop Δp must not exceed the inlet gauge budget. \
For the 1D lumped-element Hagen–Poiseuille model, Δp = Σ (128 μ L Q / π D_h⁴) across all \
channel segments (Casson-corrected blood viscosity μ ≈ 3.5 mPa·s at γ̇ > 100 s⁻¹). \
This ensures the device can be driven by clinically available extracorporeal pump heads.\n\n\
2. **Plate fits** — chip body dimensions ≤ 127.76 × 85.47 mm (ANSI/SLAS 1-2004, 96-well \
format), ensuring compatibility with standard fluorescence microscopy stages and optical \
delivery windows for sonosensitiser activation.\n\n\
3. **FDA main-channel compliance** — sustained wall shear stress τ_w = μ γ̇_w must not \
exceed 150 Pa in any non-venturi channel. This threshold is derived from the Giersiepen \
(1990) hemolysis correlation at exposure times typical of millifluidic transit (10–500 ms); \
the Taskin (2012) strain-based hemolysis model is available as an alternative predictor \
for designs where cumulative strain history dominates over instantaneous shear. Venturi \
throats are excluded from this gate and evaluated separately under the transient shear \
exception (≤ 300 Pa for ≤ 5 ms transit).\n\n\
4. **σ finite** — the Bernoulli cavitation number σ = (p∞ − pᵥ)/(½ρv²) must converge \
to a finite real value, confirming the 1D network solver produced a physically meaningful \
flow field (no zero-flow or infinite-velocity singularities).\n\n\
5. **σ < 1** — confirms entry into the incipient-cavitation window at the venturi throat. \
At 0 < σ < 1, the throat is cavitation-capable and susceptible to vapor/gas nucleus growth; \
σ < 0 is the stronger regime in which the local static pressure at the vena contracta drops \
below the vapour pressure of blood (pᵥ ≈ 6.3 kPa at 37 °C). Designs with σ ≥ 1 have no hydrodynamic \
cavitation capability and must rely on externally applied acoustic energy alone for \
sonosensitiser activation."
        .to_string()
}

pub(super) fn build_cavitation_formulas_intro(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let mut out = format!(
        "These derived metrics translate raw 1D Hagen–Poiseuille network physics (pressure drop, \
wall shear rate, flow partition fractions) into clinically interpretable quantities.\n\n\
**Cavitation number** σ = (p∞ − pᵥ) / (½ρv²_throat), where p∞ is the far-field \
(inlet) absolute pressure, pᵥ ≈ 6.3 kPa (blood vapor pressure at 37 °C), and \
v_throat is the cross-section averaged throat velocity (White 2011). For short venturi \
throats where L/D_h < 20, the Durst (2005) developing-flow entrance correction is \
available to account for the non-parabolic velocity profile in the entrance region, \
which increases the effective centreline velocity and wall shear relative to the \
fully-developed assumption. σ < 1 indicates incipient hydrodynamic cavitation; \
σ < 0 implies p_throat < pᵥ (strong cavitation). \
Cavitation potential C_p = max(0, 1 − σ).\n\n\
**Cavitation intensity** I_cav = C_p × (0.5 + 0.5 × κ), where κ is the constriction \
score — the logarithmic velocity ratio ln(v_throat/v_upstream) normalised by ln(v_ref). \
This couples bubble nucleation probability (C_p) with collapse violence (proportional \
to kinetic energy at the vena contracta).\n\n\
**Tumor-targeted cavitation index** TTCI = f_cancer,center × I_cav, where \
f_cancer,center is the cancer center-lane enrichment fraction from Zweifach–Fung \
flow partitioning at asymmetric bifurcations.\n\n\
**Therapeutic window score** TWS = clamp(TTCI / (10⁻⁶ + L_risk) / L_ref, 0, 1), \
where L_risk = HI × (1 + 5 × f_rbc,venturi × H_local) is the lysis risk index \
(HI = Giersiepen hemolysis index per pass, H_local = local hematocrit at the throat). \
TWS → 1 when cancer-targeted cavitation is high relative to blood damage.\n\n\
**Sonoluminescence proxy** S = clamp(C_p × (p_abs/p_vap)^((κ_poly−1)/κ_poly) / S_ref, 0, 1), \
modelling the adiabatic collapse temperature ratio for sonosensitiser (5-ALA, Ce6) \
photoactivation under polytropic bubble dynamics (κ_poly = 1.4).\n\n\
For the selected Option 2 design: TTCI = {:.4}, TWS = {:.4}, \
lysis_risk_index = {:.4e}, sonoluminescence_proxy = {:.4}, \
oncology_selectivity_index = {:.4}, cancer_rbc_cavitation_bias = {:.4}.",
        m.cancer_targeted_cavitation,
        m.therapeutic_window_score,
        m.lysis_risk_index,
        m.sonoluminescence_proxy,
        m.oncology_selectivity_index,
        m.cancer_rbc_cavitation_bias_index,
    );
    if m.cavitation_number < 0.0 {
        let _ = write!(
            out,
            "\n\n**Negative σ interpretation:** The selected Option 2 design has σ = {:.4}, \
confirming that the local static pressure at the vena contracta drops below the vapor \
pressure of blood (p_v ≈ 6.3 kPa at 37 °C). Negative σ is the desired operating regime \
for Option 2 — it indicates active hydrodynamic cavitation at the venturi throat without \
requiring external acoustic energy input.",
            m.cavitation_number,
        );
    }
    out
}

pub(super) fn build_limits_of_usage_intro(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let q_ml_min = option2.flow_rate_ml_min();
    let gauge_kpa = option2.inlet_gauge_kpa();
    let q_ml_s = q_ml_min / 60.0;
    let q600_text = if m.clotting_flow_compliant_10ml_s {
        format!(
            "At {q_ml_min:.1} mL/min ({q_ml_s:.2} mL/s) the operating flow rate exceeds \
the conservative 10 mL/s (600 mL/min) threshold; low-flow stasis-clotting risk is minimal."
        )
    } else {
        format!(
            "At {q_ml_min:.1} mL/min ({q_ml_s:.2} mL/s) the operating flow rate is below \
the conservative 10 mL/s (600 mL/min) threshold (`Q>=600=FAIL`). Anticoagulation protocol \
adjustments should be evaluated per institutional guidelines for moderate-flow \
extracorporeal circuits."
        )
    };

    // Use < 0.01 for sub-millikelvin heating to avoid misleading "0.00 K" or confusing scientific notation.
    let temp_str = if m.throat_temperature_rise_k < 0.01 && m.throat_temperature_rise_k > 0.0 {
        "< 0.01".to_string()
    } else {
        format!("{:.2}", m.throat_temperature_rise_k)
    };
    let thermal_verdict = if m.fda_thermal_compliant {
        "well within"
    } else {
        "**exceeding**"
    };
    let cavitation_clause = cavitation_regime_clause(m.cavitation_number);

    format!(
        "Limits are derived from the selected Option 2 operating point ({q_ml_min:.1} mL/min, \
{gauge_kpa:.0} kPa gauge). {q600_text} At this operating point, {cavitation_clause}. Viscous heating in the \
{:.0} µm venturi throat is {temp_str} K — {thermal_verdict}: the FDA 42 °C \
(37 + 5 K) temperature ceiling.",
        option2.throat_width_um().unwrap_or(0.0),
    )
}

/// Build the CRI expansion-sensitivity table for the PST topology parameter space.
///
/// # Theorem (Expansion Stasis Risk)
/// For a channel–throat pair with geometric expansion ratio `r = channel_width / throat_diameter`,
/// the per-stage stasis risk follows a log-linear model:
/// `p_stage = clamp((ln r − ln r_lo) / (ln r_hi − ln r_lo), 0, 1)`
/// where `r_lo = EXPANSION_RATIO_LOW_RISK = 4` and `r_hi = VENTURI_EXPANSION_RATIO_HIGH_RISK = 10 000`.
///
/// The cumulative risk after `n` independent stages is
/// `P(n) = 1 − (1 − p_stage)^n`,
/// and the expansion contribution to `clotting_risk_index` is `0.20 × P(n)`.
///
/// The table spans all `SELECTIVE_TREE_WIDTHS × SELECTIVE_TREE_THROATS` combinations
/// and reports vt1 (n=1) and vt2 (n=2) columns so design engineers can read the
/// hemocompatibility cost directly from the chosen width/throat pairing.
pub(super) fn build_cri_expansion_sensitivity() -> String {
    const WIDTHS_MM: [f64; 3] = [4.0, 6.0, 8.0];
    const THROATS_UM: [f64; 6] = [35.0, 45.0, 55.0, 75.0, 100.0, 120.0];
    const CRI_WEIGHT: f64 = 0.20;

    let ln_lo = EXPANSION_RATIO_LOW_RISK.ln();
    let ln_hi = VENTURI_EXPANSION_RATIO_HIGH_RISK.ln();
    let range = ln_hi - ln_lo;

    let mut out = String::new();
    out.push_str("<div align=\"center\">\n");
    out.push_str("<table>\n");
    out.push_str("<caption><strong>Table 5.7.</strong> Expansion-driven clot-risk contribution across the swept Milestone 12 channel-width / venturi-throat design space. Lower <code>CRI_exp</code> values are preferred.</caption>\n");
    out.push_str("<thead>\n<tr><th>Channel width (mm)</th><th>Throat Ø (µm)</th><th>w/d ratio</th><th>p_stage</th><th>P(vt1)</th><th>P(vt2)</th><th>CRI_exp (vt1)</th><th>CRI_exp (vt2)</th></tr>\n</thead>\n<tbody>\n");

    for &w_mm in &WIDTHS_MM {
        for &d_um in &THROATS_UM {
            let ratio = (w_mm * 1e-3) / (d_um * 1e-6);
            let p_stage = ((ratio.ln() - ln_lo) / range).clamp(0.0, 1.0);
            let p_vt1 = p_stage;
            let p_vt2 = 1.0 - (1.0 - p_stage).powi(2);
            let _ = writeln!(
                out,
                "<tr><td>{w_mm:.0}</td><td>{d_um:.0}</td><td>{ratio:.1}</td><td>{p_stage:.3}</td><td>{p_vt1:.3}</td><td>{p_vt2:.3}</td><td>{:.4}</td><td>{:.4}</td></tr>",
                CRI_WEIGHT * p_vt1,
                CRI_WEIGHT * p_vt2,
            );
        }
    }
    out.push_str("</tbody>\n</table>\n</div>");
    out
}

/// Generate narrative section describing the physics models available for
/// multi-fidelity simulation.
///
/// Lists all validated physics models and their literature references
/// in a markdown table format suitable for inclusion in the Milestone 12 report.
pub(super) fn physics_model_inventory() -> String {
    "\
### Physics Model Inventory

The following table lists all validated physics models integrated or available \
in the CFDrs multi-fidelity simulation stack, grouped by physical domain. Models \
marked **Active** are used in the current 1D scoring pipeline; models marked \
**Available** are implemented and validated but reserved for higher-fidelity \
2D/3D refinement or future integration.

| # | Domain | Model | Reference | Status | Key parameter affected |
|---:|---|---|---|---|---|
| 1 | Flow | Hagen-Poiseuille network solver | Sutera & Skalak (1993) | Active | `total_pressure_drop_pa` |
| 2 | Flow | Casson viscosity (blood) | Casson (1959), Merrill (1969) | Active | `BLOOD_VISCOSITY_PA_S` |
| 3 | Flow | Fahraeus-Lindqvist apparent viscosity | Pries et al. (1992) | Available | Effective resistance in D_h < 300 um |
| 4 | Flow | Quemada rouleaux aggregation viscosity | Quemada (1978) | Available | Low-shear viscosity, `local_hematocrit_venturi` |
| 5 | Flow | Developing-flow entrance correction | Durst et al. (2005) | Available | Wall shear in short venturi throats (L/D_h < 20) |
| 6 | Flow | Cross-section averaged throat velocity | White (2011) | Active | `cavitation_number`, v_throat |
| 7 | Separation | Inertial lift focusing (Di Carlo) | Di Carlo (2009) | Active | `cancer_center_fraction`, equilibrium positions |
| 8 | Separation | Confinement-dependent lift correction | Amini et al. (2014) | Available | Lift for a/D_h > 0.1 (CTC-scale particles) |
| 9 | Separation | Zweifach-Fung flow partitioning | Zweifach & Fung (1971) | Active | Branch flow fractions at bifurcations |
| 10 | Separation | Plasma skimming hematocrit partitioning | Pries et al. (1989) | Available | `rbc_venturi_exposure_fraction` |
| 11 | Cavitation | Bernoulli cavitation number | Brennen (1995) | Active | `cavitation_number` (sigma) |
| 12 | Cavitation | Polytropic bubble collapse | Brennen (1995), Rayleigh-Plesset | Active | `sonoluminescence_proxy` |
| 13 | Cavitation | Sonoporation lysis amplification | Ohl et al. (2006) | Active | `LYSIS_CAVITATION_AMPLIFICATION` (5x) |
| 14 | Hemolysis | Giersiepen power-law HI | Giersiepen et al. (1990) | Active | `hemolysis_index_per_pass` |
| 15 | Hemolysis | Taskin strain-based HI | Taskin et al. (2012) | Available | Alternative `hemolysis_index_per_pass` |
| 16 | Hemolysis | Hellums PAI power-law | Hellums (1994) | Active | `platelet_activation_index` |
| 17 | Hemolysis | Cavitation-amplified HI | Ohl (2006) + Giersiepen (1990) | Active | `hemolysis_index_per_pass_cavitation_amplified` |
| 18 | Secondary flow | Dean vortex secondary flow | Dean (1927) | Active | Dean number in GA serpentine scoring |
| 19 | Secondary flow | Millifluidic Dean correlation | Bayat-Rezai (2017) | Active | Dean number for D_h > 500 um channels |
| 20 | Thermal | Viscous dissipation heating | First-law energy balance | Active | `throat_temperature_rise_k` |
| 21 | Acoustic | 412 kHz resonance matching | Standing wave half-wavelength | Active | `channel_resonance_score` |
| 22 | Optical | Beer-Lambert 405 nm attenuation | Blood optical properties | Active | `blue_light_delivery_index_405nm` |

**Note:** The 1D lumped-element pipeline uses models marked Active. The cascade \
2D FVM and 3D FEM validation pipelines (`cascade_2d_3d_validation.rs`) can \
incorporate Available models for higher-fidelity spot-checks on selected designs.
"
    .to_string()
}

#[cfg(test)]
mod tests {
    use super::{build_cri_expansion_sensitivity, cavitation_regime_clause};

    #[test]
    fn cavitation_regime_clause_distinguishes_sigma_windows() {
        assert!(cavitation_regime_clause(-0.25).contains("σ = -0.2500 < 0"));
        assert!(cavitation_regime_clause(0.63).contains("(0 < σ < 1)"));
    }

    #[test]
    fn cri_expansion_sensitivity_renders_centered_html_table() {
        let html = build_cri_expansion_sensitivity();
        assert!(html.contains("<div align=\"center\">"));
        assert!(html.contains("<caption><strong>Table 5.7."));
        assert!(html.contains("<table>"));
        assert!(html.contains("</table>"));
    }
}
